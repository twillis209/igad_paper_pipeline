import requests
import json
import pandas as pd
import re

chr_col = snakemake.params.chr_col
bp_col = snakemake.params.bp_col
snp_col = snakemake.params.snp_col
ref_col = snakemake.params.ref_col
alt_col = snakemake.params.alt_col

variant_query = """
    query annotateLeadSnp($inputVariantId: String!) {
        variantInfo(variantId: $inputVariantId) {
            rsId
            gnomadAFR
            gnomadAMR
            gnomadASJ
            gnomadEAS
            gnomadFIN
            gnomadNFE
            gnomadNFEEST
            gnomadNFENWE
            gnomadNFESEU
            }
        }
"""

base_url = "https://api.genetics.opentargets.org/graphql"

daf = pd.read_csv(snakemake.input[0], delim_whitespace = True)

result_dict = {}

def query_snp(snp_id):
    r = requests.post(base_url, json={"query": variant_query, "variables": {'inputVariantId': snp_id}})

    variant_response_data = json.loads(r.text)['data']['variantInfo']

    return {'variantInfo' : variant_response_data}

d = []

for index, row in daf.iterrows():
    res = query_snp('_'.join([str(row.CHR38), str(row.BP38), row.REF, row.ALT]))

    if res['variantInfo'] is None or res['variantInfo']['rsId'] == 'null':
        res = query_snp('_'.join([str(row.CHR38), str(row.BP38), row.ALT, row.REF]))

    res_dict = {}

    res_dict['SNPID'] = row.SNPID

    for k,v in res['variantInfo'].items():
        res_dict[k] = v

    d.append(res_dict)

daf.merge(pd.DataFrame(d), left_on = 'SNPID', right_on = 'SNPID').to_csv(snakemake.output[0], sep = '\t', index = False)
