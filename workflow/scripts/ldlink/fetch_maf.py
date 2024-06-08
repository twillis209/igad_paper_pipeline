import requests
import json
import pandas as pd
import re

# NB: Based on the sample script provided by Open Targets Genetics here: https://genetics-docs.opentargets.org/data-access/graphql-api#available-endpoints

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
    res = query_snp(row.SNP2.replace(':', '_'))

    if res['variantInfo'] is None or res['variantInfo']['rsId'] == 'null':
        tokens = row.SNP2.split(':')
        alt_snp = '_'.join([tokens[0], tokens[1], tokens[3], tokens[2]])
        res = query_snp(alt_snp)

    res_dict = {}

    res_dict['SNP2_ID'] = row.SNP2

    for k,v in res['variantInfo'].items():
        res_dict[k] = v

    d.append(res_dict)

daf.merge(pd.DataFrame(d), left_on = 'SNP2', right_on = 'SNP2_ID').to_csv(snakemake.output[0], sep = '\t', index = False)
