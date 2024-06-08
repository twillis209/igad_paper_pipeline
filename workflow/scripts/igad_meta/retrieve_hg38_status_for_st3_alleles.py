import requests
import json
import pandas as pd
import re

variant_query = """
    query fetch_allele_info($rsID: String!) {
        search(queryString: $rsID) {
            variants {
            rsId
            refAllele
            altAllele
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
}
"""

base_url = "https://api.genetics.opentargets.org/graphql"

daf = pd.read_csv(snakemake.input[0], delim_whitespace = True)

result_dict = {}

def query_snp(rsid):
    variables = {"rsID": rsid}

    r = requests.post(base_url, json={"query": variant_query, "variables": variables})

    variant_response_data = json.loads(r.text)['data']['search']['variants']

    if not variant_response_data:
        return {}
    elif variant_response_data[0]['rsId'] != rsid:
        raise Exception('Response rsID does not match input')
    else:
        return {'rsID': rsid,
                'REF' : variant_response_data[0]['refAllele'],
                'ALT' : variant_response_data[0]['altAllele'],
                'gnomadNFE' : variant_response_data[0]['gnomadNFE'],
                'gnomadFIN' : variant_response_data[0]['gnomadFIN'],
                'gnomadNFENWE' : variant_response_data[0]['gnomadNFENWE'],
                'gnomadNFESEU' : variant_response_data[0]['gnomadNFESEU']
            }


d = []

for index, row in daf.iterrows():
    res_dict = query_snp(row.Variant)

    if not res_dict:
        continue
    else:
        d.append(res_dict)

otg_daf = pd.DataFrame(d)

merged = daf.merge(right = otg_daf, how = 'left', left_on = 'Variant', right_on = 'rsID')

merged = merged.astype({'allele_hg38_status': str})

columns_to_round = ["gnomadNFE", "gnomadFIN", "gnomadNFENWE", "gnomadNFESEU"]

merged[columns_to_round] = merged[columns_to_round].round(2)

merged.loc[merged.Allele == merged.REF, 'allele_hg38_status'] = 'REF'
merged.loc[merged.Allele == merged.ALT, 'allele_hg38_status'] = 'ALT'

merged.to_csv(snakemake.output[0], sep = '\t', index = False)
