import requests
import json
import pandas as pd
import re

print('Running...')

variant_query = """
    query getSnpCoordinates($inputVariantId: String!) {
        search(queryString: $inputVariantId) {
        variants {
            chromosome
            position
            refAllele
            altAllele
            }
        }
    }
"""

base_url = "https://api.genetics.opentargets.org/graphql"

daf = pd.read_csv(snakemake.input[0])

daf.drop(['id', 'pos', 'chromosome'], axis = 1, inplace = True)

def query_snp(rs_id):
    r = requests.post(base_url, json={"query": variant_query, "variables": {'inputVariantId': rs_id}})

    try:
        variant_response_data = json.loads(r.text)['data']['search']['variants'][0]
    except KeyError:
        return {'chromosome': None, 'position': None, 'refAllele': None, altAllele: None}
    except IndexError:
        return {'chromosome': None, 'position': None, 'refAllele': None, altAllele: None}

    return variant_response_data

result_dict = {}

for index, row in daf.iterrows():
    print(index)
    result_dict[row.rsid] = query_snp(row.rsid)

res_daf = pd.DataFrame(result_dict)

res_daf.to_csv(snakemake.output[0], sep = '\t', index = False)
