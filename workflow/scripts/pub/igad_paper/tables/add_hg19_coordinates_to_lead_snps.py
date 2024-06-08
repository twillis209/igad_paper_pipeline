import requests
import json
import pandas as pd
import re

variant_query = """
    query getSnpCoordinates($inputVariantId: String!) {
        search(queryString: $inputVariantId) {
        variants {
            chromosomeB37
            positionB37
            refAllele
            altAllele
            }
        }
    }
"""

base_url = "https://api.genetics.opentargets.org/graphql"

daf = pd.read_csv(snakemake.input[0], sep = '\t')

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
    result_dict[row.rsID] = query_snp(row.rsID)

res_daf = pd.DataFrame.from_dict(result_dict, orient = 'index')

res_daf.rename({'chromosomeB37': 'CHR19', 'positionB37': 'BP19', 'refAllele': 'REF19', 'altAllele': 'ALT19'}, axis = 1, inplace = True)
res_daf['rsID'] = res_daf.index

daf = daf.merge(res_daf, on = 'rsID')

daf.to_csv(snakemake.output[0], sep = '\t', index = False)
