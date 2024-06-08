import pandas as pd
import re
import concurrent.futures
import requests
import json
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry

# Request handling from: https://stackoverflow.com/a/70336440/22252122

retry_strategy = Retry(
    total=3,
    backoff_factor=1,
    status_forcelist=[500, 502, 503, 504],
    allowed_methods=["HEAD", "GET", "OPTIONS"]
)
adapter = HTTPAdapter(max_retries=retry_strategy)
http = requests.Session()
http.mount("https://", adapter)
http.mount("http://", adapter)

gene_daf = pd.read_csv(snakemake.input[0], sep = '\t', header = 0)

def query_gene_with_ensembl(ensembl_id):
    base_url = "https://rest.ensembl.org"
    endpoint = f"/lookup/id/{ensembl_id}?expand=1"

    headers = {
        "Content-Type": "application/json"
    }


    response = http.get(base_url + endpoint, headers=headers)

    if response.status_code == 200:
        gene_data = response.json()
        if "seq_region_name" in gene_data and "start" in gene_data and "end" in gene_data:
            print(ensembl_id)
            chromosome = gene_data["seq_region_name"]
            start = gene_data["start"]
            end = gene_data["end"]
            strand = gene_data['strand']
            symbol = gene_data['display_name'] if 'display_name' in gene_data else None
            tss = start if strand == 1 else end

            return pd.DataFrame([{'id': ensembl_id, 'symbol': symbol, 'chromosome': chromosome, 'start': start, 'end': end, 'tss': tss, 'fwdStrand': strand == 1}])
        else:
            return None
    else:
        print("Error:", response.status_code)
        return None

with concurrent.futures.ThreadPoolExecutor(max_workers = snakemake.threads) as executor:
    d = list(executor.map(query_gene_with_ensembl, gene_daf['Ensembl ID'].tolist()))

daf = pd.concat(d)

#daf.rename(columns = {'symbol': 'gene_name'}, inplace = True)

daf.to_csv(snakemake.output[0], sep = '\t', index = False)
