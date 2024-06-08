rule tabulate_metabolic_genes:
    input:
        "resources/metab/metabolic_genes.xlsx"
    output:
        "resources/metab/metabolic_genes.tsv"
    params:
        cols = ['ENSEMBL ID', 'GENE NAME', 'PATHWAY'],
        renamed_cols = {'ENSEMBL ID': 'Ensembl ID', 'GENE NAME': 'gene', 'PATHWAY': 'pathway'}
    localrule: True
    run:
        excel_file = pd.ExcelFile(input[0])

        dafs = []

        for sheet_name in excel_file.sheet_names:
            daf = excel_file.parse(sheet_name)
            daf = daf[params.cols]
            daf.rename(params.renamed_cols, axis = 'columns', inplace = True)
            daf['sheet'] = sheet_name
            dafs.append(daf)

        daf = pd.concat(dafs, ignore_index=True)
        daf.dropna().to_csv(output[0], sep = '\t', index = False)

use rule fetch_ndd_gene_metadata as fetch_metabolic_gene_metadata with:
    input:
        "resources/metab/metabolic_genes.tsv"
    output:
        "resources/metab/metabolic_gene_coordinates.tsv"
        
use rule create_pid_gene_intervals as create_metab_gene_intervals with:
    input:
        "resources/metab/metabolic_gene_coordinates.tsv"
    output:
        "results/metab/{window}kb_gene_intervals.tsv"
