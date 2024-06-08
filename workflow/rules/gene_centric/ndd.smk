rule download_ndd_genes:
    output:
        "resources/ndd/genetrek_genes.tsv"
    localrule: True
    shell:
        "wget -O {output} https://genetrek.pasteur.fr/downloadAllData?filetype=tsv"

rule tabulate_ndd_genes:
    input:
        "resources/ndd/genetrek_genes.tsv"
    output:
        "resources/ndd/ndd_genes.tsv"
    localrule: True
    conda: env_path("pid_cfdr_pipeline.yaml")
    script: script_path("gene_centric/ndd/tabulate_ndd_genes.R")

rule fetch_ndd_gene_metadata:
    input:
        "resources/ndd/ndd_genes.tsv"
    output:
        "resources/ndd/ndd_gene_coordinates.tsv"
    threads: 8
    resources:
        runtime = 90
    script: script_path("gene_centric/ndd/fetch_ndd_gene_metadata.py")

rule downsample_ndd_genes:
    input:
        "resources/ndd/ndd_gene_coordinates.tsv"
    output:
        "resources/ndd/{no_of_genes}_downsampled_ndd_genes_seed_{seed}.tsv"
    params:
        seed = lambda w: int(w.seed),
        n = lambda w: int(w.no_of_genes)
    localrule: True
    run:
        daf = pd.read_csv(input[0], sep = '\t', header = 0)

        exc_chrom = ['X', 'MT', 'Y']

        daf = daf.query('chromosome not in @exc_chrom')

        daf = daf.sample(params.n, random_state = params.seed)

        daf.to_csv(output[0], index = False, sep = '\t')

use rule create_pid_gene_intervals as create_ndd_gene_intervals with:
    input:
        "resources/ndd/{no_of_genes}_downsampled_ndd_genes_seed_{seed}.tsv"
    output:
        "results/ndd/{no_of_genes}_{window}kb_gene_intervals_seed_{seed}.tsv"
    localrule: True
