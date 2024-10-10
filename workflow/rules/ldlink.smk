rule fetch_ld_matrix:
    input:
        "results/processed_gwas/{trait}/{variant_set}/{clump_threshold}/{variant_id}/sum_stats.tsv.gz"
    output:
        "results/processed_gwas/{trait}/{variant_set}/{clump_threshold}/{variant_id}/ld_matrix.tsv.gz"
    params:
        chr_col = 'CHR38',
        bp_col = 'BP38',
        snp_col = 'SNPID',
        ldlink_api_token = config['ldlink_api_token'],
        population = 'EUR',
        genome_build = 'grch38_high_coverage',
    threads: 1
    resources:
        runtime = 15
    group: "gwas"
    #conda: env_path("pid_cfdr_pipeline.yaml")
    script: script_path("ldlink/fetch_ld_matrix.R")
