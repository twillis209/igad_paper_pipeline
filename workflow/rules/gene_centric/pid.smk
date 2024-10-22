# NB: In supplement to Tangye et al. (2022) IUIS paper
rule extract_gene_names_from_tangye_2022_list:
    input:
        "resources/pid_genes/tangye_table.xlsx"
    output:
        "resources/pid/pid_genes.tsv"
    #conda: env_path("pid_cfdr_pipeline.yaml")
    localrule: True
    script: script_path("gene_centric/pid/tabulate_pid_genes.R")

rule fetch_pid_gene_metadata:
    input:
        "resources/pid/pid_genes.tsv"
    output:
        "resources/pid/pid_gene_coordinates.tsv"
    localrule: True
    script: script_path("gene_centric/pid/fetch_pid_gene_metadata.py")

rule create_pid_gene_intervals:
    input:
        "resources/pid/pid_gene_coordinates.tsv"
    output:
        tsv = "results/pid/{window}kb_gene_intervals.tsv",
        rds = "results/pid/{window}kb_pid_genes.rds"
    params:
        window = lambda w: float(w.window.replace('_', '.'))*1000
    #conda: env_path("pid_cfdr_pipeline.yaml")
    localrule: True
    script: script_path("gene_centric/pid/create_pid_gene_intervals.R")

rule subset_gwas_using_pid_gene_intervals:
    input:
        gwas = "results/processed_gwas/{trait}.tsv.gz",
        gene_intervals = "results/pid/{window}kb_gene_intervals.tsv"
    output:
        "results/pid/{trait}/{window}kb/subset.tsv.gz"
    params:
        chr_col = 'CHR38',
        bp_col = 'BP38',
        ref_col = 'REF',
        alt_col = 'ALT',
        beta_col = 'BETA',
        se_col = 'SE',
        p_col = 'P'
    threads: 8
    resources:
        runtime = 30
    #conda: env_path("pid_cfdr_pipeline.yaml")
    script: script_path("gene_centric/pid/intersect_gwas_and_gene_intervals.R")

rule draw_manhattan_for_subset_gwas:
    input:
        gwas = "results/pid/{trait}/{window}kb/subset.tsv.gz"
    output:
        "results/pid/{trait}/{window}kb/manhattan.png"
    params:
        chr_col = 'CHR38',
        bp_col = 'BP38',
        p_col = 'P',
        snp_col = 'SNPID'
    threads: 8
    resources:
        runtime = 30
    #conda: env_path("pid_cfdr_pipeline.yaml")
    script: script_path("gwas/plot_gwas_manhattan.R")

rule create_all_gene_granges:
    output:
        "results/pid/{window}kb_all_genes.rds"
    params:
        window = lambda w: float(w.window.replace('_', '.'))*1000
    localrule: True
    #conda: env_path("pid_cfdr_pipeline.yaml")
    script: script_path("gene_centric/pid/create_all_gene_intervals.R")

rule join_pad_meta_on_pid_and_non_pid_granges:
    input:
        gwas = "results/pad_meta/li_left/with_ukb/with_igad/meta_prescreen.tsv.gz",
        pid_genes = "results/pid/{window}kb_pid_genes.rds",
        all_genes = "results/pid/{window}kb_all_genes.rds"
    output:
        "results/pid/{window}kb/pad_with_pid_and_non_pid.tsv.gz"
    params:
        window = lambda w: float(w.window.replace('_', '.'))*1000
    #conda: env_path("pid_cfdr_pipeline.yaml")
    script: script_path("gene_centric/pid/join_gwas_on_pid_and_non_pid_granges.R")

use rule join_pad_meta_on_pid_and_non_pid_granges as join_igad_meta_on_pid_and_non_pid_granges with:
    input:
        gwas = "results/igad_meta/meta.tsv.gz",
        pid_genes = "results/pid/{window}kb_pid_genes.rds",
        all_genes = "results/pid/{window}kb_all_genes.rds"
    output:
        "results/pid/{window}kb/bronson-finngen-igad_with_pid_and_non_pid.tsv.gz"
