rule fetch_associations_from_ebi_for_igad:
    output:
        json = "results/igad_meta/allele_investigation/ebi_associations.json",
        tsv = "results/igad_meta/allele_investigation/ebi_associations.tsv"
    localrule: True
    script: script_path("igad_meta/fetch_associations_from_ebi.py")

rule retrieve_lead_snps_from_lim_igad:
    input:
        "resources/gwas/lim-igad.tsv.gz"
    output:
        "results/igad_meta/allele_investigation/lead_snps_in_lim.tsv"
    threads: 8
    resources:
    #conda: env_path("pid_cfdr_pipeline.yaml")
    script: script_path("igad_meta/retrieve_lead_snps_from_lim_igad.R")

rule retrieve_lead_snp_or_and_alleles_for_unprocessed_igad:
    input:
        igad = "resources/gwas/igad.tsv.gz",
        lead_snps = "results/igad_meta/allele_investigation/ebi_associations.tsv",
        table_two = "resources/gwas/igad/bronson_table_2.tsv",
        supp_table_two = "resources/gwas/igad/bronson_supp_table_2.tsv",
        lim_lead_snps = "results/igad_meta/allele_investigation/lead_snps_in_lim.tsv"
    output:
        "results/igad_meta/allele_investigation/all_alleles_and_ors.tsv"
    threads: 8
    resources:
    #conda: env_path("pid_cfdr_pipeline.yaml")
    script: script_path("igad_meta/retrieve_lead_snp_or_and_alleles.R")

rule add_ref_alt_status_to_st3_allele:
    input:
        "resources/gwas/igad/bronson_supp_table_3.tsv"
    output:
        "results/igad_meta/allele_investigation/st3.tsv"
    localrule: True
    script: script_path("igad_meta/retrieve_hg38_status_for_st3_alleles.py")

rule merge_suggestive_snps_from_supp_table_3_with_unprocessed_bronson_and_lim:
    input:
        supp_table_three = "results/igad_meta/allele_investigation/st3.tsv",
        igad = "resources/gwas/igad.tsv.gz",
        lim = "results/processed_gwas/lim-igad.tsv.gz",
    output:
        "results/igad_meta/allele_investigation/st3_bronson_and_lim.tsv"
    threads: 8
    resources:
    #conda: env_path("pid_cfdr_pipeline.yaml")
    script: script_path("igad_meta/merge_st3_with_bronson_and_lim.R")

rule tabulate_st3_relationships:
    input:
        "results/igad_meta/allele_investigation/st3_bronson_and_lim.tsv"
    output:
        "results/igad_meta/allele_investigation/st3_bronson_and_lim_relationships.tsv"
    localrule: True
    #conda: env_path("pid_cfdr_pipeline.yaml")
    script: script_path("igad_meta/tabulate_st3_relationships.R")

rule handle_allele_flips_for_igad_gwas_bespokely:
    input:
        lim = "results/gwas_pipeline/reheaded/fixed_alleles/min_col_set/recalc_sumstats/liftover/realigned_alleles/lim-igad_with_allele_codes.tsv.gz",
        bronson = "results/gwas_pipeline/reheaded/fixed_alleles/min_col_set/recalc_sumstats/liftover/realigned_alleles/igad_with_allele_codes.tsv.gz"
    output:
        join = "results/igad_meta/allele_investigation/lim_bronson_conservative_join.tsv.gz",
        plot = "results/igad_meta/allele_investigation/lim_bronson_conservative_join_zscores.png"
    threads: 8
    resources:
    #conda: env_path("pid_cfdr_pipeline.yaml")
    script: script_path("igad_meta/handle_allele_flips.R")

rule plot_bronson_and_lim_zscores:
    input:
        bronson = "results/processed_gwas/bronson-uncorrected.tsv.gz",
        lim = "results/processed_gwas/lim-igad.tsv.gz"
    output:
        "results/igad_meta/allele_investigation/bronson_and_lim_zscores.png"
    params:
        suffix_one = ".bronson",
        suffix_two = ".lim"
    threads: 10
    resources:
    #conda: env_path("pid_cfdr_pipeline.yaml")
    script: script_path("igad_meta/plot_paired_zscores.R")
