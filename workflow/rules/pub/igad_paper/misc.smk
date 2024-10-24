rule compute_gif_for_igad_data:
    input:
        "results/pid/100kb/bronson-finngen-igad_with_pid_and_non_pid.tsv.gz"
    output:
        iei = "results/pub/igad_paper/misc/100kb_{mhc,with_mhc|sans_mhc}_iei_genes_bronson-finngen-igad_gif.tsv",
        non_iei = "results/pub/igad_paper/misc/100kb_{mhc,with_mhc|sans_mhc}_non_iei_genes_bronson-finngen-igad_gif.tsv",
        all_genes = "results/pub/igad_paper/misc/100kb_{mhc,with_mhc|sans_mhc}_all_genes_bronson-finngen-igad_gif.tsv"
    params:
        sans_mhc = lambda w: True if w.mhc == 'sans_mhc' else False,
        percentiles = [0.5, 0.9, 0.95, 0.99]
    threads: 8
    resources:
    localrule: True
    script: script_path("pub/igad_paper/misc/compute_gif_for_iei_genic_snps.R")

rule run_gif_permutations_and_compute_pvalue:
    input:
        iei_result = "results/pub/igad_paper/misc/{window}kb_{mhc}_iei_genes_bronson-finngen-igad_gif.tsv",
        all_genes = "results/pid/{window}kb_all_genes.rds",
        gwas = "results/igad_meta/meta.tsv.gz"
    output:
        permutations = "results/pub/igad_paper/misc/{window}kb_{mhc}_{seed}_{permutations}.tsv.gz",
        pvalue = "results/pub/igad_paper/misc/{window}kb_{mhc}_{seed}_{permutations}_pvalue.tsv"
    params:
        sans_mhc = lambda w: True if w.mhc == 'sans_mhc' else False,
        no_of_genes = 448,
        seed = lambda w: int(w.seed),
        no_of_permutations = lambda w: int(w.permutations),
        percentiles = [0.5, 0.9, 0.95, 0.99]
    threads: 20
    resources:
        runtime = 120
    script: script_path("pub/igad_paper/misc/run_gif_permutations_and_compute_pvalue.R")

rule merge_meta_cfdr_aux_for_igad_paper:
    input:
        igad = "results/igad_meta/meta.tsv.gz",
        cfdr = "results/cfdr/bronson-finngen-igad-gc_on_asthma-ex_and_ra-ishigaki_and_liu-decode-lyons-dennis-iga/left/sans_mhc/snps_only/1000kb_1_0_2/3/postscreen/cfdr.tsv.gz",
        ra = "results/processed_gwas/ra-ishigaki.tsv.gz",
        asthma = "results/processed_gwas/asthma-ex.tsv.gz",
        iga = "results/processed_gwas/liu-decode-lyons-dennis-iga.tsv.gz",
        igan = "results/processed_gwas/igan.tsv.gz",
        maf = "results/1kG/hg38/eur/snps_only/005/merged.afreq"
    output:
        "results/pub/igad_paper/misc/merged_meta_cfdr_and_aux_data.tsv.gz"
    params:
        stat_cols = ['P', 'BETA', 'SE']
    threads: 16
    resources:
        runtime = 15
    script: script_path("pub/igad_paper/misc/merge_meta_cfdr_aux.R")

rule merge_suggestive_igad_with_iga_meta:
    input:
        igad = "results/igad_meta/suggestive/lead_snps.distance_clumped.rsIDs",
        iga = "results/pub/igad_paper/tables/iga_non_mhc_lead_snps.tsv",
        lit_iga = "results/iga_meta/existing_associations.tsv"
    output:
        "results/pub/igad_paper/misc/suggestive_igad_on_iga_meta.tsv"
    params:
        window = 1e6
    localrule: True
    script: script_path("pub/igad_paper/misc/merge_suggestive_igad_with_iga_meta.R")

rule merge_suggestive_igad_with_existing_igan_associations:
    input:
        igad = "results/igad_meta/suggestive/lead_snps.distance_clumped.rsIDs",
        igan = "resources/gwas/igan/ebi_associations.tsv"
    output:
        "results/pub/igad_paper/misc/suggestive_igad_on_igan.tsv"
    params:
        window = 1e6
    localrule: True
    script: script_path("pub/igad_paper/misc/merge_suggestive_igad_with_igan.R")

rule plot_iga_lead_snp_betas:
    input:
        lead_snps = "results/pub/igad_paper/tables/iga_non_mhc_lead_snps.tsv",
    output:
        "results/pub/igad_paper/misc/iga_lead_snp_betas.png"
    threads: 10
    resources:
    script: script_path("pub/igad_paper/misc/get_iga_lead_snp_betas.R")

rule compile_stats_for_iga_igad_igan_at_gws_iga_snps:
    input:
        merged = "results/pub/igad_paper/misc/merged_meta_cfdr_and_aux_data.tsv.gz",
        iga_lead_snps = "results/pub/igad_paper/tables/iga_non_mhc_lead_snps.tsv"
    output:
        "results/pub/igad_paper/misc/stats_at_gws_iga_snps.tsv"
    threads: 8
    resources:
    localrule: True
    script: script_path("pub/igad_paper/misc/compile_stats_for_iga_igad_igan_at_gws_iga_snps.R")

rule compile_stats_for_iga_igad_igan_at_gws_igad_snps:
    input:
        igan = "results/processed_gwas/igan.tsv.gz",
        igad_lead_snps = "results/pub/igad_paper/tables/igad_lead_snps.tsv"
    output:
        "results/pub/igad_paper/misc/stats_at_gws_igad_snps.tsv"
    threads: 8
    resources:
    localrule: True
    script: script_path("pub/igad_paper/misc/compile_stats_for_iga_igad_igan_at_gws_igad_snps.R")
