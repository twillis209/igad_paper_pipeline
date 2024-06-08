rule run_bronson_and_finngen_igad_meta_analysis:
    input:
        bronson = "results/processed_gwas/igad.tsv.gz",
        finngen = "results/processed_gwas/finngen-igad.tsv.gz",
        ukb = "results/processed_gwas/ukb-igad.tsv.gz"
    output:
        without_ukb = "results/igad_meta/meta.tsv.gz",
        with_ukb = "results/igad_meta/with_ukb/meta.tsv.gz"
    params:
        chr_col = 'CHR38',
        bp_col = 'BP38',
        ref_col = 'REF',
        alt_col = 'ALT',
        beta_col = 'BETA',
        se_col = 'SE',
        p_col = 'P'
    threads: 16
    resources:
        runtime = 20
    group: "gwas"
    conda: env_path("pid_cfdr_pipeline.yaml")
    script: script_path("igad_meta/run_meta_analysis.R")

rule copy_bronson_finngen_meta_analysis_to_results:
    input:
        "results/igad_meta/meta.tsv.gz"
    output:
        "results/processed_gwas/bronson-finngen-igad.tsv.gz"
    localrule: True
    shell: "cp {input} {output}"

rule apply_genomic_control_to_igad_meta_for_cfdr:
    input:
        gwas = "results/igad_meta/meta.tsv.gz",
        gif = "results/igad_meta/gif.tsv"
    output:
        "results/igad_meta/meta_gc.tsv.gz"
    params:
        chr_col = 'CHR38',
        bp_col = 'BP38',
        ref_col = 'REF',
        alt_col = 'ALT',
        snp_col = 'SNPID',
        beta_col = 'BETA',
        se_col = 'SE',
        p_col = 'P'
    threads: 16
    resources:
        runtime = 20
    group: "gwas"
    conda: env_path("pid_cfdr_pipeline.yaml")
    script: script_path("igad_meta/apply_gc_to_meta.R")

rule copy_gc_bronson_finngen_meta_analysis_to_results:
    input:
        "results/igad_meta/meta_gc.tsv.gz"
    output:
        "results/processed_gwas/bronson-finngen-igad-gc.tsv.gz"
    localrule: True
    shell: "cp {input} {output}"

rule join_all_igad_data_sets:
    input:
        bronson = "results/processed_gwas/igad.tsv.gz",
        finngen = "results/processed_gwas/finngen-igad.tsv.gz",
        meta = "results/igad_meta/meta.tsv.gz"
    output:
        "results/igad_meta/all_sum_stats.tsv.gz"
    params:
        chr_col = 'CHR38',
        bp_col = 'BP38',
        ref_col = 'REF',
        alt_col = 'ALT',
        beta_col = 'BETA',
        se_col = 'SE',
        p_col = 'P'
    threads: 16
    resources:
        runtime = 20
    group: "gwas"
    conda: env_path("pid_cfdr_pipeline.yaml")
    script: script_path("igad_meta/join_all_igad_data.R")

rule distance_clump_igad_meta:
    input:
        gwas = "results/igad_meta/meta.tsv.gz"
    output:
        "results/igad_meta/{threshold}/lead_snps.distance_clumped"
    params:
        mhc = lambda wildcards: False, #if wildcards.snp_set == 'sans_mhc' else True,
        chr_col = 'CHR38',
        bp_col = 'BP38',
        ref_col = 'REF',
        alt_col = 'ALT',
        snp_col = 'SNPID',
        p_col = 'P',
        beta_col = 'BETA',
        se_col = 'SE',
        index_threshold = lambda wildcards: 5e-8 if wildcards.threshold == 'gws' else 1e-5,
        distance_window = 2e6,
    threads: 16
    resources:
        runtime = 5
    group: "gwas"
    conda: env_path("pid_cfdr_pipeline.yaml")
    script: script_path("gwas/distance_clump.R")

use rule annotate_lead_snps as annotate_igad_meta_lead_snps with:
    input:
        "results/igad_meta/{threshold}/lead_snps.distance_clumped"
    output:
        annotations = "results/igad_meta/{threshold}/lead_snps.distance_clumped.annotations",
        rsIDs = "results/igad_meta/{threshold}/lead_snps.distance_clumped.rsIDs"
    localrule: True

use rule draw_manhattan_with_lead_snp_annotation as draw_igad_meta_manhattan_with_lead_snp_annotation with:
    input:
        gwas = "results/igad_meta/meta.tsv.gz",
        rsIDs = "results/igad_meta/{threshold}/lead_snps.distance_clumped.rsIDs"
    output:
        "results/igad_meta/{threshold}/annotated_manhattan.distance_clumped.png"
    threads: 8

rule collate_igad_meta_lead_snp_sum_stats:
    input:
        lead_snps = "results/igad_meta/{threshold}/lead_snps.distance_clumped",
        sum_stats = "results/igad_meta/all_sum_stats.tsv.gz"
    output:
        "results/igad_meta/{threshold}/lead_snps_sum_stats.tsv.gz"
    params:
        snp_col = 'SNPID'
    threads: 16
    resources:
        runtime = 5
    group: "gwas"
    conda: env_path("pid_cfdr_pipeline.yaml")
    script: script_path("gwas/collate_lead_snp_sum_stats.R")

use rule compute_genomic_inflation_factor as compute_genomic_inflation_factors_for_igad_meta with:
    input:
        "results/igad_meta/meta.tsv.gz"
    output:
        "results/igad_meta/{variant_set}_gif.tsv"
    params:
        chr_col = 'CHR38',
        bp_col = 'BP38',
        p_col = 'P',
        percentiles = [10, 20, 50, 75],
        controls = lambda w: get_metadata_field('bronson-finngen-igad', 'N0'),
        cases = lambda w: get_metadata_field('bronson-finngen-igad', 'N1'),

use rule draw_qqplot as draw_igad_meta_qqplot with:
    input:
        gwas = "results/igad_meta/meta.tsv.gz",
        gif = "results/igad_meta/gif.tsv"
    output:
        "results/igad_meta/meta_qqplot.png"

use rule subset_summary_statistics_about_variant as subset_summary_statistics_about_variant_for_igad_meta with:
    input:
        "results/igad_meta/meta.tsv.gz"
    output:
        sum_stats = temp("results/igad_meta/{threshold}/{variant_id}/{window_size}/sum_stats.tsv.gz"),
        ids = temp("results/igad_meta/{threshold}/{variant_id}/{window_size}/ids.txt")

rule fetch_rsids_for_sum_stats_about_igad_meta_variant:
    input:
        "results/igad_meta/{threshold}/{variant_id}/{window_size}/sum_stats.tsv.gz"
    output:
        temp("results/igad_meta/{threshold}/{variant_id}/{window_size}/sum_stats_with_rsids.tsv.gz")
    params:
        chr_col = 'CHR38',
        bp_col = 'BP38',
        snp_col = 'SNPID',
        ref_col = 'REF',
        alt_col = 'ALT'
    resources:
        runtime = 15
    script: script_path("ldlink/fetch_rsids.py")

use rule annotate_lead_snps as lookup_ifih1_variants with:
    input:
        "results/igad_meta/gws/2_162267541_C_T/100kb/sum_stats.tsv.gz"
    output:
        annotations = "results/igad_meta/gws/2_162267541_C_T/100kb/snps.annotations",
        rsIDs = "results/igad_meta/gws/2_162267541_C_T/100kb/snps.rsIDs"
    localrule: True

rule plot_ahi1_tnfaip3_locus:
    input:
        "results/igad_meta/meta.tsv.gz"
    output:
        "results/igad_meta/gws/ahi1_tnfaip3.png"
    params:
        chrom = 6,
        start_bp = 133000000,
        stop_bp = 140000000
    threads: 8
    resources:
        runtime = 20
    group: "gwas"
    conda: env_path("locuszoomr.yaml")
    script: script_path("igad_meta/plot_ahi1_tnfaip3.R")

use rule distance_clump_igad_meta as distance_clumped_igad_meta_with_ukb with:
    input:
        gwas = "results/igad_meta/with_ukb/meta.tsv.gz"
    output:
        "results/igad_meta/with_ukb/{threshold}/lead_snps.distance_clumped"

use rule annotate_igad_meta_lead_snps as annotate_igad_meta_with_ukb_lead_snps with:
    input:
        "results/igad_meta/with_ukb/{threshold}/lead_snps.distance_clumped"
    output:
        annotations = "results/igad_meta/with_ukb/{threshold}/lead_snps.distance_clumped.annotations",
        rsIDs = "results/igad_meta/with_ukb/{threshold}/lead_snps.distance_clumped.rsIDs"

use rule draw_igad_meta_manhattan_with_lead_snp_annotation as draw_igad_meta_with_ukb_manhattan_with_lead_snp_annotation with:
    input:
        gwas = "results/igad_meta/with_ukb/meta.tsv.gz",
        rsIDs = "results/igad_meta/with_ukb/{threshold}/lead_snps.distance_clumped.rsIDs"
    output:
        "results/igad_meta/with_ukb/{threshold}/annotated_manhattan.distance_clumped.png"
