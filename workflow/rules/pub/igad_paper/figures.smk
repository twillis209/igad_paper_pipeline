rule plot_annotated_igad_manhattan_for_igad_paper:
    input:
        gwas = "results/igad_meta/meta.tsv.gz",
        annotations = "results/pub/igad_paper/tables/igad_lead_snps.tsv"
    output:
        "results/pub/igad_paper/figures/annotated_igad_manhattan.png"
    threads: 8
    resources:
        runtime = 10
    localrule: True
    #conda: env_path("pid_cfdr_pipeline.yaml")
    script: script_path('pub/igad_paper/figures/plot_annotated_igad_manhattan.R')

rule plot_annotated_iga_manhattan_for_igad_paper:
    input:
        gwas = "results/iga_meta/with_decode/with_dennis/meta_prescreen.tsv.gz",
        annotations = "results/pub/igad_paper/tables/iga_lead_snps.tsv"
    output:
        "results/pub/igad_paper/figures/annotated_iga_manhattan.png"
    threads: 8
    resources:
        runtime = 10
    localrule: True
    #conda: env_path("pid_cfdr_pipeline.yaml")
    script: script_path('pub/igad_paper/figures/plot_annotated_iga_manhattan.R')

rule fine_tune_annotations_for_lyons_iga_manhattan:
    input:
        "results/processed_gwas/iga/sans_mhc/gws/lead_snps.distance_clumped.rsIDs"
    output:
        "results/pub/igad_paper/figures/lyons_iga_manhattan_annotations.tsv"
    localrule: True
    #conda: env_path("pid_cfdr_pipeline.yaml")
    script: script_path('pub/igad_paper/figures/fine_tune_annotations_for_lyons_iga_manhattan.R')

rule plot_annotated_lyons_iga_manhattan_for_igad_paper:
    input:
        gwas = "results/processed_gwas/iga.tsv.gz",
        annotations = "results/pub/igad_paper/figures/lyons_iga_manhattan_annotations.tsv"
    output:
        "results/pub/igad_paper/figures/annotated_lyons_iga_manhattan.png"
    threads: 8
    resources:
        runtime = 10
    localrule: True
    #conda: env_path("pid_cfdr_pipeline.yaml")
    script: script_path('pub/igad_paper/figures/plot_annotated_lyons_iga_manhattan.R')

rule process_igad_cfdr_lead_snps_for_manhattan_plot:
    input:
        "results/cfdr/bronson-finngen-igad-gc_on_asthma-ex_and_ra-ishigaki_and_liu-decode-lyons-dennis-iga/left/sans_mhc/snps_only/1000kb_1_0_2/3/gws/lead_snps.distance_clumped.rsIDs"
    output:
        "results/pub/igad_paper/data/igad_cfdr_lead_snps.tsv"
    localrule: True
    #conda: env_path("pid_cfdr_pipeline.yaml")
    script: script_path('pub/igad_paper/figures/process_igad_cfdr_lead_snps_for_manhattan_plot.R')

rule plot_annotated_igad_cfdr_manhattan_for_igad_paper:
    input:
        gwas = "results/cfdr/bronson-finngen-igad-gc_on_asthma-ex_and_ra-ishigaki_and_liu-decode-lyons-dennis-iga/left/sans_mhc/snps_only/1000kb_1_0_2/3/postscreen/cfdr.tsv.gz",
        annotations = "results/pub/igad_paper/tables/cfdr_lead_snps.tsv"
    output:
        "results/pub/igad_paper/figures/annotated_igad_cfdr_manhattan.png"
    threads: 8
    resources:
        runtime = 10
    localrule: True
    #conda: env_path("pid_cfdr_pipeline.yaml")
    script: script_path('pub/igad_paper/figures/plot_annotated_igad_cfdr_manhattan.R')

rule plot_igad_rg_estimates_with_sample_size:
    input:
        metadata = "resources/gwas/metadata/metadata_all_fields.tsv",
        rg = "results/ldak/ldak-thin/combined/sans_mhc/snps_only/bronson-finngen-igad_and_imds.tsv"
    output:
        "results/pub/igad_paper/figures/rg_with_sample_size.png"
    localrule: True
    #conda: env_path('pid_cfdr_pipeline.yaml')
    script: script_path('pub/igad_paper/figures/plot_rg_estimates.R')

rule draw_tnfaip3_locuszoom_plot:
    input:
        gwas = "results/pub/igad_paper/misc/tnfaip3_locus_snps_with_r2_for_{variant_id}.tsv.gz"
    output:
        "results/pub/igad_paper/figures/tnfaip3_locus_with_r2_for_{variant_id,6_137833918_C_T}.png"
    params:
        xrange = [137700000, 138000000],
        lead_snp_rsID = 'rs112920346',
        gene_name = ['TNFAIP3', 'RP11-356I2.4', 'RP11-356I2.1', 'RP11-10J5.1', 'RP11-240M16.1']
    localrule: True
    container: "docker://twillis209/r-locuszoomr:latest"
    script: script_path("pub/igad_paper/figures/plot_tnfaip3_with_r2.R")

use rule draw_tnfaip3_locuszoom_plot as draw_cd86_locuszoom_plot with:
    input:
        gwas = "results/pub/igad_paper/misc/cd86_locus_snps_with_r2_for_{variant_id}.tsv.gz"
    output:
        "results/pub/igad_paper/figures/cd86_locus_with_r2_for_{variant_id,3_122081640_A_C}.png"
    params:
        xrange = [122e6, 1222e5],
        lead_snp_rsID = 'rs9831894',
        gene_name = ['CD86', 'ILDR1']
    localrule: True

# 6:135320231:T:G is the SNP with smallest P-value present in the 1kG data set
rule draw_ahi1_tnfaip3_locuszoom_plot:
    input:
        gwas = "results/pub/igad_paper/misc/ahi1_tnfaip3_locus_snps_with_r2_for_{variant_id}.tsv.gz"
    output:
        "results/pub/igad_paper/figures/ahi1_tnfaip3_locus_with_r2_for_{variant_id,6_135320231_T_G}.png"
    params:
        xrange = [134.2e6, 139.5e6],
        lead_snp_rsid = 'rs9321500',
        lead_snp_snpid = '6:135320231:T:G',
        tnfaip3_lead_snp_rsid = 'rs112920346',
        tnfaip3_lead_snp_snpid = "6:137833918:C:T",
        gene_biotype = 'protein_coding'
    #localrule: True
    container: "docker://twillis209/r-locuszoomr:latest"
    script: script_path("pub/igad_paper/figures/plot_ahi1_tnfaip3_with_r2.R")

rule draw_cd68_locuszoom_plot:
    input:
        gwas = "results/pub/igad_paper/misc/cd68_locus_snps_with_r2_for_{variant_id}.tsv.gz"
    output:
        "results/pub/igad_paper/figures/cd68_locus_with_r2_for_{variant_id,17_7581494_G_A|17_7451821_C_A}.png"
    params:
        xrange = [745e4, 77e5],
        lead_snp_rsid = 'rs9901675',
        lead_snp_snpid = '17:7581494:G:A',
        gene_biotype = 'protein_coding'
    container: "docker://twillis209/r-locuszoomr:latest"
    script: script_path("pub/igad_paper/figures/plot_cd68_with_r2.R")

use rule draw_cd68_locuszoom_plot as draw_tnfsf13_locuszoom_plot with:
    input:
        gwas = "results/pub/igad_paper/misc/tnfsf13_locus_snps_with_r2_for_{variant_id}.tsv.gz"
    output:
        "results/pub/igad_paper/figures/tnfsf13_locus_with_r2_for_{variant_id,17_7559652_A_G}.png"
    params:
        xrange = [745e4, 77e5],
        lead_snp_rsid = 'rs3803800',
        lead_snp_snpid = '17:7559652:A:G',
        gene_biotype = 'protein_coding'

use rule draw_cd68_locuszoom_plot as draw_inava_locuszoom_plot with:
    input:
        gwas = "results/pub/igad_paper/misc/inava_locus_snps_with_r2_for_{variant_id}.tsv.gz"
    output:
        "results/pub/igad_paper/figures/inava_locus_with_r2_for_{variant_id}.png"
    params:
        xrange = [200.8e6, 201.1e6],
        lead_snp_rsid = 'rs7522462',
        lead_snp_snpid = '1:200912467:G:A',
        gene_biotype = 'protein_coding'

rule draw_iei_qqplot:
    input:
        "results/pid/100kb/bronson-finngen-igad_with_pid_and_non_pid.tsv.gz"
    output:
        "results/pub/igad_paper/figures/100kb_iei_qqplot.png"
    threads: 8
    resources:
        runtime = 10
    #conda: env_path("pid_cfdr_pipeline.yaml")
    script: script_path("pub/igad_paper/figures/plot_iei_qqplot.R")

rule draw_rg_gps_plot_for_igad_paper:
    input:
        rg = "results/ldak/ldak-thin/combined/sans_mhc/snps_only/bronson-finngen-igad_and_imds.tsv",
        gps = "results/gps/combined/bronson-finngen-igad/inner/sans_mhc/snps_only/1000kb_1_0_8/3000_draws/pvalues.tsv",
        metadata = "resources/gwas/metadata/metadata_all_fields.tsv"
    output:
        "results/pub/igad_paper/figures/rg_gps.png"
    params:
        imd_traits = config.get('igad_paper').get('imd_traits')
    localrule: True
    #conda: env_path("pid_cfdr_pipeline.yaml")
    script: script_path("pub/igad_paper/figures/draw_rg_gps_plot.R")


rule draw_rg_plots_for_igad_and_iga:
    input:
        igad = "results/ldak/ldak-thin/combined/sans_mhc/snps_only/bronson-finngen-igad_and_imds.tsv",
        imds = "results/ldak/ldak-thin/combined/sans_mhc/snps_only/imds.tsv",
        metadata = "resources/gwas/metadata/metadata_all_fields.tsv"
    output:
        "results/pub/igad_paper/figures/rg_igad_iga.png"
    params:
        imd_traits = config.get('igad_paper').get('imd_traits')
    localrule: True
    #conda: env_path("pid_cfdr_pipeline.yaml")
    script: script_path("pub/igad_paper/figures/draw_rg_plots_for_igad_and_iga.R")

rule draw_imd_rg_plot_for_igad_paper:
    input:
        rg = "results/pub/igad_paper/export_data/rg_estimates.tsv",
        metadata = "resources/gwas/metadata/metadata_all_fields.tsv"
    output:
        pheatmap = "results/pub/igad_paper/figures/imd_rg_pheatmap.png",
        ggcorr = "results/pub/igad_paper/figures/imd_rg_ggcorr.png"
    params:
        no_of_clusters = 3,
        traits = ['bronson-finngen-igad']+config.get('igad_paper').get('imd_traits_rg_with_igad')
    localrule: True
    #conda: env_path("pid_cfdr_pipeline.yaml")
    script: script_path("pub/igad_paper/figures/draw_imd_rg_heatmap.R")

rule plot_igad_and_igan_zscores:
    input:
        "results/merged_gwas/bronson-finngen-igad_and_igan/inner/sans_mhc/merged.tsv.gz"
    output:
        "results/pub/igad_paper/figures/zscore_plots/bronson-finngen-igad_and_igan.png"
    params:
        xlab = "SIgAD",
        ylab = "IgAN"
    threads: 8
    resources:
        runtime = 10
    #conda: env_path("pid_cfdr_pipeline.yaml")
    script: script_path("pub/igad_paper/figures/plot_zscores.R")

use rule plot_igad_and_igan_zscores as plot_igad_and_iga_zscores with:
    input:
        "results/merged_gwas/bronson-finngen-igad_and_liu-decode-lyons-dennis-iga/inner/sans_mhc/merged.tsv.gz"
    output:
        "results/pub/igad_paper/figures/zscore_plots/bronson-finngen-igad_and_liu-decode-lyons-dennis-iga.png"
    params:
        xlab = "SIgAD",
        ylab = "serum IgA"

use rule plot_igad_and_igan_zscores as plot_igad_and_iga_zscores_with_cfdr_labels with:
    input:
        "results/merged_gwas/bronson-finngen-igad_and_liu-decode-lyons-dennis-iga/inner/sans_mhc/merged.tsv.gz",
        labels = "results/pub/igad_paper/tables/non_mhc_lead_snps.tsv"
    output:
        "results/pub/igad_paper/figures/zscore_plots/bronson-finngen-igad_and_liu-decode-lyons-dennis-iga_with_labels.png"
    params:
        xlab = "SIgAD",
        ylab = "serum IgA"

rule plot_iga_igad_igan_at_gws_iga_snps:
    input:
        "results/pub/igad_paper/misc/stats_at_gws_iga_snps.tsv"
    output:
        igan = "results/pub/igad_paper/figures/iga_vs_igan_betas.png",
        igan_no_ci = "results/pub/igad_paper/figures/iga_vs_igan_betas_no_ci.png",
        igad = "results/pub/igad_paper/figures/iga_vs_igad_betas.png",
        igad_no_ci = "results/pub/igad_paper/figures/iga_vs_igad_betas_no_ci.png",
        igad_no_ci_del_outlier = "results/pub/igad_paper/figures/iga_vs_igad_betas_no_ci_del_outlier.png"
    localrule: True
    #conda: env_path("pid_cfdr_pipeline.yaml")
    script: script_path("pub/igad_paper/figures/plot_iga_igad_igan_at_gws_iga_snps.R")

rule plot_iga_igad_igan_at_gws_igad_snps:
    input:
        "results/pub/igad_paper/misc/stats_at_gws_igad_snps.tsv"
    output:
        igan = "results/pub/igad_paper/figures/igad_vs_igan_betas.png",
        igan_no_ci = "results/pub/igad_paper/figures/igad_vs_igan_betas_no_ci.png",
        iga = "results/pub/igad_paper/figures/igad_vs_iga_betas.png",
        iga_no_ci = "results/pub/igad_paper/figures/igad_vs_iga_betas_no_ci.png",
        iga_no_ci_del_outlier = "results/pub/igad_paper/figures/igad_vs_iga_betas_no_ci_del_outlier.png"
    localrule: True
    #conda: env_path("pid_cfdr_pipeline.yaml")
    script: script_path("pub/igad_paper/figures/plot_iga_igad_igan_at_gws_igad_snps.R")

rule all_igad_paper_figures:
    input:
        # Figure 1
        "results/pub/igad_paper/figures/annotated_igad_manhattan.png",
        # Figure 2
        "results/pub/igad_paper/figures/100kb_iei_qqplot.png",
        # Figure 3
        "results/pub/igad_paper/figures/annotated_igad_cfdr_manhattan.png",
        # Supplementary Figure 1
        "results/pub/igad_paper/figures/tnfaip3_locus_with_r2_for_6_137833918_C_T.png",
        # Supplementary Figure 2
        "results/pub/igad_paper/figures/ahi1_tnfaip3_locus_with_r2_for_6_135320231_T_G.png",
        # Supplementary Figure 3
        "results/pub/igad_paper/figures/cd86_locus_with_r2_for_3_122081640_A_C.png",
        # Supplementary Figure 4
        "results/pub/igad_paper/figures/annotated_lyons_iga_manhattan.png"
        # Supplementary Figure 5
        "results/pub/igad_paper/figures/annotated_iga_manhattan.png",
        # Supplementary Figure 6
        "results/pub/igad_paper/figures/rg_gps.png",
        # Supplementary Figure 7
        "results/pub/igad_paper/misc/mr/igan_and_iga_ivw.png",
        # Supplementary Figure 8
        "results/pub/igad_paper/misc/mr/igad_and_iga_ivw.png",
        # Supplementary Figure 9
        "results/pub/igad_paper/figures/imd_rg_pheatmap.png",
        # Supplementary Figure 10
        "results/igad_meta/allele_investigation/bronson_and_lim_zscores.png"
