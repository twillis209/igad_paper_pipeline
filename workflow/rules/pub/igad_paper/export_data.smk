rule compile_rg_estimates_for_igad_paper:
    input:
        imds = "results/ldak/ldak-thin/combined/sans_mhc/snps_only/imds.tsv",
        igad = "results/ldak/ldak-thin/combined/sans_mhc/snps_only/bronson-finngen-igad_and_imds.tsv"
    output:
        "results/pub/igad_paper/export_data/rg_estimates.tsv"
    params:
        traits = ['bronson-finngen-igad']+config.get('igad_paper').get('imd_traits')
    localrule: True
    run:
        traits = params.traits

        imds = pd.read_csv(input.imds, sep = '\t')
        imds.rename(columns = {'trait.A': 'trait_A', 'trait.B': 'trait_B'}, inplace = True)

        igad = pd.read_csv(input.igad, sep = '\t')
        igad.rename(columns = {'trait.A': 'trait_A', 'trait.B': 'trait_B'}, inplace = True)

        imds = imds.query('trait_A in @traits and trait_B in @traits')
        igad = igad.query('trait_A in @traits and trait_B in @traits')

        pd.concat([imds, igad]).to_csv(output[0], sep = '\t')

rule write_out_iei_gene_list:
    input:
        "resources/pid/pid_gene_coordinates.tsv"

rule temp_rule_for_all_igad_paper_artifacts:
    input:
        # Table 1
        "results/pub/igad_paper/tables/igad_lead_snps.tsv",
        # Table 2
        "results/pub/igad_paper/tables/cfdr_lead_snps.tsv",
        # Supplementary Table 2
        "results/pub/igad_paper/tables/lyons_iga_lead_snps.tsv",
        # Supplementary Table 3
        "results/pub/igad_paper/tables/iga_non_mhc_lead_snps.tsv",
        # Supplementary Table 4
        "results/pub/igad_paper/tables/lead_snps_with_other_associations.tsv",
        # Supplementary Table 7
        "results/pub/igad_paper/misc/stats_at_gws_iga_snps.tsv",
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
        "results/pub/igad_paper/figures/annotated_lyons_iga_manhattan.png",
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
        "results/igad_meta/allele_investigation/bronson_and_lim_zscores.png",
        "results/pub/igad_paper/misc/merged_meta_cfdr_and_aux_data.tsv.gz",
        "results/igad_meta/meta.tsv.gz",
        "results/iga_meta/with_decode/with_dennis/meta_prescreen.tsv.gz",
        "results/processed_gwas/iga.tsv.gz"
