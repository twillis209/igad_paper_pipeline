rule process_igad_snp_loci_table_for_google_sheets:
    input:
        gwas_lead_snps = "results/igad_meta/gws/lead_snps.distance_clumped.rsIDs",
        cfdr_lead_snps = "results/cfdr/bronson-finngen-igad-gc_on_asthma-ex_and_ra-ishigaki_and_liu-decode-lyons-dennis-iga/left/sans_mhc/snps_only/1000kb_1_0_2/3/gws/lead_snps.distance_clumped.rsIDs",
        merged_data = "results/pub/igad_paper/misc/merged_meta_cfdr_and_aux_data.tsv.gz",
        genes = "results/pid/0kb_all_genes.rds"
    output:
        all = "results/pub/igad_paper/tables/non_mhc_lead_snps.tsv",
    params:
        known_gwas_genes = config.get('known_igad_gwas_loci'),
        novel_cfdr = config.get('igad').get('novel_cfdr'),
        iei_genes = config.get('igad').get('iei_genes'),
        top_gene_to_chosen_gene = config.get('igad').get('top_gene_to_chosen_gene'),
        window = 2e5,
        full_precision = True
    localrule: True
    script: script_path('pub/igad_paper/tables/process_snp_table_for_google_sheets.R')

rule add_hg19_to_google_sheet_igad_lead_snps:
    input:
        "results/pub/igad_paper/tables/non_mhc_lead_snps.tsv"
    output:
        "results/pub/igad_paper/tables/non_mhc_lead_snps_with_hg19.tsv"
    localrule: True
    script: script_path('pub/igad_paper/tables/add_hg19_coordinates_to_lead_snps.py')

rule merge_lead_snps_with_unprocessed_gwas:
    input:
        bronson = "resources/gwas/igad.tsv.gz",
        finngen = "resources/gwas/finngen-igad.tsv.gz",
        liu = "resources/gwas/liu-decode-iga.tsv.gz",
        dennis = "resources/gwas/dennis-iga.tsv.gz",
        lyons = "resources/gwas/iga.tsv.gz",
        ra = "resources/gwas/ra-ishigaki.tsv.gz",
        asthma = "resources/gwas/asthma-ex.tsv.gz",
        lead_snps = "results/pub/igad_paper/tables/non_mhc_lead_snps_with_hg19.tsv"
    output:
        "results/pub/igad_paper/tables/non_mhc_lead_snps_with_raw_data.tsv"
    threads: 8
    resources:
        runtime = 15
    script: script_path('pub/igad_paper/tables/merge_with_unprocessed_gwas.R')

rule compare_lead_snps_with_unprocessed_gwas:
    input:
        "results/pub/igad_paper/tables/non_mhc_lead_snps_with_raw_data.tsv"
    output:
        "results/pub/igad_paper/tables/non_mhc_lead_snp_allele_sign_concordance.tsv"
    localrule: True
    script: script_path('pub/igad_paper/tables/compare_with_unprocessed_gwas.R')

rule print_meta_analysis_gwas_metadata_table:
    input:
        "resources/gwas/metadata/metadata_all_fields.tsv"
    output:
        tex = "results/pub/igad_paper/tables/meta_analysis_gwas_metadata.tex",
        tsv = "results/pub/igad_paper/tables/meta_analysis_gwas_metadata.tsv"
    params:
        abbrvs = ['asthma-ex', 'ra-ishigaki', 'igad', 'finngen-igad', 'liu-decode-iga', 'dennis-iga', 'iga'],
        caption = "The GWAS data sets included in the SIgAD and IgA meta-analyses, and SIgAD cFDR analysis. Pan-UKB is the Pan-UK Biobank study.",
        label = "meta_analysis_gwas_metadata"
    localrule: True
    script: script_path('pub/igad_paper/tables/print_gwas_metadata_table.R')

use rule add_hg19_to_google_sheet_igad_lead_snps as add_hg19_to_iga_lead_snps with:
    input:
        "results/iga_meta/with_decode/with_dennis/prescreen/gws/lead_snps.distance_clumped.rsIDs"
    output:
        "results/pub/igad_paper/tables/non_mhc_iga_lead_snps_with_hg19.tsv"

rule merge_iga_lead_snps_with_ieis:
    input:
        gwas = "results/pub/igad_paper/tables/non_mhc_iga_lead_snps_with_hg19.tsv",
        genes = "resources/pid/pid_gene_coordinates.tsv"
    output:
        "results/pub/igad_paper/tables/non_mhc_iga_lead_snps_with_hg19_and_ieis.tsv"
    params:
        window = 1e5
    localrule: True
    script: script_path('pub/igad_paper/tables/merge_lead_snps_with_pid_genes.R')

rule process_iga_snp_loci_table_for_google_sheets:
    input:
        gwas_lead_snps = "results/pub/igad_paper/tables/non_mhc_iga_lead_snps_with_hg19_and_ieis.tsv",
        existing_associations = "results/iga_meta/existing_associations.tsv",
        merged_data = "results/iga_meta/with_decode/with_dennis/all_sum_stats.tsv.gz",
    output:
        "results/pub/igad_paper/tables/iga_non_mhc_lead_snps.tsv",
    params:
        known_gwas_genes = config.get('iga').get('known_iga_gwas_loci'),
        iei_genes = config.get('iga').get('iei_genes'),
        window = 2e5,
    threads: 8
    resources:
    localrule: True
    script: script_path('pub/igad_paper/tables/process_iga_snp_table_for_google_sheets.R')

rule print_rg_gwas_metadata_table:
    input:
        "resources/gwas/metadata/metadata_all_fields.tsv"
    output:
        tex = "results/pub/igad_paper/tables/rg_gwas_metadata.tex",
        tsv = "results/pub/igad_paper/tables/rg_gwas_metadata.tsv"
    params:
        abbrvs = config.get('igad_paper').get('imd_traits'),
        caption = "The GWAS data sets with which we measured SIgAD's genetic similarity. IMSGC is the International Multiple Sclerosis Genetics Consortium. Pan-UKB is the Pan-UK Biobank study.",
        label = "rg_gwas_metadata"
    localrule: True
    script: script_path('pub/igad_paper/tables/print_rg_gwas_metadata_table.R')

rule print_lead_snp_table:
    input:
        "results/pub/igad_paper/tables/non_mhc_lead_snps.tsv",
    output:
        "results/pub/igad_paper/tables/lead_snps.tex"
    params:
        caption = "Lead SNPs from genome-wide significant associations in the GWAS meta-analysis and cFDR analysis. `Novel' indicates whether an association with SIgAD has previously been reported for a SNP. `IEI gene' indicates whether the SNP is located in, near, or is otherwise associated with a gene known to harbour variants causal for IEIs. `Analysis' indicates whether the variant was first identified in the GWAS meta-analysis or the cFDR analysis. `OR' is odds ratio.",
        label = "lead_snp_table"
    localrule: True
    script: script_path('pub/igad_paper/tables/print_lead_snp_table.R')

rule print_igad_lead_snp_table:
    input:
        "results/pub/igad_paper/tables/non_mhc_lead_snps.tsv",
    output:
        tex = "results/pub/igad_paper/tables/igad_lead_snps.tex",
        tsv = "results/pub/igad_paper/tables/igad_lead_snps.tsv"
    params:
        caption = "Lead SNPs from genome-wide significant associations in the SIgAD GWAS meta-analysis. The `Variant' column gives the rsID of each SNP, and the reference and effect alleles separated by `$>$'. `Gene' gives the gene(s) with the most evidence linking it/them to the association signal. `Novel' indicates whether an association with SIgAD has previously been reported for a SNP. `IEI gene' indicates whether the SNP is located in, near, or is otherwise associated with a gene known to harbour variants causal for IEIs. `OR' is odds ratio.",
        label = "igad_lead_snp_table"
    localrule: True
    script: script_path('pub/igad_paper/tables/print_igad_lead_snp_table.R')

rule print_cfdr_lead_snp_table:
    input:
        "results/pub/igad_paper/tables/non_mhc_lead_snps.tsv",
    output:
        tex = "results/pub/igad_paper/tables/cfdr_lead_snps.tex",
        tsv = "results/pub/igad_paper/tables/cfdr_lead_snps.tsv"
    params:
        caption = r"Lead SNPs from genome-wide significant associations in the cFDR analysis. The `Variant' column gives the rsID of each SNP, and the reference and effect alleles separated by `$>$'. `Gene' gives the gene(s) with the most evidence linking it/them to the association signal. `Novel' indicates whether an association with SIgAD has previously been reported for a SNP. `IEI gene' indicates whether the SNP is located in, near, or is otherwise associated with a gene known to harbour variants causal for IEIs. We list the p-values for SIgAD and the three auxiliary traits; `RA' is rheumatoid arthritis. `Auxiliary trait significance' indicates whether each auxiliary trait reached genome-wide significance: `$+$' indicates a significant risk effect, `$-$' a significant protective effect, and `$\cdot$' no significant effect. The `v-value' is the output of the cFDR procedure and a p-value against a null hypothesis of no association of the SNP with SIgAD after conditioning on the auxiliary traits.",
        label = "cfdr_lead_snp_table"
    localrule: True
    script: script_path('pub/igad_paper/tables/print_cfdr_lead_snp_table.R')

rule print_iga_lead_snp_table:
    input:
        "results/pub/igad_paper/tables/iga_non_mhc_lead_snps.tsv",
    output:
        tex = "results/pub/igad_paper/tables/iga_lead_snps.tex",
        tsv = "results/pub/igad_paper/tables/iga_lead_snps.tsv"
    params:
        caption = r"Lead SNPs from genome-wide significant associations in the serum IgA GWAS meta-analysis. The `Variant' column gives the rsID of each SNP, and the reference and effect alleles separated by `$>$'. `Gene' gives the gene(s) with the most evidence linking it/them to the association signal. `Novel' indicates whether an association with SIgAD has previously been reported for a SNP. `IEI gene' indicates whether the SNP is located in, near, or is otherwise associated with a gene known to harbour variants causal for IEIs. `GWAS p-value' gives the meta-analytic p-value. `Effect direction' indicates whether the effect allele is associated with an IgA-increasing (`$+$') or decreasing (`$-$') effect. `Study effects' indicates whether a significant effect was found in the component GWAS of our meta-analysis: Liu, our own GWAS, and Dennis, respectively. `$\cdot$' indicates no significant effect.",
        label = "iga_lead_snps",
        # Missense variant so we know it has nothing to do with TNFSF12-TNFSF13
        iei_false_positives = ['CD68']
    localrule: True
    script: script_path('pub/igad_paper/tables/print_iga_lead_snp_table.R')

rule print_lyons_iga_lead_snp_table:
    input:
        "results/processed_gwas/iga/sans_mhc/gws/lead_snps.distance_clumped.rsIDs",
    output:
        tex = "results/pub/igad_paper/tables/lyons_iga_lead_snps.tex",
        tsv = "results/pub/igad_paper/tables/lyons_iga_lead_snps.tsv"
    params:
        caption = "Lead SNPs from genome-wide significant associations in our GWAS of serum IgA. The `Variant' column gives the rsID of each SNP, and the reference and effect alleles separated by `$>$'. `Gene' gives the gene(s) with the most evidence linking it/them to the association signal. `Novel' indicates whether an association with serum IgA had previously been reported for a SNP. `Effect direction' indicates whether the effect allele is associated with an IgA-increasing (`$+$') or decreasing (`$-$') effect.",
        label = "lyons_iga_lead_snps"
    localrule: True
    script: script_path('pub/igad_paper/tables/print_lyons_iga_lead_snp_table.R')

rule locate_lead_snps:
    input:
        "results/pub/igad_paper/tables/non_mhc_lead_snps.tsv",
    output:
        "results/pub/igad_paper/tables/gene_exon_intron_lists.RData"
    localrule: True
    conda: env_path('gene_annotations.yaml')
    script: script_path('pub/igad_paper/tables/locate_lead_snps.R')


rule print_tnfaip3_gws_table:
    input:
        "results/pub/igad_paper/resources/gwas-association-downloaded_2024-02-13-ensemblMappedGenes_TNFAIP3.tsv"
    output:
        "results/pub/igad_paper/tables/tnfaip3_gws_table.tex"
    params:
        cited_study_accessions = {"yin": "GCST011956",
                                   "ishigaki": "GCST90132222",
                                   "beecham": "GCST005531",
                                   "li": "GCST002217",
                                   "tsoi": "GCST005527"}
    localrule: True
    script: script_path('pub/igad_paper/tables/print_tnfaip3_gws_table.R')

rule print_iga_gwas_metadata_table:
    input:
        "resources/gwas/metadata/metadata_all_fields.tsv"
    output:
        tex = "results/pub/igad_paper/tables/iga_gwas_metadata.tex",
        tsv = "results/pub/igad_paper/tables/iga_gwas_metadata.tsv"
    params:
        abbrvs = config.get('igad_paper').get('iga_gwas'),
        caption = "The serum IgA GWAS data sets we subjected to meta-analysis. The 8,000-sample GWAS is as yet unpublished.",
        label = "iga_gwas_metadata"
    localrule: True
    script: script_path('pub/igad_paper/tables/print_iga_gwas_metadata_table.R')

rule collate_existing_igan_associations:
    input:
        ebi = "resources/gwas/igan/ebi_associations.tsv",
        kiryluk = "resources/gwas/igan/kiryluk_table_1.tsv"
    output:
        "results/iga_meta/existing_igan_associations.tsv"
    localrule: True
    script: script_path("iga_meta/collate_existing_igan_associations.R")

rule print_gwas_lead_snps_with_other_associations_table:
    input:
        igad = "results/pub/igad_paper/tables/non_mhc_lead_snps.tsv",
        others = "resources/igad_meta/gws/lead_snps_and_non_igad_associations.tsv"
    output:
        tex = "results/pub/igad_paper/tables/lead_snps_with_other_associations.tex",
        tsv = "results/pub/igad_paper/tables/lead_snps_with_other_associations.tsv"
    params:
        caption = "The association of the lead SNPs from the SIgAD meta-analysis with other immune-mediated diseases (IMDs) as identified by GWAS. The `Variant' column gives the rsID of each SNP, and the reference and effect alleles separated by `$>$'. `Gene' gives the gene(s) with the most evidence linking it/them to the association signal. `Novel' indicates whether an association with SIgAD has previously been reported for a SNP. `IMD associations in LD' lists IMDs which have significantly associated variants in linkage disequilibrium with the SIgAD-associated variant given in `Variant'.",
        label = "lead_snp_table_with_others"
    localrule: True
    script: script_path('pub/igad_paper/tables/print_lead_snps_with_other_associations_table.R')

rule outer_join_of_gwas_and_cfdr_lead_snps:
    input:
        gwas = "results/igad_meta/gws/lead_snps.distance_clumped.rsIDs",
        cfdr = "results/cfdr/bronson-finngen-igad-gc_on_asthma-ex_and_ra-ishigaki_and_liu-decode-lyons-dennis-iga/left/sans_mhc/snps_only/1000kb_1_0_2/3/gws/lead_snps.distance_clumped.rsIDs",
    output:
        "results/pub/igad_paper/tables/gwas_and_cfdr_lead_snps.tsv",
    params:
        window = 2e5,
        top_gene_to_chosen_gene = config.get('igad').get('top_gene_to_chosen_gene'),
    localrule: True
    script: script_path('pub/igad_paper/tables/outer_join_of_gwas_and_cfdr_lead_snps.R')

rule gwas_and_cfdr_gene_coordinates:
    input:
        "results/pub/igad_paper/tables/gwas_and_cfdr_lead_snps.tsv"
    output:
        gwas = "results/pub/igad_paper/tables/gwas_gene_coordinates.tsv",
        cfdr = "results/pub/igad_paper/tables/cfdr_gene_coordinates.tsv"
    params:
        window = 1e5,
    localrule: True
    script: script_path('pub/igad_paper/tables/fetch_gwas_and_cfdr_gene_coordinates.R')

rule tabulate_iei_genes_across_analyses:
    input:
        igad = "results/pub/igad_paper/tables/non_mhc_lead_snps.tsv",
        iga = "results/pub/igad_paper/tables/iga_non_mhc_lead_snps.tsv"
    output:
        "results/pub/igad_paper/tables/iei_genes.tsv"
    localrule: True
    script: script_path('pub/igad_paper/tables/tabulate_iei_genes_across_analyses.R')

rule write_out_mr_table_for_supplement:
    input:
        "results/pub/igad_paper/misc/stats_at_gws_iga_snps.tsv"
    output:
        "results/pub/igad_paper/tables/mr_stats.tsv"
    localrule: True
    run:
        daf = pd.read_csv(input[0], sep = '\t')

        daf['Variant'] = daf['rsID'] + ':' + daf['REF'] + '>' + daf['ALT']

        daf = daf.rename(columns = {'BP38': 'Position',
                              'CHR38': 'Chromosome',
                              'BETA.meta': 'Effect size (IgA)',
                              'SE.meta': 'Standard error (IgA)',
                              'BETA.igad': 'Effect size (SIgAD)',
                              'SE.igad': 'Standard error (SIgAD)',
                              'BETA.igan': 'Effect size (IgAN)',
                              'SE.igan': 'Standard error (IgAN)',
                              })

        daf[['Variant',
            'Position',
            'Chromosome',
            'Effect size (IgA)',
            'Standard error (IgA)',
            'Effect size (SIgAD)',
            'Standard error (SIgAD)',
            'Effect size (IgAN)',
            'Standard error (IgAN)']].to_csv(output[0], sep = '\t', index = False)

rule all_igad_paper_tables:
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
        "results/pub/igad_paper/tables/lead_snps_with_other_associations.tsv"
        # Supplementary Table 7
        "results/pub/igad_paper/misc/stats_at_gws_iga_snps.tsv"
