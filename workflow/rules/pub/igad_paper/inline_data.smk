rule make_igad_paper_data_file:
    input:
        sans_mhc_enrichment_gif = "results/pub/igad_paper/misc/100kb_sans_mhc_iei_genes_bronson-finngen-igad_gif.tsv",
        with_mhc_enrichment_gif = "results/pub/igad_paper/misc/100kb_with_mhc_iei_genes_bronson-finngen-igad_gif.tsv",
        sans_mhc_pan_genes_enrichment_gif = "results/pub/igad_paper/misc/100kb_sans_mhc_all_genes_bronson-finngen-igad_gif.tsv",
        with_mhc_pan_genes_enrichment_gif = "results/pub/igad_paper/misc/100kb_with_mhc_all_genes_bronson-finngen-igad_gif.tsv",
        sans_mhc_enrichment_pvalue = "results/pub/igad_paper/misc/100kb_sans_mhc_100_5000_pvalue.tsv",
        with_mhc_enrichment_pvalue = "results/pub/igad_paper/misc/100kb_with_mhc_101_5000_pvalue.tsv",
        gwas_gif = "results/pub/igad_paper/misc/compiled_gif.tsv",
        igad_rg = "results/ldak/ldak-thin/combined/sans_mhc/snps_only/bronson-finngen-igad_and_imds.tsv",
        imd_rg = "results/ldak/ldak-thin/combined/sans_mhc/snps_only/imds.tsv"
    output:
        "results/pub/igad_paper/data/inline_data.csv"
    params:
        imd_traits = config.get('igad_paper').get('imd_traits')
    localrule: True
    script: script_path('pub/igad_paper/data/make_igad_paper_data_file.py')

rule gif_to_update:
    input:
        sans_mhc_enrichment_gif = "results/pub/igad_paper/misc/100kb_sans_mhc_iei_genes_bronson-finngen-igad_gif.tsv",
        with_mhc_enrichment_gif = "results/pub/igad_paper/misc/100kb_with_mhc_iei_genes_bronson-finngen-igad_gif.tsv",
        sans_mhc_enrichment_pvalue = "results/pub/igad_paper/misc/100kb_sans_mhc_100_5000_pvalue.tsv",
        with_mhc_enrichment_pvalue = "results/pub/igad_paper/misc/100kb_with_mhc_101_5000_pvalue.tsv",
        gwas_gif = "results/pub/igad_paper/misc/compiled_gif.tsv"
