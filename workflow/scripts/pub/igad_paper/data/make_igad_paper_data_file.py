import pandas as pd
import statistics as stat

sans_mhc_enrichment_gif = pd.read_csv(snakemake.input.sans_mhc_enrichment_gif, sep = '\t')
with_mhc_enrichment_gif = pd.read_csv(snakemake.input.with_mhc_enrichment_gif, sep = '\t')
sans_mhc_pan_genes_enrichment_gif = pd.read_csv(snakemake.input.sans_mhc_pan_genes_enrichment_gif, sep = '\t')
with_mhc_pan_genes_enrichment_gif = pd.read_csv(snakemake.input.with_mhc_pan_genes_enrichment_gif, sep = '\t')
sans_mhc_enrichment_pvalue = pd.read_csv(snakemake.input.sans_mhc_enrichment_pvalue, sep = '\t')
with_mhc_enrichment_pvalue = pd.read_csv(snakemake.input.with_mhc_enrichment_pvalue, sep = '\t')

gifs = pd.read_csv(snakemake.input.gwas_gif, sep = '\t')

imd_rg = pd.read_csv(snakemake.input.imd_rg, sep = '\t')
imd_rg = imd_rg.query("`trait.A` in @snakemake.params.imd_traits and `trait.B` in @snakemake.params.imd_traits")
igad_rg = pd.read_csv(snakemake.input.igad_rg, sep = '\t')
igad_rg = igad_rg.query("`trait.B` in @snakemake.params.imd_traits")

daf = pd.DataFrame(columns = ['key', 'value'])

l = []

l.append({'key': 'sans_mhc_enrichment_gif', 'value': sans_mhc_enrichment_gif['gif'].values[0]})
l.append({'key': 'with_mhc_enrichment_gif', 'value': with_mhc_enrichment_gif['gif'].values[0]})
l.append({'key': 'sans_mhc_pan_genes_enrichment_gif', 'value': sans_mhc_pan_genes_enrichment_gif['gif'].values[0]})
l.append({'key': 'with_mhc_pan_genes_enrichment_gif', 'value': with_mhc_pan_genes_enrichment_gif['gif'].values[0]})
l.append({'key': 'sans_mhc_enrichment_pvalue', 'value': sans_mhc_enrichment_pvalue['pvalue'].values[0]})
l.append({'key': 'with_mhc_enrichment_pvalue', 'value': with_mhc_enrichment_pvalue['pvalue'].values[0]})
l.append({'key': 'with_mhc_dennis_iga_gif', 'value': gifs.query('variant_set == \'with_mhc\' and dataset == \'dennis-iga\'')['lambda_0_50'].values[0]})
l.append({'key': 'sans_mhc_dennis_iga_gif', 'value': gifs.query('variant_set == \'sans_mhc\' and dataset == \'dennis-iga\'')['lambda_0_50'].values[0]})
l.append({'key': 'with_mhc_liu_decode_iga_gif', 'value': gifs.query('variant_set == \'with_mhc\' and dataset == \'liu-decode-iga\'')['lambda_0_50'].values[0]})
l.append({'key': 'sans_mhc_liu_decode_iga_gif', 'value': gifs.query('variant_set == \'sans_mhc\' and dataset == \'liu-decode-iga\'')['lambda_0_50'].values[0]})
l.append({'key': 'with_mhc_lyons_iga_gif', 'value': gifs.query('variant_set == \'with_mhc\' and dataset == \'lyons-iga\'')['lambda_0_50'].values[0]})
l.append({'key': 'sans_mhc_lyons_iga_gif', 'value': gifs.query('variant_set == \'sans_mhc\' and dataset == \'lyons-iga\'')['lambda_0_50'].values[0]})
l.append({'key': 'with_mhc_iga_meta_gif', 'value': gifs.query('variant_set == \'with_mhc\' and dataset == \'iga-meta\'')['lambda_0_50'].values[0]})
l.append({'key': 'sans_mhc_iga_meta_gif', 'value': gifs.query('variant_set == \'sans_mhc\' and dataset == \'iga-meta\'')['lambda_0_50'].values[0]})
l.append({'key': 'with_mhc_bronson_igad_gif', 'value': gifs.query('variant_set == \'with_mhc\' and dataset == \'bronson-igad\'')['lambda_0_50'].values[0]})
l.append({'key': 'sans_mhc_bronson_igad_gif', 'value': gifs.query('variant_set == \'sans_mhc\' and dataset == \'bronson-igad\'')['lambda_0_50'].values[0]})
l.append({'key': 'with_mhc_finngen_igad_gif', 'value': gifs.query('variant_set == \'with_mhc\' and dataset == \'finngen-igad\'')['lambda_0_50'].values[0]})
l.append({'key': 'sans_mhc_finngen_igad_gif', 'value': gifs.query('variant_set == \'sans_mhc\' and dataset == \'finngen-igad\'')['lambda_0_50'].values[0]})
l.append({'key': 'with_mhc_igad_meta_gif', 'value': gifs.query('variant_set == \'with_mhc\' and dataset == \'igad-meta\'')['lambda_0_50'].values[0]})
l.append({'key': 'sans_mhc_igad_meta_gif', 'value': gifs.query('variant_set == \'sans_mhc\' and dataset == \'igad-meta\'')['lambda_0_50'].values[0]})

l.append({'key': 'asthma_ra_rg', 'value': imd_rg.query('`trait.A` == \'asthma-ex\' and `trait.B` == \'ra\'')['rg.sr'].values[0]})
l.append({'key': 'asthma_ra_rg_p', 'value': imd_rg.query('`trait.A` == \'asthma-ex\' and `trait.B` == \'ra\'')['rg.p.sr'].values[0]})
l.append({'key': 'asthma_iga_rg', 'value': imd_rg.query('`trait.A` == \'asthma-ex\' and `trait.B` == \'liu-decode-lyons-dennis-iga\'')['rg.sr'].values[0]})
l.append({'key': 'asthma_iga_rg_p', 'value': imd_rg.query('`trait.A` == \'asthma-ex\' and `trait.B` == \'liu-decode-lyons-dennis-iga\'')['rg.p.sr'].values[0]})
l.append({'key': 'ra_iga_rg', 'value': imd_rg.query('`trait.A` == \'ra\' and `trait.B` == \'liu-decode-lyons-dennis-iga\'')['rg.sr'].values[0]})
l.append({'key': 'ra_iga_rg_p', 'value': imd_rg.query('`trait.A` == \'ra\' and `trait.B` == \'liu-decode-lyons-dennis-iga\'')['rg.p.sr'].values[0]})
l.append({'key': 'igad_ra_rg', 'value': igad_rg.query('`trait.A` == \'bronson-finngen-igad\' and `trait.B` == \'ra\'')['rg.sr'].values[0]})
l.append({'key': 'igad_ra_rg_p', 'value': igad_rg.query('`trait.A` == \'bronson-finngen-igad\' and `trait.B` == \'ra\'')['rg.p.sr'].values[0]})
l.append({'key': 'igad_asthma_rg', 'value': igad_rg.query('`trait.A` == \'bronson-finngen-igad\' and `trait.B` == \'asthma-ex\'')['rg.sr'].values[0]})
l.append({'key': 'igad_asthma_rg_p', 'value': igad_rg.query('`trait.A` == \'bronson-finngen-igad\' and `trait.B` == \'asthma-ex\'')['rg.p.sr'].values[0]})
l.append({'key': 'igad_iga_rg', 'value': igad_rg.query('`trait.A` == \'bronson-finngen-igad\' and `trait.B` == \'liu-decode-lyons-dennis-iga\'')['rg.sr'].values[0]})
l.append({'key': 'igad_iga_rg_p', 'value': igad_rg.query('`trait.A` == \'bronson-finngen-igad\' and `trait.B` == \'liu-decode-lyons-dennis-iga\'')['rg.p.sr'].values[0]})

pd.DataFrame(l).to_csv(snakemake.output[0], sep = ',', index = False)
