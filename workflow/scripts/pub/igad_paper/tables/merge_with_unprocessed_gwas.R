library(data.table)
setDTthreads(snakemake@threads)

lead_snps <- fread(snakemake@input$lead_snps)

# TODO I think this is Bronson with incorrect effect direction
bronson_cols <- c('variant_id', 'chromosome', 'base_pair_location', 'effect_allele', 'other_allele', 'odds_ratio', 'p_value')
bronson <- fread(snakemake@input$bronson, select = bronson_cols)
setnames(bronson, bronson_cols, paste0(bronson_cols, '.bronson'))

lead_snps <- merge(lead_snps, bronson, all.x = T, by.x = 'rsID', by.y = 'variant_id.bronson')
rm(bronson)

finngen_cols <- c('#chrom', 'pos', 'ref', 'alt', 'rsids', 'pval', 'beta', 'sebeta')
finngen <- fread(snakemake@input$finngen, select = finngen_cols)
setnames(finngen, finngen_cols, paste0(finngen_cols, '.finngen'))

lead_snps <- merge(lead_snps, finngen, all.x = T, by.x = 'rsID', by.y = 'rsids.finngen')
rm(finngen)

liu_cols <- c('CHR', 'BP_hg19', 'A1', 'A2', 'BETA', 'SE', 'P')
liu <- fread(snakemake@input$liu, select = liu_cols)
setnames(liu, liu_cols, paste0(liu_cols, '.liu'))

lead_snps <- merge(lead_snps, liu, all.x = T, by.x = c('CHR19', 'BP19'), by.y = c('CHR.liu', 'BP_hg19.liu'))
lead_snps[, `:=` (A1.liu = toupper(A1.liu), A2.liu = toupper(A2.liu))]
rm(liu)

dennis_cols <- c('variant_id', 'A1', 'A2', 'BETA', 'SE', 'p_value')
dennis <- fread(snakemake@input$dennis, select = dennis_cols)
setnames(dennis, dennis_cols, paste0(dennis_cols, '.dennis'))

lead_snps <- merge(lead_snps, dennis, all.x = T, by.x = 'rsID', by.y = 'variant_id.dennis')
rm(dennis)

lyons_cols <- c('rsid', 'REF', 'ALT', 'p_value', 'beta', 'se')
lyons <- fread(snakemake@input$lyons, select = lyons_cols)
setnames(lyons, lyons_cols, paste0(lyons_cols, '.lyons'))

lead_snps <- merge(lead_snps, lyons, by.x = 'rsID', by.y = 'rsid.lyons', all.x = T)
rm(lyons)

ra_cols <- c('variant_id', 'effect_allele', 'other_allele', 'beta', 'p_value')
ra <- fread(snakemake@input$ra, select = ra_cols)
setnames(ra, ra_cols, paste0(ra_cols, '.ra'))

lead_snps <- merge(lead_snps, ra, by.x = 'rsID', by.y = 'variant_id.ra', all.x = T)
rm(ra)

asthma_cols <- c('chr', 'pos', 'ref', 'alt', 'beta_EUR')
asthma <- fread(snakemake@input$asthma, select = asthma_cols)
setnames(asthma, asthma_cols, paste0(asthma_cols, '.asthma'))
asthma[, chr.asthma := as.integer(chr.asthma)]

lead_snps <- merge(lead_snps, asthma, by.x = c('CHR19', 'BP19'), by.y = c('chr.asthma', 'pos.asthma'), all.x = T)
rm(asthma)

fwrite(lead_snps, file = snakemake@output[[1]], sep = '\t')
