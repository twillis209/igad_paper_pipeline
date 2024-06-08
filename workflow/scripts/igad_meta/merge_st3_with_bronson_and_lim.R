library(data.table)
setDTthreads(snakemake@threads)
library(magrittr)

st3 <- fread(snakemake@input[['supp_table_three']], select = c("Variant", "Allele", "allele_hg38_status", "OR", "gnomadNFE"))
setnames(st3, c('Allele', 'allele_hg38_status', 'OR'), c('allele.ST3', 'allele_status.ST3', 'OR.ST3'))

bronson <- fread(snakemake@input[['igad']], select = c('variant_id', 'chromosome', 'base_pair_location', 'effect_allele', 'other_allele', 'odds_ratio'))
setnames(bronson, c('effect_allele', 'other_allele', 'odds_ratio'), c('effect_allele.bronson', 'other_allele.bronson', 'OR.bronson'))

lim <- fread(snakemake@input[['lim']], select = c('CHR19', 'BP19', 'REF', 'ALT', 'BETA'))
lim[, OR.lim := exp(BETA)]
setnames(lim, c('REF', 'ALT'), c('other_allele.lim', 'effect_allele.lim'))

merge(bronson, st3, by.x = 'variant_id', by.y = 'Variant', all.y = T) %>%
  merge(lim, by.x = c('chromosome', 'base_pair_location'), by.y = c('CHR19', 'BP19'), all.x = T) -> merged

or_cols <- c('OR.bronson', 'OR.lim', 'OR.ST3')

merged <- merged[, .(variant_id, effect_allele.bronson, effect_allele.lim, other_allele.bronson, other_allele.lim, allele.ST3, allele_status.ST3, gnomadNFE, OR.bronson, OR.lim, OR.ST3)]
merged[, (or_cols) := lapply(.SD, signif, 2), .SDcols = or_cols]

fwrite(merged, file = snakemake@output[[1]], sep = '\t')
