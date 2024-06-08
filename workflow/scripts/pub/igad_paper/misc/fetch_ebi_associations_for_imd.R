library(gwasrapidd)
library(dplyr)

print(snakemake@params$efo_trait_code)

assocs <- get_associations(efo_id = snakemake@params$efo_trait_code)
variants <- get_variants(efo_id = snakemake@params$efo_trait_code)

subset(assocs@associations, pvalue <= 5e-8) %>%
  select(c('association_id', 'pvalue', 'or_per_copy_number', 'beta_unit', 'beta_direction')) %>%
  left_join(assocs@risk_alleles[, c('association_id', 'variant_id', 'risk_allele')], by = 'association_id') %>%
  left_join(variants@variants[, c('variant_id', 'chromosome_name', 'chromosome_position')], by = 'variant_id') %>%
  distinct(variant_id, .keep_all = T) %>%
  arrange(chromosome_name, chromosome_position) %>%
  write.table(sep = '\t', file = snakemake@output[[1]], row.names = F)

