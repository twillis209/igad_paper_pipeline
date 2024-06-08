library(data.table)
library(kableExtra)
library(magrittr)

dat <- fread(snakemake@input[[1]])

dat[topGene == 'ASCC2', topGene := 'HORMAD2']

dat[, gnomadNFE := round(as.numeric(gnomadNFE), 2)]
dat[, BP := format(BP, big.mark = ',')]


cols_to_keep <- c('Variant' = 'SNPID',
                  'Chromosome' = 'CHR',
                  'Position' = 'BP',
                  'Effect allele frequency' = 'gnomadNFE',
                  'Gene' = 'topGene',
                  'Novel' = 'novel',
                  'Effect size' = 'BETA',
                  'GWAS p-value' = 'P',
                  'Effect direction' = 'symbol'
                  )

dat[, already_reported := T]
dat[rsID == 'rs184702468', already_reported := F]
dat[, iei := '']
dat[, novel := ifelse(already_reported == T, '', "true"), by = 1:nrow(dat)]

dat[, symbol := ifelse(BETA > 0, '+', '-')]

dat[, SNPID := paste0(rsID, ':', REF, '>', ALT)]

dat[, mostSevereConsequence := gsub('upstream_gene_variant|intergenic_variant', 'intergenic', mostSevereConsequence)]
dat[, mostSevereConsequence := gsub('intron_variant', 'intron', mostSevereConsequence)]
dat[, mostSevereConsequence := gsub('missense_variant', 'missense', mostSevereConsequence)]

dat[order(CHR, BP)] %>%
  .[, ..cols_to_keep] %>%
  setnames(cols_to_keep, names(cols_to_keep)) %>%
fwrite(sep = '\t', file = snakemake@output[['tsv']])

#format for tex
dat[, SNPID := paste0(rsID, ':', REF, '$>$', ALT)]
dat[, topGene := paste0('\\textit{', topGene, '}')]
dat[, symbol := paste0("$", symbol, "$")]
dat[, novel := ifelse(already_reported == T, '', r"(\checkmark)"), by = 1:nrow(dat)]

scientific_cols <- c('P', 'BETA')
dat[, (scientific_cols) := lapply(.SD, format, scientific = T, digits = 2), .SDcols = scientific_cols]

dat[order(CHR, BP)] %>%
  .[, ..cols_to_keep] %>%
  kbl(
    format = 'latex',
    col.names = names(cols_to_keep),
    caption = snakemake@params$caption,
    label = snakemake@params$label,
    escape = F,
    align = c('l', rep('r', 3), rep('l', 1), rep('r', 3)),
    booktabs = T
  ) %>%
  column_spec(4, width = '1.5cm' ) %>%
  column_spec(7, width = '1.5cm' ) %>%
  column_spec(8, width = '1.5cm' ) %>%
  kable_styling(latex_options = 'scale_down') %>%
  save_kable(snakemake@output[['tex']])
