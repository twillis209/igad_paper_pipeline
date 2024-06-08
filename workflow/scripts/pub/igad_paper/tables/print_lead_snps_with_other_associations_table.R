library(data.table)
library(kableExtra)
library(magrittr)

dat <- fread(snakemake@input$igad)

others <- fread(snakemake@input$others)

dat <- merge(dat, others, by = 'rsID')


cols_to_keep <- c('Variant' = 'SNPID', 'Chromosome' = 'CHR', 'Position' = 'BP', 'Gene' = 'chosen_gene', 'Novel' = 'novel', 'IMD associations' = 'imds', 'IMD associations in LD' = 'ld_imds')

dat[, BP := format(BP, big.mark = ',')]
dat[, novel := ifelse(already_reported == T, '', "true")]
dat[, SNPID := paste0(rsID, ':', REF, '>', ALT)]

dat[order(CHR, BP)] %>%
  .[, ..cols_to_keep] %>%
  setnames(cols_to_keep, names(cols_to_keep)) %>%
  fwrite(sep = '\t', file = snakemake@output[['tsv']])

# format for tex
dat[, novel := ifelse(already_reported == T, '', r"(\checkmark)"), by = 1:nrow(dat)]
dat[, SNPID := paste0(rsID, ':', REF, '$>$', ALT)]

dat[order(CHR, BP)] %>%
  .[, ..cols_to_keep] %>%
  kbl(
    format = 'latex',
    col.names = names(cols_to_keep),
    caption = snakemake@params$caption,
    label = snakemake@params$label,
    escape = F,
    align = c('l', rep('r', 2), 'l', 'r', rep('l', 2)),
    booktabs = T
  ) %>%
  column_spec(6:7, width = '4cm' ) %>%
  kable_styling(latex_options = 'scale_down') %>%
  save_kable(snakemake@output[['tex']])
