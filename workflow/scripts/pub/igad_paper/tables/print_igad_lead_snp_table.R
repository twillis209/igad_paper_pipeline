library(data.table)
library(kableExtra)
library(magrittr)

dat <- fread(snakemake@input[[1]])

dat <- dat[origin == 'GWAS']

dat[, gnomadNFE := format(gnomadNFE, digits = 2)]
dat[, BP := format(BP, big.mark = ',')]
dat[, novel := ifelse(already_reported == T, '', "true"), by = 1:nrow(dat)]
dat[, iei := ifelse(in_or_near_iei_gene == T, "true", ''), by = 1:nrow(dat)]
dat[, SNPID := paste0(rsID, ':', REF, '>', ALT)]

dat[mostSevereConsequence == '3 prime UTR', mostSevereConsequence := r'($3\text{'}$ UTR)']

cols_to_keep <- c('Variant' = 'SNPID', 'Chromosome' = 'CHR', 'Position' = 'BP', 'Effect allele frequency' = 'gnomadNFE', 'Gene' = 'chosen_gene', 'Novel' = 'novel', 'IEI gene' = 'iei', 'OR' = 'OR.igad', 'p-value' = 'P.igad')

dat[order(CHR, BP)] %>%
  .[, ..cols_to_keep] %>%
  setnames(cols_to_keep, names(cols_to_keep)) %>%
  fwrite(sep = '\t', file = snakemake@output[['tsv']])

# Format for tex
dat[, SNPID := paste0(rsID, ':', REF, '$>$', ALT)]
dat[, novel := ifelse(already_reported == T, '', r"(\checkmark)"), by = 1:nrow(dat)]
dat[, iei := ifelse(in_or_near_iei_gene == T, r"(\checkmark)", ''), by = 1:nrow(dat)]
dat[, chosen_gene := paste0("\\textit{", chosen_gene, '}')]

scientific_cols <- c('P.igad', 'P.cfdr')
dat[, (scientific_cols) := lapply(.SD, format, scientific = T, digits = 2), .SDcols = scientific_cols]
dat[, OR.igad := format(OR.igad, digits = 2)]
dat[, gnomadNFE := round(as.numeric(gnomadNFE), 2)]

dat[order(CHR, BP)] %>%
  .[, ..cols_to_keep] %>%
  kbl(
    format = 'latex',
    col.names = names(cols_to_keep),
    caption = snakemake@params$caption,
    label = snakemake@params$label,
    escape = F,
    align = c('l', rep('r', 3), rep('l', 3), rep('r', 2)),
    booktabs = T
  ) %>%
  column_spec(4, width = '1.5cm' ) %>%
  column_spec(10, width = '1.5cm' ) %>%
  save_kable(snakemake@output[['tex']])
