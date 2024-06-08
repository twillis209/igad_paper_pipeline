library(data.table)
library(kableExtra)
library(magrittr)

dat <- fread(snakemake@input[[1]])

beta_cols <- c('BETA.meta', 'BETA.liu_decode', 'BETA.lyons', 'BETA.dennis')
dat[, (gsub('BETA', 'sign', beta_cols)) := lapply(.SD, sign), .SDcols = beta_cols]

dat[topGene %in% snakemake@params$iei_false_positives, iei_gene := ""]

dat[topGene == 'NT5M', topGene := 'TNFRSF13B']

cols_to_keep <- c('Variant' = 'SNPID',
                  'Chromosome' = 'CHR38',
                  'Position' = 'BP38',
                  'Effect allele frequency' = 'gnomadNFE',
                  'Gene' = 'topGene',
                  'Novel' = 'novel',
                  'IEI' = 'iei',
                  'GWAS p-value' = 'P.meta',
                  'Effect size' = 'BETA.meta',
                  'Effect direction' = 'SYMBOL.meta',
                  'Study effects' = 'SYMBOL')

dat[, gnomadNFE := round(as.numeric(gnomadNFE), 2)]
dat[, BP38 := format(BP38, big.mark = ',')]
dat[, novel := ifelse(already_reported == T, '', "true")]
dat[, iei := ifelse(iei_gene != "", "true", '')]
dat[, SNPID := paste0(rsID, ':', REF, '>', ALT)]

symbol_cols <- c('SYMBOL.meta', 'SYMBOL.liu_decode', 'SYMBOL.lyons', 'SYMBOL.dennis')
dat[, (symbol_cols) := "."]
dat[P.liu_decode < 5e-8, SYMBOL.liu_decode := ifelse(sign.liu_decode == 1, '+', '-'), by = 1:nrow(dat)]
dat[P.lyons < 5e-8, SYMBOL.lyons := ifelse(sign.lyons == 1, '+', '-'), by = 1:nrow(dat)]
dat[P.dennis < 5e-8, SYMBOL.dennis := ifelse(sign.dennis == 1, '+', '-'), by = 1:nrow(dat)]
dat[P.meta < 5e-8, SYMBOL.meta := ifelse(sign.meta == 1, '+', '-'), by = 1:nrow(dat)]
dat[, SYMBOL := paste0(SYMBOL.liu_decode, SYMBOL.lyons, SYMBOL.dennis)]

dat[order(CHR38, BP38)] %>%
  .[, ..cols_to_keep] %>%
  setnames(cols_to_keep, names(cols_to_keep)) %>%
  fwrite(sep = '\t', file = snakemake@output[['tsv']])

#format for tex
dat[, novel := ifelse(already_reported == T, '', r"(\checkmark)"), by = 1:nrow(dat)]
dat[, iei := ifelse(iei_gene != "", r"(\checkmark)", ''), by = 1:nrow(dat)]
dat[, topGene := paste0("\\textit{", topGene, '}')]

dat[, (symbol_cols) := r"(\cdot)"]
dat[P.liu_decode < 5e-8, SYMBOL.liu_decode := ifelse(sign.liu_decode == 1, '+', '-'), by = 1:nrow(dat)]
dat[P.lyons < 5e-8, SYMBOL.lyons := ifelse(sign.lyons == 1, '+', '-'), by = 1:nrow(dat)]
dat[P.dennis < 5e-8, SYMBOL.dennis := ifelse(sign.dennis == 1, '+', '-'), by = 1:nrow(dat)]
dat[P.meta < 5e-8, SYMBOL.meta := ifelse(sign.meta == 1, '+', '-'), by = 1:nrow(dat)]
dat[, SYMBOL := paste0("$", SYMBOL.liu_decode, SYMBOL.lyons, SYMBOL.dennis, "$")]
dat[, SYMBOL.meta := paste0("$", SYMBOL.meta, "$")]
dat[, SNPID := paste0(rsID, ':', REF, '$>$', ALT)]

dat[mostSevereConsequence == '3 prime UTR', mostSevereConsequence := r'($3\text{'}$ UTR)']
dat[mostSevereConsequence == '5 prime UTR', mostSevereConsequence := r'($5\text{'}$ UTR)']

scientific_cols <- c('P.meta', beta_cols)
dat[, (scientific_cols) := lapply(.SD, format, scientific = T, digits = 2), .SDcols = scientific_cols]

dat[order(CHR38, BP38)] %>%
  .[, ..cols_to_keep] %>%
  kbl(
    format = 'latex',
    col.names = names(cols_to_keep),
    caption = snakemake@params$caption,
    label = snakemake@params$label,
    escape = F,
    align = c('l', rep('r', 3), rep('l', 3), rep('r', 3)),
    booktabs = T
  ) %>%
  column_spec(4, width = '1.5cm' ) %>%
  column_spec(10, width = '1.5cm' ) %>%
  column_spec(11, width = '1cm' ) %>%
  kable_styling(latex_options = 'scale_down') %>%
  save_kable(snakemake@output[['tex']])
