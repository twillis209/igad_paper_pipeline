library(data.table)
library(kableExtra)
library(magrittr)

dat <- fread(snakemake@input[[1]])

dat <- dat[origin == 'cFDR']

dat[, gnomadNFE := format(gnomadNFE, digits = 2)]
dat[, BP := format(BP, big.mark = ',')]
dat[, novel := ifelse(already_reported == T, '', "true"), by = 1:nrow(dat)]
dat[, iei := ifelse(in_or_near_iei_gene == T, "true", ''), by = 1:nrow(dat)]
dat[, SNPID := paste0(rsID, ':', REF, '>', ALT)]

dat[mostSevereConsequence == '3 prime UTR', mostSevereConsequence := r'($3\text{'}$ UTR)']

symbol_cols <- c('SYMBOL.iga', 'SYMBOL.ra', 'SYMBOL.asthma')
dat[, (symbol_cols) := "."]
dat[P.iga < 5e-8, SYMBOL.iga := ifelse(BETA.iga > 0, '+', '-')]
dat[P.ra < 5e-8, SYMBOL.ra := ifelse(BETA.ra > 0, '+', '-')]
dat[P.asthma < 5e-8, SYMBOL.asthma := ifelse(BETA.asthma > 0, '+', '-')]
dat[, SYMBOL := paste0(SYMBOL.iga, SYMBOL.asthma, SYMBOL.ra)]

cols_to_keep <- c('Variant' = 'SNPID', 'Chromosome' = 'CHR', 'Position' = 'BP', 'Effect allele frequency' = 'gnomadNFE', 'Gene' = 'chosen_gene', 'IEI gene' = 'iei', 'SIgAD p-value' = 'P.igad', 'SIgAD OR' = 'OR.igad', 'IgA p-value' = 'P.iga', 'Asthma p-value' = 'P.asthma', 'RA p-value' = 'P.ra', 'Auxiliary trait significance' = 'SYMBOL', 'SIgAD v-value' = 'P.cfdr')

dat[chosen_gene == 'RUNX3', chosen_gene := 'RUNX3, RP11-84D1.2']

dat[order(CHR, BP)] %>%
  .[, ..cols_to_keep] %>%
  setnames(cols_to_keep, names(cols_to_keep)) %>%
  fwrite(sep = '\t', file = snakemake@output[['tsv']])

# Format for tex
scientific_cols <- c('P.igad', 'P.cfdr', 'P.ra', 'P.asthma', 'P.iga')
dat[, (scientific_cols) := lapply(.SD, format, scientific = T, digits = 2), .SDcols = scientific_cols]
dat[, OR.igad := format(OR.igad, digits = 2)]
dat[, gnomadNFE := round(as.numeric(gnomadNFE), 2)]
dat[, novel := ifelse(already_reported == T, '', r"(\checkmark)"), by = 1:nrow(dat)]
dat[, iei := ifelse(in_or_near_iei_gene == T, r"(\checkmark)", ''), by = 1:nrow(dat)]
dat[, SNPID := paste0(rsID, ':', REF, '$>$', ALT)]

dat[, (symbol_cols) := r"(\cdot)"]
dat[P.iga < 5e-8, SYMBOL.iga := ifelse(BETA.iga > 0, '+', '-')]
dat[P.ra < 5e-8, SYMBOL.ra := ifelse(BETA.ra > 0, '+', '-')]
dat[P.asthma < 5e-8, SYMBOL.asthma := ifelse(BETA.asthma > 0, '+', '-')]
dat[, SYMBOL := paste0("$", SYMBOL.iga, SYMBOL.asthma, SYMBOL.ra, "$")]

dat[, chosen_gene := paste0("\\textit{", chosen_gene, '}')]

dat[order(CHR, BP)] %>%
  .[, ..cols_to_keep] %>%
  kbl(
    format = 'latex',
    col.names = names(cols_to_keep),
    caption = snakemake@params$caption,
    label = snakemake@params$label,
    escape = F,
    align = c('l', rep('r', 3), rep('l', 1), rep('r', 7)),
    booktabs = T
  ) %>%
  # gnomadNFE
  column_spec(4, width = '1.5cm' ) %>%
  # v-value
  column_spec(7:11, width = '1.5cm' ) %>%
  kable_styling(latex_options = 'scale_down') %>%
  save_kable(snakemake@output[[1]])
