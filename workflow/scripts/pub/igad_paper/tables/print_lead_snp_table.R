library(data.table)
library(kableExtra)
library(magrittr)

save.image('igad_lead_snps.RData')

dat <- fread(snakemake@input[[1]])

dat[, gnomadNFE := format(gnomadNFE, digits = 2)]
dat[, BP := format(BP, big.mark = ',')]
dat[, novel := ifelse(already_reported == T, '', r"(\checkmark)"), by = 1:nrow(dat)]
dat[, iei := ifelse(in_or_near_iei_gene == T, r"(\checkmark)", ''), by = 1:nrow(dat)]
dat[, chosen_gene := paste0("\textit{", chosen_gene, '}')]

dat[, SNPID := paste0(rsID, ':', REF, '$>$', ALT)]

dat[mostSevereConsequence == '3 prime UTR', mostSevereConsequence := r'($3\text{'}$ UTR)']

scientific_cols <- c('P.igad_meta', 'P.igad_cfdr')
dat[, (scientific_cols) := lapply(.SD, format, scientific = T, digits = 2), .SDcols = scientific_cols]
dat[, OR.igad_meta := format(OR.igad_meta, digits = 2)]
dat[, gnomadNFE := round(as.numeric(gnomadNFE), 2)]

cols_to_keep <- c('Variant' = 'SNPID', 'Chr.' = 'CHR', 'Position' = 'BP', 'Effect allele frequency' = 'gnomadNFE', 'Gene' = 'chosen_gene', 'Location/type' = 'mostSevereConsequence', 'Novel' = 'novel', 'IEI gene' = 'iei', 'Analysis' = 'origin', 'OR' = 'OR.igad_meta', 'GWAS p-value' = 'P.igad_meta',  'cFDR v-value' = 'P.igad_cfdr')

dat[order(CHR, BP)] %>%
  .[, ..cols_to_keep] %>%
  kbl(
    format = 'latex',
    col.names = names(cols_to_keep),
    caption = snakemake@params$caption,
    label = snakemake@params$label,
    escape = F,
    align = c('l', rep('r', 3), rep('l', 5), rep('r', 3)),
    booktabs = T
  ) %>%
  column_spec(4, width = '1.5cm' ) %>%
  column_spec(11, width = '1.5cm' ) %>%
  column_spec(12, width = '1.5cm' ) %>%
  collapse_rows(c(2)) %>%
  kable_styling(latex_options = 'scale_down') %>%
  save_kable(snakemake@output[[1]])
