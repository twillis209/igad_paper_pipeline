library(data.table)
library(kableExtra)
library(magrittr)

dat <- fread(snakemake@input[[1]], select = c('abbrv', 'pretty_name', 'N0', 'Ancestry', 'citation', 'First_Author'))

dat <- dat[abbrv %in% snakemake@params$abbrvs]

dat[!(First_Author %in% c('FinnGen', 'Pan-UKB', 'IMSGC', 'UKB')), First_Author := paste(First_Author, 'et al.', sep = ' ')]

fwrite(dat[order(Ancestry), .(Phenotype = pretty_name, `Sample size` = N0, Ancestry, Citation = First_Author)], file = snakemake@output[['tsv']], sep = '\t')

sample_size_cols <- c('N0')
dat[, (sample_size_cols) := lapply(.SD, format, big.mark = ','), .SDcols = sample_size_cols]

dat[citation != '', citation := paste0('\\cite{', citation, '}')]

dat[order(Ancestry), .(pretty_name, N0, Ancestry, citation)] %>%
  kbl(
    format = 'latex',
    col.names = c('Phenotype', 'Sample size', 'Ancestry', 'Citation'),
    caption = snakemake@params$caption,
    label = snakemake@params$label,
    escape = F,
    align = c('l', rep('r', 1), rep('l', 2)),
    booktabs = T
  ) %>%
  collapse_rows(c(1,3)) %>%
  save_kable(snakemake@output[['tex']])
