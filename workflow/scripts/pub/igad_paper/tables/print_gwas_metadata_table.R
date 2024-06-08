library(data.table)
library(kableExtra)
library(magrittr)

dat <- fread(snakemake@input[[1]], select = c('abbrv', 'pretty_name', 'N0', 'N1', 'Ancestry', 'citation', 'First_Author', 'URL'))

dat <- dat[abbrv %in% snakemake@params$abbrvs]

dat[!is.null(N1), Neff := round(4/(1/N0 + 1/N1))]
dat[N1 == 0, Neff := N0]

sample_size_cols <- c('N0', 'N1', 'Neff')
dat[, (sample_size_cols) := lapply(.SD, format, big.mark = ','), .SDcols = sample_size_cols]

dat[citation != '', citation := paste0('\\cite{', citation, '}')]

dat[!(First_Author %in% c('FinnGen', 'Pan-UKB', 'IMSGC', 'UKB')), First_Author := paste(First_Author, 'et al.', sep = ' ')]
dat[First_Author == 'Lyons et al.', First_Author := 'This work']

fwrite(dat[order(pretty_name, Ancestry), .(Phenotype = pretty_name, Cases = N1, Controls = N0, `Effective sample size` = Neff, Ancestry, Citation = First_Author, URL)], file = snakemake@output[['tsv']], sep = '\t')

dat[First_Author %in% c('FinnGen', 'Pan-UKB', 'IMSGC', 'UKB'), citation := paste(First_Author, citation, sep = ' ')]

dat[order(pretty_name, Ancestry), .(pretty_name, N1, N0, Neff, Ancestry, citation)] %>%
  kbl(
    format = 'latex',
    col.names = c('Phenotype', 'Cases', 'Controls', 'Effective sample size', 'Ancestry', 'Citation'),
    caption = snakemake@params$caption,
    label = snakemake@params$label,
    escape = F,
    align = c('l', rep('r', 3), rep('l', 2)),
    booktabs = T
  ) %>%
  collapse_rows(c(1,5)) %>%
  kable_styling(latex_options = 'scale_down') %>%
  save_kable(snakemake@output[['tex']])
