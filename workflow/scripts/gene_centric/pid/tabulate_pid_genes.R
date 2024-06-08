library(data.table)
library(xlsx)

dat <- data.table(read.xlsx(snakemake@input[[1]], sheetIndex = 1))

dat <- na.omit(unique(dat[, .(Genetic.defect)]))

dat[, gene := trimws(Genetic.defect)]

# NB: This was the regex used to identify non-standard gene field contents
#dat[grepl('^[A-Z0-9]*$', gene)]

gene_fields_to_exclude <- c('11q23del',
'Del10p13-p14',
'Large (3Mb) deletion of 22q11.2',
'Unknown / environment',
'Unknown',
'Mutation or chromosomal deletion at 14q32')

dat <- dat[!(gene %in% gene_fields_to_exclude)]

dat[gene == "CD40 (TNFRSF5)", gene := 'CD40']
dat[gene == "CD40LG (TNFSF5)", gene := 'CD40LG']
dat[gene %like% 'BLM', gene := 'BLM']
dat[gene %like% 'KMT2D', gene := 'KMT2D']
dat[gene %like% 'MOGS', gene := 'MOGS']
dat[gene %like% 'PIK3CD', gene := 'PIK3CD']
dat[gene %like% 'IKBKG', gene := 'IKBKG']
dat <- dat[!(gene %like% 'C4A')]
dat <- dat[!(gene %like% 'CFHR1')]

out_dat <- rbindlist(list
(dat[, .(gene)],
  data.table(gene = c('CFHR1', 'CFHR2', 'CFHR3', 'CFHR4', 'CFHR5', 'C4A', 'C4B'))
)
)

out_dat <- unique(out_dat, by = 'gene')

fwrite(out_dat[, .(gene)], file = snakemake@output[[1]], sep = '\t')
