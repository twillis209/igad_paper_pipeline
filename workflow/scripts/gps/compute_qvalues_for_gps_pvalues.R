library(data.table)
library(qvalue)

dat <- fread(snakemake@input[[1]], sep = '\t', header = T)

dat[, qvalue := qval(pvalue)$qvalue]

fwrite(dat, file = snakemake@output[[1]], sep = '\t')
