library(data.table)

dat <- fread(snakemake@input[[1]], select = c('rsID', 'CHR', 'BP', 'chosen_gene', 'P.igad'))

dat <- dat[P.igad < 5e-8]

dat <- unique(dat)

dat[, manhattan_label := chosen_gene]

fwrite(dat, file = snakemake@output[[1]], sep = '\t')
