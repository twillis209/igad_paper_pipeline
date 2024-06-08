library(data.table)

dat <- fread(snakemake@input[[1]], select = c('rsID', 'CHR', 'BP', 'topGene', 'nearestGene', 'P', 'SNPID'))

dat <- dat[P < 5e-8]

dat <- unique(dat)

dat[, gene_label := topGene]

fwrite(dat, file = snakemake@output[[1]], sep = '\t')
