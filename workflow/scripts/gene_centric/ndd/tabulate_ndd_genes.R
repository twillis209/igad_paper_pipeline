library(data.table)

dat <- fread(snakemake@input[[1]])

dat <- dat[`Candidate NDD genes v3` == T]

dat <- dat[, .(`Ensembl ID`, `Gene`, `Gene type`, chr = `hg38 - chromosome`, start  = `hg38 - start`, stop = `hg38 - end`, `pLI`)]

fwrite(dat, file = snakemake@output[[1]], sep = '\t')
