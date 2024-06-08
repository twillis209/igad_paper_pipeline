library(data.table)
setDTthreads(snakemake@threads)

ichip <- fread(snakemake@input[['ichip']], sep = '\t', header = T)
reading_snps <- fread(snakemake@input[['reading_snps']], sep =',', header = T, select = c(snakemake@params[['reading_snpid']]))

ichip <- merge(ichip, reading_snps, by.x = snakemake@params[['snpid']], by.y = snakemake@params[['reading_snpid']])

fwrite(ichip, file = snakemake@output[[1]], sep = '\t')
