library(data.table)
setDTthreads(snakemake@threads)

anns <- fread(snakemake@input[[1]])

fwrite(anns[id == T & imd == F & iga == F, .(SNPID)], sep = ' ', col.names = F, file = snakemake@output[['in_snps']])
fwrite(anns[!(id == T & imd == F & iga == F), .(SNPID)], sep = ' ', col.names = F, file = snakemake@output[['out_snps']])
