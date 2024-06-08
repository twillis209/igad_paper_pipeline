library(data.table)
setDTthreads(snakemake@threads)

snp_col <- snakemake@params[['snp_col']]

dat <- fread(snakemake@input[[1]], sep = '\t', header = T)

dat <- dat[!(get(snp_col) %in% snakemake@params[['snps_to_drop']])]

fwrite(dat, file = snakemake@output[[1]], sep = '\t')
