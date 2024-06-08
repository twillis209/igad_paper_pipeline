library(data.table)
setDTthreads(snakemake@threads)

snp_col <- snakemake@params[['snp_col']]

lead_snps <- fread(snakemake@input[['lead_snps']])

sum_stats <- fread(snakemake@input[['sum_stats']])

sum_stats <- sum_stats[get(snp_col) %in% lead_snps[, get(snp_col)]]

fwrite(sum_stats, file = snakemake@output[[1]], sep = '\t')
