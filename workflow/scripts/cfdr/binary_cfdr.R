library(data.table)
setDTthreads(snakemake@threads)
library(fcfdr)

snp_col <- snakemake@params[['snp_col']]
chr_col <- snakemake@params[['chr_col']]
bp_col <- snakemake@params[['bp_col']]
ref_col <- snakemake@params[['ref_col']]
alt_col <- snakemake@params[['alt_col']]
prin_col <- snakemake@params[['prin_col']]
aux_cols <- snakemake@params[['aux_cols']]
tmpdir <- snakemake@resources[['tmpdir']]

gwas_dat <- fread(snakemake@input[[1]], sep = '\t', header = T, tmpdir = tmpdir)

gwas_dat[, (chr_col) := as.character(get(chr_col))]

gwas_dat <- na.omit(gwas_dat, cols = c(chr_col, bp_col, prin_col))

gwas_dat <- unique(gwas_dat, by = c(chr_col, bp_col))

# Transforming p-values of 0 to Z-scores yields Inf values
gwas_dat[get(prin_col) < 1e-300, (prin_col) := 1e-300]

res <- binary_cfdr(p = gwas_dat[[prin_col]],
                   q = gwas_dat[[aux_cols[1]]],
                   group = gwas_dat[[chr_col]],
                   threads = snakemake@threads
                   )

gwas_dat[, v.1 := res$v]

fwrite(gwas_dat, file = snakemake@output[[1]], sep = '\t')
