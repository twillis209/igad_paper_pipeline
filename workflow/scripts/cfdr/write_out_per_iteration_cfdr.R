library(data.table)
library(stringr)
library(magrittr)

setDTthreads(snakemake@threads)

chr_col <- snakemake@params[['chr_col']]
bp_col <- snakemake@params[['bp_col']]
ref_col <- snakemake@params[['ref_col']]
alt_col <- snakemake@params[['alt_col']]
snp_col <- snakemake@params[['snp_col']]
prin_col <- snakemake@params[['prin_col']]
aux_cols <- snakemake@params[['aux_cols']]
prune_cols <- snakemake@params[['prune_cols']]
v_cols <- snakemake@params[['v_cols']]

dat <- fread(snakemake@input[[1]], sep = '\t', header = T)

dir.create(path = snakemake@output[[1]])

for(i in seq_along(snakemake@params[['aux_traits']])) {
  colnames <- c(snp_col, chr_col, bp_col, ref_col, alt_col, prin_col, aux_cols[i], prune_cols[i], v_cols[i])

  print(i)
  print(colnames)

  sub_dat <- dat[, ..colnames]

  fwrite(sub_dat, file = sprintf(snakemake@params[['out_format']], i), sep = '\t')
}
