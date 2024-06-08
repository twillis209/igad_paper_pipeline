library(data.table)
setDTthreads(snakemake@threads)


chr_col <- snakemake@params[['chr_col']]
bp_col <- snakemake@params[['bp_col']]
ref_col <- snakemake@params[['ref_col']]
alt_col <- snakemake@params[['alt_col']]
p_col <- snakemake@params[['p_col']]
beta_col <- snakemake@params[['beta_col']]
se_col <- snakemake@params[['se_col']]
id_col <- snakemake@params[['id_col']]

gwas <- fread(snakemake@input[['gwas']], sep = '\t', header = T, select = c(chr_col, bp_col, ref_col, alt_col, p_col, beta_col, se_col, id_col))

ld <- fread(snakemake@input[['ld']], sep = '\t', header = T)
