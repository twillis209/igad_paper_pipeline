library(data.table)
setDTthreads(snakemake@threads)

chr_col <- snakemake@params[['chr_col']]
bp_col <- snakemake@params[['bp_col']]
ref_col <- snakemake@params[['ref_col']]
alt_col <- snakemake@params[['alt_col']]
snp_col <- snakemake@params[['snp_col']]
se_col <- snakemake@params[['se_col']]
beta_col <- snakemake@params[['beta_col']]
p_col <- snakemake@params[['p_col']]
li_p_col <- snakemake@params[['li_p_col']]

cols <- c(chr_col, bp_col, ref_col, alt_col, se_col, beta_col, p_col, snp_col, li_p_col)

gwas <- fread(snakemake@input[['gwas']], sep = '\t', header = T, select = cols)

gif_dat <- fread(snakemake@input[['gif']], sep = '\t', header = T)

gif <- gif_dat[, lambda_0_50]

gwas[, (p_col) := pchisq(qchisq(get(p_col), df = 1, lower.tail = F)/gif, df = 1, lower.tail = F)]

fwrite(gwas, file = snakemake@output[[1]], sep = '\t')
