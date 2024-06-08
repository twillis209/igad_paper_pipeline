library(data.table)
setDTthreads(snakemake@threads)

chr_col <- snakemake@params[['chr_col']]
bp_col <- snakemake@params[['bp_col']]
ref_col <- snakemake@params[['ref_col']]
alt_col <- snakemake@params[['alt_col']]
beta_col <- snakemake@params[['beta_col']]
se_col <- snakemake@params[['se_col']]
p_col <- snakemake@params[['p_col']]

cols <- c(chr_col, bp_col, ref_col, alt_col, beta_col, se_col, p_col)

gwas <- fread(snakemake@input[['gwas']], sep = '\t', header = T, select = cols)
gwas[, end := BP38+1]

gene_intervals <- fread(snakemake@input[['gene_intervals']], sep = '\t', header = T)

setkey(gwas, CHR38, BP38, end)
setkey(gene_intervals, chr, start, end)

overlap_res <- foverlaps(gwas, gene_intervals)

overlap_res <- overlap_res[!is.na(start)]

overlap_res[, c('chr', 'start', 'end', 'width', 'strand', 'i.end') := NULL]

fwrite(overlap_res, file = snakemake@output[[1]], sep = '\t')
