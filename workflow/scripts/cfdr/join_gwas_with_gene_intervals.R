library(data.table)
setDTthreads(snakemake@threads)

snp_col <- snakemake@params[['snp_col']]
chr_col <- snakemake@params[['chr_col']]
bp_col <- snakemake@params[['bp_col']]
ref_col <- snakemake@params[['ref_col']]
alt_col <- snakemake@params[['alt_col']]
beta_col <- snakemake@params[['beta_col']]
se_col <- snakemake@params[['se_col']]
p_col <- snakemake@params[['p_col']]

cols <- c(snp_col, chr_col, bp_col, ref_col, alt_col, beta_col, se_col, p_col)

gwas <- fread(snakemake@input[['gwas']], sep = '\t', header = T, select = cols)
gwas[, CHR38 := as.character(CHR38)]
gwas[, end := BP38+1]

if(!snakemake@params[['mhc']]) {
  gwas <- gwas[!(get(chr_col) == 6 & get(bp_col) %between% c(24e6, 45e6))]
}

gene_intervals <- fread(snakemake@input[['gene_intervals']], sep = '\t', header = T)
gene_intervals[, chr := as.character(chr)]

setkey(gwas, CHR38, BP38, end)
setkey(gene_intervals, chr, start, end)

overlap_res <- foverlaps(gwas, gene_intervals)

overlap_res <- overlap_res[!is.na(start)]

overlap_res[, c('chr', 'start', 'end', 'width', 'strand', 'i.end') := NULL]

fwrite(overlap_res, file = snakemake@output[['inner']], sep = '\t')

gwas[, end := NULL]

gwas[overlap_res, q.1 := 1, on = c('CHR38', 'BP38', 'REF', 'ALT')]
gwas[is.na(q.1), q.1 := 0]

fwrite(gwas, file = snakemake@output[['gwas_with_binary_covariate']], sep = '\t')
