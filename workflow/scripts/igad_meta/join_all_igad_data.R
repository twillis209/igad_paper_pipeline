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

bronson <- fread(snakemake@input[['bronson']], sep = '\t', header = T, select = cols)
setnames(bronson, c('BETA', 'SE', 'P'), c('BETA.bronson', 'SE.bronson', 'P.bronson'))
bronson <- unique(bronson, by = c(chr_col, bp_col, ref_col, alt_col))

finngen <- fread(snakemake@input[['finngen']], sep = '\t', header = T, select = cols)
setnames(finngen, c('BETA', 'SE', 'P'), c('BETA.finngen', 'SE.finngen', 'P.finngen'))
finngen <- unique(finngen, by = c(chr_col, bp_col, ref_col, alt_col))

meta <- fread(snakemake@input[['meta']], sep = '\t', header = T, select = cols)
setnames(meta, c('BETA', 'SE', 'P'), c('BETA.meta', 'SE.meta', 'P.meta'))
meta <- unique(meta, by = c(chr_col, bp_col, ref_col, alt_col))

bronson[, (chr_col) := as.character(get(chr_col))]
finngen[, (chr_col) := as.character(get(chr_col))]
meta[, (chr_col) := as.character(get(chr_col))]

merged <- merge(bronson, finngen, by = c(chr_col, bp_col, ref_col, alt_col), all.x = T)
merged <- merge(merged, meta, by = c(chr_col, bp_col, ref_col, alt_col))

merged[, SNPID := paste(get(chr_col), get(bp_col), get(ref_col), get(alt_col), sep = ':')]

fwrite(merged, file = snakemake@output[[1]], sep = '\t')
