library(data.table)
setDTthreads(snakemake@threads)

chr_col <- snakemake@params[['chr_col']]
bp_col <- snakemake@params[['bp_col']]
ref_col <- snakemake@params[['ref_col']]
alt_col <- snakemake@params[['alt_col']]
p_col <- snakemake@params[['p_col']]
snp_col <- snakemake@params[['snp_col']]
beta_col <- snakemake@params$beta_cols
se_col <- snakemake@params$se_cols
p_col <- snakemake@params$p_cols
maf_col <- snakemake@params$maf_col
window <- as.numeric(snakemake@params[['window']])
variant_id <- snakemake@params[['variant_id']]

variant <- strsplit(variant_id, split = ':')[[1]]
names(variant) <- c('chr', 'pos', 'ref', 'alt')

cols <- c(snp_col, chr_col, bp_col, ref_col, alt_col, beta_col, se_col, p_col, maf_col)

dat <- fread(snakemake@input[[1]], sep = '\t', select = cols)

if(dat[get(snp_col) == variant_id, .N] == 0) {
  warning("Couldn't find variant in summary statistics")
  variant_id <- dat[get(chr_col) == as.integer(variant[['chr']])][order(abs(get(bp_col) - as.integer(variant[['pos']])))][1, get(snp_col)]
}

variant_bp <- dat[get(snp_col) == variant_id, get(bp_col)]
variant_chr <- dat[get(snp_col) == variant_id, get(chr_col)]

dat <- dat[get(chr_col) == variant_chr]

dat <- dat[get(bp_col) %between% c(variant_bp-window/2, variant_bp+window/2)]

fwrite(dat[, ..cols], file = snakemake@output[['sum_stats']], sep = '\t')

fwrite(dat[, .(get(snp_col))], file = snakemake@output[['ids']], sep = ' ')
