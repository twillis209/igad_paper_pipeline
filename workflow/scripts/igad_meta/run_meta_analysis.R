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

ukb <- fread(snakemake@input[['ukb']], sep = '\t', header = T , select = cols)
setnames(ukb, c('BETA', 'SE', 'P'), c('BETA.ukb', 'SE.ukb', 'P.ukb'))
ukb <- unique(ukb, by = c(chr_col, bp_col, ref_col, alt_col))
ukb[, (chr_col) := as.character(get(chr_col))]

bronson[, (chr_col) := as.character(get(chr_col))]
finngen[, (chr_col) := as.character(get(chr_col))]
ukb[, (chr_col) := as.character(get(chr_col))]
bronson[, (bp_col) := as.character(get(bp_col))]
finngen[, (bp_col) := as.character(get(bp_col))]
ukb[, (bp_col) := as.character(get(bp_col))]
bronson[, (ref_col) := as.character(get(ref_col))]
finngen[, (ref_col) := as.character(get(ref_col))]
ukb[, (ref_col) := as.character(get(ref_col))]
bronson[, (alt_col) := as.character(get(alt_col))]
finngen[, (alt_col) := as.character(get(alt_col))]
ukb[, (alt_col) := as.character(get(alt_col))]

merged <- merge(bronson, finngen, by = c(chr_col, bp_col, ref_col, alt_col), all.x = T)
merged <- merge(merged, ukb, by = c(chr_col, bp_col, ref_col, alt_col), all.x = T)

merged[, `:=` (wt.finngen = SE.finngen^-2, wt.bronson = SE.bronson^-2)]

merged[!is.na(BETA.bronson) & !is.na(BETA.finngen), `:=` (BETA = (BETA.finngen * wt.finngen + BETA.bronson * wt.bronson)/(wt.finngen + wt.bronson), SE = (wt.finngen + wt.bronson)^-0.5)]

merged[is.na(BETA.finngen) & !is.na(BETA.bronson), `:=` (BETA = BETA.bronson, SE = SE.bronson)]

merged[, Z2 := (BETA/SE)^2]

merged[, P := pchisq(Z2, df = 1, lower.tail = F)]

merged[, SNPID := paste(get(chr_col), get(bp_col), get(ref_col), get(alt_col), sep = ':')]

fwrite(merged, sep = '\t', file = snakemake@output[['without_ukb']])

merged[, `:=` (wt.meta = SE^-2, wt.ukb = SE.ukb^-2)]

merged[!is.na(BETA) & !is.na(BETA.ukb), `:=` (BETA = (BETA.ukb * wt.ukb + BETA * wt.meta)/(wt.ukb + wt.meta), SE = (wt.meta + wt.ukb)^-0.5)]

merged[, Z2 := (BETA/SE)^2]

merged[, P := pchisq(Z2, df = 1, lower.tail = F)]

fwrite(merged, sep = '\t', file = snakemake@output[['with_ukb']])
