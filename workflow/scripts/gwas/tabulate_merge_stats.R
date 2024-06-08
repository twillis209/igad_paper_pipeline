library(data.table)
setDTthreads(snakemake@threads)
library(pidPipelineCode)

gwas_file_a <- snakemake@input[['A']]
gwas_file_b <- snakemake@input[['B']]
chr_col <- snakemake@params[['chr_col']]
bp_col <- snakemake@params[['bp_col']]
ref_col <- snakemake@params[['ref_col']]
alt_col <- snakemake@params[['alt_col']]

dat_a <- fread(gwas_file_a, sep = '\t', header = T, select = c(chr_col, bp_col, ref_col, alt_col))
dat_b <- fread(gwas_file_b, sep = '\t', header = T, select = c(chr_col, bp_col, ref_col, alt_col))

dat_a[ , c(ref_col, alt_col) := list(toupper(get(ref_col)), toupper(get(alt_col)))]
dat_a <- dat_a[get(ref_col) %in% c('A','T','C','G') & get(alt_col) %in% c('A','T','C','G')]
dat_a[, (chr_col) := as.character(get(chr_col))]
dat_a <- na.omit(dat_a)

dat_b[ , c(ref_col, alt_col) := list(toupper(get(ref_col)), toupper(get(alt_col)))]
dat_b <- dat_b[get(ref_col) %in% c('A','T','C','G') & get(alt_col) %in% c('A','T','C','G')]
dat_b[, (chr_col) := as.character(get(chr_col))]
dat_b <- na.omit(dat_b)

if(snakemake@params[['join']] == 'inner') {
  merged_dat <- merge(dat_a, dat_b, by = c(chr_col, bp_col), suffixes = c('.A', '.B'))
} else {
  stop(sprintf("Unrecognised join param: %s", join))
}

# Removes the MHC
if(!snakemake@params[['mhc']]) {
  merged_dat <- merged_dat[!(get(chr_col) == 6 & get(bp_col) %between% c(24e6, 45e6))]
}

ref_a <- paste0(ref_col, '.A')
ref_b <- paste0(ref_col, '.B')
alt_a <- paste0(alt_col, '.A')
alt_b <- paste0(alt_col, '.B')

merged_dat[, A.CODE := paste(get(ref_a), get(alt_a), sep = '/')]
merged_dat[, B.CODE := paste(get(ref_b), get(alt_b), sep = '/')]
merged_dat[, code := g.class(A.CODE, B.CODE)]

fwrite(merged_dat, file = snakemake@output[['AB_with_codes']], sep = '\t')

res_dat <- data.table(code = g.class(merged_dat[,A.CODE], merged_dat[, B.CODE]))[, .N, by = code]

fwrite(res_dat, file = snakemake@output[['merge_stats']], sep = '\t')
