library(data.table)

chr_col <- snakemake@params[['chr_col']]
bp_col <- snakemake@params[['bp_col']]
ref_col <- snakemake@params[['ref_col']]
alt_col <- snakemake@params[['alt_col']]
p_col <- snakemake@params[['p_col']]
snp_col <- snakemake@params[['snp_col']]

setDTthreads(snakemake@threads)

prin_dat <- fread(snakemake@input[['principal_trait_gwas_file']], sep = '\t', header = T, select = c(chr_col, bp_col, ref_col, alt_col, p_col, snp_col))

prin_dat <- prin_dat[!is.na(get(p_col))]

prin_dat[, (chr_col) := as.character(get(chr_col))]

for(i in seq_along(snakemake@input[['auxiliary_trait_gwas_files']])) {
  aux_file <- snakemake@input[['auxiliary_trait_gwas_files']][[i]]
  pruned_aux_file <- snakemake@input[['pruned_auxiliary_trait_gwas_files']][[i]]

  aux_dat <- fread(aux_file, sep = '\t', header = T, select = c(chr_col, bp_col, ref_col, alt_col, p_col))

  aux_p_col <- sprintf('%s.%d', p_col, i)

  setnames(aux_dat, p_col, aux_p_col)

  aux_dat[, (aux_p_col) := as.numeric(get(aux_p_col))]

  aux_dat <- aux_dat[!is.na(get(aux_p_col))]

  aux_dat[, (chr_col) := as.character(get(chr_col))]

  aux_dat[, (ref_col) := toupper(get(ref_col))]
  aux_dat[, (alt_col) := toupper(get(alt_col))]

  pruned_dat <- fread(pruned_aux_file, sep = '\t', header = T, select = c(chr_col, bp_col, ref_col, alt_col))

  pruned_dat[, (chr_col) := as.character(get(chr_col))]

  pruned_dat[, prune_in := T]

  aux_dat <- merge(aux_dat, pruned_dat, all.x = T, by = c(chr_col, bp_col, ref_col, alt_col))

  aux_dat[is.na(prune_in), prune_in := F]

  setnames(aux_dat, 'prune_in', sprintf('prune_in.%d', i))

  aux_dat <- na.omit(aux_dat, cols = c(chr_col, bp_col))

  prin_dat <- merge(prin_dat, aux_dat, all.x = T, by = c(chr_col, bp_col))

  aux_p <- paste0('P.', i)
  aux_prune_in <- sprintf('prune_in.%d', i)

  ref_a <- paste0(ref_col, '.x')
  alt_a <- paste0(alt_col, '.x')
  ref_b <- paste0(ref_col, '.y')
  alt_b <- paste0(alt_col, '.y')

  # NB: Working with p-values so no need to flip beta signs
  prin_dat[!(get(ref_a) == get(ref_b) & get(alt_a) == get(alt_b)) &
  !(get(ref_a) == get(alt_b) & get(alt_a) == get(ref_b)), (aux_p) := NA]
  prin_dat[!(get(ref_a) == get(ref_b) & get(alt_a) == get(alt_b)) &
           !(get(ref_a) == get(alt_b) & get(alt_a) == get(ref_b)), (aux_prune_in) := NA]

  prin_dat[, c(paste0(ref_col, '.y'), paste0(alt_col, '.y')) := NULL]

  setnames(prin_dat, c(paste0(ref_col, '.x'), paste0(alt_col, '.x')), c(ref_col, alt_col))

  prin_dat <- na.omit(prin_dat, cols = c(chr_col, bp_col))
}

if(!snakemake@params[['mhc']]) {
  prin_dat <- prin_dat[!(get(chr_col) == 6 & get(bp_col) %between% c(24e6, 45e6))]
}

prin_dat <- unique(prin_dat, by = c(chr_col, bp_col, ref_col, alt_col))

fwrite(prin_dat, file = snakemake@output[[1]], sep = '\t')
