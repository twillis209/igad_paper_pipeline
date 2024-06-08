library(cfdr)
library(data.table)
library(parallel)

gwas_file <- snakemake@input[['gwas_file']]
chr_col <- snakemake@params[['chr_col']]
bp_col <- snakemake@params[['bp_col']]
ref_col <- snakemake@params[['ref_col']]
alt_col <- snakemake@params[['alt_col']]
prin_col <- snakemake@params[['prin_col']]
aux_cols <- snakemake@params[['aux_cols']]
prune_cols <- snakemake@params[['prune_cols']]
v_cols <- snakemake@params[['v_cols']]
v_cols <- c(prin_col, v_cols)
no_of_threads <- snakemake@threads
p_threshold <- snakemake@params[['p_threshold']]
output_file <- snakemake@output[[1]]
use_ldetect_blocks <- snakemake@params[['use_ldetect_blocks']]
tmpdir <- snakemake@resources[['tmpdir']]

setDTthreads(no_of_threads)

gwas_dat <- fread(gwas_file, sep = '\t', header = T, tmpdir = tmpdir)

gwas_dat[, (chr_col) := as.integer(get(chr_col))]

gwas_dat <- na.omit(gwas_dat, cols = c(chr_col, bp_col, prin_col))

gwas_dat <- unique(gwas_dat, by = c(chr_col, bp_col))

# Transforming p-values of 0 to Z-scores yields Inf values
gwas_dat[get(prin_col) < 1e-300, (prin_col) := 1e-300]

# v_cols is p and initial v_cols
for(i in seq_along(aux_cols)) {
  print(sprintf("Iteration %d", i))

  gwas_dat[get(aux_cols[i]) < 1e-300, (aux_cols[i]) := 1e-300]

  if(i == 1) {
    sub_cols <- c(chr_col, bp_col, ref_col, alt_col, prin_col, aux_cols[i], prune_cols[i])
  } else {
    sub_cols <- c(chr_col, bp_col, ref_col, alt_col, prin_col, v_cols[i], aux_cols[i], prune_cols[i])
  }

  sub_dat <- gwas_dat[!is.na(get(aux_cols[i])), ..sub_cols]

  # Estimate joint null distribution
  est_q0_pars <- fit.2g(P = sub_dat[!is.na(get(aux_cols[i])) & get(prune_cols[i]) == T & get(v_cols[i]) > 0.5, get(aux_cols[i])])$pars

  sub_dat[, is_candidate := ifelse(get(prin_col) <= p_threshold, T, F)]

  non_empty_folds <- sub_dat[is_candidate == T, .N, by = .(chrom_col = get(chr_col))][N > 0, chrom_col]

  # Compute L-regions
  v <- mclapply(non_empty_folds, function(x) vl(
                                sub_dat[[v_cols[i]]],
                                sub_dat[[aux_cols[i]]],
                                indices = sub_dat[is_candidate == T & get(chr_col) == x, which = T],
                                mode = 2,
                                fold = sub_dat[get(chr_col) == x, which = T]),
                mc.cores = no_of_threads)

  saveRDS(list(v, sub_dat, i), file = snakemake@log[['v']])

  # Integrate over L-regions to obtain v-values
  for(j in seq_along(non_empty_folds)) {
    tryCatch( {
      sub_dat[is_candidate == T & get(chr_col) == non_empty_folds[j], (v_cols[i+1]) := il(v[[j]], pi0_null = est_q0_pars[1], sigma_null = est_q0_pars[2], distx = "norm")]; },
      error = function (e) {save(list = c(v, j, sub_dat, est_q0_pars, file = snakemake@log[['v']])); stop("tryCatch")})
  }

  sub_cols <- c(chr_col, bp_col, ref_col, alt_col, v_cols[i+1])

  gwas_dat <- merge(gwas_dat, sub_dat[, ..sub_cols], all.x = T, suffixes = c('', ''), by = c(chr_col, bp_col, ref_col, alt_col))

  gwas_dat[is.na(get(v_cols[i+1])), (v_cols[i+1]) := get(v_cols[i])]

  fwrite(gwas_dat, file = output_file, sep = '\t')
}

fwrite(gwas_dat, file = output_file, sep = '\t')
