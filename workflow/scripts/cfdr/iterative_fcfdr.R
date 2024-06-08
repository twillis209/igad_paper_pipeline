library(data.table)
setDTthreads(snakemake@threads)
library(fcfdr)

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
output_file <- snakemake@output[[1]]
tmpdir <- snakemake@resources[['tmpdir']]

gwas_dat <- fread(gwas_file, sep = '\t', header = T, tmpdir = tmpdir)

gwas_dat[, (chr_col) := as.integer(get(chr_col))]

gwas_dat <- na.omit(gwas_dat, cols = c(chr_col, bp_col, prin_col))

gwas_dat <- unique(gwas_dat, by = c(chr_col, bp_col))

# Transforming p-values of 0 to Z-scores yields Inf values
gwas_dat[get(prin_col) < 1e-300, (prin_col) := 1e-300]

# v_cols is p and initial v_cols
for(i in seq_along(aux_cols)) {

  print(sprintf('Iteration %d', i))

  gwas_dat[get(aux_cols[i]) < 1e-300, (aux_cols[i]) := 1e-300]

  if(i == 1) {
    sub_cols <- c(chr_col, bp_col, ref_col, alt_col, prin_col, aux_cols[i], prune_cols[i])
  } else {
    sub_cols <- c(chr_col, bp_col, ref_col, alt_col, prin_col, v_cols[i], aux_cols[i], prune_cols[i])
  }

  sub_dat <- gwas_dat[!is.na(get(aux_cols[i])), ..sub_cols]

  ind_indices <- sub_dat[!is.na(get(aux_cols[i])) & get(prune_cols[i]) == T, which = T]

  fcfdr_result <- flexible_cfdr(p = sub_dat[!is.na(get(aux_cols[i])), get(v_cols[i])],
                                q = sub_dat[!is.na(get(aux_cols[i])), get(aux_cols[i])],
                                indep_index = ind_indices,
                                maf = NULL,
                                check_indep_cor = F,
                                enforce_p_q_cor = F)

  gwas_dat[, (v_cols[i+1]) := Inf]

  gwas_dat[!is.na(get(aux_cols[i])), (v_cols[i+1]) := fcfdr_result[[1]]$v]

  gwas_dat[is.infinite(get(v_cols[i+1])), (v_cols[i+1]) := get(v_cols[i])]

  if(gwas_dat[is.infinite(get(v_cols[i+1])), .N] > 0) {
    stop(sprintf("%d infinite-valued rows left in %s", gwas_dat[is.infinite(get(v_cols[i+1])), .N], v_cols[i+1]))
  }

  fwrite(gwas_dat, file = output_file, sep = '\t')
}

fwrite(gwas_dat, file = output_file, sep = '\t')
