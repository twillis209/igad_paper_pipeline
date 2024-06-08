library(data.table)

coloc_results <- list()

for(x in snakemake@input) {
  res <- readRDS(x)

  window <- regmatches(x, regexpr('[0-9]+kb', x))
  index_snp <- gsub('_', ':', regmatches(x, regexpr('[0-9]+_[0-9]+_[ACTG]+_[ACTG]+', x)))

 dat <- data.table(t(res$summary))
  dat[, `:=` ('window' = window, 'index_snp' = index_snp)]
  coloc_results[[index_snp]] <- dat
}

dat <- rbindlist(coloc_results)
dat$gene <- names(snakemake@params$genes)

hypothesis_pp_cols <- c('PP.H0.abf', 'PP.H1.abf', 'PP.H2.abf', 'PP.H3.abf', 'PP.H4.abf')

dat$PP.largest <- hypothesis_pp_cols[dat[, max.col(.SD), .SDcols = hypothesis_pp_cols]]

fwrite(dat, file = snakemake@output[[1]], sep = '\t')
