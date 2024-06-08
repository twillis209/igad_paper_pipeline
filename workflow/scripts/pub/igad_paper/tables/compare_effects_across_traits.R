library(data.table)

dat <- fread(snakemake@input[[1]])
beta_cols <- c("BETA.igad_meta", "BETA.ra", "BETA.asthma", "BETA.iga")

dat[, gsub('BETA', 'sgn', beta_cols) := lapply(.SD, sign), .SDcols = beta_cols]
