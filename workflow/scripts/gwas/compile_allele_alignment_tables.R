library(data.table)
setDTthreads(snakemake@threads)

results <- data.table(trait = character(), nochange = integer(), ambig = integer(), impossible = integer(), rev = integer(), revcomp = integer(), comp = integer())

for(i in seq_along(snakemake@input)) {
  dat <- fread(snakemake@input[[i]])
  dat <- dcast(dat, . ~ code, value.var = 'N')

  dat[, trait := snakemake@params[['traits']][i]]

  results <- rbind(results, dat, fill = T)
}

fwrite(results, file = snakemake@output[[1]], sep = '\t')
