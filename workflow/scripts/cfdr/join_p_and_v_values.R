library(data.table)
setDTthreads(snakemake@threads)

inputs <- names(snakemake@input)[names(snakemake@input) != ""]

dats <- list()

stat_cols <- snakemake@params[['stat_cols']]

all_cols <- c('CHR38', 'BP38', 'REF', 'ALT', stat_cols)

for(x in inputs) {
  dat <-  fread(snakemake@input[[x]], sep = '\t')

  if(!(x %like% 'cfdr')) {
    dat <- dat[, ..all_cols]
    setnames(dat, stat_cols, paste(stat_cols, x, sep = '.'))
  } else {
    dat <- dat[, .(CHR38, BP38, REF, ALT, v = v.3)]
    setnames(dat, 'v', paste('P', x, sep = '.'))
  }

  dats[[x]] <- dat
}

fwrite(Reduce(function(x,y) merge(x, y, by = c('CHR38', 'BP38', 'REF', 'ALT'), all = T), dats), file = snakemake@output[[1]], sep = '\t')
