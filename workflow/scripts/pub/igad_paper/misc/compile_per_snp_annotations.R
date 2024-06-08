library(data.table)
setDTthreads(snakemake@threads)

snps <- fread(snakemake@input$snps, header = F)
names(snps) <- 'SNPID'

for(x in names(snakemake@input)[!(names(snakemake@input) %in% c('snps', ""))]) {
  ann <- fread(snakemake@input[[x]], header = F)

  names(ann) <- 'SNPID'

  ann[, (x) := T]

  snps <- merge(snps, ann, all.x = T, by= 'SNPID')

  snps[is.na(get(x)), (x) := F]
}

fwrite(snps, sep = '\t', file = snakemake@output[[1]])
