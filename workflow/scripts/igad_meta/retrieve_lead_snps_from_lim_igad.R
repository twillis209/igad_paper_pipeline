library(data.table)
setDTthreads(snakemake@threads)

lim <- fread(snakemake@input[[1]])

# In hg19
snps <- data.table(rsID = c("rs1990760", "rs34069391", "rs4565870", "rs7773987"),
                   chr = c(2, 16, 11, 6),
                   pos = c(163124051, 11161215, 46349869, 135707486)
                   )

merged <- merge(lim, snps, by.x = c('CHR', 'POS'), by.y = c('chr', 'pos'))

fwrite(merged[, .(rsID, REF, ALT, BETA, AF1, P)], file = snakemake@output[[1]], sep = '\t')
