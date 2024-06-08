library(GenomicRanges)
library(data.table)

# NB: We halve the window first
window <- snakemake@params[['window']]/2.

dat <- fread(snakemake@input[[1]], sep = '\t', header = T)

dat[, chromosome := paste0('chr', chromosome)]
dat[, strand := ifelse(fwdStrand == T, '+', '-')]

dat[, start := start - window]
dat[, end := end + window]

granges <- reduce(
  GRanges(dat$chromosome,
                   IRanges(dat$start, dat$end),
                   strand = dat$strand,
                   metadata = data.frame(gene_name = dat$gene_name)
          )
)

saveRDS(granges, file = snakemake@output[['rds']])

out_dat <- data.table(data.frame(granges))

setnames(out_dat, 'seqnames', 'chr')

out_dat[, chr := stringr::str_replace(chr, 'chr', '')]
out_dat <- out_dat[order(chr)]

fwrite(out_dat, file = snakemake@output[['tsv']], sep = '\t')
