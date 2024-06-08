library(data.table)
library(dplyr)
library(GenomicRanges)

flank <- snakemake@params$window/2

snakemake@input %>%
  lapply(fread) %>%
  rbindlist %>%
  na.omit(c('chromosome_name', 'chromosome_position')) %>%
  unique(by = 'variant_id') -> rbound

rbound[, `:=` (start = chromosome_position - flank, stop = chromosome_position + flank)]

granges <- GRanges(seqnames = rbound$chromosome_name,
                   ranges = IRanges(start = rbound$start, end = rbound$stop),
                   variant_id = rbound$variant_id)

dat <- data.table(data.frame(reduce(granges)))

fwrite(dat[, .(CHR38 = seqnames, start, end)], sep = '\t', file = snakemake@output[[1]])


