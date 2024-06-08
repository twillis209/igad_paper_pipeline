library(data.table)
library(GenomicRanges)

flank <- snakemake@params$window/2

dat <- fread(snakemake@input[[1]])

dat[, `:=` (start = as.numeric(BP38) - flank, end = as.numeric(BP38) + flank)]
dat <- dat[!is.na(start) & !is.na(end)]

granges <- GRanges(seqnames = dat$CHR38,
                   ranges = IRanges(start = dat$start, end = dat$end),
                   variant_id = dat$rsID)

dat <- data.table(data.frame(reduce(granges)))

fwrite(dat[, .(CHR38 = seqnames, start, end)], sep = '\t', file = snakemake@output[[1]])
