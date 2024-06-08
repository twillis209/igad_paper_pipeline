library(data.table)
library(dplyr)
library(GenomicRanges)

flank <- snakemake@params$window/2

fread(snakemake@input[[1]]) %>%
  .[!((`MAPPED_TRAIT` %like% ',')|(`DISEASE/TRAIT` %like% 'MTAG'))] %>%
  .[`P-VALUE` <= 5e-8] %>%
  unique(by = 'SNPS') %>%
  .[, .(SNPS, CHR_ID, CHR_POS)] -> dat

dat[, `:=` (start = as.numeric(CHR_POS) - flank, stop = as.numeric(CHR_POS) + flank)]
dat <- dat[!is.na(start) & !is.na(stop)]

granges <- GRanges(seqnames = dat$CHR_ID,
          ranges = IRanges(start = dat$start, end = dat$stop),
          variant_id = dat$SNPS)

dat <- data.table(data.frame(reduce(granges)))

fwrite(dat[, .(CHR38 = seqnames, start, end)], sep = '\t', file = snakemake@output[[1]])
