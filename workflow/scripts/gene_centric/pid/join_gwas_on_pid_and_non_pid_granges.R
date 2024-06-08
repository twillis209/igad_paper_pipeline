library(data.table)
setDTthreads(snakemake@threads)
library(GenomicRanges)

gwas <- fread(snakemake@input[['gwas']], select = c('SNPID', 'CHR38', 'BP38', 'REF', 'ALT', 'P'))
gwas[, CHR38 := as.character(CHR38)]
gwas[, end2 := BP38+1]
setkey(gwas, CHR38, BP38, end2)

pid_granges <- readRDS(snakemake@input[['pid_genes']])
all_granges <- readRDS(snakemake@input[['all_genes']])

pid_granges <- reduce(pid_granges)
all_granges <- reduce(all_granges)

non_pid_gene_granges <- setdiff(all_granges, pid_granges)

non_pid_dat <- data.table(data.frame(non_pid_gene_granges))
# Remove non-canonical autosome assemblies
non_pid_dat <- non_pid_dat[seqnames %like% '^chr\\d+$']
non_pid_dat[, 'chr' := tstrsplit(seqnames, 'chr', keep = 2)]
setkey(non_pid_dat, chr, start, end)

pid_dat <- data.table(data.frame(pid_granges))
pid_dat <- pid_dat[seqnames %like% '^chr\\d+$']
pid_dat[, 'chr' := tstrsplit(seqnames, 'chr', keep = 2)]
setkey(pid_dat, chr, start, end)

pid_overlap <- foverlaps(gwas, pid_dat, mult = 'first', nomatch = NA)
pid_overlap[, in_pid := ifelse(!is.na(start), T, F)]
pid_overlap[, c('seqnames', 'start', 'end', 'width', 'strand') := NULL]

pid_non_pid_overlap <- foverlaps(pid_overlap, non_pid_dat, mult = 'first', nomatch = NA)
pid_non_pid_overlap[, in_non_pid := ifelse(!is.na(start), T, F)]
pid_non_pid_overlap[, c('seqnames', 'start', 'end', 'width', 'strand', 'end2') := NULL]

fwrite(pid_non_pid_overlap, file = snakemake@output[[1]], sep = '\t')
