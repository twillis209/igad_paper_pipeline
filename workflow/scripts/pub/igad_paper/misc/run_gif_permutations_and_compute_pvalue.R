library(data.table)
setDTthreads(snakemake@threads)
library(parallel)
library(magrittr)
library(GenomicRanges)

gwas <- fread(snakemake@input$gwas, sep = '\t', select = c('CHR38', 'BP38', 'P'))
gwas <- gwas[!is.na(P)]

if(snakemake@params$sans_mhc == T) {
  gwas <- gwas[!(CHR38 == 6 & BP38 %between% c(24e6, 45e6))]
}

gwas[, CHR38 := as.character(CHR38)]
gwas[, end2 := BP38+1]
setkey(gwas, CHR38, BP38, end2)

all_granges <- readRDS(snakemake@input[['all_genes']])

gene_dat <- data.table(data.frame(all_granges))
gene_dat <- gene_dat[seqnames %like% '^chr\\d+$']
gene_dat[, 'chr' := tstrsplit(seqnames, 'chr', keep = 2)]

sample_genes_and_compute_gif <- function(no_of_genes) {
  sub_gene_dat <- gene_dat[sample(.N, no_of_genes)]
  setkey(sub_gene_dat, chr, start, end)
  sub_gene_overlap <- foverlaps(gwas, sub_gene_dat, mult = 'first', nomatch = NULL)

  sapply(snakemake@params$percentiles, function(x)  quantile(sub_gene_overlap[, qchisq(P, lower.tail = F, df = 1)], probs = x, na.rm = T)/qchisq(x, df = 1))
}

RNGkind("L'Ecuyer-CMRG")

set.seed(snakemake@params$seed)

res <- data.table(gif = Reduce(rbind, mclapply(1:snakemake@params$no_of_permutations, FUN = function(i) sample_genes_and_compute_gif(snakemake@params$no_of_genes), mc.cores = snakemake@threads)))

names(res) <- gsub('\\.', '_', snakemake@params$percentiles)

fwrite(res, sep = '\t', file = snakemake@output[['permutations']])

iei_stats <- fread(snakemake@input[['iei_result']])$stat

pvals <- numeric()

for(i in 1:length(snakemake@params$percentiles)) {
  pvals[i] <- sum(res[[i]] >= iei_stats[i])/length(res[[i]])
}

fwrite(data.table(percentile = snakemake@params$percentiles, pvalue = pvals), file = snakemake@output[['pvalue']])
