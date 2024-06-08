library(data.table)
setDTthreads(snakemake@threads)
library(pracma)
library(parallel)
library(magrittr)
library(GenomicRanges)

RNGkind("L'Ecuyer-CMRG")

compute_auc <- function(x, y) -polyarea(c(0, x, tail(x, 1)), c(0, y, 0))
compute_aul <- function(x) 0.5 * tail(x,1)^2

gwas <- fread(snakemake@input$gwas, sep = '\t', select = c('CHR38', 'BP38', 'P'))
gwas <- gwas[!is.na(P)]
gwas[, CHR38 := as.character(CHR38)]
gwas[, end2 := BP38+1]
setkey(gwas, CHR38, BP38, end2)

all_granges <- readRDS(snakemake@input[['all_genes']])

gene_dat <- data.table(data.frame(all_granges))
gene_dat <- gene_dat[seqnames %like% '^chr\\d+$']
gene_dat[, 'chr' := tstrsplit(seqnames, 'chr', keep = 2)]

sample_genes_and_compute_auc <- function(no_of_genes) {
  sub_gene_dat <- gene_dat[sample(.N, no_of_genes)]
  setkey(sub_gene_dat, chr, start, end)
  sub_gene_overlap <- foverlaps(gwas, sub_gene_dat, mult = 'first', nomatch = NULL)

  x <- rev(sub_gene_overlap[, -log10(ppoints(.N))])
  y <- sub_gene_overlap[order(-log10(P), decreasing = F), -log10(P)]

  gif <- sub_gene_overlap[, median(qchisq(P, lower.tail = F, df = 1), na.rm = T)]/qchisq(0.5, df = 1)

  data.frame(auc = compute_auc(x, y), aul = compute_aul(x), gif = gif)
}

set.seed(snakemake@params$seed)

res <- Reduce(rbind, mclapply(1:1000, FUN = function(i) sample_genes_and_compute_auc(snakemake@params$no_of_genes), mc.cores = snakemake@threads))

fwrite(res, sep = '\t', file = snakemake@output[[1]])
