library(data.table)
setDTthreads(snakemake@threads)
library(pracma)
library(parallel)
library(magrittr)

dat <- fread(snakemake@input[[1]])

dat <- dat[!is.na(P)]

x <- rev(dat[in_pid == T, -log10(ppoints(.N))])
y <- dat[in_pid == T][order(-log10(P), decreasing = F), -log10(P)]

compute_auc <- function(x, y) -polyarea(c(0, x, tail(x, 1)), c(0, y, 0))
compute_aul <- function(x) 0.5 * tail(x,1)^2

in_pid.auc <- compute_auc(x, y)
in_pid.aul <- compute_aul(x)

resample_and_compute_auc <- function(no_of_snps) {
  sub_dat <- dat[!is.na(P) & (in_pid == T | in_non_pid == T)][sample(.N, no_of_snps)][, .(P)]

  x <- rev(sub_dat[, -log10(ppoints(.N))])
  y <- sub_dat[order(-log10(P), decreasing = F), -log10(P)]

  compute_auc(x, y)
}

set.seed(snakemake@params$seed)

res <- unlist(mclapply(1:1000, FUN = function(i) resample_and_compute_auc(no_of_snps = dat[in_pid == T, .N]), mc.cores = snakemake@threads))

save(res, in_pid.auc, file = snakemake@output[[1]])
