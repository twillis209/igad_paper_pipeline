library(fcfdr)
library(data.table)
setDTthreads(snakemake@threads)

set.seed(snakemake@params[['seed']])
n <- snakemake@params[['no_snps']]

n1p <- n*snakemake@params[['n1_prop']]
zp <- c(rnorm(n1p, sd=5), rnorm(n-n1p, sd=1))
p <- 2*pnorm(-abs(zp))

# generate q
q <- rbinom(n, 1, 0.1)

group <- c(rep("A", n/2), rep("B", n/2)) 

res <- binary_cfdr(p, q, group)

fwrite(data.table(res), file = snakemake@output[[1]], sep = '\t')
