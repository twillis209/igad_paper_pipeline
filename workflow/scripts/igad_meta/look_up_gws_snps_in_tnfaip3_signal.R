library(data.table)
setDTthreads(snakemake@threads)

dat <- fread(snakemake@input[[1]])

dat <- dat[CHR38 == 6][BP38 %between% c(137833918-1e4, 137833918+1e5)][order(P)][1:9]

dat[, in_wakmar2 := BP38 %between% snakemake@params$wakmar2_interval]
dat[, in_tnfaip3 := BP38 %between% snakemake@params$tnfaip3_interval]

fwrite(dat[order(CHR38, BP38)], file = snakemake@output[[1]], sep = '\t')
