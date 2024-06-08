library(data.table)
setDTthreads(snakemake@threads)
library(qvalue)
library(boot)

dat <- fread(snakemake@input[[1]])

bootobj <- boot(dat[in_pid == T, P], statistic = function(x, y) pi0est(x[y])$pi0, R = 100, parallel = 'multicore', ncpus = snakemake@threads)
bootobj <- boot(dat[in_pid == T & !is.na(P), P], function(x, y) pi0est(x[y])$pi0, R = 2, parallel = 'multicore', ncpus = snakemake@threads)
bootobj_non_pid <- boot(dat[in_pid == F & in_non_pid == T & !is.na(P), P], function(x, y) pi0est(x[y])$pi0, R = 1000, parallel = 'multicore', ncpus = snakemake@threads)

saveRDS(bootobj, file = snakemake@output[[1]])
