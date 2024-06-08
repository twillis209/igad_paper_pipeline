library(data.table)
setDTthreads(snakemake@threads)

dat <- fread(snakemake@input[[1]])

if(snakemake@params$sans_mhc == T) {
  dat <- dat[!(CHR38 == 6 & BP38 %between% c(24e6, 45e6))]
}

dat[in_pid == T, set := 'pid_only']
dat[in_pid == F & in_non_pid == T, set := 'non_pid']
dat[is.na(set), set := 'intergenic']

gif <- sapply(snakemake@params$percentiles, function(x)  quantile(dat[in_pid == T & !is.na(P), qchisq(P, lower.tail = F, df = 1)], probs = x, na.rm = T)/qchisq(x, df = 1))
fwrite(data.table(stat = gif, percentile = snakemake@params$percentiles), sep = '\t', file = snakemake@output[['iei']])

gif <- sapply(snakemake@params$percentiles, function(x)  quantile(dat[in_pid == F & in_non_pid == T & !is.na(P), qchisq(P, lower.tail = F, df = 1)], probs = x, na.rm = T)/qchisq(x, df = 1))
fwrite(data.table(stat = gif, percentile = snakemake@params$percentiles), sep = '\t', file = snakemake@output[['non_iei']])

gif <- sapply(snakemake@params$percentiles, function(x)  quantile(dat[(in_pid == T | in_non_pid == T) & !is.na(P), qchisq(P, lower.tail = F, df = 1)], probs = x, na.rm = T)/qchisq(x, df = 1))
fwrite(data.table(stat = gif, percentile = snakemake@params$percentiles), sep = '\t', file = snakemake@output[['all_genes']])
