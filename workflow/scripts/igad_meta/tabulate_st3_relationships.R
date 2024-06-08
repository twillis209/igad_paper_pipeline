library(data.table)

dat <- fread(snakemake@input[[1]], sep = '\t', header = T)

dat[, `:=` (risk.bronson = OR.bronson >= 1, risk.lim = OR.lim >= 1, risk.ST3 = OR.ST3 >= 1)]

fwrite(dat[, .N, by = .(effect_allele.bronson == allele.ST3, other_allele.bronson == allele.ST3, risk.bronson == risk.ST3)], file = snakemake@output[[1]], sep = '\t')
