library(data.table)

dat <- fread(snakemake@input[[1]], sep = '\t')

dat <- dat[, .(SNPID, CHR, BP, nearestGene)]

dat <- dat[nearestGene != 'CFAP91']

dat[, manhattan_label := nearestGene]

dat[nearestGene == 'FAP', manhattan_label := 'IFIH1']
dat[nearestGene == 'HARBI1', manhattan_label := 'ATG13']
dat[nearestGene == 'SUOX', manhattan_label := 'RAB5B']
dat[nearestGene == 'RPL6', manhattan_label := 'ALDH2']

fwrite(dat[, .(SNPID, CHR, BP, manhattan_label)], file = snakemake@output[[1]], sep = '\t')
