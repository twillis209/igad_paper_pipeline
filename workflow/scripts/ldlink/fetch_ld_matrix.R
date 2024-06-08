library(LDlinkR)
library(data.table)

chr_col <- snakemake@params[['chr_col']]
bp_col <- snakemake@params[['bp_col']]
snp_col <- snakemake@params[['snp_col']]

dat <- fread(snakemake@input[[1]], sep = '\t', header = T)

lead_snp_id <- dat[SNPID == snakemake@params[['variant_id']], rsId]

dat[, ldlink_id := paste0('chr', get(chr_col), ':', get(bp_col))]

ld <- LDmatrix(
  snps = dat$ldlink_id,
  pop = snakemake@params[['population']],
  token = snakemake@params[['ldlink_api_token']],
  genome_build = snakemake@params[['genome_build']]
)

molten_ld <- data.table(melt(ld, id.var = 'RS_number'))

names(molten_ld) <- c('SNP1', 'SNP2', 'r2')

molten_ld <- molten_ld[SNP1 == lead_snp_id]

ld_friends <- merge(dat, molten_ld[, .(SNP2, r2)], by.x = 'rsId', by.y = 'SNP2')

#ld_friends <- merged[r2 > snakemake@params[['ld_friend_threshold']]]

ld_friends[, ldlink_id := NULL]

fwrite(ld_friends, file = snakemake@output[['ld_friends']], sep = '\t')
