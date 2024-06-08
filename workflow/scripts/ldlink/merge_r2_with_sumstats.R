library(data.table)
setDTthreads(snakemake@threads)

var_id_sub <- stringr::str_replace_all(snakemake@wildcards['variant_id'], '_', ':')

gwas <- fread(snakemake@input[['gwas']])

ld <- fread(snakemake@input[['ld']], header = F)

ld_vars <- fread(snakemake@input[['ld_vars']], header = F, col.names = 'ID')

names(ld) <- ld_vars$ID

ld[, ID := ld_vars$ID]

if(ld[ID == var_id_sub, .N] == 0) {
  stop("Lead SNP not present in LD matrix")
}

molten_ld <- melt(ld[ID == var_id_sub], id.vars = 'ID')
molten_ld[, ID := NULL]
names(molten_ld) <- c('SNPID', 'r2')

merged <- merge(gwas, molten_ld, by = 'SNPID', all.x = T)

fwrite(merged, file = snakemake@output[[1]], sep = '\t')
