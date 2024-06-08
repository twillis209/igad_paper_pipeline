library(data.table)
setDTthreads(snakemake@threads)

gwas <- fread(snakemake@input$gwas)
setkey(gwas, CHR, start_with_flank, end_with_flank)

cfdr <- fread(snakemake@input$cfdr)
setkey(cfdr, CHR, start_with_flank, end_with_flank)

snps <- fread(snakemake@input$snps, header = F)
names(snps) <- 'SNPID'
snps[, c('CHR', 'BP') := tstrsplit(SNPID, split = ':', keep = 1:2)]
snps[, `:=` (CHR = as.integer(CHR), BP = as.integer(BP))]
snps[, `:=` (start = BP, end = BP+1)]
setkey(snps, CHR, start, end)

gwas_snps <- foverlaps(snps, gwas, mult = 'first', nomatch = NA)
gwas_snps[, gwas_gene := ifelse(!is.na(chosen_gene), 1, 0)]
fwrite(gwas_snps[gwas_gene == 0, .(SNPID)], sep = ' ', col.names = F, file = snakemake@output[['out_gwas']])
fwrite(gwas_snps[gwas_gene == 1, .(SNPID)], sep = ' ', col.names = F, file = snakemake@output[['in_gwas']])

cfdr_snps <- foverlaps(snps, cfdr, mult = 'first', nomatch = NA)
cfdr_snps[, cfdr_gene := ifelse(!is.na(chosen_gene), 1, 0)]
fwrite(cfdr_snps[cfdr_gene == 0, .(SNPID)], sep = ' ', col.names = F, file = snakemake@output[['out_cfdr']])
fwrite(cfdr_snps[cfdr_gene == 1, .(SNPID)], sep = ' ', col.names = F, file = snakemake@output[['in_cfdr']])
