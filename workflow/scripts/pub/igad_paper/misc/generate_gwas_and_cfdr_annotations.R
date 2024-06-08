library(data.table)
setDTthreads(snakemake@threads)

flank <- snakemake@params$window/2

snps <- fread(snakemake@input$snps, header = F)
snps[, c('CHR', 'BP') := tstrsplit(V1, split = ':', keep = 1:2)]
snps[, `:=` (CHR = as.character(CHR), start = as.integer(BP), end = as.integer(BP)+1)]
setkey(snps, CHR, start, end)

lead_snps <- fread(snakemake@input$lead_snps)

gwas <- lead_snps[dataset %in% c('both', 'gwas')]
gwas[, c('CHR', 'BP') := tstrsplit(SNPID.gwas, split = ':', keep = 1:2)]
gwas[, `:=` (CHR = as.character(CHR), BP = as.integer(BP))]
gwas[, `:=` (start = BP-flank, end = BP+flank)]
setkey(gwas, CHR, start, end)

cfdr <- lead_snps[dataset %in% c('both', 'cfdr')]
cfdr[, c('CHR', 'BP') := tstrsplit(SNPID.cfdr, split = ':', keep = 1:2)]
cfdr[, `:=` (CHR = as.character(CHR), BP = as.integer(BP))]
cfdr[, `:=` (start = BP-flank, end = BP+flank)]
setkey(cfdr, CHR, start, end)

gwas_snps <- foverlaps(snps, gwas, mult = 'first', nomatch = NA)
gwas_snps[, gwas_snp := ifelse(!is.na(chosen_gene), 1, 0)]
fwrite(gwas_snps[gwas_snp == 0, .(V1)], sep = ' ', col.names = F, file = snakemake@output[['out_gwas']])
fwrite(gwas_snps[gwas_snp == 1, .(V1)], sep = ' ', col.names = F, file = snakemake@output[['in_gwas']])

cfdr_snps <- foverlaps(snps, cfdr, mult = 'first', nomatch = NA)
cfdr_snps[, cfdr_snp := ifelse(!is.na(chosen_gene), 1, 0)]
fwrite(cfdr_snps[cfdr_snp == 0, .(V1)], sep = ' ', col.names = F, file = snakemake@output[['out_cfdr']])
fwrite(cfdr_snps[cfdr_snp == 1, .(V1)], sep = ' ', col.names = F, file = snakemake@output[['in_cfdr']])
