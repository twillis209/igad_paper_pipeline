library(data.table)
setDTthreads(snakemake@threads)

imd <- fread(snakemake@input$imd)
imd[, CHR38 := as.character(CHR38)]
setkey(imd, CHR38, start, end)

snps <- fread(snakemake@input$snps, header = F)
names(snps) <- 'SNPID'
snps[, c('CHR', 'BP') := tstrsplit(SNPID, split = ':', keep = 1:2)]
snps[, `:=` (BP = as.integer(BP))]
snps[, `:=` (start = BP, end = BP+1)]
setkey(snps, CHR, start, end)

imd_snps <- foverlaps(snps, imd, mult = 'first', nomatch = NA)
imd_snps[, imd := ifelse(!is.na(start), 1, 0)]
fwrite(imd_snps[imd == 0, .(SNPID)], sep = ' ', col.names = F, file = snakemake@output[['out_snps']])
fwrite(imd_snps[imd == 1, .(SNPID)], sep = ' ', col.names = F, file = snakemake@output[['in_snps']])
