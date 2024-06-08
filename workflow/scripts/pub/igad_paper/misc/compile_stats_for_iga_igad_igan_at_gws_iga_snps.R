library(data.table)
setDTthreads(snakemake@threads)

merged <- fread(snakemake@input$merged)

lead_snps <- fread(snakemake@input$iga_lead_snps)

merged2 <- merge(lead_snps, merged[, .(CHR38, BP38, REF, ALT, P.igan, BETA.igan, SE.igan, P.igad, BETA.igad, SE.igad, P.cfdr)], all.x = T, by = c('CHR38', 'BP38', 'REF', 'ALT'))

fwrite(merged2, file = snakemake@output[[1]], sep = '\t')
