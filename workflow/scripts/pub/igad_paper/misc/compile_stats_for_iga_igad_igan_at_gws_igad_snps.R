library(data.table)
setDTthreads(snakemake@threads)

lead_snps <- fread(snakemake@input$igad_lead_snps)

igan <- fread(snakemake@input$igan)

merged <- merge(lead_snps, igan[, .(CHR = CHR38, BP = BP38, REF, ALT, P.igan = P, BETA.igan = BETA, SE.igan = SE)], all.x = T, by = c('CHR', 'BP', 'REF', 'ALT'))

fwrite(merged, file = snakemake@output[[1]], sep = '\t')
