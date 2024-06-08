library(data.table)

igad <- fread(snakemake@input$igad)
iga <- fread(snakemake@input$iga)

igad <- igad[in_or_near_iei_gene == T, .(Gene = chosen_gene, Analysis = origin)]
igad[, Analysis := ifelse(Analysis == 'GWAS', 'SIgAD meta-analysis', 'SIgAD cFDR')]
igad <- igad[order(-Analysis)]

iga <- iga[iei_gene != '', .(Gene = iei_gene, Analysis = 'IgA meta-analysis')]

fwrite(rbindlist(list(igad, iga)), sep = '\t', file = snakemake@output[[1]])
