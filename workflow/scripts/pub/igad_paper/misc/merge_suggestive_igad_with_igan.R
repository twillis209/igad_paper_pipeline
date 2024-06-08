library(data.table)
library(httr)

igad <- fread(snakemake@input$igad)
igan <- fread(snakemake@input$igan)

bp_flank <- snakemake@params$window/2

igad <- igad[P < 1e-5]
igad <- igad[!(CHR == 6 & BP %between% c(24e6, 45e6))]
igad <- igad[, .(rsID.igad = rsID, CHR = CHR, BP.igad = BP, REF.igad = REF, ALT.igad = ALT, topGene.igad = topGene, nearestGene.igad = nearestGene, P.igad = P)]
igad[, `:=` (BP.left.igad = BP.igad - bp_flank, BP.right.igad = BP.igad + bp_flank, CHR = as.character(CHR))]
setkey(igad, CHR, BP.left.igad, BP.right.igad)

igan <- igan[`P-VALUE` <= 5e-8]
igan <- igan[, .(P.igan = `P-VALUE`, CHR = CHR_ID, BP = as.integer(CHR_POS), reported_gene.igan = `REPORTED GENE(S)`, mapped_gene.igan = MAPPED_GENE, rsID.igan = SNPS, OR_beta = `OR or BETA`)]
igan[, `:=` (BP.igan = BP, BP.igan.succ = BP + 1, CHR = as.character(CHR))]
setkey(igan, CHR, BP.igan, BP.igan.succ)

# i.x are IgAD columns
overlaps <- foverlaps(igad, igan[!is.na(BP.igan)])[!is.na(reported_gene.igan), .(rsID.igad, rsID.igan, CHR, BP.igan, BP.igad, topGene.igad, nearestGene.igad, reported_gene.igan, mapped_gene.igan, P.igad, P.igan)]

overlaps <- overlaps[!duplicated(overlaps, by = c('CHR', 'BP.igad'))]

fwrite(overlaps[order(as.integer(CHR), BP.igan)], file = snakemake@output[[1]], sep = '\t')
