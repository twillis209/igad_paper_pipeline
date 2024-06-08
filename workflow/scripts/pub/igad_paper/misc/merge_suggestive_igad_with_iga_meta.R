library(data.table)

igad <- fread(snakemake@input$igad)
iga <- fread(snakemake@input$iga)
lit_iga <- fread(snakemake@input$lit_iga)

bp_flank <- snakemake@params$window/2

igad <- igad[P < 1e-5]
igad[, `:=` (BP.igad = BP, BP.left.igad = BP - bp_flank, BP.right.igad = BP + bp_flank, CHR = as.character(CHR), P.igad = P, BETA.igad = BETA, REF.igad = REF, ALT.igad = ALT)]
setkey(igad, CHR, BP.left.igad, BP.right.igad)

iga <- iga[already_reported == F, .(rsID, CHR, BP, REF, ALT, topGene, nearestGene, BETA.meta, P.meta)]
lit_iga <- lit_iga[, .(rsID, CHR = CHR38, BP = BP38, REF = REF38, ALT = ALT38, risk_allele.ebi, effect_allele.liu, topGene = genes, P.ebi = p.ebi, BETA.ebi = beta.ebi, BETA.liu = beta.liu, P.liu = p.liu)]

combined_iga <- rbindlist(list(
  iga,
  lit_iga
), fill = T
)[order(CHR, BP)]

combined_iga <- combined_iga[!(CHR == 6 & BP %between% c(24e6, 45e6))]

combined_iga[, `:=` (BP.iga = BP, BP.iga.succ = BP + 1, CHR = as.character(CHR))]
setkey(combined_iga, CHR, BP.iga, BP.iga.succ)

# i.x are IgAD columns
overlaps <- foverlaps(igad, combined_iga)[!is.na(topGene), .(CHR, BP.iga, BP.igad, topGene.iga = topGene, topGene.igad = i.topGene, nearestGene.igad = i.nearestGene, P.igad, P.meta, P.ebi, P.liu, BETA.igad, BETA.meta, BETA.ebi, BETA.liu, REF.igad, ALT.igad, risk_allele.ebi, effect_allele.liu)]

fwrite(overlaps, file = snakemake@output[[1]], sep = '\t')
