library(data.table)

flank <- snakemake@params$window/2

gwas <- fread(snakemake@input$gwas)
gwas[, `:=` (start = BP-flank, stop = BP+flank)]
setkey(gwas, CHR, start, stop)

cfdr <- fread(snakemake@input$cfdr)
cfdr[, `:=` (start = BP-flank, stop = BP+flank)]
setkey(cfdr, CHR, start, stop)

cfdr <- cfdr[!(topGene %like% 'CFAP91')]

overlap <- foverlaps(gwas[, .(SNPID, CHR, start, stop, topGene, nearestGene)],
                     cfdr[, .(SNPID, CHR, start, stop, topGene, nearestGene)],
                     mult = 'first', nomatch = NA)

merged <- overlap[!is.na(SNPID), .(SNPID.cfdr = SNPID, SNPID.gwas = i.SNPID, CHR, topGene, nearestGene)]
merged[, dataset := 'both']

cfdr_only <- cfdr[!(SNPID %in% merged$SNPID.cfdr)]
cfdr_only[, dataset := 'cfdr']
cfdr_only <- cfdr_only[, .(SNPID.cfdr = SNPID, CHR, topGene, nearestGene, dataset)]

gwas_only <- gwas[!(SNPID %in% merged$SNPID.gwas)]
gwas_only[, dataset := 'gwas']
gwas_only <- gwas_only[, .(SNPID.gwas = SNPID, CHR, topGene, nearestGene, dataset)]

rbound <- rbindlist(
  list(
    merged,
    gwas_only,
    cfdr_only
  ),
  fill = T
)

rbound[, chosen_gene := snakemake@params$top_gene_to_chosen_gene[[topGene]], by = 1:nrow(rbound)]

fwrite(rbound, sep = '\t', file = snakemake@output[[1]])
