library(data.table)
setDTthreads(snakemake@threads)

gwas_cols <- c("rsID", "CHR", "BP", "REF", "ALT", "gnomadNFE", "topGene", "nearestGene", "mostSevereConsequence", "P", "iei_gene")
gwas <- fread(snakemake@input$gwas_lead_snps, select = gwas_cols)
setnames(gwas, 'P', 'P.meta')

gwas <- gwas[!(CHR == 6 & BP %between% c(24e6, 45e6))]
gwas[, `:=` (start = BP, stop = BP+1, CHR = as.character(CHR))]
setnames(gwas, c('CHR', 'BP'), c('CHR38', 'BP38'))
setkey(gwas, CHR38, start, stop)

existing <- fread(snakemake@input$existing_associations, select = c('rsID', 'CHR38', 'BP38', 'genes'))
#existing <- fread(existing_associations, select = c('rsID', 'CHR38', 'BP38', 'genes'))

window <- snakemake@params$window
flank <- window/2
existing[, `:=` (start = BP38 - flank, stop = BP38 + flank, CHR38 = as.character(CHR38))]
setkey(existing, CHR38, start, stop)

overlap <- foverlaps(gwas, existing[, .(rsID, genes, start, stop, CHR38)], mult = 'first', nomatch = NA)
overlap[, already_reported := ifelse(!is.na(start), T, F)]
overlap[, c('rsID', 'genes', 'start', 'stop', 'i.start', 'i.stop') := NULL]
setnames(overlap, 'i.rsID', 'rsID')

# Signal is present in Liu, just lead SNP changed
if(overlap[rsID == 'rs9625935', .N] == 0) stop("Missing lead SNP for LIF/HORMAD2 association which is *not* novel")
overlap[rsID == 'rs9625935', already_reported := T]

merged_data_cols <- c("CHR38", "BP38", "REF", "ALT", "BETA.meta", "SE.meta", "P.meta", "BETA.liu_decode", "SE.liu_decode", "P.liu_decode", "BETA.lyons", "SE.lyons", "P.lyons", "BETA.dennis", "SE.dennis", "P.dennis")

merged_data <- fread(snakemake@input$merged_data, select = merged_data_cols)
#merged_data <- fread(merged_data, select = merged_data_cols)
merged_data[, CHR38 := as.character(CHR38)]

merged <- merge(overlap, merged_data[, -c('P.meta')], by = c('CHR38', 'BP38', 'REF', 'ALT'), all.x = T)
merged[, maf := min(gnomadNFE, 1-gnomadNFE), by = 1:nrow(merged)]

merged[, mostSevereConsequence := gsub('_', ' ', mostSevereConsequence)]
merged[, mostSevereConsequence := gsub('variant', '', mostSevereConsequence)]
merged[, mostSevereConsequence := gsub('non coding transcript exon', 'ncRNA exon', mostSevereConsequence)]
merged[mostSevereConsequence %like% 'downstream|upstream', mostSevereConsequence := 'intergenic']

fwrite(merged, file = snakemake@output[[1]], sep = '\t')
