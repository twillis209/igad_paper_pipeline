library(data.table)
setDTthreads(snakemake@threads)
library(GenomicRanges)
library(magrittr)

gwas_cols <- c("rsID", "CHR", "BP", "REF", "ALT", "gnomadNFE", "topGene", "nearestGene", "mostSevereConsequence", "P")
gwas <- fread(snakemake@input$gwas_lead_snps, select = gwas_cols)

gwas <- gwas[!(CHR == 6 & BP %between% c(24e6, 45e6))]

gwas[, already_reported := topGene %like% paste(snakemake@params$known_gwas_genes, collapse = '|'), by = 1:nrow(gwas)]
gwas[, origin := 'gwas']

setnames(gwas, 'P', 'P.igad')

cfdr_cols <- c("rsID", "CHR", "BP", "REF", "ALT", "gnomadNFE", "topGene", "nearestGene", "mostSevereConsequence", "v.3")
cfdr <- fread(snakemake@input$cfdr_lead_snps, select = cfdr_cols)

setnames(cfdr, 'v.3', 'P.cfdr')

cfdr <- cfdr[topGene %like% paste(snakemake@params$novel_cfdr, collapse = '|')]
cfdr[, `:=` (origin = 'cfdr', already_reported = FALSE)]

merged_data_cols <- c("CHR38", "BP38", "REF", "ALT","P.igad", "BETA.igad", "SE.igad", "P.cfdr", "P.ra", "BETA.ra", "SE.ra", "P.asthma", "BETA.asthma", "SE.asthma", "P.iga", "BETA.iga", "SE.iga")

merged_data <- fread(snakemake@input$merged_data, select = merged_data_cols)
setnames(merged_data, c('CHR38', 'BP38'), c('CHR', 'BP'))

merged_on_gwas <- merge(gwas, merged_data[, -c('P.igad')], by = c('CHR', 'BP', 'REF', 'ALT'), all.x = T)

merged_on_cfdr <- merge(cfdr, merged_data[, -c('P.cfdr')], by = c('CHR', 'BP', 'REF', 'ALT'), all.x = T)

merged <- rbindlist(list(
  merged_on_gwas,
  merged_on_cfdr
),
use.names = T
)

merged[, in_or_near_iei_gene := (topGene %like% paste(snakemake@params$iei_genes, collapse = '|') | nearestGene %like% paste(snakemake@params$iei_genes, collapse = '|'))]

merged[, chosen_gene := snakemake@params$top_gene_to_chosen_gene[[topGene]], by = 1:nrow(merged)]

merged[, maf := min(gnomadNFE, 1-gnomadNFE), by = 1:nrow(merged)]

if(!snakemake@params$full_precision) {
  cols_to_round <- c("maf", "gnomadNFE", "P.igad", "BETA.igad", "SE.igad", "P.cfdr", "P.ra", "BETA.ra", "SE.ra", "P.asthma", "BETA.asthma", "SE.asthma", "P.iga", "BETA.iga", "SE.iga")

  merged[, (cols_to_round) := lapply(.SD, format, digits = 2), .SDcols = cols_to_round]
}

gene_granges <- readRDS(snakemake@input$genes)

window <- snakemake@params$window

neighbour_dat <- GRanges(merged[, .(seqnames = paste0('chr', CHR), start = BP-(window/2), end = BP+(window/2), strand = '*', rsID = rsID)]) %>%
    split(., .$rsID) %>%
    lapply(., function(x) subsetByOverlaps(gene_granges, x)) %>%
    lapply(., function(x) paste(x$symbol, collapse = ', ')) %>%
    data.table(rsID = names(.), neighbours = .)

merged_with_neighbours <- merge(merged, neighbour_dat, by = 'rsID', all.x = T)

merged_with_neighbours[, mostSevereConsequence := gsub('_', ' ', mostSevereConsequence)]
merged_with_neighbours[, mostSevereConsequence := gsub('variant', '', mostSevereConsequence)]
merged_with_neighbours[, mostSevereConsequence := gsub('non coding transcript exon', 'ncRNA exon', mostSevereConsequence)]
merged_with_neighbours[mostSevereConsequence %like% 'downstream|upstream', mostSevereConsequence := 'intergenic']

merged_with_neighbours[origin == 'gwas', origin := 'GWAS']
merged_with_neighbours[origin == 'cfdr', origin := 'cFDR']

merged_with_neighbours[, OR.igad := exp(BETA.igad)]

cols <- c("rsID", "CHR", "BP", "REF", "ALT", "maf", "gnomadNFE", "chosen_gene", "mostSevereConsequence", "neighbours", "already_reported", "origin", "in_or_near_iei_gene", "P.igad", "OR.igad", "BETA.igad", "SE.igad", "P.cfdr", "P.ra", "BETA.ra", "SE.ra", "P.asthma", "BETA.asthma", "SE.asthma", "P.iga", "BETA.iga", "SE.iga")

fwrite(merged_with_neighbours[, ..cols], file = snakemake@output[['all']], sep = '\t')
#fwrite(merged_with_neighbours[, ..cols][origin == 'GWAS'], file = snakemake@output[['gwas_only']], sep = '\t')

# Get min. auxiliary p-value in 200kb neighbourhood
cfdr_dat <- merged_with_neighbours[, ..cols][origin == 'cFDR']
aux_cols <- c('P.ra', 'P.asthma', 'P.iga')

for(i in 1:nrow(cfdr_dat)) {
  chr <- cfdr_dat[i, CHR]
  bp <- cfdr_dat[i, BP]

  window_minima <- merged_data[CHR == chr & BP %between% c(bp - window/2, bp + window/2), lapply(.SD, min, na.rm = T), .SDcols = aux_cols]

  cfdr_dat[i, `:=` (P.ra.min = window_minima$P.ra, P.asthma.min = window_minima$P.asthma, P.iga.min = window_minima$P.iga)]
}

#fwrite(cfdr_dat, file = snakemake@output[['cfdr_only']], sep = '\t')
