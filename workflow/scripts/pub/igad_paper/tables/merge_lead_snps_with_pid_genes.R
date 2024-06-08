library(data.table)
setDTthreads(snakemake@threads)
#library(GenomicRanges)
#library(magrittr)

gwas <- fread(snakemake@input[['gwas']])

gwas <- gwas[!(CHR == 6 & BP %between% c(24e6, 45e6))]
gwas[, `:=` (start = BP, stop = BP+1, CHR = as.character(CHR))]
setkey(gwas, CHR, start, stop)

genes <- fread(snakemake@input$genes)

flank <- snakemake@params$window/2

genes[, `:=` (start = start-flank, stop = end+flank, CHR = as.character(chromosome))]
setkey(genes, CHR, start, stop)

overlap <- foverlaps(gwas, genes[, .(CHR, start, stop, iei_gene = gene_name)], mult = 'first', nomatch = NA)
overlap[, c('start', 'stop', 'i.start', 'i.stop') := NULL]

cols <- c("CHR", "SNPID", "rsID", "nearestGene", "nearestGeneDistance", "topGene", "mostSevereConsequence", "gnomadFIN", "gnomadNFE", "gnomadAFR", "gnomadAMR", "gnomadASJ", "gnomadEAS", "gnomadNFEEST", "gnomadNFENWE", "gnomadNFESEU", "associatedTraits", "BP", "REF", "ALT", "P", "SNP", "CHR19", "BP19", "REF19", "ALT19", "iei_gene")

fwrite(overlap[, ..cols], file = snakemake@output[[1]], sep = '\t')
