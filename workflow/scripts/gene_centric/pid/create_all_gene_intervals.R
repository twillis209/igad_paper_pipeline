library(org.Hs.eg.db)
library(data.table)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)

id_col_daf <- select(org.Hs.eg.db,
                    keys=keys(org.Hs.eg.db),
                    columns=c("ENTREZID", "ENSEMBL", "GENENAME", "SYMBOL"))

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# I believe the GENEID column is the Entrez ID
# I checked that strand and TXSTRAND always match
strand_daf <- data.frame(genes(txdb, columns=c("GENEID")))

merged <- data.table(merge(id_col_daf, strand_daf, by.x="ENTREZID", by.y="GENEID"))

granges <- GRanges(seqnames = merged$seqnames,
                   ranges = IRanges(start = merged$start, end = merged$end),
                   strand = merged$strand,
                   ensembl_id = merged$ENSEMBL,
                   gene_name = merged$GENENAME,
                   symbol = merged$SYMBOL)

# This adds window/2 up- and downstream
granges <- resize(granges, width = width(granges) + snakemake@params[['window']], fix = 'center')

saveRDS(granges, file = snakemake@output[[1]])
