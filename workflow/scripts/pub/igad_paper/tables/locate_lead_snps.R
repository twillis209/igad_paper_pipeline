library(data.table)
library(GenomicRanges)
library(GenomicFeatures)
library(EnsDb.Hsapiens.v86)
library(magrittr)

edb <- EnsDb.Hsapiens.v86

dat <- fread(snakemake@input[[1]])

dat[, granges_in := paste0('chr', CHR, ':', BP, '-', BP)]

gr <- GRanges(dat$granges_in)
seqlevelsStyle(gr)  <- 'NCBI'
names(gr) <- dat$rsID

gr_list <- split(gr, names(gr))

gr_list_genes <- lapply(gr_list, function(x) genes(edb, filter = AnnotationFilterList(
                                                          GRangesFilter(x),
                                                          GeneIdFilter("ENSG", "startsWith")
                                                        )
                                                   ))

gr_list_exons <- lapply(gr_list, function(x) exons(edb, filter = GRangesFilter(x)))

# NB: unlist and unique to combine introns across transcripts, note also that introns have same name according to transcript so will appear to be duplicated
gr_list_introns <- lapply(gr_list, function(x) intronsByTranscript(edb, filter = GRangesFilter(x)) %>% unlist %>% unique %>% countOverlaps(., x) %>% .[. > 0])

save(gr_list_exons, gr_list_introns, gr_list_genes, file = snakemake@output[[1]])
