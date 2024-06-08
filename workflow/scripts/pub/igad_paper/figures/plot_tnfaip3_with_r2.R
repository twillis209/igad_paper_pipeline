library(ggplot2)
library(locuszoomr)
library(EnsDb.Hsapiens.v86)
library(data.table)

gwas  <- fread(snakemake@input[['gwas']])

variant_id <- snakemake@wildcards[['variant_id']]

variant <- strsplit(variant_id, split = '_')[[1]]
names(variant) <- c('chr', 'pos', 'ref', 'alt')
colon_variant_id <- paste(variant, collapse = ':')

gwas[SNPID == colon_variant_id, SNPID := snakemake@params$lead_snp_rsID]

loc <- locus(data = gwas, chrom = 'CHR38', pos = 'BP38', p = 'P', ens_db = "EnsDb.Hsapiens.v86", seqname = as.integer(variant[['chr']]), LD = 'r2', xrange = snakemake@params$xrange, index_snp = snakemake@params$lead_snp_rsID, labs = 'SNPID')

ggsave(locus_ggplot(loc, showLD = T, labels = 'index', filter_gene_biotype = snakemake@params$gene_biotype, filter_gene_name = snakemake@params$gene_name, min.segment.length = 0, nudge_y = 1), file = snakemake@output[[1]])
