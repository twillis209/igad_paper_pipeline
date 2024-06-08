library(ggplot2)
library(locuszoomr)
library(EnsDb.Hsapiens.v86)
library(data.table)

gwas  <- fread(snakemake@input[['gwas']])

gwas <- gwas[!is.na(P.iga)]

variant_id <- snakemake@wildcards[['variant_id']]

variant <- strsplit(variant_id, split = '_')[[1]]
names(variant) <- c('chr', 'pos', 'ref', 'alt')
colon_variant_id <- paste(variant, collapse = ':')

gwas[SNPID == snakemake@params$lead_snp_snpid, SNPID := snakemake@params$lead_snp_rsid]

loc <- locus(data = gwas, chrom = 'CHR38', pos = 'BP38', p = 'P.iga', labs = 'SNPID', ens_db = "EnsDb.Hsapiens.v86", seqname = as.integer(variant[['chr']]), LD = 'r2', xrange = snakemake@params$xrange)

pl <- locus_ggplot(loc, showLD = T, labels = c(snakemake@params$lead_snp_rsid),
                   filter_gene_biotype = snakemake@params$gene_biotype,
                   filter_gene_name = snakemake@params$gene_name, min.segment.length = 0, nudge_y = 1.5)

ggsave(pl, file = snakemake@output[[1]])
