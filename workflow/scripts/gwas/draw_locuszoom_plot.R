library(ggplot2)
library(locuszoomr)
library(EnsDb.Hsapiens.v86)
library(data.table)

gwas  <- fread(snakemake@input[['gwas']])

variant_id <- snakemake@wildcards[['variant_id']]

variant <- strsplit(variant_id, split = '_')[[1]]
names(variant) <- c('chr', 'pos', 'ref', 'alt')

has_r2 <- 'r2' %in% names(gwas)

if(has_r2) {
  loc <- locus(data = gwas, chrom = 'CHR38', pos = 'BP38', p = 'P', labs = 'SNPID', ens_db = "EnsDb.Hsapiens.v86", seqname = variant[['chr']], xrange = c(as.integer(variant[['pos']]) - (snakemake@params[['window']]/2), as.integer(variant[['pos']]) + (snakemake@params[['window']]/2)), LD = 'r2', index_snp = paste(variant, ':'))
} else {
  loc <- locus(data = gwas, chrom = 'CHR38', pos = 'BP38', p = 'P', labs = 'SNPID', ens_db = "EnsDb.Hsapiens.v86", seqname = variant[['chr']], xrange = c(as.integer(variant[['pos']]) - (snakemake@params[['window']]/2), as.integer(variant[['pos']]) + (snakemake@params[['window']]/2)), index_snp = paste(variant, ':'))
}


if(snakemake@params[['with_genes']]) {
  pl <- locus_ggplot(loc, showLD = has_r2)
} else {
  pl <- gg_scatter(loc, showLD = has_r2)
}

ggsave(pl, file = snakemake@output[[1]])
