library(data.table)
setDTthreads(snakemake@threads)
library(EnsDb.Hsapiens.v86)
library(locuszoomr)
library(ggplot2)
library(gridExtra)

a4_mm_width <- 210
a4_mm_height <- 297

dat <- fread(snakemake@input[[1]])

topGenes <- fread(snakemake@input[['top_genes']])

dat[, SNPID := paste(CHR38, BP38, REF, ALT, sep = ':')]

setnames(dat, c('CHR38', 'BP38'), c('chrom', 'pos'))

p_meta_cols <- snakemake@params$p_meta_cols
meta_subtitles <- snakemake@params$meta_subtitles
p_cfdr_cols <- snakemake@params$p_cfdr_cols
cfdr_subtitles <- snakemake@params$cfdr_subtitles
p_aux_cols <- snakemake@params$p_aux_cols
aux_subtitles <- snakemake@params$aux_subtitles

p_cols <- p_meta_cols
subtitles <- meta_subtitles

if(snakemake@params[['include_cfdr']] == T) {
  p_cols <- c(p_cols, p_cfdr_cols)
  subtitles <- c(subtitles, cfdr_subtitles)
}

if(snakemake@params[['include_aux']] == T) {
  p_cols <- c(p_cols, p_aux_cols)
  subtitles <- c(subtitles, aux_subtitles)
}

gene_id <- topGenes[topGene == snakemake@wildcards[['locus']], topGeneId]

if(length(gene_id) != 1) {
  if(length(unique(topGenes[topGene == snakemake@wildcards[['locus']], topGeneId])) == 1) {
    gene_id <- unique(topGenes[topGene == snakemake@wildcards[['locus']], topGeneId])
  } else {
    stop(sprintf("Unexpected number of gene IDs matching name of locus: %d", length(gene_id)))
  }
}

gene_name <- genes(EnsDb.Hsapiens.v86, filter = GeneIdFilter(gene_id))$gene_name

if(!is.null(snakemake@params$seqname)) {
  loci <- lapply(p_cols, function(x) locus(data = dat, p = x, labs = 'SNPID', seqname = snakemake@params$seqname, xrange = as.integer(snakemake@params$xrange), ens_db = 'EnsDb.Hsapiens.v86'))
} else if(!is.null(snakemake@params$index_snp_pos)) {
  xrange <- c(snakemake@params$index_snp_pos - snakemake@params$flank, snakemake@params$index_snp_pos + snakemake@params$flank)
  loci <- lapply(p_cols, function(x) locus(data = dat, p = x, labs = 'SNPID', seqname = snakemake@params$index_snp_seqname, xrange = xrange, ens_db = 'EnsDb.Hsapiens.v86'))
} else {
  loci <- lapply(p_cols, function(x) locus(data = dat, p = x, labs = 'SNPID', gene = gene_name, flank = snakemake@params[['flank']], ens_db = 'EnsDb.Hsapiens.v86'))
}

names(loci) <- subtitles

if(snakemake@params[['include_cfdr']] == T) {
  ymax <- max(sapply(loci[1:6], function(x) max(x$data$logP, na.rm = T)))
} else {
  ymax <- max(sapply(loci[1:3], function(x) max(x$data$logP, na.rm = T)))
}

gene_track_pl <- gg_genetracks(loci[[1]])

pls <- lapply(seq_along(loci), function(i) gg_scatter(loci[[i]])+ggtitle(subtitles[i])+ylim(0, ymax)+geom_hline(yintercept = -log10(5e-8), linetype = 'dashed'))

out_pls <- list()

layout_mat <- matrix(1:4)
out_pls[1:3] <- pls[1:3]
out_pls[[4]] <- gene_track_pl

if(snakemake@params[['include_cfdr']] == T) {
  layout_mat <- cbind(layout_mat, 5:8)
  out_pls[5:7] <- pls[4:6]
  out_pls[[8]] <- gene_track_pl
}

if(snakemake@params[['include_aux']] == T) {
  layout_mat <- cbind(layout_mat, 9:12)
  out_pls[9:11] <- pls[7:9]
  out_pls[[12]] <- gene_track_pl
  a4_mm_width <- a4_mm_width * 1.5
}

title_chr <- loci[[1]]$seqname
title_range_left <- format(loci[[1]]$xrange[1], big.mark = ',')
title_range_right <- format(loci[[1]]$xrange[2], big.mark = ',')

if(!snakemake@params[['gene_track']]) {
  layout_mat <- layout_mat[1:3, ]
}

ggsave(arrangeGrob(grobs = out_pls, layout_matrix = layout_mat, top = sprintf('%s; %s:%s-%s', snakemake@wildcards[['locus']], title_chr, title_range_left, title_range_right)), file = snakemake@output[[1]], width = a4_mm_width, height = a4_mm_height, unit = 'mm')
