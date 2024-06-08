library(data.table)
setDTthreads(snakemake@threads)
library(EnsDb.Hsapiens.v86)
library(locuszoomr)
library(ggplot2)
library(gridExtra)

a4_mm_width <- 210
a4_mm_height <- 297

dat <- fread(snakemake@input[[1]])

dat[, SNPID := paste(CHR38, BP38, REF, ALT, sep = ':')]

setnames(dat, c('CHR38', 'BP38'), c('chrom', 'pos'))

p_cols <- c('P.igad', 'P.cfdr', 'P.asthma', 'P.ra', 'P.iga')
subtitles <- c('IgAD meta-analysis', 'IgAD cFDR analysis', 'asthma', 'rheumatoid arthritis', 'serum IgA')

if(!is.null(snakemake@params$yfloor)) {
  dat[, (p_cols) := lapply(.SD, function(x) ifelse(x < snakemake@params$yfloor, snakemake@params$yfloor, x)), .SDcols = p_cols]
}

xrange <- c(snakemake@params$index_snp_pos - snakemake@params$flank, snakemake@params$index_snp_pos + snakemake@params$flank)

loci <- lapply(p_cols, function(x) locus(data = dat[!is.na(get(x))], p = x, labs = 'SNPID', seqname = snakemake@params$index_snp_seqname, xrange = xrange, ens_db = 'EnsDb.Hsapiens.v86'))

names(loci) <- subtitles

if(is.null(snakemake@params$ymax)) {
  ymax <- max(sapply(loci[1:5], function(x) max(x$data$logP, na.rm = T)))
} else {
  ymax <- snakemake@params$ymax
}

gene_track_pl <- gg_genetracks(loci[[1]], filter_gene_biotype = 'protein_coding')

pls <- lapply(seq_along(loci), function(i) gg_scatter(loci[[i]], labels = 'index')+ggtitle(subtitles[i])+ylim(0, ymax)+geom_hline(yintercept = -log10(5e-8), linetype = 'dashed'))

out_pls <- list()

layout_mat <- cbind(1:4, 5:8)
out_pls[2:3] <- pls[1:2]
out_pls[[4]] <- gene_track_pl

out_pls[5:7] <- pls[3:5]
out_pls[[8]] <- gene_track_pl

title_chr <- loci[[1]]$seqname
title_range_left <- format(loci[[1]]$xrange[1], big.mark = ',')
title_range_right <- format(loci[[1]]$xrange[2], big.mark = ',')

ggsave(arrangeGrob(grobs = out_pls, layout_matrix = layout_mat, top = sprintf('%s; %s:%s-%s', snakemake@wildcards[['locus']], title_chr, title_range_left, title_range_right), padding = unit(2.5, 'lines')), file = snakemake@output[[1]], width = a4_mm_width, height = a4_mm_height, unit = 'mm')
