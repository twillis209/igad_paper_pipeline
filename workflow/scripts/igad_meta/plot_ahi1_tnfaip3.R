library(data.table)
setDTthreads(snakemake@threads)
library(EnsDb.Hsapiens.v86)
library(locuszoomr)
library(ggplot2)
library(gridExtra)

dat <- fread(snakemake@input[[1]])

dat <- dat[CHR38 == 6 & BP38 %between% c(snakemake@params[['start_bp']], snakemake@params[['stop_bp']])]
setnames(dat, c('CHR38', 'BP38'), c('chrom', 'pos'))

igad_locus <- locus(data = dat, p = 'P', labs = 'SNPID', seqname = snakemake@params[['chrom']], xrange = c(snakemake@params[['start_bp']], snakemake@params[['stop_bp']]), ens_db = 'EnsDb.Hsapiens.v86')

scatter_pl <- gg_scatter(igad_locus)+
  geom_hline(yintercept = -log10(5e-8), linetype = 'dashed')

gene_track_pl <- gg_genetracks(igad_locus)

ggsave(arrangeGrob(grobs = list(scatter_pl, gene_track_pl), top = sprintf('AHI1 and TNFAIP3; %s:%s-%s', snakemake@params[['chrom']], format(snakemake@params[['start_bp']], big.mark = ',', scientific = F), format(snakemake@params[['stop_bp']], big.mark = ',', scientific = F))), file = snakemake@output[[1]])
