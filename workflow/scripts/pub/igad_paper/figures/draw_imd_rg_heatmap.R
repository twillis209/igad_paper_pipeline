library(data.table)
library(pheatmap)
library(ggplot2)
library(ggcorrplot)
library(magrittr)

theme_set(theme_bw()+
          theme(
            axis.text = element_text(size= 20, color="black"),
            axis.title = element_text(size = 24, color = 'black'),
            legend.text = element_text(size = 20),
            legend.title = element_text(size = 24)
          ))

rg <- fread(snakemake@input[['rg']])
meta <- fread(snakemake@input[['metadata']], sep = '\t', header = T, select = c('abbrv', 'pretty_name', 'N1'))

meta <- meta[abbrv %in% snakemake@params$traits]

rg <- rg[trait_A %in% snakemake@params$traits & trait_B %in% snakemake@params$traits]

rg <- merge(rg, meta[, .(abbrv, pretty_A = pretty_name)], by.x = 'trait_A', by.y = 'abbrv')
rg <- merge(rg, meta[, .(abbrv, pretty_B = pretty_name)], by.x = 'trait_B', by.y = 'abbrv')

pretty_names <- unique(c(rg$pretty_A, rg$pretty_B))

rg[is.na(rg.sr), rg.sr := 0]

rbind(rg[, .(pretty_A, pretty_B, rg.sr)],
      rg[, .(pretty_A = pretty_B, pretty_B = pretty_A, rg.sr)]) %>%
  dcast('pretty_A ~ pretty_B', fill = 1, value.var = 'rg.sr') %>%
  .[order(match(pretty_A, pretty_names))] %>%
  .[, ..pretty_names] %>%
  as.matrix %>%
  set_rownames(pretty_names) -> rg_mat

rbind(rg[, .(pretty_A, pretty_B, rg.p.sr)],
      rg[, .(pretty_A = pretty_B, pretty_B = pretty_A, rg.p.sr)]) %>%
  dcast('pretty_A ~ pretty_B', fill = 1, value.var = 'rg.p.sr') %>%
  .[order(match(pretty_A, pretty_names))] %>%
  .[, ..pretty_names] %>%
  as.matrix %>%
  set_rownames(pretty_names) -> rg_p_mat

rg_p_mat[is.na(rg_p_mat)] <- 1

phmap <- pheatmap::pheatmap(rg_mat, silent = T, cutree_rows = snakemake@params$no_of_clusters, cutree_cols = snakemake@params$no_of_clusters)

clustered_row_order <- phmap$tree_row$order

#sample_size_dat <- meta[pretty_name %in% pretty_names][, .(abbrv, pretty_name, N1)][order(match(pretty_name, pretty_names))][clustered_row_order]
#
#annotation_row_daf <- sample_size_dat[, .(N1)]
#rownames(annotation_row_daf) <- sample_size_dat[, pretty_name]

png(filename = snakemake@output[['pheatmap']], width = 8, height = 6, res = 400, units = 'in', type = 'cairo')
phmap <- pheatmap::pheatmap(rg_mat, cellwidth = 20, cellheight = 20, breaks = seq(-1, 1, length.out = 100), cutree_rows = snakemake@params$no_of_clusters, cutree_cols = snakemake@params$no_of_clusters, display_numbers = T)
dev.off()

pl <- ggcorrplot(rg_mat, p.mat = rg_p_mat, insig = 'blank', lab = T, hc.order = T, ggtheme = theme_get())

ggsave(pl, file = snakemake@output[['ggcorr']], width = 12, height = 12)
