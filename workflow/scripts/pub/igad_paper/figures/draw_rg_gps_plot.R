library(data.table)
library(ggplot2)
library(patchwork)
library(magrittr)
library(pidPipelineCode)

theme_set(theme_bw()+
          theme(
            axis.text = element_text(size= 20, color="black"),
            axis.title = element_text(size = 24, color = 'black'),
            legend.text = element_text(size = 20),
            legend.title = element_text(size = 24)
            ))

update_geom_defaults("point", list(size = 5))

if(!is.null(snakemake@params$plot_background)) {
  if(snakemake@params$plot_background == F) {
    theme_update(plot.background = element_rect(fill = 'transparent', color = NA),
                 legend.key = element_rect(fill = 'transparent', color = NA)
                 )
  }
}

rg <- fread(snakemake@input[['rg']])

metadat <- fread(snakemake@input[['metadata']])

gps <- fread(snakemake@input[['gps']])

gps <- gps[trait_B %in% snakemake@params$imd_traits]

rg <- rg[trait.B %in% snakemake@params$imd_traits]

merged <- merge(rg, metadat[, .(abbrv, pretty_name)], by.x = 'trait.B', by.y = 'abbrv', all.x = T)

merged <- merge(merged, gps[, .(trait_B, gps.p = pvalue)], by.x = 'trait.B', by.y = 'trait_B')

merged[, `:=` (rg.fdr = p.adjust(rg.p.sr, method = 'BH'), gps.fdr = p.adjust(gps.p, method = 'BH'))]
merged[, `:=` (is.rg.fdr = rg.fdr <= 0.05, is.gps.fdr = gps.fdr <= 0.05)]

merged[, is.rg.fdr := factor(is.rg.fdr, levels = c(T, F), labels = c('True', 'False'))]

merged[, .(pretty_name, is.rg.fdr, rg.p.sr, rg.sr, rg.se.sr)] %>%
ggplot()+
  geom_point(aes(y = reorder(pretty_name, rg.p.sr), x = rg.sr, col = is.rg.fdr))+
  geom_errorbarh(aes(y = reorder(pretty_name, rg.p.sr), xmin = rg.sr-1.96*rg.se.sr, xmax = rg.sr+1.96*rg.se.sr, col = is.rg.fdr), linewidth = 2)+
  geom_vline(xintercept = 0, linetype = 'dashed')+
  ylab('Trait')+
  xlab('Genetic correlation')+
  xlim(-2, 2)+
  labs(colour = 'FDR < 0.05') -> rg_pl

merged[, .(pretty_name, fdr = gps.fdr, is.gps.fdr, rg.p.sr, gps.p)] %>%
ggplot()+
  geom_point(aes(y = reorder(pretty_name, rg.p.sr), x = gps.p, col = is.gps.fdr))+
  ylab('Trait')+
  xlab('-log10 GPS p-value ')+
  scale_x_neglog10()+
  labs(colour = 'FDR < 0.05')+
  guides(colour = F) -> gps_pl

ggsave(rg_pl+gps_pl+plot_layout(axes = 'collect')+plot_annotation(tag_levels = 'A') & theme(legend.position = 'bottom', plot.tag = element_text(size = 28)), file = snakemake@output[[1]], width = 370, height = 210, unit = 'mm')
