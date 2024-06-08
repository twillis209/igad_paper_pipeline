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

igad_rg <- fread(snakemake@input$igad)
imd_rg <- fread(snakemake@input$imds)
metadat <- fread(snakemake@input$metadata)

igad_rg <- igad_rg[trait.B %in% snakemake@params$imd_traits]

imd_rg <- imd_rg[(trait.A %in% snakemake@params$imd_traits) | (trait.B %in% snakemake@params$imd_traits)][trait.A == 'liu-decode-lyons-dennis-iga' | trait.B == 'liu-decode-lyons-dennis-iga']

igad_rg <- igad_rg[, .(trait.A = 'SIgAD', imd = trait.B, rg.sr, rg.se.sr, rg.p.sr)]

imd_rg[trait.B == 'liu-decode-lyons-dennis-iga', trait.B := trait.A]
imd_rg <- imd_rg[, .(trait.A = 'Serum IgA', imd = trait.B, rg.sr, rg.se.sr, rg.p.sr)]

imd_rg <- imd_rg[imd %in% igad_rg$imd]

imd_rg[, rg.fdr := p.adjust(rg.p.sr, method = 'BH')]
igad_rg[, rg.fdr := p.adjust(rg.p.sr, method = 'BH')]
imd_rg[, is.rg.fdr := rg.fdr <= 0.05]
imd_rg[, is.rg.fdr := factor(is.rg.fdr, levels = c(T, F), labels = c('True', 'False'))]
igad_rg[, is.rg.fdr := rg.fdr <= 0.05]
igad_rg[, is.rg.fdr := factor(is.rg.fdr, levels = c(T, F), labels = c('True', 'False'))]

igad_rg <- igad_rg[imd != 'liu-decode-lyons-dennis-iga']

rbound <- rbindlist(list(imd_rg,
                         igad_rg
                         )
                    )

merged <- merge(rbound, metadat[, .(abbrv, pretty_name)], by.x = 'imd', by.y = 'abbrv', all.x = T)

pretty_name_order <- merged[trait.A == 'SIgAD'][order(rg.p.sr), pretty_name]
merged[, pretty_name := factor(pretty_name, levels = pretty_name_order)]

merged[trait.A == 'SIgAD'] %>%
ggplot()+
  geom_point(aes(y = pretty_name, x = rg.sr, col = is.rg.fdr))+
  geom_errorbarh(aes(y = pretty_name, xmin = rg.sr-1.96*rg.se.sr, xmax = rg.sr+1.96*rg.se.sr, col = is.rg.fdr), linewidth = 2)+
  geom_vline(xintercept = 0, linetype = 'dashed')+
  ylab('Trait')+
  xlab('Genetic correlation')+
  xlim(-2, 2)+
  labs(colour = 'FDR < 0.05') -> igad_pl

merged[trait.A == 'Serum IgA'] %>%
ggplot()+
  geom_point(aes(y = pretty_name, x = rg.sr, col = is.rg.fdr))+
  geom_errorbarh(aes(y = pretty_name, xmin = rg.sr-1.96*rg.se.sr, xmax = rg.sr+1.96*rg.se.sr, col = is.rg.fdr), linewidth = 2)+
  geom_vline(xintercept = 0, linetype = 'dashed')+
  ylab('Trait')+
  xlab('Genetic correlation')+
  xlim(-2, 2)+
  labs(colour = 'FDR < 0.05') -> iga_pl

ggsave(igad_pl + iga_pl +plot_layout(axes = 'collect')+plot_annotation(tag_levels = 'A') & theme(legend.position = 'bottom', plot.tag = element_text(size = 28)), file = snakemake@output[[1]], width = 370, height = 210, unit = 'mm')
