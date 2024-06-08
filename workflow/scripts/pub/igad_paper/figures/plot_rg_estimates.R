library(data.table)
library(ggplot2)
library(pidPipelineCode)

theme_set(theme_bw()+
          theme(
            axis.text.x=element_text(size= 16, color="black"),
            axis.text.y=element_text(size= 16, color="black"),
            axis.title.y = element_text(size = 20, color = 'black'),
            axis.title.x = element_text(size = 20, color = 'black')
          ))

rg <- fread(snakemake@input$rg)

rg <- rg[!(trait.B %in% c('igan-sans-ic', 'spondylo'))]
metadat <- fread(snakemake@input$metadata)
merged <- merge(rg, metadat[, .(abbrv, pretty_name, N0, N1)], by.x = 'trait.B', by.y = 'abbrv')
merged[N1 > 0, Neff := 2/(1/N0 + 1/N1)]
merged[N1 == 0, Neff := N0]

merged[, `:=` (rg.fdr = p.adjust(rg.p.sr, method = 'BH'))]

rg_pl <- ggplot(merged)+
  geom_point(aes(y = reorder(pretty_name, rg.p.sr), x = rg.sr, col = rg.fdr <= 0.05, size = Neff))+
  geom_errorbarh(aes(y = reorder(pretty_name, rg.p.sr), xmin = rg.sr-1.96*rg.se.sr, xmax = rg.sr+1.96*rg.se.sr, col = rg.fdr <= 0.05))+
  geom_vline(xintercept = 0, linetype = 'dashed')+
  xlim(-3, 3)+
  ylab('Trait')+
  xlab('rg')

ggsave(rg_pl, file = snakemake@output[[1]])
