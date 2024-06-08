library(data.table)
setDTthreads(snakemake@threads)
library(ggplot2)
library(ggtext)
library(ggrepel)
library(tidyverse)
library(pidPipelineCode)

theme_set(theme_bw()+
          theme(
            axis.title = element_text(size=12),
            plot.title=element_text(hjust=0.5, size=12),
            strip.text=element_text(size=10),
            axis.text.x=element_text(size=6, angle=90, color="black"),
            axis.text.y=element_text(size=10, color="black"),
            legend.title=element_text(size=10),
            legend.text=element_text(size=10)
          )
          )

colorblind_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

chr_col <- 'CHR38'
bp_col <- 'BP38'
p_col <- 'P'
snp_col <- 'SNPID'

gwas_dat <- fread(snakemake@input[['gwas']], sep = '\t', select = c(chr_col, bp_col, p_col, snp_col))

setnames(gwas_dat, c(chr_col, bp_col, p_col), c('chr', 'bp', 'p'))

gwas_dat <- gwas_dat[chr %in% seq(1, 22)]

gwas_dat[, chr := as.integer(chr)]

# From https://danielroelfs.com/blog/how-i-create-manhattan-plots-using-ggplot/

data_cum <- gwas_dat %>% 
  group_by(chr) %>% 
  summarise(max_bp = max(bp)) %>% 
  mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) %>% 
  select(chr, bp_add)

gwas_data <- gwas_dat %>% 
  inner_join(data_cum, by = "chr") %>% 
  mutate(bp_cum = as.numeric(bp + bp_add))

axis_set <- gwas_data %>% 
  group_by(chr) %>% 
  summarize(center = mean(bp_cum))

manhplot <- ggplot(gwas_data, aes(x = bp_cum, y = p, color = as_factor(chr))) +
  geom_hline(yintercept = 5e-8, color = "grey40", linetype = "dashed") + 
  geom_point(size = 0.3) +
  scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
  scale_color_manual(values = rep(colorblind_palette, unique(length(axis_set$chr)))) +
  scale_size_continuous(range = c(0.5,3)) +
  scale_y_neglog10(limits = c(1, 1e-25))+
  labs(x = NULL, 
       y = "-log<sub>10</sub>(p)") + 
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_markdown(),
    axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
  )+
  ggtitle(snakemake@params[['title']])

annot_dat <- fread(snakemake@input[['annotations']], sep = '\t')

merged_dat <- merge(annot_dat, gwas_data, by = 'SNPID')

merged_dat <- merged_dat[!(chr == '6' & bp %between% c(24e6, 45e6))]

manhplot <- manhplot+
  geom_point(size = 0.9, pch = 21, colour = 'black', data = merged_dat)+
  geom_text_repel(aes(label = gene_label), size = 3, data = merged_dat[gene_label == 'RUNX3'], colour = 'black', min.segment.length = 0, nudge_y = -3, angle = 90)+
  geom_text_repel(aes(label = gene_label), size = 3, data = merged_dat[gene_label != 'RUNX3'], colour = 'black', min.segment.length = 0, nudge_y = 10, angle = 90)

ggsave(manhplot, file = snakemake@output[[1]], width = 6, height = 3)
