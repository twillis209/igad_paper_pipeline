library(data.table)
setDTthreads(snakemake@threads)
library(ggplot2)
library(ggtext)
library(ggrepel)
library(ggnewscale)
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

if(!is.null(snakemake@params$plot_background)) {
  if(snakemake@params$plot_background == F) {
    theme_update(plot.background = element_rect(fill = 'transparent', color = NA))
  }
}

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
  summarise(max_bp = max(bp) ) %>% 
  mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) %>% 
  select(chr, bp_add)

gwas_data <- gwas_dat %>% 
  inner_join(data_cum, by = "chr") %>% 
  mutate(bp_cum = as.numeric(bp + bp_add))

axis_set <- gwas_data %>% 
  group_by(chr) %>% 
  summarize(center = mean(bp_cum))

# Ignore typically large peak from MHC when setting y scale
ylim <- gwas_data[chr != 6] %>% 
  filter(p == min(p)) %>% 
  mutate(ylim = abs(floor(log10(p))) + 2) %>% 
  pull(ylim)

manhplot <- ggplot(gwas_data, aes(x = bp_cum, y = p, color = as_factor(chr))) +
  geom_hline(yintercept = 5e-8, color = "grey40", linetype = "dashed") + 
  geom_point(size = 0.3) +
  scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
  scale_color_manual(values = rep(colorblind_palette, unique(length(axis_set$chr)))) +
  scale_size_continuous(range = c(0.5,3)) +
  scale_y_neglog10(limits = c(1, 1e-20))+
  labs(x = NULL, 
       y = "-log<sub>10</sub>(p)") + 
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_markdown(angle = 0, vjust = 0.5),
    axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
  )+
  ggtitle(snakemake@params[['title']])

if(!is.null(snakemake@input$annotations)) {
  annot_dat <- fread(snakemake@input[['annotations']], sep = '\t')
  annot_dat[, `:=` (CHR = as.integer(Chromosome), BP = as.integer(gsub(',', '', Position)))]

  merged_dat <- merge(y = annot_dat, x = gwas_data, by.y = c('CHR', 'BP'), by.x = c('chr', 'bp'))

  merged_dat <- merged_dat[!(chr == '6' & bp %between% c(24e6, 45e6))]

  merged_dat[, text_colour := 'black']

  if(!is.null(snakemake@params$green_highlights)) {
    merged_dat[Gene %in% snakemake@params$green_highlights, text_colour := 'green']
  }

  manhplot <- manhplot+
    geom_point(size = 0.9, pch = 21, col = 'black', data = merged_dat)+
    new_scale_color()+
    geom_label_repel(aes(label = Gene, col = text_colour), size = 3, data = merged_dat, min.segment.length = 0, nudge_y = 5, angle = 90)+
    scale_color_manual(values =  c('green' = "#009E73", 'black' = '#000000'))
}

ggsave(manhplot, file = snakemake@output[[1]], width = 6, height = 3)
