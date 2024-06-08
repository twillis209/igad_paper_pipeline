library(data.table)
setDTthreads(snakemake@threads)
library(pidPipelineCode)

library(ggplot2)
theme_set(theme_bw()+
          theme(
            plot.title=element_text(hjust=0.5, size=12),
            strip.text=element_text(size=10),
            legend.title=element_text(size=16),
            legend.text=element_text(size=16),
            axis.text.x=element_text(size= 16, color="black"),
            axis.text.y=element_text(size= 16, color="black"),
            axis.title = element_text(size = 20, color = 'black'),
            plot.background = element_rect(color = 'transparent', fill = 'transparent'),
          )
          )

snp_col <- snakemake@params[['snp_col']]
prin_col <- snakemake@params[['prin_col']]
aux_col <- snakemake@params[['aux_col']]

set_one <- fread(snakemake@input[['set_one']], sep = '\t', header = T, select = c(snp_col, prin_col, aux_col))
set_two <- fread(snakemake@input[['set_two']], sep = '\t', header = T, select = c(snp_col, prin_col, aux_col))

set_one[, set := 'One']
set_two[, set := 'Two']

dat <- rbindlist(list(set_one, set_two))

ggplot(dat)+
  geom_qq(aes(sample = get(prin_col), group = as.factor(get(aux_col)), col = as.factor(get(aux_col))), distribution = qunif, size = 0.2)+
  scale_x_neglog10()+
  scale_y_neglog10()+
  geom_abline(intercept = 0, slope = 1, linetype = 2)+
  xlab('Expected')+
  ylab('Observed')+
  scale_color_discrete(name = 'SNP location', labels = c(`0` = 'extragenic', `1` = 'intragenic'))+
  facet_wrap(. ~ set)+
  theme(strip.text = element_blank()) -> pl

ggsave(pl, width = 12, height = 8, file = snakemake@output[[1]])
