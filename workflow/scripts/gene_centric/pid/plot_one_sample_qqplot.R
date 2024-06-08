library(data.table)
library(ggplot2)
library(pidPipelineCode)

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
            aspect.ratio = 1
            )
          )

dat <- fread(snakemake@input[[1]])

if(snakemake@wildcards[['snp_set']] == 'sans_mhc') {
  dat <- dat[!(CHR38 == 6 & BP38 %between% c(24e6, 45e6))]
}

dat[in_pid == T, set := 'pid_only']
dat[in_pid == F & in_non_pid == T, set := 'non_pid']
dat[is.na(set), set := 'intergenic']
dat[, set := as.factor(set)]

axis_max <- 10^-(ceiling(-log10(min(dat[, P]))))

ggplot(dat)+
  geom_qq(aes(sample = P, group = set, col = set), distribution = qunif, size = 0.2)+
  scale_x_neglog10(limits = c(1, 1e-15))+
  scale_y_neglog10(limits = c(1, 1e-15))+
  geom_abline(intercept = 0, slope = 1, linetype = 2)+
  xlab('Expected')+
  ylab('Observed')+
  coord_fixed()+
  scale_color_discrete(name = 'Proximate gene', labels = c('pid_only' = 'PID', 'non_pid' = 'Non-PID', 'extragenic' = 'Extragenic')) -> pl


ggsave(pl, width = 8, height = 8, file = snakemake@output[['one_sample']])
