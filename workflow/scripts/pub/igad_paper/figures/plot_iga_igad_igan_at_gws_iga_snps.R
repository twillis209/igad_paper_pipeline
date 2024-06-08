library(data.table)
library(ggplot2)

theme_set(theme_bw()+
          theme(
            axis.title = element_text(size=12),
            axis.text = element_text(size = 10, color = "black"),
            legend.title = element_text(size = 10),
            legend.text = element_text(size = 10)
          )
          )


dat <- fread(snakemake@input[[1]])

pl <- ggplot(dat)+
  geom_point(aes(x = BETA.meta, y = BETA.igan), alpha = 0.5)+
  geom_errorbar(aes(x = BETA.meta, ymin = BETA.igan-1.96*SE.igan, ymax = BETA.igan+1.96*SE.igan))+
  geom_errorbarh(aes(y = BETA.igan, xmin = BETA.meta-1.96*SE.meta, xmax = BETA.meta+1.96*SE.meta))+
  geom_smooth(aes(x = BETA.meta, y = BETA.igan), method = 'lm', formula = y ~ x)+
  coord_fixed()+
  xlim(-0.5, 0.5)+
  ylim(-0.5, 0.5)

ggsave(pl, file = snakemake@output[['igan']], width = 10, height = 10)

pl <- ggplot(dat)+
  geom_point(aes(x = BETA.meta, y = BETA.igan), alpha = 0.5)+
  geom_smooth(aes(x = BETA.meta, y = BETA.igan), method = 'lm', formula = y ~ x)+
  coord_fixed()+
  xlim(-0.5, 0.5)+
  ylim(-0.5, 0.5)

ggsave(pl, file = snakemake@output[['igan_no_ci']], width = 10, height = 10)

pl <- ggplot(dat)+
  geom_point(aes(x = BETA.meta, y = BETA.igad), alpha = 0.5)+
  geom_errorbar(aes(x = BETA.meta, ymin = BETA.igad-1.96*SE.igad, ymax = BETA.igad+1.96*SE.igad))+
  geom_errorbarh(aes(y = BETA.igad, xmin = BETA.meta-1.96*SE.meta, xmax = BETA.meta+1.96*SE.meta))+
  geom_smooth(aes(x = BETA.meta, y = BETA.igad), method = 'lm', formula = y ~ x)+
  coord_fixed()+
  xlim(-0.5, 0.5)+
  ylim(-0.5, 0.5)

ggsave(pl, file = snakemake@output[['igad']])

pl <- ggplot(dat)+
  geom_point(aes(x = BETA.meta, y = BETA.igad), alpha = 0.5)+
  geom_smooth(aes(x = BETA.meta, y = BETA.igad), method = 'lm', formula = y ~ x)+
  coord_fixed()+
  xlim(-0.5, 0.5)+
  ylim(-0.5, 0.5)

ggsave(pl, file = snakemake@output[['igad_no_ci']])

pl <- ggplot(dat[BETA.meta > -0.25])+
  geom_point(aes(x = BETA.meta, y = BETA.igad), alpha = 0.5)+
  geom_smooth(aes(x = BETA.meta, y = BETA.igad), method = 'lm', formula = y ~ x)+
  coord_fixed()+
  xlim(-0.5, 0.5)+
  ylim(-0.5, 0.5)

ggsave(pl, file = snakemake@output[['igad_no_ci_del_outlier']])
