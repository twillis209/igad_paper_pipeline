library(data.table)
setDTthreads(snakemake@threads)
library(ggplot2)

theme_set(theme_bw()+
          theme(
            axis.title = element_text(size=12),
            axis.title.y = element_text(angle = 0, vjust = 0.5),
            plot.title=element_text(hjust=0.5, size=12),
            strip.text=element_text(size=10),
            axis.text.x=element_text(size=10, color="black"),
            axis.text.y=element_text(size=10, color="black"),
            legend.title=element_text(size=10),
            legend.text=element_text(size=10)
          )
          )

bronson <- fread(snakemake@input[['bronson']], sep = '\t', select = c('CHR38', 'BP38', 'REF', 'ALT', 'BETA', 'SE'))
bronson[, Z.bronson := BETA/SE]
cvid <- fread(snakemake@input[['cvid']], sep = '\t', select = c('CHR38', 'BP38', 'REF', 'ALT', 'BETA', 'SE'))
cvid[, Z.cvid := BETA/SE]
pad <- fread(snakemake@input[['pad']], sep = '\t', select = c('CHR38', 'BP38', 'REF', 'ALT', 'BETA', 'SE'))
pad[, Z.pad := BETA/SE]

merged <- merge(bronson[, .(CHR38, BP38, REF, ALT, Z.bronson)],
                cvid[, .(CHR38, BP38, REF, ALT, Z.cvid)], by = c('CHR38', 'BP38', 'REF', 'ALT'))
merged <- merge(merged,
                pad[, .(CHR38, BP38, REF, ALT, Z.pad)], by = c('CHR38', 'BP38', 'REF', 'ALT'))

molten <- melt(merged, id.vars = c('CHR38', 'BP38', 'REF', 'ALT', 'Z.bronson'))

pl <- ggplot(molten)+
  geom_point(aes(x = Z.bronson, y = value), alpha = 0.2)+
  xlim(c(-20, 20))+
  ylim(c(-20, 20))+
  xlab('Bronson Z-score')+
  ylab('Non-Bronson Z-score')+
  facet_wrap(. ~ variable)

ggsave(pl, width = 8, height = 4, file = snakemake@output[[1]])
