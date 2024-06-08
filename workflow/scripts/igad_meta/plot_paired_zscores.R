library(data.table)
setDTthreads(snakemake@threads)
library(ggplot2)
library(stringr)

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

dat_one <- fread(snakemake@input[[1]])

dat_two <- fread(snakemake@input[[2]])

dat_one[, Z := BETA/SE]
dat_two[, Z := BETA/SE]

dat_one[, SNPID := paste(CHR38, BP38, REF, ALT, sep = ':')]
dat_two[, SNPID := paste(CHR38, BP38, REF, ALT, sep = ':')]

merged <- merge(dat_one, dat_two, by = 'SNPID', suffixes = c(snakemake@params[['suffix_one']], snakemake@params[['suffix_two']]))

pl <- ggplot(merged)+
  geom_point(aes_string(x = paste0("Z", snakemake@params[['suffix_one']]), y = paste0("Z", snakemake@params[['suffix_two']])), alpha = 0.1)+
  xlim(c(-25, 25))+
  ylim(c(-25, 25))+
  xlab(sprintf("%s Z-score", str_to_title(str_replace(snakemake@params[['suffix_one']], '\\.', ''))))+
  ylab(sprintf("%s Z-score", str_to_title(str_replace(snakemake@params[['suffix_two']], '\\.', ''))))+
  geom_hline(yintercept = qnorm(5e-8), linetype = 'dashed')+
  geom_hline(yintercept = -qnorm(5e-8), linetype = 'dashed')+
  geom_vline(xintercept = qnorm(5e-8), linetype = 'dashed')+
  geom_vline(xintercept = -qnorm(5e-8), linetype = 'dashed')+
  coord_fixed()

ggsave(pl, width = 6, height = 4, file = snakemake@output[[1]])
