library(data.table)
setDTthreads(snakemake@threads)
library(ggplot2)
library(ggrepel)

theme_set(theme_bw()+
          theme(
            axis.title = element_text(size=12),
            axis.text = element_text(size = 10, color = "black"),
            legend.title = element_text(size = 10),
            legend.text = element_text(size = 10)
          )
          )

dat <- fread(snakemake@input[[1]])

dat <- dat[!(CHR38 == 6 & BP38 %between% c(24e6, 45e6))]

dat[, `:=` (Z.A = BETA.A/SE.A, Z.B = BETA.B/SE.B)]

pl <- ggplot(dat)+
  geom_point(aes(x = Z.A, y = Z.B), alpha = 0.05)+
  geom_hline(yintercept = qnorm(5e-8), linetype = 'dashed')+
  geom_hline(yintercept = -qnorm(5e-8), linetype = 'dashed')+
  geom_vline(xintercept = qnorm(5e-8), linetype = 'dashed')+
  geom_vline(xintercept = -qnorm(5e-8), linetype = 'dashed')+
  ylim(c(-10, 10))+
  xlim(c(-10, 10))+
  xlab(snakemake@params$xlab)+
  ylab(snakemake@params$ylab)+
  coord_fixed()

if(!is.null(snakemake@input$labels)) {
  labels <- fread(snakemake@input$labels)

  labels <- labels[origin == 'cFDR']

  labels[, SNPID := paste(CHR, BP, REF, ALT, sep = ':')]

  merged <- merge(dat, labels, by.x = 'SNPID.A', by.y = 'SNPID')

  pl <- pl+
    geom_text_repel(aes(x = Z.A, y = Z.B, label = chosen_gene), col = 'blue', data = merged, min.segment.length = 0)+
    geom_point(aes(x = Z.A, y = Z.B), col = 'red', data = merged)
}

ggsave(pl, width = 4, height = 3, file = snakemake@output[[1]])
