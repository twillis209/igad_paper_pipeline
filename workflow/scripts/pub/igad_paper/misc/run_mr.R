library(data.table)
library(MendelianRandomization)

dat <- fread(snakemake@input[[1]])

dat <- dat[!(get(snakemake@params$snp_col) %in% snakemake@params$snps_to_exclude)]

mr_input_obj <- mr_input(bx = dat[[snakemake@params$x_beta]],
                     bxse = dat[[snakemake@params$x_se]],
                     by = dat[[snakemake@params$y_beta]],
                     byse = dat[[snakemake@params$y_se]],
                     snps = dat[[snakemake@params$snp_col]])

res <- mr_allmethods(mr_input_obj, method = 'all')

fwrite(res@Values, file = snakemake@output$tsv, sep = '\t')

png(snakemake@output$png, units = 'in', width = 6, height = 4, res = 300)
mr_plot(res)
dev.off()

png(snakemake@output$png_ivw, units = 'in', width = 6, height = 4, res = 300)
mr_plot(mr_input_obj, line = 'ivw', interactive = F)
dev.off()
