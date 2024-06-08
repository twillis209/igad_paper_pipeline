library(coloc)
library(data.table)

N0.first_trait <- as.integer(snakemake@params$first_controls)
N1.first_trait <- as.integer(snakemake@params$first_cases)
N0.second_trait <- as.integer(snakemake@params$second_controls)
N1.second_trait <- as.integer(snakemake@params$second_cases)
maf_col <- snakemake@params$maf_col

dat <- fread(snakemake@input[[1]])

dat <- na.omit(dat,
               c(snakemake@params$first_beta, snakemake@params$first_se, snakemake@params$second_beta, snakemake@params$second_se))

if(snakemake@wildcards$first_trait == 'iga') {
  first_dataset <- list(snp = dat$SNPID,
                        beta = dat[[snakemake@params$first_beta]],
                        varbeta = (dat[[snakemake@params$first_se]])^2,
                        type = 'quant',
                        N = max(N0.first_trait,
                                N1.first_trait),
                        sdY = 1)
} else {
  first_dataset <- list(snp = dat$SNPID,
                        beta = dat[[snakemake@params$first_beta]],
                        varbeta = (dat[[snakemake@params$first_se]])^2,
                        type = 'cc',
                        N = N0.first_trait + N1.first_trait,
                        s = N1.first_trait / (N0.first_trait + N1.first_trait))

}

if(snakemake@wildcards$second_trait == 'iga') {
  second_dataset <- list(snp = dat$SNPID,
                          beta = dat[[snakemake@params$second_beta]],
                          varbeta = (dat[[snakemake@params$second_se]])^2,
                          type = 'quant',
                          N = max(N0.second_trait,
                          N1.second_trait),
                          sdY = 1)
} else {
  second_dataset <- list(snp = dat$SNPID,
                          beta = dat[[snakemake@params$second_beta]],
                          varbeta = (dat[[snakemake@params$second_se]])^2,
                          type = 'cc',
                          N = N0.second_trait + N1.second_trait ,
                          s = N1.second_trait / (N0.second_trait + N1.second_trait))
}

saveRDS(coloc.abf(first_dataset, second_dataset, MAF = dat[[snakemake@params$maf_col]]), file = snakemake@output[[1]])
