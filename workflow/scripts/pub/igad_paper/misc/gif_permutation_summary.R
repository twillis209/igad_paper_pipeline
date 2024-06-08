library(data.table)
library(magrittr)

pattern <- "results/pid/100kb_1000_(\\d+)/(\\d+)_gif_permutations\\.tsv\\.gz"

matches <- lapply(snakemake@input, function(x) regmatches(x, regexec(pattern, x))[[1]]) %>%
  lapply(., '[', 2:3) %>%
  lapply(., function(x) data.frame(gene_ratio = x[1], seed = x[2]))

matches_daf <- Reduce(rbind, matches)

get_p_values <- function(x) {
  dat <- fread(x)
  dat[gif.V1 > snakemake@params$pid_gif, .N]/dat[, .N]
}

p_value_daf <- lapply(snakemake@input, get_p_values) %>%
  Reduce(rbind, .)

fwrite(cbind(matches_daf, p = p_value_daf), file = snakemake@output[[1]], sep = '\t')
