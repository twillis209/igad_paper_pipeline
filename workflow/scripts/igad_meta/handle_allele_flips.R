library(data.table)
setDTthreads(snakemake@threads)
library(ggplot2)

lim <- fread(snakemake@input[['lim']])

lim[, Z := BETA/SE]

lim <- lim[!(code %in% c('impossible', 'revcomp', 'comp'))]

lim <- lim[!(code == 'ambig' & !(REF == ALT.1kG & ALT == REF.1kG))]

lim[code == 'rev', `:=` (REF = REF.1kG, ALT = ALT.1kG, Z = -1 * Z)]
lim[code == 'ambig' & REF == ALT.1kG & ALT == REF.1kG, `:=` (REF = REF.1kG, ALT = ALT.1kG, Z = -1 * Z)]

bronson <- fread(snakemake@input[['bronson']])

bronson[, Z := BETA/SE]

bronson <- bronson[!(code %in% c('ambig', 'impossible'))]

bronson[code %in% c('rev', 'revcomp'), `:=` (REF = REF.1kG, ALT = ALT.1kG, Z = -1 * Z)]
#bronson[code %in% c('rev', 'revcomp'), `:=` (REF = REF.1kG, ALT = ALT.1kG)]

merged <- merge(lim[, .(CHR38, BP38, REF, ALT, Z, code)],
                bronson[, .(CHR38, BP38, REF, ALT, Z, code)],
                by = c('CHR38', 'BP38', 'REF', 'ALT'),
                suffixes = c('.lim', '.bronson')
                )

pl <- ggplot(merged)+
  geom_point(aes(x = Z.bronson, y = Z.lim), alpha = 0.2)+
  xlim(-20, 20)+
  ylim(-20, 20)

ggsave(pl, file = snakemake@output[[1]], width = 7, height = 7)
