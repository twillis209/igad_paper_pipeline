library(data.table)
library(magrittr)

inputs <- names(snakemake@input)[names(snakemake@input) != ""]

dats <- list()

for(x in inputs[inputs != 'merged_data']) {
  if(x %like% 'cfdr') {
    dat <-  fread(snakemake@input[[x]], sep = '\t', select = c('SNPID', 'rsID', 'CHR', 'BP', 'REF', 'ALT', 'gnomadNFE', 'topGene', 'topGeneId', 'nearestGene', 'nearestGeneDistance', 'mostSevereConsequence', 'v.3'))
    setnames(dat, 'v.3', 'P')
  } else if(x %like% 'meta') {
    dat <-  fread(snakemake@input[[x]], sep = '\t', select = c('SNPID', 'rsID','CHR', 'BP', 'REF', 'ALT', 'gnomadNFE', 'topGene', 'topGeneId', 'nearestGene', 'nearestGeneDistance', 'mostSevereConsequence', 'P', 'BETA', 'SE'))
  }
  dat[, dataset := x]
  dats[[x]] <- dat
}

# We want p-values for each row in this from all traits
by_snp <- rbindlist(dats, fill = T)[order(CHR, BP, topGene)]

all_stats_dat <- fread(snakemake@input[['merged_data']], sep = '\t', header = T)

by_snp <- merge(by_snp, all_stats_dat, all.x = T,
                by.x = c('CHR', 'BP', 'REF', 'ALT'),
                by.y = c('CHR38', 'BP38', 'REF', 'ALT'))
by_snp[, c('BETA', 'P', 'SE') := NULL]

fwrite(by_snp, sep = '\t', file = snakemake@output[['by_snp']])

relabelled_dats <- list()

for(x in inputs[inputs != 'merged_data']) {
  if(x %like% 'cfdr') {
    dat <-  fread(snakemake@input[[x]], sep = '\t', select = c('SNPID', 'rsID', 'CHR', 'BP', 'REF', 'ALT', 'topGene', 'topGeneId', 'nearestGene', 'nearestGeneDistance', 'mostSevereConsequence', 'v.3'))
    setnames(dat, 'v.3', sprintf('P.%s', x))
  } else if(x %like% 'meta') {
    dat <-  fread(snakemake@input[[x]], sep = '\t', select = c('SNPID', 'rsID', 'CHR', 'BP', 'REF', 'ALT', 'topGene', 'topGeneId', 'nearestGene', 'nearestGeneDistance',  'mostSevereConsequence', 'P', 'BETA', 'SE'))
    setnames(dat, 'P', sprintf('P.%s', x))
  }

  dat[, `:=` (BP.left = BP - snakemake@params[['window']]/2, BP.right = BP + snakemake@params[['window']]/2, is_mhc = CHR == 6 & BP %between% c(24e6, 45e6))]

  cols_to_rename <- c('SNPID', 'rsID', 'CHR', 'BP', 'REF', 'ALT', 'BP.left', 'BP.right', 'mostSevereConsequence', 'is_mhc')
  setnames(dat, cols_to_rename, paste(cols_to_rename, x, sep = '.'))
  setkeyv(dat, paste(c('CHR', 'BP.left', 'BP.right'), x, sep = '.'))
  relabelled_dats[[x]] <- dat
}

lapply(relabelled_dats, function(x) x[, .(topGene, topGeneId)]) %>%
  rbindlist %>%
  unique -> topGene_dat

by_topGene <- Reduce(function(x, y) merge(x[topGene != '' & topGeneId != '', -c('nearestGene', 'nearestGeneDistance')], y[topGene != '' & topGeneId != '', -c('nearestGene', 'nearestGeneDistance')], by = c('topGene', 'topGeneId'), all.x = T), relabelled_dats, init = topGene_dat)

by_topGene[!is.na(SNPID.cvid_meta) | !is.na(SNPID.igad_meta) | !is.na(SNPID.pad_meta), is_gwas_gws := T]
by_topGene[is.na(is_gwas_gws), is_gwas_gws := F]

by_topGene[!is.na(SNPID.cvid_cfdr) | !is.na(SNPID.igad_cfdr) | !is.na(SNPID.pad_cfdr), is_cfdr_gws := T]
by_topGene[is.na(is_cfdr_gws), is_cfdr_gws := F]

by_topGene[P.cvid_meta %between% c(5e-8, 1e-5), is_cvid_suggestive := T]
by_topGene[P.igad_meta %between% c(5e-8, 1e-5), is_igad_suggestive := T]
by_topGene[P.pad_meta %between% c(5e-8, 1e-5), is_pad_suggestive := T]

by_topGene[, is_mhc := any(.SD, na.rm = T), .SDcols = names(by_topGene) %like% 'mhc', by = topGene]

by_topGene <- by_topGene[topGene != '' & is_mhc == F]

fwrite(by_topGene[, .SD, .SDcols = names(by_topGene) %like% 'SNPID|rsID|mostSevereConsequence|topGene|is_gwas_gws|is_cfdr_gws|^is_.+_suggestive$|^is_mhc$|^P\\.'], sep = '\t', file = snakemake@output[['by_gene']])
