library(data.table)

dat <- fread(snakemake@input[[1]])

dat_sub <- dat[, .(rsID, REF.meta = REF, ALT.meta = ALT, BETA.meta = BETA.igad_meta, P.meta = P.igad_meta, REF.bronson = effect_allele.bronson, ALT.bronson = other_allele.bronson, BETA.bronson = exp(odds_ratio.bronson), P.bronson = p_value.bronson, REF.finngen = ref.finngen, ALT.finngen = alt.finngen, P.finngen = pval.finngen, BETA.finngen = beta.finngen, ALT.liu = A1.liu, REF.liu = A2.liu, BETA.liu, P.liu, ALT.dennis = A1.dennis, REF.dennis = A2.dennis, BETA.dennis, P.dennis = p_value.dennis, REF.lyons, ALT.lyons, BETA.lyons = beta.lyons, P.lyons = p_value.lyons, ALT.ra= effect_allele.ra, REF.ra =  other_allele.ra, BETA.ra = beta.ra, P.ra = p_value.ra, REF.asthma = ref.asthma, ALT.asthma = alt.asthma, BETA.asthma = beta_EUR.asthma, P.asthma)]

beta_cols <- names(dat_sub)[names(dat_sub) %like% 'BETA']
dat_sub[, gsub('BETA', 'SIGN', beta_cols) := lapply(.SD, sign), .SDcols = beta_cols]

# TODO where's asthma
molten <- dat_sub[, .SD, .SDcols = names(dat_sub) %like% 'rsID|REF|ALT|BETA|SIGN|^P'] %>%
  melt(measure.vars = patterns('^REF')) %>%
  .[, variable := gsub('^REF\\.', '', variable)] %>%
  setnames(c('variable', 'value'), c('dataset', 'REF'))

molten2 <- melt(molten, measure.vars = patterns('^ALT')) %>%
  .[, variable := gsub('^ALT\\.', '', variable)] %>%
  .[dataset == variable] %>%
  .[, variable := NULL] %>%
  setnames(c('value'), c('ALT'))

molten3 <- melt(molten2, measure.vars = patterns('^SIGN')) %>%
  .[, variable := gsub('^SIGN\\.', '', variable)] %>%
  .[dataset == variable] %>%
  .[, variable := NULL] %>%
  setnames(c('value'), c('SIGN'))

molten4 <- melt(molten3, measure.vars = patterns('^BETA')) %>%
  .[, variable := gsub('^BETA\\.', '', variable)] %>%
  .[dataset == variable] %>%
  .[, variable := NULL] %>%
  setnames(c('value'), c('BETA'))

molten5 <- melt(molten4, measure.vars = patterns('^P')) %>%
  .[, variable := gsub('^P\\.', '', variable)] %>%
  .[dataset == variable] %>%
  .[, variable := NULL] %>%
  setnames(c('value'), c('P'))

dataset_and_pheno <- data.table(dataset = c('meta', 'bronson', 'finngen', 'liu', 'dennis', 'lyons', 'ra', 'asthma'), phenotype = c('IgAD', 'IgAD', 'IgAD', 'IgA', 'IgA', 'IgA', 'RA', 'asthma'))

merged <- merge(molten5, dataset_and_pheno, by = 'dataset')

fwrite(merged[, .(rsID, dataset, phenotype, REF, ALT, SIGN, BETA, P)][order(rsID, phenotype)], file = snakemake@output[[1]], sep = '\t')

