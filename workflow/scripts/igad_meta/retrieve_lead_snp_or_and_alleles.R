library(data.table)
setDTthreads(snakemake@threads)
library(magrittr)

lead_snps <- fread(snakemake@input[['lead_snps']], select = c("rsID", "risk_allele", "orPerCopyNum"))

supp_table <- fread(snakemake@input[['supp_table_two']], select = c('rsID', 'Cohort', 'OR', 'Total.controls', 'Total.cases', 'AA.cases', 'BB.cases', 'AA.controls', 'BB.controls', 'AB.controls', 'AB.cases', 'MAF.cases', 'MAF.controls', 'A', 'B'))

ebi_sum_stats <- fread(snakemake@input[['igad']], select = c("hm_rsid", "effect_allele", "hm_effect_allele", "other_allele", "hm_other_allele", "odds_ratio", "hm_odds_ratio"))

lim <- fread(snakemake@input[['lim_lead_snps']])

table_2 <- fread(snakemake@input[['table_two']], sep = '\t', select = c('Variant', 'minor_allele', 'OR'))
setnames(table_2, c('minor_allele', 'OR'), c('minor_allele.T2', 'OR.T2'))

supp_table <- supp_table[Cohort == 'Swedish']
controls_props <- c('AA.controls', 'AB.controls', 'BB.controls')
cases_props <- c('AA.cases', 'AB.cases', 'BB.cases')
supp_table[, (paste0(controls_props, '.n')) := lapply(.SD, function(x) Total.controls*x*0.01), .SDcols = controls_props]
supp_table[, (paste0(cases_props, '.n')) := lapply(.SD, function(x) Total.cases*x*0.01), .SDcols = cases_props]

supp_table[, `:=` (A.cases.n = 2 * AA.cases.n + AB.cases.n,
            B.cases.n = 2 * BB.cases.n + AB.cases.n,
            A.controls.n = 2 * AA.controls.n + AB.controls.n,
            B.controls.n = 2 * BB.controls.n + AB.controls.n
            )]

supp_table[, al.or.A := (A.cases.n*B.controls.n) / (A.controls.n*B.cases.n)]
supp_table[, al.or.B := (B.cases.n*A.controls.n) / (B.controls.n*A.cases.n)]

merge(lead_snps, ebi_sum_stats, by.x = 'rsID', by.y = 'hm_rsid', all.x = T) %>%
merge(supp_table, by.x = 'rsID', by.y = 'rsID', all.x = T) %>%
merge(lim, by.x = 'rsID', by.y = 'rsID', all.x = T) %>%
merge(table_2, by.x = 'rsID', by.y = 'Variant', all.x = T) -> merged

or_cols <- c('orPerCopyNum', 'odds_ratio', 'OR', 'OR.T2', 'al.or.A', 'al.or.B')
merged[, (or_cols) := lapply(.SD, signif, 2), .SDcols = or_cols]

merged <- merged[, .(rsID, risk_allele.ebi_lead = risk_allele, effect_allele.ebi_sum = effect_allele, other_allele.ebi_sum = other_allele, A.supp = A, B.supp = B, effect_allele.lim = ALT, other_allele.lim = REF, minor_allele.T2, OR.ebi_lead = orPerCopyNum, OR.ebi_sum = odds_ratio, OR.supp = OR, OR.A = al.or.A, OR.B = al.or.B, OR.lim = exp(BETA), OR.T2)]

merged <- merged[!is.na(A.supp)]

fwrite(merged, file = snakemake@output[[1]], sep = '\t')
