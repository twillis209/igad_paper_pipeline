library(data.table)

dat <- fread(snakemake@input)

cited_study_accessions <- c("GCST011956",
"GCST90132222",
"GCST005531",
"GCST002217",
"GCST005527")

dat <- dat[`P-VALUE` <= 5e-8][order(`P-VALUE`)][, .(rsID = SNPS, P = `P-VALUE`, BP = CHR_POS, OR = `OR or BETA`, trait = `DISEASE/TRAIT`, risk_allele = `STRONGEST SNP-RISK ALLELE`, risk_allele_frequency = `RISK ALLELE FREQUENCY`)]
