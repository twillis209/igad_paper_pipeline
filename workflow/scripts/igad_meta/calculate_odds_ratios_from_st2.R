library(data.table)

dat <- fread(snakemake@input[[1]])
dat <- dat[Cohort == 'Swedish']

controls_props <- c('AA.controls', 'AB.controls', 'BB.controls')
cases_props <- c('AA.cases', 'AB.cases', 'BB.cases')
dat[, (paste0(controls_props, '.n')) := lapply(.SD, function(x) Total.controls*x*0.01), .SDcols = controls_props]
dat[, (paste0(cases_props, '.n')) := lapply(.SD, function(x) Total.cases*x*0.01), .SDcols = cases_props]

dat[, `:=` (A.cases.n = 2 * AA.cases.n + AB.cases.n,
            B.cases.n = 2 * BB.cases.n + AB.cases.n,
            A.controls.n = 2 * AA.controls.n + AB.controls.n,
            B.controls.n = 2 * BB.controls.n + AB.controls.n
            )]

dat[, al.or.B := (B.cases.n*A.controls.n) / (B.controls.n*A.cases.n)]
dat[, al.or.A := (A.cases.n*B.controls.n) / (A.controls.n*B.cases.n)]
