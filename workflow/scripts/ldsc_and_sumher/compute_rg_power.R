library(data.table)

dat <- fread(snakemake@input[[1]], sep = '\t')

##' NB: code from the GCTA website: https://yanglab.westlake.edu.cn/software/gcta/#GREMLpowercalculator
var_vg_func <- function(N, var_pi=2e-5){
  return(2/(N^2*var_pi))
}

var_rg_func <- function(N1, N2, hsq1, hsq2, rg, rp, overlap=TRUE, var_pi=2e-5){
  if(overlap==T) var_rg=((1-rg*rp)^2+(rg-rp)^2)/(hsq1*hsq2*N1^2*var_pi)
  if(overlap==F) var_rg=(rg^2*(N1^2*hsq1^2+N2^2*hsq2^2)+2*hsq1*hsq2*N1*N2)/(2*hsq1^2*hsq2^2*N1^2*N2^2*var_pi)
  return(var_rg)
}

power_func <- function(ncp, alpha){
  pchisq(qchisq(alpha, df=1,lower.tail=F), ncp=ncp, df=1, lower.tail=F)
}

h2O_func <- function(ncase, ncontrol, K, h2L, var_pi=2e-5){
  n=ncase+ncontrol
  v=ncase/(ncase+ncontrol)
  z=dnorm(qnorm(K))
  c=(K*(1-K))^2/(v*(1-v)*z^2)
  h2O=h2L/c
  var_h2O=var_vg_func(n, var_pi)
  var_h2L=c^2*var_h2O
  return(list(h2L=h2L, var_h2L=var_h2L, h2O=h2O, var_h2O=var_h2O))
}

calcUniQt <- function(
                      n     =1000, 
                      hsq   =0.5, 
                      alpha =0.05,
                      var_pi=2e-5
                      ){
  l <- list()
  var_vg <- var_vg_func(n, var_pi)
  l$se <- sqrt(var_vg)
  l$ncp <- hsq^2/var_vg;
  l$power <- power_func(l$ncp, alpha)
  return(l)
}

calcUniCc <- function(
    ncase    = 1000, 
    ncontrol = 1000, 
    hsq      = 0.5, 
    K        = 0.1, 
    alpha    = 0.05,
    var_pi=2e-5
){
    h <- h2O_func(ncase, ncontrol, K, hsq, var_pi)
    l <- list()
    l$se <- sqrt(h$var_h2L)
    l$ncp <- h$h2L^2/h$var_h2L
    l$power <- power_func(l$ncp, alpha)
    return(l)
}

calcBiQt <- function(
    n1      = 1000, 
    n2      = 1000, 
    hsq1    = 0.5, 
    hsq2    = 0.5, 
    rg      = 0.5, 
    rp      = 0.5, 
    overlap = FALSE, 
    alpha   = 0.05,
    var_pi=2e-5
){
    var_rg <- var_rg_func(n1, n2, hsq1, hsq2, rg, rp, overlap, var_pi)
    l <- list()
    l$se <- sqrt(var_rg)
    l$ncp <- rg^2/var_rg;
    l$power <- power_func(l$ncp, alpha)
    return(l)
}

calcBiCc <- function(
    ncase1    = 1000, 
    ncase2    = 1000, 
    ncontrol1 = 1000, 
    ncontrol2 = 1000, 
    hsq1      = 0.5, 
    hsq2      = 0.5, 
    K1        = 0.1, 
    K2        = 0.1, 
    rg        = 0.5, 
    overlap   = FALSE, 
    alpha     = 0.05,
    var_pi=2e-5
){
    h1 <- h2O_func(ncase1, ncontrol1, K1, hsq1, var_pi)
    h2 <- h2O_func(ncase2, ncontrol2, K2, hsq2, var_pi)
    n1 <- ncase1+ncontrol1
    n2 <- ncase2+ncontrol2
    var_rg <- var_rg_func(n1, n2, h1$h2O, h2$h2O, rg, rg, overlap, var_pi)
    l <- list()
    l$se <- sqrt(var_rg)
    l$ncp <- rg^2/var_rg;
    l$power <- power_func(l$ncp, alpha)
    return(l)
}

calcBiQtCc <- function(
    n        = 1000, 
    ncase    = 1000, 
    ncontrol = 1000, 
    hsq1     = 0.5, 
    hsq2     = 0.5, 
    K        = 0.1, 
    rg       = 0.5, 
    overlap  = FALSE, 
    alpha    = 0.05,
    var_pi=2e-5
){
    h2=h2O_func(ncase, ncontrol, K, hsq2, var_pi)
    n2=ncase+ncontrol
    var_rg=var_rg_func(n, n2, hsq1, h2$h2O, rg, rg, overlap, var_pi)
    l <- list()
    l$se <- sqrt(var_rg)
    l$ncp <- rg^2/var_rg;
    l$power <- power_func(l$ncp, alpha)
    return(l)
}

dat[, calcBiCc(ncase1 = ncase1, ncase2 = ncase2, ncontrol1 = ncontrol1, ncontrol2 = ncontrol2, hsq1 = hsq1, hsq2 = hsq2, K1 = K1, K2 = K2, rg = rg)$se]
dat[, se := calcBiCc(ncase1, ncase2, ncontrol1, ncontrol2, hsq1, hsq2, K1, K2, rg)$se, by = 1:nrow(dat)]
dat[, power := calcBiCc(ncase1, ncase2, ncontrol1, ncontrol2, hsq1, hsq2, K1, K2, rg)$power, by = 1:nrow(dat)]

fwrite(dat, file = snakemake@output[[1]], sep = '\t')
