rm(list=ls())
setwd("/Users/abarbour/nute.processing/development/rlpSpec")
source('rsrc/.sourceloads.R')
Xo <- 1e3
Xseq <- 1:Xo
X <- Xseq + rnorm(Xo)
##
load("data/bsm/bsm.rda")
X <- bsm$pinyon_dat$CH0[1:5e3]
##
##
doit1 <- function(){
  ## base test
  initEnv()
  psd0 <- .devpsdcore(X, ntaper=10, plotpsd=T)
  Xseq <- psd0$freq
  #   summary(psd0)
  #
  kopt.slope <- riedsid(psd0$spec, ntaper=10, restrict.deriv="slope")
  kopt.loess <- riedsid(psd0$spec, ntaper=10, restrict.deriv="loess")
  kopt.super <- riedsid(psd0$spec, ntaper=10, restrict.deriv="friedman.super")
  kopt       <- riedsid(psd0$spec, ntaper=10, restrict.deriv="none")
  #   plot(Xseq, 1/sqrt(kopt), type="l", col="red") #,ylim=c(0,100))
  #   lines(Xseq, 1/sqrt(kopt.super), col="blue")
  #   lines(Xseq, 1/sqrt(kopt.loess), col="dark green")
  #   lines(Xseq, 1/sqrt(kopt.slope), col="black")
  plot(Xseq, log10(kopt), type="l", col="red") #,ylim=c(0,100))
  lines(Xseq, log10(kopt.super), col="blue", lwd=3)
  lines(Xseq, log10(kopt.loess), col="dark green", lwd=3)
  lines(Xseq, log10(kopt.slope), col="black")
}
# doit1()
##
##
psd1 <- pspectrum(X, fsamp=20, niter=5, devmode=T)
class(psd1) <- 'spec'

##
##