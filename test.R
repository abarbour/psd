rm(list=ls())
setwd("/Users/abarbour/nute.processing/development/rlpSpec")
source('rsrc/.sourceloads.R')
Xo <- 1e3
X <- 1:Xo + rnorm(Xo)
##
## base test
initEnv()
psd0 <- .devpsdcore(X, plotpsd=F)
summary(psd0)
kopt <- riedsid(psd0$spec, ntaper=10)
##
pspectrum(X, devmode=T)
#