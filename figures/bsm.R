##
##  Creates figure(s) in text for CO2 data
##
# setwd("~/kook.processing/developmen/rlpSpec/figures")
source('funcload.R')
load("../data/bsm/bsm.rda")

library(zoo)
library(RColorBrewer)
library(ggplot2)

bsmnm <- as.data.frame(na.approx(read.table('../data/bsm/barbour_bsmnoise_2011.txt', h=T)))
spec <- 10^((bsmnm$P50-5)/10)
psdf <- data.frame(f=bsmnm$freq, psd=spec, src="rlpSpec", niter=0, sta="PBO", cha="")
psdf <- rbind(psdf, data.frame(f=bsmnm$freq, psd=spec, src="rlpSpec", niter=1, sta="PBO", cha=""))
psdf <- rbind(psdf, data.frame(f=bsmnm$freq, psd=spec, src="rlpSpec", niter=2, sta="PBO", cha=""))
psdf <- rbind(psdf, data.frame(f=bsmnm$freq, psd=spec, src="spec.pgram", niter=0, sta="PBO", cha=""))

doBSMspec <- function(dat, sps=20, stacha="", span=NULL, niter=0, doRspec=TRUE){
  ## return a data frame with spectral estimates
  # multitaper spectra
  rlpps <- pspectrum(dat, niter=niter, fsamp=sps, ndec=1, plotpsd=F)
  dfspec <- data.frame(f=rlpps$f, psd=rlpps$psd, 
                       src="rlpSpec", niter=niter, sta=stacha[1], cha=stacha[2])
  # r built-in (optional)
  if (doRspec || niter==0){
    cat("\t>>>> Doing spec.pgram\n")
    pad<-1
    tap<-0.2
    datts <- ts(dat, frequency=sps)
    rps <- spectrum(datts, spans=rep(span,2), pad=pad, taper=tap, plot=F)
    dfspec <- rbind(dfspec, data.frame(f=rps$freq, psd=rps$spec, 
                                       src="spec.pgram", niter=0, sta=stacha[1], cha=stacha[2]))
  }
  return(dfspec)
}
## do something with the data
sc <- 1e-9
nt <- 2e4
attach(bsm)
# varian
# varian_dat
dat <- sc*varian_dat$CH0[1:nt]
stacha <- c("Varian","CH0")
psdf <- rbind(psdf, doBSMspec(dat, stacha=stacha))
psdf <- rbind(psdf, doBSMspec(dat, niter=1, doRspec=FALSE, stacha=stacha))
psdf <- rbind(psdf, doBSMspec(dat, niter=2, doRspec=FALSE, stacha=stacha))
# and
# pinyon
# pinyon_dat
dat <- sc*pinyon_dat$CH0[1:nt]
stacha <- c("Pinyon","CH0")
psdf <- rbind(psdf, doBSMspec(dat, stacha=stacha))
psdf <- rbind(psdf, doBSMspec(dat, niter=1, doRspec=FALSE, stacha=stacha))
psdf <- rbind(psdf, doBSMspec(dat, niter=2, doRspec=FALSE, stacha=stacha))
#
detach(bsm)

## Now start plotting                                                      
pal <- brewer.pal(n=8,"Paired")
g <- ggplot(psdf, aes(x=log10(f), y=10*log10(psd), colour=factor(sta)))
(p <- g+geom_path(alpha=0.8)+facet_grid(factor(niter)~src)+
#   scale_x_log10(limits=c(0.01,10),expand=c(0,0))+
#   coord_trans(x = "log10", y="identity")+
  scale_y_continuous(limits=c(-230, -160))+
  scale_size_manual(values=c(0.5, 0.5))+
  scale_color_discrete()
 )
# ggsave("./co2.pdf", height=3.2, width=7)
##
##

