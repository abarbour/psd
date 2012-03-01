##
##  Creates figure(s) in text for CO2 data
##
setwd("~/kook.processing/R/dev/packages/rlpSpec/figures")
source('funcload.R')
load("../data/bsm/bsm.rda")

library(zoo)
library(RColorBrewer)
library(ggplot2)

bsmnm <- as.data.frame(na.approx(read.table('../data/bsm/barbour_bsmnoise_2011.txt', h=T)))
psdf <- data.frame(f=bsmnm$freq, psd=10^bsmnm$P50, src="", niter=0, sta="PBO", cha="")
psdf <- rbind(psdf, data.frame(f=bsmnm$freq, psd=bsmnm$P50, src="", niter=2, sta="PBO", cha="")
psdf <- rbind(psdf, data.frame(f=bsmnm$freq, psd=bsmnm$P50, src="", niter=2, sta="PBO", cha="")

doBSMspec <- function(dat, sps=20, stacha="", span=NULL, niter=0, doRspec=TRUE){
  ## return a data frame with spectral estimates
  # multitaper spectra
  rlpps <- pspectrum(dat, niter=niter, fsamp=sps, ndec=sps, plotpsd=F)
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
attach(bsm)
# varian
# varian_dat
dat <- 1e-9*varian_dat$CH0[1:5e3]
stacha <- c("Varian","CH0")
psdf <- rbind(psdf, doBSMspec(dat, stacha=stacha))
psdf <- rbind(psdf, doBSMspec(dat, niter=1, doRspec=FALSE, stacha=stacha))
psdf <- rbind(psdf, doBSMspec(dat, niter=2, doRspec=FALSE, stacha=stacha))
# and
# pinyon
# pinyon_dat
dat <- 1e-9*pinyon_dat$CH0[1:5e3]
stacha <- c("Pinyon","CH0")
psdf <- rbind(psdf, doBSMspec(dat, stacha=stacha))
psdf <- rbind(psdf, doBSMspec(dat, niter=1, doRspec=FALSE, stacha=stacha))
psdf <- rbind(psdf, doBSMspec(dat, niter=2, doRspec=FALSE, stacha=stacha))
#
detach(bsm)

## Now start plotting                                                      
pal <- brewer.pal(n=8,"Paired")
g <- ggplot(psdf, aes(x=f, y=10*log10(psd), colour=factor(sta)))
(p <- g+geom_path(alpha=0.8)+facet_grid(niter~src)+
  scale_x_log10(limits=c(0.01,10),expand=c(0,0))+
  scale_y_continuous(limits=c(-240, -160))+
  scale_size_manual(values=c(0.5, 0.5), legend=F)+
  scale_color_discrete()
 )
# ggsave("./co2.pdf", height=3.2, width=7)
##
##

