##
##  Creates figure(s) in text for CO2 data
##
setwd("~/kook.processing/R/dev/packages/rlpSpec/figures")
source('funcload.R')
load("../data/bsm/bsm.rda")

library(zoo)
library(RColorBrewer)
library(ggplot2)

rm(bsmnm, psdf, psdf.a, psdf.b)
nit <- 5

## Noise models from Barbour and Agnew 2011, BSSA
bsmnm <- as.data.frame(na.approx(read.table('../data/bsm/barbour_bsmnoise_2011.txt', h=T)))
lsmnm <- as.data.frame(na.approx(read.table('../data/bsm/barbour_lsmnoise_2011.txt', h=T)))

# begin assembling model dataframe
# LSM: last frame only
spec <- c(lsmnm$P10, rev(lsmnm$P50))
freq <- c(lsmnm$freq, rev(lsmnm$freq))
psdf.a <- data.frame(f=freq, psd=spec, ntap="", src="rlpSpec", niter=factor(5), sta="B&A 11: LSM", cha="LSM", polyf=factor(1))
# BSM: all frames
spec <- c(bsmnm$P10, rev(bsmnm$P50))
freq <- c(bsmnm$freq, rev(bsmnm$freq))
for (n in 0:nit){
  # -5 is a correction, as justified in B&A 2011
  psdf.a <- rbind(psdf.a, data.frame(f=freq, psd=spec-5, ntap="", src="rlpSpec", niter=factor(n), sta="B&A 11: BSM", cha="BSM", polyf=factor(2)))
}

# psdf.a <- rbind(psdf.a, data.frame(f=bsmnm$freq, psd=spec, src="spec.pgram", niter=0, sta="PBO", cha=""))

# wrapper for spectum estimators
doBSMspec <- function(dat, sps=20, stacha="", span=NULL, niter=0, doRspec=TRUE){
  ## return a data frame with spectral estimates
  # multitaper spectra
  rlpps <- pspectrum(dat, niter=niter, fsamp=sps, ndec=1, plotpsd=F)
  dfspec <- data.frame(f=rlpps$f, psd=rlpps$psd, ntap=rlpps$ntaper,
                       src="rlpSpec", niter=niter, sta=stacha[1], cha=stacha[2], polyf=factor(0))
  # r built-in (optional)
  if (doRspec){
    cat("\t>>>> Doing spec.pgram\n")
    pad<-1
    tap<-0.2
    datts <- ts(dat, frequency=sps)
    rps <- spectrum(datts, spans=rep(span,2), pad=pad, taper=tap, plot=F)
    dfspec <- rbind(dfspec, data.frame(f=rps$freq, psd=rps$spec, ntap=1,
                                       src="spec.pgram", niter=0, sta=stacha[1], cha=stacha[2], polyf=factor(0)))
  }
  return(dfspec)
}

linInvert <- function(linbsm, gap=0.0001, constant=0){
  # inverts linearized gauge strain to
  # raw counts
  D <- .087
  d <- gap
  c <- constant
  e <- linbsm
  dR <- 1e8
  toret <- dR*D*(e+c)/(d + D*(e+c))
  return(toret)
}

dbGain <- function(gap=0.0001, db=TRUE){
  ## returns least count gain for gap
  # see barbour agnew 2011 appendix a
  D <- .087
  d <- gap
  dR <- 1e-8
  toret <- dR*(d/D)/(1-(0.5+dR))**2
  if (db){
    toret <- 20*log10(toret)
  }
  return(toret)
}

## do something with them raw BSM data
gap <- 0.0001
sc <- 1e-9
nt <- 1e5
attach(bsm)
# varian
# varian_dat
dat.b <- linInvert(sc*varian_dat$CH0[1:nt], gap=gap)
stacha.b <- c("Varian","CH0")
# and
# pinyon
# pinyon_dat
dat.c <- linInvert(sc*pinyon_dat$CH0[1:nt], gap=gap)
stacha.c <- c("Pinyon","CH0")
#
detach(bsm)

## assemble spectrum dataframes
psdf.b <- doBSMspec(dat.b, doRspec=FALSE, stacha=stacha.b)
for (n in 1:nit){
  psdf.b <- rbind(psdf.b, doBSMspec(dat.b, niter=n, doRspec=FALSE, stacha=stacha.b))
}
# for Pinyon, just do last for plotting comparison (clutter reduction)
psdf.c <- doBSMspec(dat.c, niter=5, doRspec=FALSE, stacha=stacha.c)
# for (n in 1:nit){
#   psdf.c <- rbind(psdf.c, doBSMspec(dat.c, niter=n, doRspec=FALSE, stacha=stacha.c))
# }

## bind all the data frames
gain <- dbGain(gap=gap)
psdf.a$psdc <- psdf.a$psd
psdf.b$psdc <- 10*log10(psdf.b$psd) + gain
psdf.c$psdc <- 10*log10(psdf.c$psd) + gain
psdf <- rbind(psdf.c,
              psdf.b,
              psdf.a)
## Now start plotting                                
# pal <- brewer.pal(n=8,"Paired")
psdf$niter<-factor(psdf$niter, levels=c(0, 3, 1, 4, 2, 5) )
g <- ggplot(psdf, aes(x=log10(f), y=psdc, colour=factor(sta)))
(p <- g+geom_path()+facet_wrap(~niter,, as.table=T, ncol=2)+
  #aes(fill=polyf)
#   scale_x_log10(limits=c(0.01,10),expand=c(0,0))+
  scale_x_continuous("Frequency, log Hz", limits=c(-2.5, 1), breaks=-2:1, labels=-2:1, expand=c(0,0))+
  scale_y_continuous("Power, dB (strain^2 * sec)", limits=c(-230, -140), expand=c(0,0))+
#   scale_size_manual(values=c(0.5, 0.5))+
#   scale_color_discrete("Station")+
#   scale_fill_manual(breaks=c(0,1,2), values=c(NA,"gray","gray"), labels=c("","LSM","BSM"), legend=F)+
  scale_color_brewer("PSD: Station or Model", palette="Paired")+
  opts(title="BSM Spectra by Iteration")
 )
## save it
ggsave("./bsm.png")#, height=4.2, width=7)
##
psdf <- psdf.b[psdf.b$f>0,c(1,3,5)]
g <- ggplot(psdf, aes(x=log10(f), y=1/sqrt(ntap), colour=factor(niter)))
p <- g + geom_line(size=0.2)+
  scale_x_continuous("Frequency, log Hz", limits=c(-2.5, 1), breaks=-2:1, labels=-2:1, expand=c(0,0))+
  opts(title="BSM Spectra, Variance Reduction by Iteration")

g <- ggplot(psdf, aes(x=niter, y=1/sqrt(ntap), colour=log10(f), group=niter))
p <- g + geom_jitter(size=0.2)+
  scale_x_discrete("Iteration")+
  opts(title="BSM Spectra, Variance Reduction by Frequency")
##

