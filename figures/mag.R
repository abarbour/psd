##
##  Creates figure(s) in text for MAGSAT data
##
##
# We expect the true spectrum of the magnetic field to behave as
# 
#	S(k) = A exp(-kB)					(1)
#
# Each spectral estimate is distributed as X^2 (chi-squared); hence, the 
# pdf (phi) of the single taper estimate S' at wavenumber k will be
#
#	phi(k) = exp(-k/S)/S, k>0				(2)
#
# where S is the true PSD. The variance of the estimation is
#
#	var[S']=S^2						(3)
#
# Since we have an exact expression for phi and S we can use the Maximum Likelihood
# Estimate for A and B, where the likelihood function is the product of pdfs with the
# observed spectrum in place of k, or:
#
#	L(A,B)	= ln[ phi(S1') * phi(S2') * ... * phi(Sn') ]
#		= sum_n{ ln[phi(Sn'(k))] }				
#		= sum_n{ -ln[S(k)] - Sn'(k)/S(k) }			(5)
#
# where Sn'(k) = Sx'(k_n) is the n-th measurement, the spectral estimate at k_n, and S(k)
# is (1).  Thus, the log likelihood function to maximize :
#
#	L(A,B)  = - sum_n{ ln[A] - k_n*B + Sn'*exp(k_n*B)/A  }
#
##
setwd("~/kook.processing/R/dev/packages/rlpSpec/figures")
# setwd("~/nute.processing/development/rlpSpec/figures")
source('funcload.R')
load("../data/mag/mag.rda")

library(psych)
library(ggplot2)
library(plyr)
library(reshape)
library(RColorBrewer)

#
# 
skip.c <-c(1,1) # the number of spectrum estimates to skip (head and tail)
skip.r <- c(1,600)

doLinMag <- function(dat, src, psdskip, tapinit, rspec=F){
  ## 
  ##   Use lm to estimate spectrum params
  ##
  if (!rspec){
    est <- "rlpSpec"
    psd <- pspectrum(dat, niter=0, ntapinit=tapinit, plot=F)
    psd.c <- psd[(psdskip[1]+1):(length(psd$psd)-psdskip[2]),]
    psd.lm <- lm(log(psd)~f, psd.c)
  } else {
    tapinit <- 1
    est <- "spc.pgram"
    psd <- spectrum(dat, taper=0.2, pad=T, plot=F)
    psd.s <- psd$spec[(psdskip[1]+1):(length(psd$spec)-psdskip[2])]
    psd.f <- psd$freq[(psdskip[1]+1):(length(psd$freq)-psdskip[2])]
    psd.lm <- lm(log(psd.s)~psd.f)
  }
  #
  lmpar <- data.frame(spec=src, est=est,
                      A=exp(coefficients(psd.lm)[1]), 
                      seA=exp(summary(psd.lm)$coefficients[1,2]),
                      B=coefficients(psd.lm)[2], 
                      seB=summary(psd.lm)$coefficients[2,2],
                      tap=tapinit)
  return(lmpar)
}
# spec.pgram
tmpdf <- doLinMag(mag$raw, "raw", skip.r, 0, rspec=TRUE)
tmpdf <- rbind(tmpdf, doLinMag(mag$clean, "clean", skip.c, 0, rspec=TRUE))
# rownames(lmOptR) <- NULL
# rlpSpec
tap <- c(seq(1,14,1),seq(15,99,15),seq(100,1000,100))
taps <- 1:length(tap)
for (nt in taps){
  tmpdf <- rbind(tmpdf, doLinMag(mag$clean, "clean", skip.c, tap[nt]))
  tmpdf <- rbind(tmpdf, doLinMag(mag$raw, "raw", skip.r, tap[nt]))
}

lmOpt <- tmpdf
rownames(lmOpt) <- NULL

# LM parameter A
g <- ggplot(data=lmOpt, aes(x=log10(tap), y=(A), colour=est))
# g+ geom_hline()+
g + geom_path() + 
  geom_pointrange(size=0.5, aes(ymin=A*(1-(seA-1)), ymax=A*(1+(seA-1)), shape=spec, fill=est)) +
  scale_shape_manual(guide="none", values=c(21,22)) +
  scale_x_continuous("Tapers, log10", limits=c(0,3),breaks=0:3,labels=10^(0:3))+
#   scale_size_continuous("Standard\ndeviation")+ #, limits=c(1,1.5))+
  scale_color_discrete("Estimator")+
  scale_linetype(guide="none")+
  scale_fill_discrete(guide="none")+
  ylab("nT^2 km") + 
  opts(title="MAGSAT Spectrum: Linearized model (A)") + theme_bw() + facet_grid(.~spec)
ggsave("mag_lmA.pdf",height=3.5)

# LM parameter B
g <- ggplot(data=lmOpt, aes(x=log10(tap), y=(-B), colour=est, fill=est))
# g+ geom_hline()+
g + geom_path() + 
  geom_pointrange(size=0.5, aes(ymin=-B*(1-seB), ymax=-B*(1+seB), shape=spec, fill=est)) +
  scale_shape_manual(guide="none", values=c(21,22)) +
  scale_x_continuous("Tapers, log10", limits=c(0,3),breaks=0:3,labels=10^(0:3))+
#   scale_size_continuous("Standard\ndeviation", limits=c(0.5,1.5))+
  scale_color_discrete("Estimator")+
  scale_linetype(guide="none")+
  scale_fill_discrete(guide="none")+
  ylab("Wavenumber, km") + 
  opts(title="MAGSAT Spectrum: Linearized model (B)") + theme_bw() + facet_grid(.~spec)
ggsave("mag_lmB.pdf",height=3.5)

##
##  Maximum Likelihood Estimate of the spectrum parameters for the MAGSAT data
##

# REDO for single taper periodograms (for statistical independence)
psmag.cln <- pspectrum(mag$clean, niter=0, ntapinit=1, plot=F)
psmag.raw <- pspectrum(mag$raw, niter=0, ntapinit=1, plot=F)
# cull
psmag.c <- psmag.cln[(skip.c[1]+1):(length(psmag.cln$psd)-skip.c[2]),]
psmag.r <- psmag.raw[(skip.r[1]+1):(length(psmag.raw$psd)-skip.r[2]),]

# setup likelihood function (2D: A and B to optimize)
mag.lik <- function(AB, magpsd){
  A <- AB[1]
  B <- AB[2]
  S <- magpsd$psd
  k <- magpsd$f
  lik <- log(A) - k*B + S*exp(k*B)/A
  L <- sum(lik)
  # dont return -L!
  return(L)
}

# do the suite optimization methods
doMagOpt <- function(magpsd, lik, src,
                     ABinit=c(500, 15), 
                     optmeth=c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN")){
  # optmeth: gradient/finite-diff methods built into 'optim': go through them all
  # ("Brent" only for 1D opt)
  meths <- 1:length(optmeth)
  ABopt <- data.frame(spec=NA, Ao=NA, A=NA, Bo=NA, B=NA, method=NA, convergence=NA, lik_calls=NA, grad_calls=NA, log2_norm_sys_time=NA)
  for (nm in meths){
    meth <- optmeth[nm]
    sysT <- system.time(mag.opt <<- optim(ABinit, lik, magpsd=magpsd, method=meth, hessian=TRUE))[1]
    ABopt <- rbind(ABopt, c(src, ABinit[1], mag.opt$par[1], ABinit[2], mag.opt$par[2], meth, mag.opt$convergence, mag.opt$counts[1], mag.opt$counts[2], log2(sysT/mag.opt$counts[1])))
  }
  return(ABopt[meths+1,])
}

# A and B optimization
MagSpecOpt <- as.data.frame(rbind(doMagOpt(psmag.c, mag.lik, "clean"), doMagOpt(psmag.r, mag.lik, "raw")))
rownames(MagSpecOpt) <- NULL
MagSpecOpt
##