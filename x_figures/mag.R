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
#source('funcload.R')
source('../rsrc/.sourceloads.R')
load("../data/mag/mag.rda")

library(psych)
library(ggplot2)
library(plyr)
library(reshape)
library(RColorBrewer)
library(Hmisc)
library(xtable)

###
#
mag$km <- 1:length(mag$raw)
pdf("./magsat.pdf",height=3.5)
plot(raw ~ km, mag, type="s", 
     main="MAGSAT Airborne-magnetometer data",
     xlab="Along-path distance, km", yaxs="i", ylim=280*c(-1,1),
     ylab="nanoTesla")
mtext("Horizontal-component field")
text(200, 100, "raw series")
lines(clean-50 ~ km, mag, type="s", col="red")
text(350, -150, "cleaned series (-50)", col="red")
lines(edit-100 ~ km, mag)
dev.off()
#
psmag1 <- pspectrum(mag$clean, niter=5, plot=F)
psmag1$src <- "clean"
psmag2 <- pspectrum(mag$raw, niter=5, plot=F)
psmag2$src <- "raw"
nfrqfull <- length(psmag2$src)
psmagFULL <- rbind(psmag1[1:(nfrqfull-1),],psmag2[2:nfrqfull,])

xax <- seq(0,0.5,0.1)
plot.psd(psmagFULL, logx=F, ylabel="PSD, nT**2 km")
xaxl <- c("0",xax[2:6])
last_plot()+
  facet_grid(.~src)+
  scale_fill_manual("",breaks=rev(c("sig","tap")),
                    values=c("dark gray","light grey"),
                    labels=rev(c("Uncertainty: +- Sigma", "Optim. tapers: Rel. pilot spec.")))+
  scale_x_continuous("Wavenumber, 1/km", limits=range(xax), breaks=xax, labels=xaxl, expand=c(0,0))+
  opts(title=sprintf("MAGSAT PSD Estimation after %i iterations",5),
       legend.position=c(0.3, 0.90),
       panel.margin = unit(1, "lines"))
ggsave("./magsatPSD.pdf",height=4.5)

###
###  Spectum fitting
###

##
##  Maximum Likelihood Estimate of the spectrum parameters for the MAGSAT data
##
skip.c <-c(10,10) # the number of spectrum estimates to skip (head and tail)
skip.r <- c(10,500)
# Do for single taper periodograms (for statistical independence)
psmag.cln <- pspectrum(mag$clean, niter=0, ntapinit=1, plot=F)
psmag.raw <- pspectrum(mag$raw, niter=0, ntapinit=1, plot=F)
# cull
psmag.c <- psmag.cln[(skip.c[1]+1):(length(psmag.cln$psd)-skip.c[2]),]
psmag.r <- psmag.raw[(skip.r[1]+1):(length(psmag.raw$psd)-skip.r[2]),]

# setup likelihood function (2D: A and B to optimize)
mag.lik <- function(AB, magpsd, uncert=FALSE){
  A <- AB[1]
  B <- AB[2]
  S <- magpsd$psd
  k <- magpsd$f
  if (!uncert){
    # return the likelihood for A,B
    lik <- log(A) - k*B + S*exp(k*B)/A
    L <- sum(lik)
    # but dont return -L!
    return(L)
  } else {
    # return the standard deviation
    # approx (sigmaA**2 ~ 1/sum(d2L/dA2) and similarly for sigmaB**2)
    ekB <- exp(k*B)
    var.A <- -1/sum( 1/(A**2) - (2*S*ekB)/(A**3) )
    var.B <- 1/sum( (k**2)*S*ekB/A )
    return(list(sigmaA=sqrt(var.A), sigmaB=sqrt(var.B)))
  }
}
mag.bias <- function(S){
  return(log(S*exp(-0.577721155))) # eulers constant
}

# do the suite optimization methods
doMagOpt <- function(magpsd, LikFunc, src,
                     ABinit=c(700, 14), 
                     optmeth=c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN")){
  # optmeth: gradient/finite-diff methods built into 'optim': go through them all
  # ("Brent" only for 1D opt)
  meths <- 1:length(optmeth)
  ABopt <- data.frame(spec=NA, Ao=NA, A=NA, sigA=NA, Bo=NA, B=NA, sigB=NA, method=NA, convergence=NA, lik_calls=NA, grad_calls=NA, log2_norm_sys_time=NA)
  for (nm in meths){
    meth <- optmeth[nm]
    sysT <- system.time(mag.opt <<- optim(ABinit, LikFunc, magpsd=magpsd, method=meth, hessian=TRUE))[1]
    Aopt <- mag.opt$par[1]
    Bopt <- mag.opt$par[2]
    Sigmas <- LikFunc(c(Aopt,Bopt), magpsd, uncert=TRUE)
    # now get uncertainty
    ABopt <- rbind(ABopt, c(src, ABinit[1], Aopt, Sigmas$sigmaA, ABinit[2], Bopt, Sigmas$sigmaB, meth, mag.opt$convergence, mag.opt$counts[1], mag.opt$counts[2], log2(sysT/mag.opt$counts[1])))
  }
  return(ABopt[meths+1,])
}

# A and B optimization
MagSpecOpt <- as.data.frame(rbind(doMagOpt(psmag.c, mag.lik, "clean"), doMagOpt(psmag.r, mag.lik, "raw")))
rownames(MagSpecOpt) <- NULL
MLEconv <- MagSpecOpt[MagSpecOpt$spec == "clean" & MagSpecOpt$convergence==0 & as.numeric(MagSpecOpt$lik_calls) < 10e3, 
                      # method, lik_calls, grad_calls, A, sigA, B, sigB
                      c(8,10,11,3,4,6,7)]
rownames(MLEconv) <- NULL
# culling
MLE <- MLEconv[c(4:7)]
doMed <- function(vals){return(mean(as.numeric(vals)))}
# calc residuals (lapply??)
MLE.res <- data.frame(A=doMed(MLE$A), sigA=doMed(MLE$sigA), B=doMed(MLE$B), sigB=doMed(MLE$sigB))
# temp df for latex output
meddf <- as.data.frame(c(c("","","means"), MLE.res))
names(meddf) <- attr(MLEconv,'names')
MLEtolatex <- rbind(MLEconv, meddf)
# out to latex
xtable::print.xtable(xtable::xtable(MLEtolatex), file="../tex/MLE.tex")
# add some more info for plotting
MLE.res$AU <- MLE.res$A + MLE.res$sigA
MLE.res$AL <- MLE.res$A - MLE.res$sigA
MLE.res$BU <- MLE.res$B + MLE.res$sigB
MLE.res$BL <- MLE.res$B - MLE.res$sigB
##
## Comparing good MLE results to linearized model
##

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
    est <- "spec.pgram"
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
tap <- c(seq(1,14,1),seq(17,99,15),seq(100,1000,100))
taps <- 1:length(tap)
for (nt in taps){
  tmpdf <- rbind(tmpdf, doLinMag(mag$clean, "clean", skip.c, tap[nt]))
  tmpdf <- rbind(tmpdf, doLinMag(mag$raw, "raw", skip.r, tap[nt]))
}

lmOpt <- tmpdf
rownames(lmOpt) <- NULL
lmOpt$est <- factor(lmOpt$est, levels=rev(c("rlpSpec","spec.pgram")))
lims <- c(.63,1250)

# LM parameter A
g <- ggplot(data=lmOpt, aes(x=log10(tap), y=(A), group=est))
txt <- sprintf("Gray band is MLE:\n%.1f +/- %.1f", MLE.res$A, MLE.res$sigA)
g + 
  geom_ribbon(fill="light gray", colour="gray", size=0.4, aes(x=log10(lims), ymin=c(MLE.res$AL), ymax=c(MLE.res$AU)))+
  geom_hline(linetype=2, size=0.3, aes(yintercept=c(MLE.res$A) ))+
  geom_path(colour="dark gray", size=0.5) + 
  geom_pointrange(colour='black', aes(ymin=A*(1-(seA-1)), ymax=A*(1+(seA-1)), fill=est, shape=est, size=est)) +
  geom_text(data=data.frame(x=2,y=250,label=txt,est=NA,spec="clean"), colour="dark gray", size=3.5, aes(x=x, y=y, label=label))+
  scale_shape_manual("PSD estimator", values=c(22,21)) +
  scale_x_continuous("Tapers, log10", limits=log10(lims), expand=c(0,0), breaks=0:3,labels=10^(0:3))+
  scale_fill_manual(guide="none", values=c("gray", "light gray"))+
  scale_colour_manual(guide="none", values=c("black","black"))+
  scale_size_manual(guide="none", values=c(0.9,0.6))+
  scale_linetype(guide="none")+
  ylab("Scaling parameter, nT^2 km") + 
  opts(title="MAGSAT Spectrum: Linearized model vs MLE (A)") + theme_bw() + facet_grid(.~spec)
ggsave("mag_lmA.pdf",height=3.5)

# LM parameter B
g <- ggplot(data=lmOpt, aes(x=log10(tap), y=(-B), group=est, colour=est))
txt <- sprintf("Gray band is MLE:\n%.1f +/- %.1f", MLE.res$B, MLE.res$sigB)
g +
  geom_ribbon(fill="light gray", colour="gray", size=0.4, aes(x=log10(lims), ymin=c(MLE.res$BL), ymax=c(MLE.res$BU)))+
  geom_hline(linetype=2, size=0.3, aes(yintercept=c(MLE.res$B)))+
  geom_path(colour="dark gray", size=0.4) + 
  geom_pointrange(aes(ymin=-B*(1-seB), ymax=-B*(1+seB), fill=est, shape=est, colour=est, size=est)) +
  geom_text(data=data.frame(x=2,y=0,label=txt,est=NA,spec="clean"), colour="dark gray", size=3.5, aes(x=x, y=y, label=label))+
  scale_shape_manual("PSD estimator", values=c(22,21)) +
  scale_x_continuous("Tapers, log10", limits=log10(lims), expand=c(0,0), breaks=0:3,labels=10^(0:3))+
  scale_fill_manual(guide="none", values=c("gray","light gray"))+
  scale_colour_manual(guide="none", values=c("black","black"))+
  scale_size_manual(guide="none", values=c(0.9,0.6))+
  scale_linetype(guide="none")+
  ylab("Decay parameter, km") + 
  opts(title="MAGSAT Spectrum: Linearized model vs MLE (B)") + theme_bw() + facet_grid(.~spec)
ggsave("mag_lmB.pdf",height=3.5)
##
