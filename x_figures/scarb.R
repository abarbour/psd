##
##  Creates figure(s) in text for Scarborough data
##
setwd("~/kook.processing/R/dev/packages/rlpSpec/figures")
# setwd("~/nute.processing/development/rlpSpec/figures")
#source('funcload.R')
source('../rsrc/.sourceloads.R')
load("../data/scarb/scarb.rda")

library(psych)
library(ggplot2)
library(plyr)
library(reshape)
library(RColorBrewer)

scarb <- scarbl$ts
nt <- length(scarb)
sps <- scarbl$sps
cha <- scarbl$component

# wrapper for spectum estimation
doSCARBspec <- function(sdat, nit=0, tst=NULL, sps=1, cha=NULL, niter=5, tapcap=1e3){
  ## return a couple data frame with spectral estimates
  secst <- (tst-1)/sps
  sxnlen <- length(sdat)/sps
  rlpps <- pspectrum(sdat, niter=niter, fsamp=sps, plotpsd=F, tapcap=tapcap)
  dfspec <- data.frame(f=rlpps$freq, psd=rlpps$spec, ntap=rlpps$ntaper, 
                       secst=secst, sxnlen=sxnlen, nit=factor(nit),
                       niter=niter, tapcap=tapcap)
  return(invisible(list(psd=rlpps, df=dfspec)))
}

nsxn <- 8
slen <- nt/nsxn
# Loop around sections of the full data vector
for (nit in 0:(nsxn-1)){
  tst <- nit*slen+1
  ten <- tst+slen-1
  dat <- scarb[tst:ten]
  scarbspec <- doSCARBspec(dat, nit=nit, tst=tst, sps=sps, cha=cha)
  if (nit==0){ 
    dfscarb1 <- scarbspec$df 
  } else {
    dfscarb1 <- rbind(dfscarb1, scarbspec$df)
  }
}
# Reshape (melt -> cast) to isolate by iteration
dfm <- reshape::melt.data.frame(dfscarb1, id.vars=c("nit","f"), measure.vars=c("psd"))
dfc <- reshape::cast(dfm, f~nit)
# f it1  it2  ...
# 0 p0/1 p0/2 ...
# ...
# cast using quantiles [ it_N replaced with CDF contours ]
dfc.quants <- reshape::cast(dfm, f~., quantile, probs = seq(0, 1, 0.10))

##
## Setup for plotting
##
frqs<-log10(dfc.quants$f)
f.p <- c(frqs,rev(frqs))
f.pmed <- 10^(frqs) #,rev(frqs))
q.frqs <- f.pmed
q.frqsp <- c(q.frqs,rev(q.frqs))
q.pmed <- log10(c(dfc.quants$X50.)) #, rev(dfc.quants$X0.)))
q.relmed <- log10(c(dfc.quants$X10., rev(dfc.quants$X90.))) - c(q.pmed, rev(q.pmed))

##
##  Figure 1: Show spectra and distribution contours
##
pal <- RColorBrewer::brewer.pal(n=10,"Paired")
pdf("./scarb.pdf", height=5)
plot(10^c(-2,2,2,-2),c(-2,-2,11,11), 
     ylim=c(0,10), xlim=10^log10(c(.1,31.25)), type="l", yaxs="i", xaxs="i",
     main="Scarborough Ez spectrum",
     xlab="Frequency, log10 Hz", ylab="PSD, log10 counts^2 / Hz",
     log="x")
lines(c(0.1,3.8),0.3*c(1,1),lty=2, lwd=2)
lines(c(4,25),0.3*c(1,1),lty=3, lwd=2)
text(c(0.3,10),0.5*c(1,1),c("spectrogram (A)","spectrogram (B)"),cex=0.7,adj=c(0.5,0))
text(5, 8.5, sprintf("PSD estimation:\n5 adaptive iterations\n%i x %.1f minute sections\n%s Hz sampling",nsxn,slen/sps/60,sps), cex=0.9, adj=c(0,0.5))
text(10^-0.5, 7, "Median PSD levels")
text(10^-0.9, 3.4, sprintf("80%s variation relative to median","%"),col=pal[2], adj=c(0,0.5))
mtext(sprintf("%s through %s", scarbl$tsst, scarbl$tsst+nt/sps))
# CDF contours relative to median
polygon(10^f.p, (q.relmed)+2, col=pal[1], border=pal[2], lwd=0.4)
lines(f.pmed, log10(dfc.quants$X50.), lwd=0.8)
dev.off()

##
## Figure 2A: Spectrogram for low frequencies
##

# now be selective
lims <- c(0,4)
fInd <- dfc$f>lims[1] & dfc$f<=lims[2]
frqs <- dfc$f[fInd]
dfc.med <- data.frame(f=dfc.quants$f, med=dfc.quants$X50.)
meds <- dfc.med$med[fInd]
iter <- 1:nsxn
dfmat <- as.matrix(dfc[fInd,(iter+1)])
#
pdf("./scarbSpecgram.pdf", height=5)
pal <- RColorBrewer::brewer.pal(n=8,"RdGy")
pal <- c(pal[5:8],pal[4:1])
filled.contour(y=(frqs), x=(iter-1)*slen/sps/60, log10(t(dfmat)), nlevels=5, 
               col=pal,
               ylim=lims,
               main="Scarborough Ez spectrogram (A)",
               ylab="Frequency, Hz", xlab=sprintf("Time, minutes from %s",scarbl$tsst), 
               plot.axes={ axis(1); axis(2); 
                           text(22.8,3.84,"PSD\nlog10 counts^2 / Hz", adj=c(1,1))
                           mtext(sprintf("%.1f minute sections",slen/sps/60))})
#
dev.off()

##
## Figure 2B: Spectrogram for high frequencies
##

nsxn <- 25
slen <- nt/nsxn
for (nit in 0:(nsxn-1)){
  tst <- nit*slen+1
  ten <- tst+slen-1
  dat <- scarb[tst:ten]
  scarbspec <- doSCARBspec(dat, nit=nit, tst=tst, sps=sps, cha=cha)
  if (nit==0){ 
    dfscarb2 <- scarbspec$df 
  } else {
    dfscarb2 <- rbind(dfscarb2, scarbspec$df)
  }
}
## Again, reshape and quantile
dfm <- reshape::melt.data.frame(dfscarb2, id.vars=c("nit","f"), measure.vars=c("psd"))
dfc <- reshape::cast(dfm, f~nit)
dfc.quants <- reshape::cast(dfm, f~., quantile, probs = seq(0, 1, 0.10))
#
lims <- c(4,25)
fInd <- dfc$f>lims[1] & dfc$f<=lims[2]
frqs <- dfc$f[fInd]
dfc.med <- data.frame(f=dfc.quants$f, med=dfc.quants$X50.)
meds <- dfc.med$med[fInd]
iter <- 1:nsxn
dfmat <- as.matrix(dfc[fInd,(iter+1)])
#
pdf("./scarbSpecgram2.pdf", height=5)
pal <- brewer.pal(n=9,"RdGy")
pal <- c(pal[5:9],pal[4:1])
filled.contour(y=(frqs), x=(iter-1)*slen/sps/60, log10(t(dfmat)), nlevels=8, 
               col=pal,
               ylim=lims,
               main="Scarborough Ez spectrogram (B)",
               ylab="Frequency, Hz", xlab=sprintf("Time, minutes from %s",scarbl$tsst), 
               plot.axes={ axis(1); axis(2); 
                           text(25,24,"PSD\nlog10 counts^2 / Hz", adj=c(1,1))
                           mtext(sprintf("%.1f minute sections",slen/sps/60))})
dev.off()

##
## Figure 3: Significant covariant correlations
##

lims <- c(0,25)
fInd<-dfc$f>lims[1] & dfc$f<=lims[2]
frqs <- dfc$f[fInd]
# length(frqs)
dfc.med <- data.frame(f=dfc.quants$f, med=dfc.quants$X50.)
meds <- dfc.med$med[fInd]
iter <- 1:nsxn
dfmat <- as.matrix(dfc[fInd,(iter+1)])

f.pmed <- dfc.quants$f
q.frqs <- f.pmed
q.frqsp <- c(q.frqs,rev(q.frqs))
q.pmed <- log10(c(dfc.quants$X50.)) #, rev(dfc.quants$X0.)))
q.relmed <- log10(c(dfc.quants$X10., rev(dfc.quants$X90.))) - c(q.pmed, rev(q.pmed))

## Pearson correlation (default) of frequencies, for iter combinations
dfmat.cor <- cor(t(dfmat)) # change method here if desired
## at N samples, correct correl. for significance (t-dist) and Bonferroni corr.
ctp <- psych::corr.p(dfmat.cor, nsxn, adjust="holm")
# output: p (Probabilities based on T), t (Students t-values), r (corrected correlation)
# decide what's significant: students t, 95% signif for 25-2=23 dof
dof <- nsxn-2
prc <- .99
sig <- qt(prc, dof)
sigt <- ctp$t > sig  
sigr <- ctp$r
# mask out insignificant correlations
sigcor <- ifelse(sigt, sigr, NA)
# Plot it
#
col <- RColorBrewer::brewer.pal("Spectral",n=11)
ylims <- lims; ylims[1] <- -10
labs <- seq(0, lims[2], by=2)
ylabs <- labs
ylabs[mod(labs,4)>0] <- ""
#
asp <- 1.1
pxht <- 1280
#
png("./scarbCorrel.png", res=200, height=asp*pxht, width=pxht)
#
filled.contour(x=(frqs), y=(frqs), z=(sigcor), col=rev(col), 
               xlim=lims, ylim=ylims, 
#                asp=1,
               zlim=c(0.45,1), nlevels=11, 
               main="Scarborough Ez spectrum correlations",
               xlab="Frequency, Hz", ylab="Frequency, Hz",
               plot.axes={ axis(1, at = labs, labels=ylabs)
                           axis(2, at = labs, labels=ylabs)
                           mtext(sprintf("Holm-Bonferroni %.0f%s Sig. at %i DoF",100*prc,'%',dof))
                           abline(h=0)
                           text(c(4.1,16),c(-3,-8),c("PSD","Variation"))
                           polygon(q.frqsp, 2*(q.relmed)-6.5, col="gray", border="black", lwd=0.4)
                           lines(q.frqs, log10(dfc.quants$X50.)-9, lwd=1)
                           })
#
dev.off()
##
##
##
