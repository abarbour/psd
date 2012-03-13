##
##  Creates figure(s) in text for Scarborough data
##
setwd("~/kook.processing/R/dev/packages/rlpSpec/figures")
# setwd("~/nute.processing/development/rlpSpec/figures")
source('funcload.R')
load("../data/scarb/scarb.rda")

library(psych)
library(fields)
library(ggplot2)
library(plyr)
library(reshape)
library(RColorBrewer)

scarb <- scarbl$ts
nt <- length(scarb)

# wrapper for spectum estimators
doSCARBspec <- function(sdat, nit=0, tst=NULL, sps=scarbl$sps, cha=scarbl$component, niter=5, tapcap=1e3){
  ## return a couple data frame with spectral estimates
  secst <- (tst-1)/sps
  sxnlen <- length(sdat)/sps
  rlpps <- pspectrum(sdat, niter=niter, fsamp=sps, plotpsd=F, tapcap=tapcap)
  dfspec <- data.frame(f=rlpps$f, psd=rlpps$psd, ntap=rlpps$ntaper, 
                       secst=secst, sxnlen=sxnlen, nit=factor(nit),
                       niter=niter, tapcap=tapcap)
  return(invisible(list(psd=rlpps, df=dfspec)))
}

nsxn <- 8
slen <- nt/nsxn
for (nit in 0:(nsxn-1)){
  tst <- nit*slen+1
  ten <- tst+slen-1
  dat <- scarb[tst:ten]
  scarbspec <- doSCARBspec(dat, nit=nit, tst=tst)
  if (nit==0){ 
    dfscarb <- scarbspec$df 
  } else {
    dfscarb <- rbind(dfscarb, scarbspec$df)
  }
}
## now lets do some reshaping
dfm <- melt.data.frame(dfscarb, id.vars=c("nit","f"), measure.vars=c("psd"))
dfc <- cast(dfm, f~nit)
# f it1  it2  ...
# 0 p0/1 p0/2 ...
# ...
# cast using quantiles (itN replaced with CDF contours)
dfc.quants <- cast(dfm, f~., quantile, probs = seq(0, 1, 0.10))

##
#  Let's start plotting
##
# g <- ggplot(dfscarb,aes(x=(f), y=log10(psd)-2*as.numeric(nit), colour=(secst/60)))
# g+geom_path()+
#   #facet_grid(secst~.)+
#   #   scale_colour_gradient(guide="none")+
#   #   scale_color_brewer()+
#   scale_y_continuous(limits=c(-18,3), expand=c(0,0))+
#   #   scale_x_continuous(limits=c(-1,log10(33)), expand=c(0,0))+
#   scale_x_continuous(limits=(c(1.8,25)), expand=c(0,0))+
#   theme_bw()

frqs<-log10(dfc.quants$f)
f.p <- c(frqs,rev(frqs))
f.pmed <- 10^(frqs) #,rev(frqs))
q.pmed <- log10(c(dfc.quants$X50.)) #, rev(dfc.quants$X0.)))
q.relmed <- log10(c(dfc.quants$X10., rev(dfc.quants$X90.))) - 
  c(q.pmed, rev(q.pmed))
# log10(c(dfc.quants$X50., rev(dfc.quants$X50.)))
##
##  Show spectra and distribution contours
##
pal <- brewer.pal(n=10,"Paired")
# some dummy data
# Ac(0,4)Bc(4,25)
pdf("./scarb.pdf", height=5)
plot(10^c(-2,2,2,-2),c(-2,-2,11,11), 
     ylim=c(0,10), xlim=10^log10(c(.1,31.25)), type="l", yaxs="i", xaxs="i",
     main="Scarborough Ez spectrum",
     xlab="Frequency, log10 Hz", ylab="PSD, log10 counts^2 / Hz",
     log="x")
lines(c(0.1,3.8),0.3*c(1,1),lty=2, lwd=2)
lines(c(4,25),0.3*c(1,1),lty=3, lwd=2)
text(c(0.3,10),0.5*c(1,1),c("spectrogram (A)","spectrogram (B)"),cex=0.7,adj=c(0.5,0))
#c(0,4),c(4,25))
text(5, 8.5, sprintf("PSD estimation:\n5 adaptive iterations\n%i x %.1f minute sections\n%s Hz sampling",nsxn,slen/sps/60,sps), cex=0.9, adj=c(0,0.5))
text(10^-0.5, 7, "Median PSD levels")
text(10^-0.9, 3.3, "Variation relative to median",col=pal[2], adj=c(0,0.5))
mtext(sprintf("%s through %s", scarbl$tsst, scarbl$tsst+nt/sps))
# CDF contours
polygon(10^f.p, (q.relmed)+2, col=pal[1], border=pal[2], lwd=0.4)
lines(10^frqs,log10(dfc.quants$X50.), lwd=0.8)
dev.off()
##
## contours relative to median
##
# pal <- brewer.pal(n=10,"Paired")
# plot(2*c(-1,1,1,-1),11*c(0,0,1,1), ylim=c(-2,2), xlim=(c(5,15)), type="l")
# # CDF contours
# q.pr <- log10(c(dfc.quants$X0.)) - q.pmed
# lines(f.pmed, q.pr, col=pal[1], lwd=0.8)
# q.pr <- log10(c(dfc.quants$X10.)) - q.pmed
# lines(f.pmed, q.pr, col=pal[2], lwd=0.8)
# q.pr <- log10(c(dfc.quants$X20.)) - q.pmed
# lines(f.pmed, q.pr, col=pal[3], lwd=0.8)
# q.pr <- log10(c(dfc.quants$X30.)) - q.pmed
# lines(f.pmed, q.pr, col=pal[4], lwd=0.8)
# q.pr <- log10(c(dfc.quants$X40.)) - q.pmed
# lines(f.pmed, q.pr, col=pal[5], lwd=0.8)
# q.pr <- log10(c(dfc.quants$X50.)) - q.pmed
# lines(f.pmed, q.pr, col=pal[6], lwd=0.8)
# q.pr <- log10(c(dfc.quants$X60.)) - q.pmed
# lines(f.pmed, q.pr, col=pal[5], lwd=0.8)
# q.pr <- log10(c(dfc.quants$X70.)) - q.pmed
# lines(f.pmed, q.pr, col=pal[4], lwd=0.8)
# q.pr <- log10(c(dfc.quants$X80.)) - q.pmed
# lines(f.pmed, q.pr, col=pal[3], lwd=0.8)
# q.pr <- log10(c(dfc.quants$X90.)) - q.pmed
# lines(f.pmed, q.pr, col=pal[2], lwd=0.8)
# q.pr <- log10(c(dfc.quants$X100.)) - q.pmed
# lines(f.pmed, q.pr, col=pal[1], lwd=0.8)

# now be selective
fInd<-dfc$f>0 & dfc$f<=4
frqs <- dfc$f[fInd]
dfc.med <- data.frame(f=dfc.quants$f, med=dfc.quants$X50.)
meds <- dfc.med$med[fInd]
iter <- 1:nsxn
dfmat <- as.matrix(dfc[fInd,(iter+1)])
#
pdf("./scarbSpecgram.pdf", height=5)
#
pal <- brewer.pal(n=8,"RdGy")
pal <- c(pal[5:8],pal[4:1])
filled.contour(y=(frqs), x=(iter-1)*slen/sps/60, log10(t(dfmat)), nlevels=5, 
               col=pal,
               ylim=c(0,4),
               main="Scarborough Ez spectrogram (A)",
               ylab="Frequency, Hz", xlab=sprintf("Time, minutes from %s",scarbl$tsst), 
               plot.axes={ axis(1); axis(2); 
                           text(22.8,3.84,"PSD\nlog10 counts^2 / Hz", adj=c(1,1))
                           mtext(sprintf("%.1f minute sections",slen/sps/60))})
#
dev.off()

nsxn <- 25
slen <- nt/nsxn
for (nit in 0:(nsxn-1)){
  tst <- nit*slen+1
  ten <- tst+slen-1
  dat <- scarb[tst:ten]
  scarbspec <- doSCARBspec(dat, nit=nit, tst=tst)
  if (nit==0){ 
    dfscarb <- scarbspec$df 
  } else {
    dfscarb <- rbind(dfscarb, scarbspec$df)
  }
}
## now lets do some reshaping
dfm <- melt.data.frame(dfscarb, id.vars=c("nit","f"), measure.vars=c("psd"))
dfc <- cast(dfm, f~nit)
# f it1  it2  ...
# 0 p0/1 p0/2 ...
# ...
# cast using quantiles (itN replaced with CDF contours)
dfc.quants <- cast(dfm, f~., quantile, probs = seq(0, 1, 0.10))
fInd<-dfc$f>4 & dfc$f<=25
frqs <- dfc$f[fInd]
dfc.med <- data.frame(f=dfc.quants$f, med=dfc.quants$X50.)
meds <- dfc.med$med[fInd]
iter <- 1:nsxn
dfmat <- as.matrix(dfc[fInd,(iter+1)])
#
#
pdf("./scarbSpecgram2.pdf", height=5)
pal <- brewer.pal(n=9,"RdGy")
pal <- c(pal[5:9],pal[4:1])
filled.contour(y=(frqs), x=(iter-1)*slen/sps/60, log10(t(dfmat)), nlevels=8, 
               col=pal,
               ylim=c(4,25),
               main="Scarborough Ez spectrogram (B)",
               ylab="Frequency, Hz", xlab=sprintf("Time, minutes from %s",scarbl$tsst), 
               plot.axes={ axis(1); axis(2); 
                           text(25,24,"PSD\nlog10 counts^2 / Hz", adj=c(1,1))
                           mtext(sprintf("%.1f minute sections",slen/sps/60))})
dev.off()
#
## Find significant covariant correlations
#
lims <- c(0,25)
fInd<-dfc$f>lims[1] & dfc$f<=lims[2]
frqs <- dfc$f[fInd]
length(frqs)
dfc.med <- data.frame(f=dfc.quants$f, med=dfc.quants$X50.)
meds <- dfc.med$med[fInd]
iter <- 1:nsxn
dfmat <- as.matrix(dfc[fInd,(iter+1)])

## Pearson correlation of frequencies, for iter combinations
dfmat.cor <- cor(t(dfmat)) #, method="spearman")
## at N samples, correct correl. for significance (t-dist) and Bonferroni corr.
ctp <- psych::corr.p(dfmat.cor, nsxn, adjust="holm")
# decide what's significant: students t, 95% signif for 25-2=23 dof
dof <- nsxn-2
prc <- .99
sig <- rev(qt(prc, dof))[1]
sigt <- ctp$t > sig  
sigr <- ctp$r
# mask out insignificant correlations
sigcor <- ifelse(sigt, sigr, NA)
# Plot it
# pdf("./scarbSpecgram2Cor.pdf")
png("./scarbSpecgram2Cor.png", res=156, height=780, width=900)
col <- RColorBrewer::brewer.pal("Spectral",n=11)
filled.contour(x=(frqs), y=(frqs), z=(sigcor), col=rev(col), 
               zlim=c(0.45,1), nlevels=11, asp=1,
               main="Scarborough Ez spectrum correlations",
               xlab="Frequency, Hz", ylab="Frequency, Hz",
               plot.axes={ axis(1); axis(2);
                           mtext(sprintf("Holm-Bonferroni %.0f%s Sig. at %i DoF",100*prc,'%',dof))
                           })
dev.off()
# fields::image.plot(x=frqs, y=frqs, z=(sigcor), add=F, 
# #                    asp=1,
#                    col=rev(col),
#                    graphics.reset=T,
#                    main="Scarborough Ez spectrum correlations",
#                    xlab="Frequency, Hz",
#                    ylab="Frequency, Hz",
#                    ylim=lims, xlim=lims)#, zlim=c(0.3,1), )
##