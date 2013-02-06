
## @knitr eval=TRUE, eval=TRUE, label="Load library."
library(rlpSpec)


## @knitr eval=TRUE, eval=TRUE, label="Load MAGSAT data."
data(magsat)


## @knitr eval=TRUE, eval=TRUE, label="Show contents of MAGSAT."
names(magsat)


## @knitr eval=TRUE, echo=TRUE, label="Outliers."
subset(magsat, abs(mdiff)>0)


## @knitr eval=TRUE, echo=TRUE, label=MAGPSDS
psdr <- pspectrum(magsat$raw)
psdc <- pspectrum(magsat$clean)


## @knitr eval=TRUE, fig=TRUE, echo=TRUE, label=RAWvCLEAN
plot(psdc, log="dB", main="Raw and Clean MAGSAT power spectral density", 
       lwd=3, ci.col=NA, ylim=c(0,32), yaxs="i")
plot(psdr, log="dB", add=TRUE, lwd=3, lty=5)
text(c(0.25,0.34), c(11,24), c("Clean","Raw"), cex=1)


## @knitr eval=FALSE, echo=TRUE, label="Naive spectrum estimation."
spec.pgram(X, pad=1, taper=0.2, detrend=FALSE, demean=FALSE, plot=F)


## @knitr eval=TRUE, echo=TRUE, fig=TRUE,  label=MAGSATNAIVE
ntap <- psdc$taper
psdcore(magsat$clean, ntaper=ntap, refresh=TRUE, plotpsd=TRUE)


## @knitr eval=TRUE, echo=TRUE, label="Load RSEIS package."
require(RSEIS)
dt=1 # km
ats <- prewhiten(ts(magsat$clean, deltat=dt), plot=FALSE) 


## @knitr eval=TRUE, echo=TRUE, fig=TRUE, fig.height=3.5, fig.width=6.5, label=ARFIT
atsar <- prewhiten(ats, AR.max=100, verbose=FALSE)


## @knitr eval=TRUE, echo=TRUE, fig=TRUE, fig.height=3, fig.width=4.5, label=ARPSD
plot(psdcore(atsar,ntaper=10), log="dB", main="PSD of MAGSAT innovations")


## @knitr eval=TRUE, echo=TRUE, label="Sampling rate versus interval."
a <- rnorm(32)
all.equal(psdcore(a,1)$spec, psdcore(a,-1)$spec)


## @knitr eval=TRUE, echo=TRUE, label="Compute PSD with mtapspec."
tapinit <- 10
Mspec <- mtapspec(ats, deltat(ats), MTP=list(kind=2, inorm=3, nwin=tapinit, npi=0))


## @knitr eval=TRUE, echo=TRUE
str(Mspec)


## @knitr eval=TRUE, echo=TRUE, label="Comparative spectra."
Xspec <- spec.pgram(ats, pad=1, taper=0.2, detr=TRUE, dem=TRUE, plot=FALSE)
Pspec <- psdcore(ats, dt, tapinit)
Aspec <- pspectrum(ats, dt, tapinit, niter=3, plot=FALSE)
# Correct for double-sidedness of spectrum and mtapspec results
class(Mspec)
Mspec <- normalize(Mspec, dt, "spectrum")
nt <- 1:Mspec$numfreqs
mspec <- Mspec$spec[nt]
class(Xspec)
Xspec <- normalize(Xspec, dt, "spectrum")


## @knitr eval=TRUE, echo=TRUE, fig=TRUE, fig.width=6.0, fig.height=5.4, label=RSEIS
require(RColorBrewer)
cols <- c("dark grey", brewer.pal(8, "Set1")[c(5:4,2)])
lwds <- c(1,2,2,5)
plot(Xspec, log="dB", ylim=40*c(-0.4,1), ci.col=NA, 
       col=cols[1], lwd=lwds[1], main="PSD Comparisons") 
pltf <- Mspec$freq
lines(pltf, pltp <- dB(mspec), col=cols[2], lwd=lwds[2]) 
plot(Pspec, log="dB",  add=TRUE, col=cols[3], lwd=lwds[3]) 
plot(Aspec, log="dB", add=TRUE, col=cols[4], lwd=lwds[4]) 
legend("topright", 
  c("spec.pgram","RSEIS::mtapspec","psdcore","pspectrum"), 
  title="Estimator", lwd=3, cex=1.1, col=cols)


## @knitr eval=TRUE, echo=TRUE, label="Interpolate results."
require(signal)
pltpi <- interp1(pltf, pltp, Pspec$freq)


## @knitr eval=TRUE, echo=TRUE, label="Summarize regression statistics."
df <- data.frame(x=dB(Pspec$spec), y=pltpi, tap=unclass(Aspec$taper))
summary(dflm <- lm(y ~ x + 0, df))
df$res <- residuals(dflm)


## @knitr eval=TRUE, echo=TRUE, label="Create ggplot objects."
require(ggplot2)
# data and regression
g <- ggplot(df, aes(x=x, y=y))
g1 <- g + geom_abline(intercept=0, slope=1, size=2, color="salmon")+
    geom_point(aes(color=dB(y/x))) + 
    geom_smooth(colour="black", formula = y ~ x + 0, method="lm", se=TRUE, fullrange=TRUE)
# and the residuals
g <- ggplot(df, aes(x=x, y=res))
g2 <- g + geom_abline(intercept=0, slope=0, size=2, color="salmon") +
    geom_point(aes(color=tap))


## @knitr eval=TRUE, echo=FALSE, fig=TRUE, fig.width=6, fig.height=5, label=RSEISvsRLP1
print(g1 + scale_colour_gradient2(mid="light grey") + theme_bw() + 
	ggtitle("Regression of mtapspec against psdcore"))


## @knitr eval=TRUE, echo=FALSE, fig=TRUE, fig.width=6, fig.height=2.5, label=RSEISvsRLP2
print(g2 + theme_bw() +
	ggtitle("Regression residuals, colored by optimized tapers"))


## @knitr eval=TRUE, echo=TRUE, fig=TRUE, fig.width=5, fig.height=4.5, label=SPECERR
sp <- spectral_properties(as.tapers(1:50), p=0.95, db.ci=TRUE)
plot(stderr.chi.upper ~ taper, sp, type="s", 
       ylim=c(-10,20), yaxs="i", xaxs="i",
       xlab=expression("number of tapers ("* nu/2 *")"), ylab="dB",
       main="Spectral uncertainties")
mtext("(additive factor)", line=.3)
lines(stderr.chi.lower ~ taper, sp, type="s")
lines(stderr.chi.median ~ taper, sp, type="s", lwd=2)
lines(stderr.chi.approx ~ taper, sp, type="s", col="red",lwd=2)
# to reach 3 db width confidence interval at p=.95
abline(v=33, lty=3)
legend("topright",
	c(expression("Based on "* chi^2 *"(p,"*nu*") and (1-p,"*nu*")"),
	  expression(""* chi^2 *"(p=0.5,"*nu*")"), 
	  "approximation"),
lwd=c(1,3,3), col=c("black","black","red"), bg="white")


## @knitr eval=TRUE, echo=TRUE, label="Compute spectral properties."
spp <- spectral_properties(Pspec$taper, db.ci=TRUE)
spa <- spectral_properties(Aspec$taper, db.ci=TRUE)
str(spa)
create_poly <- function(x, y, dy){
  xx <- c(x, rev(x))
  yy <- c(y+dy, rev(y-dy))
  return(data.frame(xx=xx, yy=yy))
}
pspp <- create_poly(Pspec$freq, dB(Pspec$spec), spp$stderr.chi.approx)
psppu <- create_poly(Pspec$freq, dB(Pspec$spec), spp$stderr.chi.upper)
pspa <- create_poly(Aspec$freq, dB(Aspec$spec), spa$stderr.chi.approx)
pspau <- create_poly(Aspec$freq, dB(Aspec$spec), spa$stderr.chi.upper)


## @knitr eval=TRUE, echo=TRUE, fig=TRUE, fig.width=6, fig.height=5.5, label=MAGERR
plot(c(0,0.5),c(-8,35),col="white", 
       main="MAGSAT Spectral Uncertainty (p > 0.95)",
       ylab="", xlab="spatial frequency, 1/km", yaxt="n", frame.plot=FALSE)
lines(c(2,1,1,2)*0.01,c(5,5,8.01,8.01)-8)
text(.05, -1.4, "3.01 dB")
polygon(psppu$xx, (psppu$yy), col="light grey", border="black", lwd=0.5)
polygon(pspp$xx, (pspp$yy), col="dark grey", border=NA)
text(0.15, 6, "With adaptive\ntaper refinement", cex=1.2)
polygon(pspau$xx, (pspau$yy)-10, col="light grey", border="black", lwd=0.5)
polygon(pspa$xx, (pspa$yy)-10, col="dark grey", border=NA)
text(0.35, 22, "Uniform tapering", cex=1.2)


## @knitr eval=TRUE, echo=TRUE, fig=TRUE, fig.width=7, fig.height=4, label=MAGRES
frq <- Aspec$freq
psd <- dB(Aspec$spec)
kopt <- unclass(Aspec$taper)
spres <- spa$resolution
pltdf <- data.frame(frq=frq, p=psd, relp=psd + dB(spres * frq), col="light grey")
pltdf$col[pltdf$relp<0] <- "light blue"
plot(relp ~ frq, pltdf, type="h",ylim=18*c(-1,1), col=pltdf$col, 
       xaxs="i", yaxs="i", main="MAGSAT Spectral Resolution",
       ylab="dB", xlab="spatial frequency, 1/km")
abline(h=0, lwd=1)
lines(relp ~ frq, pltdf, lwd=2)


## @knitr eval=TRUE, echo=TRUE, label="Get adaptive history."
pspectrum(ats, niter=6, plot=FALSE)
str(AH <- get_adapt_history())


## @knitr eval=TRUE, echo=TRUE, label="Some manipulation."
Freqs <- (AH$freq)
Dat <- AH$stg_psd
numd <- length(Freqs)
numit <- length(Dat)
StgPsd <- dB(matrix(unlist(Dat), ncol=numit))
Dat <- AH$stg_kopt
StgTap <- matrix(unlist(Dat), ncol=numit)
rm(Dat, AH)


## @knitr eval=TRUE, echo=FALSE, fig=TRUE, fig.width=6, fig.height=4, label=HIST1
seqcols <- 1:numit
itseq <- seqcols - 1
toadd <- matrix(rep(itseq, numd), ncol=numit, byrow=T)
par(xpd=TRUE)
matplot(Freqs, StgPsd + (sc<-9)*toadd, type="l", lty=1, lwd=2, col="black",
             main="PSD estimation history", ylab="", xlab="Spatial frequency",
             yaxt="n", frame.plot=FALSE)
text(0.52, 1.05*sc*itseq, itseq)
text(0.49, 1.1*sc*numit, "Stage:")


## @knitr eval=TRUE, echo=FALSE, fig=TRUE, fig.width=6, fig.height=3.5, label=HIST2
par(xpd=TRUE)
Cols <- rev(rev(brewer.pal(9, "PuBuGn"))[seqcols])
invisible(lapply(rev(seqcols), FUN=function(mcol, niter=numit, Frq=Freqs, Dat=StgTap, cols=Cols){
  iter <- (niter+1)-mcol
  y <- Dat[,mcol]
  icol <- Cols[mcol]
  if (iter==1){
    plot(Frq, y, type="h", col=icol, 
           main="Taper optimization history", ylab="", xlab="Spatial frequency",
           ylim=c(-50,650), frame.plot=FALSE)
  } else {
    lines(Frq, y, type="h", col=icol)
  }
  lines(Frq, y, type="l",  lwd=1.2)
  x <- (c(0,1)+iter-1)*.05+0.075
  y <- c(595,595,650,650,595)+10
  text(mean(x),max(y)+1.0*diff(range(y)), mcol-1)
  polygon(c(x,rev(x),x[1]),y,border="black",col=icol)
}))


## @knitr eval=TRUE, echo=FALSE, fig=TRUE, fig.width=3, fig.height=3, label=CORT
plot(cor(StgTap))


## @knitr eval=TRUE, echo=FALSE, fig=TRUE, fig.width=3, fig.height=3, label=CORP
plot(cor(StgPsd))


