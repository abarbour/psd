#RDEX#\dontrun{
require(psd)
##
## Adaptive multitaper PSD estimation
## (portions extracted from overview vignette)
##
require(RColorBrewer)
##
## adaptive estimation for the Project MAGNET dataset
data(magnet)
# adaptive psd estimation (turn off diagnostic plot)
PSDr <- pspectrum(Xr <- magnet$raw, plot=FALSE)
PSDc <- pspectrum(Xc <- magnet$clean, plot=FALSE)
# plot them on the same scale
plot(PSDc, log="dB", main="Raw and Clean Project MAGNET power spectral density",
     lwd=3, ci.col=NA, ylim=c(0,32), yaxs="i")
plot(PSDr, log="dB", add=TRUE, lwd=3, lty=5)
text(c(0.25,0.34), c(11,24), c("Clean","Raw"), cex=1)

## Change sampling, and inspect the diagnostic plot
pspectrum(Xc, niter=1, x.frqsamp=10)

## Say we forgot to assign the results: we can recover from the environment with:
PSDc_recovered <- psd:::psd_envGet("final_psd")
plot(PSDc_recovered)

##
## Visualize adaptive history
##
## Previous adaptive estimation history
pspectrum(Xc, niter=6, plot=FALSE)
AH <- get_adapt_history()
Freqs <- (AH$freq)
Dat <- AH$stg_psd
numd <- length(Freqs)
numit <- length(Dat)
StgPsd <- dB(matrix(unlist(Dat), ncol=numit))
Dat <- AH$stg_kopt
StgTap <- matrix(unlist(Dat), ncol=numit)
rm(Dat, AH)

## plot psd history
seqcols <- 1:numit
itseq <- seqcols - 1
toadd <- matrix(rep(itseq, numd), ncol=numit, byrow=TRUE)
par(xpd=TRUE)
matplot(Freqs, StgPsd + (sc<-9)*toadd, type="l", lty=1, lwd=2, col="black",
             main="PSD estimation history", ylab="", xlab="Spatial frequency",
             yaxt="n", frame.plot=FALSE)
text(.52, 1.05*sc*itseq, itseq)
text(.49, 1.1*sc*numit, "Stage:")

## plot taper history "mountain range" silhouettes
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
}
)) # end of invisible lapply
#RDEX#}
