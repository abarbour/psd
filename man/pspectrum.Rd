\name{pspectrum}
\alias{pspectrum}
\alias{pspectrum.default}
\alias{pspectrum.spec}
\title{Adaptive sine multitaper power spectral density estimation.}
\usage{
pspectrum(x, x.frqsamp = 1, ntap.init = 7, niter = 3, AR = FALSE,
  Nyquist.normalize = TRUE, verbose = TRUE, no.history = FALSE,
  plot = FALSE, ...)

\method{pspectrum}{default}(x, x.frqsamp = 1, ntap.init = 7, niter = 3,
  AR = FALSE, Nyquist.normalize = TRUE, verbose = TRUE,
  no.history = FALSE, plot = FALSE, ...)

\method{pspectrum}{spec}(x, x.frqsamp = 1, ntap.init = 7, niter = 3,
  AR = FALSE, Nyquist.normalize = TRUE, verbose = TRUE,
  no.history = FALSE, plot = FALSE, ...)
}
\arguments{
  \item{x}{vector; series to estimate PSD for.}

  \item{x.frqsamp}{scalar; the sampling rate (e.g. Hz) of
  the series \code{x}.}

  \item{ntap.init}{scalar; the number of sine tapers to use
  in the pilot spectrum estimation.}

  \item{niter}{scalar; the number of adaptive iterations to
  execute after the pilot spectrum.}

  \item{AR}{logical; should the effects of an AR model be
  removed from the pilot spectrum?}

  \item{Nyquist.normalize}{logical; should the units be
  returned on Hz, rather than Nyquist?}

  \item{verbose}{logical; Should messages be given?}

  \item{no.history}{logical; Should the adaptive history
  \emph{not} be saved?}

  \item{plot}{logical; Should the results be plotted?}

  \item{...}{Optional parameters passed to
  \code{\link{riedsid}}}
}
\value{
Object with class 'spec', invisibly. It also assigns the
object to \code{"final_psd"} in the working environment.
}
\description{
This is the primary function to be used in this package,
and returns power spectral density estimates where the
number of tapers at each frequency has been iteratively
optimized (\code{niter} times).
}
\details{
See the \strong{Adaptive estimation} section in the
description of the \code{\link{psd-package}} for details
regarding adaptive estimation.
}
\examples{
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
}
\author{
A.J. Barbour <andy.barbour@gmail.com> adapted original by
R.L. Parker.
}
\seealso{
\code{\link{psdcore}}, \code{\link{pilot_spec}},
\code{\link{riedsid}}, \code{\link{prewhiten}}
}
\keyword{numerical-derivative}
\keyword{riedel-sidorenko}
\keyword{spectrum-estimation}
\keyword{tapers}
\keyword{tapers-constraints}
\keyword{tapers-weighting}

