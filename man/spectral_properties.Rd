\name{spectral_properties}
\alias{spectral_properties}
\alias{spectral_properties.spec}
\alias{spectral_properties.tapers}
\title{Calculate spectral properties.}
\usage{
  spectral_properties(tapvec, f.samp = 1, n.freq = NULL,
    p = 0.95, db.ci = FALSE, ...)

  \method{spectral_properties}{spec} (tapvec, ...)

  \method{spectral_properties}{tapers} (tapvec, f.samp = 1,
    n.freq = NULL, p = 0.95, db.ci = FALSE, ...)
}
\arguments{
  \item{tapvec}{object with class \code{tapers} or
  \code{spec}}

  \item{f.samp}{scalar; the sampling frequency (e.g. Hz) of
  the series the tapers are for}

  \item{n.freq}{scalar; the number of frequencies of the
  original spectrum (if \code{NULL} the length of the
  tapers object is assumed to be the number)}

  \item{p}{numeric; the coverage probability, bound within
  \eqn{[0,1)}}

  \item{db.ci}{logical; should the uncertainty confidence
  intervals be returned as decibels?}

  \item{...}{additional arguments (unused)}
}
\value{
  A list with the following properties (and names):
  \itemize{ \item{\code{taper}: The original taper vector.}
  \item{\code{stderr.chi .upper, .lower, .median}: results
  returned from \code{\link{spec_confint}}.}
  \item{\code{resolution}: The effective spectral
  resolution.} \item{\code{dof}: The number of degrees of
  freedom.} \item{\code{bw}: The effective bandwidth of the
  spectrum.} }
}
\description{
  Various spectral properties may be computed from the
  vector of tapers, and if necessary the sampling
  frequency.
}
\section{Parameter Details}{
  \subsection{Uncertainty}{ See \code{\link{spec_confint}}
  for details. }

  \subsection{Resolution}{ The frequency resolution depends
  on the number of tapers (\eqn{K}), and is found from
  \deqn{\frac{K \cdot f_N}{N_f}} where \eqn{f_N} is the
  Nyquist frequency and \eqn{N_f} is the number of
  frequencies estimated. }

  \subsection{Degrees of Freedom}{ There are two degrees of
  freedom for each taper \eqn{K}: \deqn{\nu = 2 K} }

  \subsection{Bandwidth}{ The bandwidth of a multitaper
  estimate depends on the number of tapers. Following
  Walden et al (1995) the effective bandwidth is
  \eqn{\approx 2W} where \deqn{W = \frac{K + 1}{2 N}} and
  \eqn{N} is the number of terms in the series, which makes
  \eqn{N \cdot W} the approximate time-bandwidth product. }
}
\examples{
\dontrun{
require(psd)
##
## Spectral properties from the number of tapers used
## (portions extracted from overview vignette)
##

##
## Theoretical uncertainties from Chi^2 distribution
##
sp <- spectral_properties(as.tapers(1:50), p=0.95, db.ci=TRUE)
par(las=1)
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

##
## An example using the Project MAGNET dataset
##
data(magnet)
tapinit <- 15 # tapers
dt <- 1 # 1/km

# remove mean/trend (not really necessary but good practice; also, done internally)
ats <- prewhiten(ts(magnet$clean, deltat=dt), plot=FALSE)$prew_lm

# normal and adaptive multitaper spectra
Pspec <- psdcore(ats, dt, tapinit)
Aspec <- pspectrum(ats, dt, tapinit, niter=3, plot=FALSE)

# calculate spectral properties
spp <- spectral_properties(Pspec$taper, db.ci=TRUE)
spa <- spectral_properties(Aspec$taper, db.ci=TRUE)

# function to create polygon data, and create them
create_poly <- function(x, y, dy){
  xx <- c(x, rev(x))
  yy <- c(y+dy, rev(y-dy))
  return(data.frame(xx=xx, yy=yy))
}
pspp <- create_poly(Pspec$freq, dB(Pspec$spec), spp$stderr.chi.approx)
psppu <- create_poly(Pspec$freq, dB(Pspec$spec), spp$stderr.chi.upper)
pspa <- create_poly(Aspec$freq, dB(Aspec$spec), spa$stderr.chi.approx)
pspau <- create_poly(Aspec$freq, dB(Aspec$spec), spa$stderr.chi.upper)

##
## Project MAGNET uncertainties
##
plot(c(0,0.5),c(-8,35),col="white",
       main="Project MAGNET Spectral Uncertainty (p > 0.95)",
       ylab="", xlab="spatial frequency, 1/km", yaxt="n", frame.plot=FALSE)
lines(c(2,1,1,2)*0.01,c(5,5,8.01,8.01)-8)
text(.05, -1.4, "3.01 dB")
polygon(psppu$xx, (psppu$yy), col="light grey", border="black", lwd=0.5)
polygon(pspp$xx, (pspp$yy), col="dark grey", border=NA)
text(0.15, 6, "With adaptive\\ntaper refinement", cex=1.2)
polygon(pspau$xx, (pspau$yy)-10, col="light grey", border="black", lwd=0.5)
polygon(pspa$xx, (pspa$yy)-10, col="dark grey", border=NA)
text(0.35, 22, "Uniform tapering", cex=1.2)

##
## Project MAGNET resolution
##
frq <- Aspec$freq
relp <- dB(1/spa$resolution)
par(las=1)
plot(frq, relp,
     col="light grey",
     ylim=dB(c(1,5)),
     type="h", xaxs="i", yaxs="i",
     ylab="dB", xlab="frequency, 1/km",
     main="Project MAGNET Spectral Resolution and Uncertainty")
lines(frq, relp)
lines(frq, spp$stderr.chi.upper+relp, lwd=1.5, lty=3)
lines(frq, spa$stderr.chi.upper+relp, lwd=3, lty=2)
abline(h=dB(sqrt(vardiff(Aspec$spec))), lwd=1.5, lty=2, col="red")

##
}
}
\author{
  A.J. Barbour <andy.barbour@gmail.com>
}
\references{
  Prieto, G. A., R. L. Parker, D. J. Thomson, F. L. Vernon,
  and R. L. Graham (2007), Reducing the bias of multitaper
  spectrum estimates, \emph{Geophysical Journal
  International}, \strong{171}, 1269--1281, doi:
  10.1111/j.1365-246X.2007.03592.x
}
\seealso{
  \code{\link{spec_confint}}, \code{\link{psd-package}}
}
\keyword{bandwidth}
\keyword{degrees-of-freedom}
\keyword{properties}
\keyword{resolution}
\keyword{tapers}
\keyword{uncertainty}

