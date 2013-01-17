\name{psdcore}
\alias{psdcore}
\alias{psdcore.default}
\title{Multitaper power spectral density of a series.}
\usage{
  psdcore(X.d, X.frq = 1, ntaper = as.taper(1),
    ndecimate = 1L, demean = TRUE, detrend = TRUE,
    na.action = stats::na.fail, first.last = TRUE,
    Nyquist.normalize = TRUE, plotpsd = FALSE,
    as.spec = TRUE, force_calc = FALSE, ...)

  psdcore.default(X.d, X.frq = 1, ntaper = as.taper(1),
    ndecimate = 1L, demean = TRUE, detrend = TRUE,
    na.action = stats::na.fail, first.last = TRUE,
    Nyquist.normalize = TRUE, plotpsd = FALSE,
    as.spec = TRUE, force_calc = FALSE, ...)
}
\arguments{
  \item{X.d}{the series to estimate a spectrum for}

  \item{X.frq}{scalar; the sampling frequency (e.g. Hz)}

  \item{ntaper}{scalar, or vector; the number of tapers}

  \item{ndecimate}{scalar; decimation factor}

  \item{demean}{logical; should \code{X.d} be centered
  about the mean}

  \item{detrend}{logical; should \code{X.d} have a linear
  trend removed}

  \item{na.action}{function dealing with \code{NA} values}

  \item{first.last}{the extrapolates to give the zeroth and
  Nyquist frequency estimates}

  \item{Nyquist.normalize}{logical; should the units be
  returned in Hz, rather than Nyquist?}

  \item{plotpsd}{logical; should the estimate be shown
  compared to the \code{spec.pgram} estimate}

  \item{as.spec}{logical; should the object returned be of
  class 'spec'}

  \item{force_calc}{logical; force spectrum (used for
  development purposes)}

  \item{...}{(unused); parameters passed to [ NULL ]}
}
\description{
  Compute a power spectral density (PSD) estimate for the
  input series using sine multitapers.
}
\details{
  The parameter \code{ntaper} specifies the number of sine
  tapers to be used at each frequency: if it's a scalar,
  the same number of tapers will be used at every
  frequency; otherwise, use \code{ntaper(j)} sine tapers at
  \code{frequency(j)}.

  The series length \code{N} is truncated, if necessary, so
  that \code{1+N/2} evenly spaced frequencies are returned.

  The parameter \code{ndecimate} determines the PSDs
  actually computed, defined as \code{(1+n/2)/ndecimate};
  other values are found via linear interpolation.
}
\section{Warning}{
  Decimation is not well tested as of this point.
}
\examples{
X.d <- rnorm(1e3)
plot(psdcore(X.d, ntaper=10), log="dB", ylim=10*c(-1,1))
psd.n <- psdcore(X.d, ntaper=10, Nyquist.normalize=FALSE)
lines(psd.n$freq, 10*log10(psd.n$spec), col="red") # note normalization
abline(h=c(0, 3), col=c("black","red"), lwd=2)

# 10Hz sampling
plot(psdcore(X.d, X.frq=10, ntaper=10), log="dB", ylim=10*c(-0.3,1.7))
psd.n <- psdcore(X.d, X.frq=10, ntaper=10, Nyquist.normalize=FALSE)
lines(10*psd.n$freq, 10*log10(psd.n$spec), col="red") # note normalization
abline(h=c(10, 3), col=c("black","red"), lwd=2)

# if ntaper is a vector:
psdcore(X.d, ntaper=rep(10,length(X.d))
}
\author{
  Andrew Barbour <andy.barbour@gmail.com> ported original
  by R.L.Parker.
}
\seealso{
  \code{\link{pspectrum}}, \code{\link{riedsid}}
}
\keyword{normalization}
\keyword{prewhiten}
\keyword{spectrum-estimation}

