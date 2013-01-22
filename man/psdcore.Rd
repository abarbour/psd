\name{psdcore}
\alias{psdcore}
\alias{psdcore.default}
\title{Multitaper power spectral density estimates of a series.}
\usage{
  psdcore(X.d, X.frq = 1, ntaper = as.taper(1),
    ndecimate = 1L, demean = TRUE, detrend = TRUE,
    na.action = stats::na.fail, first.last = TRUE,
    Nyquist.normalize = TRUE, plotpsd = FALSE,
    as.spec = TRUE, ...)

  psdcore.default(X.d, X.frq = 1, ntaper = as.taper(1),
    ndecimate = 1L, demean = TRUE, detrend = TRUE,
    na.action = stats::na.fail, first.last = TRUE,
    Nyquist.normalize = TRUE, plotpsd = FALSE,
    as.spec = TRUE, ...)
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

  \item{...}{(unused) Optional parameters}
}
\description{
  Compute power spectral density (PSD) estimates for the
  input series using sine multitapers.
}
\details{
  The parameter \code{ntaper} specifies the number of sine
  tapers to be used at each frequency: equal tapers at each
  frequency for a scalar; otherwise, use \code{ntaper(j)}
  sine tapers at \code{frequency(j)}.

  The series length \code{N} is truncated, if necessary, so
  that \code{1+N/2} evenly spaced frequencies are returned.

  The parameter \code{ndecimate} determines the PSDs
  actually computed, defined as \code{(1+n/2)/ndecimate};
  other values are found via linear interpolation.
}
\section{Warning}{
  Decimation is not well tested as of this point.
}
\author{
  A.J. Barbour <andy.barbour@gmail.com> adapted original by
  R.L.Parker.
}
\seealso{
  \code{\link{pspectrum}}, \code{\link{riedsid}}
}
\keyword{normalization}
\keyword{prewhiten}
\keyword{spectrum-estimation}

