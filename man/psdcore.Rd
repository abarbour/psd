\name{psdcore}
\alias{psdcore}
\alias{psdcore.default}
\title{Multitaper power spectral density estimates of a series.}
\usage{
  psdcore(X.d, X.frq = 1, ntaper = as.tapers(1),
    ndecimate = 1L, demean = TRUE, detrend = TRUE,
    na.action = stats::na.fail, first.last = TRUE,
    plotpsd = FALSE, as.spec = TRUE, refresh = FALSE,
    verbose = FALSE, ...)

  \method{psdcore}{default} (X.d, X.frq = 1,
    ntaper = as.tapers(1), ndecimate = 1L, demean = TRUE,
    detrend = TRUE, na.action = stats::na.fail,
    first.last = TRUE, plotpsd = FALSE, as.spec = TRUE,
    refresh = FALSE, verbose = FALSE, ...)
}
\arguments{
  \item{X.d}{the series to estimate a spectrum for}

  \item{X.frq}{scalar; the sampling information (see
  section Sampling)}

  \item{ntaper}{scalar, or vector; the number of tapers}

  \item{ndecimate}{scalar; decimation factor}

  \item{demean}{logical; should \code{X.d} be centered
  about the mean}

  \item{detrend}{logical; should \code{X.d} have a linear
  trend removed}

  \item{na.action}{function dealing with \code{NA} values}

  \item{first.last}{the extrapolates to give the zeroth and
  Nyquist frequency estimates}

  \item{plotpsd}{logical; should the estimate be shown
  compared to the \code{spec.pgram} estimate}

  \item{as.spec}{logical; should the object returned be of
  class 'spec'?}

  \item{refresh}{logical; ensure a free environment prior
  to execution}

  \item{verbose}{logical; should messages be given?}

  \item{...}{(unused) Optional parameters}
}
\value{
  An list object, invibibly.  If \code{as.spec=TRUE} then
  an object with class 'spec'.
}
\description{
  Compute power spectral density (PSD) estimates for the
  input series using sine multitapers.
}
\details{
  \subsection{Tapering}{ The parameter \code{ntaper}
  specifies the number of sine tapers to be used at each
  frequency: equal tapers at each frequency for a scalar;
  otherwise, use \code{ntaper(j)} sine tapers at
  \code{frequency(j)}. }

  \subsection{Truncation}{ The series length \code{N} is
  truncated, if necessary, so that \code{1+N/2} evenly
  spaced frequencies are returned. }

  \subsection{Decimation}{ The parameter \code{ndecimate}
  determines the PSDs actually computed, defined as
  \code{(1+n/2)/ndecimate}; other values are found via
  linear interpolation. }

  \subsection{Sampling}{ If \code{X.frq > 0} it's assumed
  the value represents \emph{frequency} (e.g. Hz).  If
  \code{X.frq < 0} it's assumed the value represents
  \emph{interval} (e.g. seconds). }
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

