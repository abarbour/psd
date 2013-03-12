\name{psdcore}
\alias{psdcore}
\alias{psdcore.default}
\title{Multitaper power spectral density estimates of a series.}
\usage{
  psdcore(X.d, X.frq = NULL, ntaper = as.tapers(1),
    ndecimate = 1L, preproc = TRUE,
    na.action = stats::na.fail, first.last = TRUE,
    plotpsd = FALSE, as.spec = TRUE, refresh = FALSE,
    verbose = FALSE, ...)

  \method{psdcore}{default} (X.d, X.frq = NULL,
    ntaper = as.tapers(1), ndecimate = 1L, preproc = TRUE,
    na.action = stats::na.fail, first.last = TRUE,
    plotpsd = FALSE, as.spec = TRUE, refresh = FALSE,
    verbose = FALSE, ...)
}
\arguments{
  \item{X.d}{the series to estimate a spectrum for}

  \item{X.frq}{scalar; the sampling information (see
  section Sampling)}

  \item{ntaper}{scalar, or vector; the number of tapers}

  \item{ndecimate}{scalar; decimation factor}

  \item{preproc}{logical; should \code{X.d} have a linear
  trend removed?}

  \item{na.action}{the function to deal with \code{NA}
  values}

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
  An list object, invisibly.  If \code{as.spec=TRUE} then
  an object with class \code{spec}; otherwise the list
  object will have information similar to a \code{spec}
  object, but with a few additional fields.
}
\description{
  Compute power spectral density (PSD) estimates for the
  input series using sine multitapers.
}
\details{
  \subsection{Tapering}{ The parameter \code{ntaper}
  specifies the number of sine tapers to be used at each
  frequency: equal tapers at each frequency for a scalar;
  otherwise, use \code{ntaper[j]} sine tapers at
  \code{frequency[j]}. }

  \subsection{Truncation}{ The series, with length \eqn{N},
  is necessarily truncated so that \eqn{1+N/2} evenly
  spaced frequencies are returned.  This truncation makes
  the series length ``highly composite", which the discrete
  Fourier transform (DFT) is most efficient. The vignette
  "fftw" (accessed with
  \code{vignette("fftw",package="psd")}) shows how the
  performance of a DFT can be affected by series length. }

  \subsection{Decimation}{ The parameter \code{ndecimate}
  determines the number of PSD estimates actually computed.
  This number is defined as a fraction of the truncated
  length, \eqn{(1+N/2)/n_d}. Linear interpolation is used.
  }

  \subsection{Sampling}{ If \code{X.frq} is NULL, the value
  is assumed to be 1, unless \code{X.d} is a 'ts' object.
  If \code{X.frq > 0} it's assumed the value represents
  \emph{frequency} (e.g. Hz).  If \code{X.frq < 0} it's
  assumed the value represents \emph{interval} (e.g.
  seconds). }
}
\section{Warning}{
  Decimation is not well tested as of this point.

  The \code{first.last} parameter is a workaround for
  potential bug (under investigation), which causes the
  power at the zero and Nyquist frequencies to have
  anomalously low values. This argument enables using
  linear \emph{extrapolation} to correct these values.
  \strong{The feature will be deprecated if the supposed
  bug is both identified and fixed.}
}
\examples{
\dontrun{
require(psd)
##
## Multitaper PSD estimation
##
set.seed(1234)
X <- rnorm(1e3)
#
# use the defaults, and appeal to plot.spec
plot(psdcore(X))
#
# use more tapers, compare to stats::spectrum, and clear 
# env data from the previous calculation
psdcore(X, ntapers=10, plotpsd=TRUE, refresh=TRUE)
#
# change the sampling frequency to 20
psdcore(X, 20, 10, plotpsd=TRUE, refresh=TRUE) 
}
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

