\name{pspectrum}
\alias{pspectrum}
\alias{pspectrum.default}
\title{Adaptive sine multitaper power spectral density estimation.}
\usage{
  pspectrum(x, x.frqsamp = 1, ntap_pilot = 10, niter = 4,
    verbose = TRUE, ...)

  pspectrum.default(x, x.frqsamp = 1, ntap_pilot = 10,
    niter = 2, verbose = TRUE, ...)
}
\arguments{
  \item{x}{vector; series to estimate PSD for.}

  \item{x.frqsamp}{scalar; the sampling rate (e.g. Hz) of
  the series \code{x}.}

  \item{ntap_pilot}{scalar; the number of sine tapers to
  use in the pilot spectrum estimation.}

  \item{niter}{scalar; the number of adaptive iterations to
  execute after the pilot spectrum.}

  \item{verbose}{logical; Should messages be given?}

  \item{...}{optional arguments (unused)}
}
\value{
  Object with class 'spec'.
}
\description{
  Adaptive sine multitaper power spectral density
  estimation.
}
\author{
  Andrew Barbour <andy.barbour@gmail.com> based on R.L.
  Parker's algorithm.
}
\seealso{
  \code{\link{psdcore}}, \code{\link{riedsid}}
}

