\name{pilot_spec}
\alias{pilot_spec}
\alias{pilot_spec.default}
\alias{pilot_spectrum}
\alias{spec.pilot}
\title{Calculate the pilot power spectral densities.}
\usage{
  pilot_spec(...)

  pilot_spec.default(x, x.frequency = 1, ntap = 10, ...)
}
\arguments{
  \item{x}{vetor; the data series to find a pilot spectrum
  for}

  \item{x.frequency}{scalar; the sampling frequency (e.g.
  Hz) of the series}

  \item{ntap}{scalar; the number of tapers to apply during
  spectum estimation}

  \item{...}{additional parameters}
}
\value{
  An object with class 'spec'
}
\description{
  This PSD -- the pilot spectrum -- is used as the
  figurative starting point for the adaptive estimation
  routine. A fixed number of tapers is applied across all
  frequencies using \code{\link{psdcore}}, and subsequent
  taper-refinements are based on the spectral derivatives
  of this spectrum; hence, changes in the number of tapers
  can affect how many adaptive stages may be needed (though
  there are no formal convergence criteria to speak of).
}
\note{
  A mean value and linear are removed from the series prior
  to spectrum estimation, and no decimation is performed on
  the taper optimization. The taper series of the returned
  spectrum is constrained using \code{\link{as.taper(...,
  minspan=TRUE)}}.
}
\author{
  Andrew Barbour <andy.barbour@gmail.com>
}
\seealso{
  \code{\link{psdcore}}
}

