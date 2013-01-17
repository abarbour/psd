\name{spectral_properties}
\alias{spectral_properties}
\alias{spectral_properties.spec}
\alias{spectral_properties.taper}
\title{Calculate spectral properties such as standard error and resolution.}
\usage{
  spectral_properties(x, f.samp = 1, n.freq = NULL, ...)

  spectral_properties.spec(...)

  spectral_properties.taper(x, f.samp = 1, n.freq = NULL,
    ...)
}
\arguments{
  \item{x}{object with class taper}

  \item{f.samp}{scalar; the sampling frequency (e.g. Hz) of
  the series the tapers are for}

  \item{n.freq}{scalar; the number of frequencies of the
  original spectrum (if \code{NULL} the length of the taper
  object is assumed to be the number)}

  \item{...}{additional arguments (unused)}
}
\value{
  A list with the following properties (and names):
  \itemize{ \item{\code{taper}: The original taper vector.}
  \item{\code{stderr}: The standard error of the spectrum.
  Multiply by PSD for spectral uncertainty.}
  \item{\code{resolution}: The effective spectral
  resolution.} \item{\code{dof}: The number of degrees of
  freedom. May be used for jacknifing confidence
  intervals.} \item{\code{bw}: The effective bandwidth of
  the spectrum.} }
}
\description{
  Various spectral properties may be computed from the
  vector of tapers, and if necessary the sampling
  frequency.
}
\section{Uncertainty}{

}

\section{Resolution}{

}

\section{DoF}{

}

\section{Bandwidth}{

}

