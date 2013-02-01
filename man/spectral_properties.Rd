\name{spectral_properties}
\alias{spectral_properties}
\alias{spectral_properties.spec}
\alias{spectral_properties.tapers}
\title{Calculate spectral properties such as standard error and resolution.}
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
  \item{\code{stderr}: The standard error of the spectrum.}
  \item{\code{stderr.chi, .upper, .lower, .median}:A data
  frame with the results from \code{\link{spec_confint}}.}
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
  \subsection{Uncertainty}{ The errors are estimated in the
  simplest way, from the number of degrees of freedom; a
  more sophisticated (and complicated) approach is to
  estimate via jack-knifing (Prieto et al 2007) which is
  not yet available.

  Here the standard error \eqn{\delta S} is returned so
  \eqn{\delta S \cdot S} represents spectral uncertainty. }

  \subsection{Resolution}{ The frequency resolution depends
  on the number of tapers (\eqn{K}), and is found from
  \deqn{\frac{K \cdot f_N}{N_f}} where \eqn{f_N} is the
  Nyquist frequency and \eqn{N_f} is the number of
  frequencies estimated. }

  \subsection{Degrees of Freedom}{ There are two degrees of
  freedom for each taper. }

  \subsection{Bandwidth}{ The bandwidth of a multitaper
  estimate depends on the number of tapers. Following
  Walden et al (1995) the effective bandwidth is
  \eqn{\approx 2W} where \deqn{W = \frac{K + 1}{2N}} and
  \eqn{N} is the number of terms in the series, which makes
  \eqn{N \cdot W} the approximate time-bandwidth product. }
}
\author{
  A.J. Barbour <andy.barbour@gmail.com>
}
\seealso{
  \code{\link{spec_confint}}, \code{\link{rlpSpec-package}}
}
\keyword{bandwidth}
\keyword{degrees-of-freedom}
\keyword{properties}
\keyword{resolution}
\keyword{tapers}
\keyword{uncertainty}

