\name{spec_confint}
\alias{spec_confint}
\alias{spec_confint.default}
\alias{spec_confint.spec}
\alias{spec_confint.tapers}
\title{Calculate confidence intervals for a given number of degrees of freedom.}
\usage{
  spec_confint(dof, p = 0.95, as.db = FALSE)

  \method{spec_confint}{spec} (dof, p = 0.95,
    as.db = FALSE)

  \method{spec_confint}{tapers} (dof, p = 0.95,
    as.db = FALSE)

  \method{spec_confint}{default} (dof, p = 0.95,
    as.db = FALSE)
}
\arguments{
  \item{dof}{numeric; the degrees of freedom}

  \item{p}{numeric; the coverage probability, bound within
  \eqn{[0,1)}}

  \item{as.db}{logical; should the values be returned as
  decibels?}
}
\value{
  A \code{data.frame} with the following properties (and
  names): \itemize{ \item{\code{upper}: The original taper
  vector.} \item{\code{lower}: The standard error of the
  spectrum.} \item{\code{median}: The standard error of the
  spectrum from \code{\link{spec_confint}}.} }
}
\description{
  In a multitaper scheme, the degrees of freedom is two per
  taper.
}
\author{
  A.J. Barbour <andy.barbour@gmail.com> modified from the
  \code{spec.ci} function inside \code{stats::plot.spec}.
}
\seealso{
  \code{\link{spectral_properties}},
  \code{\link{rlpSpec-package}}, \code{plot.spec}
}
\keyword{degrees-of-freedom}
\keyword{properties}
\keyword{tapers}
\keyword{uncertainty}

