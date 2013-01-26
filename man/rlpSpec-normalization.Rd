\name{rlpSpec-normalization}
\alias{normalization}
\alias{nyqnorm}
\alias{nyqnorm.default}
\alias{nyqnorm.spec}
\alias{psdnorm}
\alias{psdnorm.default}
\alias{psdnorm.spec}
\alias{rlpSpec-normalization}
\title{Normalization of the power spectral densities.}
\usage{
  psdnorm(Spec, Fnyq, verbose = TRUE, ...)

  \method{psdnorm}{default} (Spec, Fnyq, verbose = TRUE,
    ...)

  \method{psdnorm}{spec} (Spec, Fnyq, verbose = TRUE, ...)

  nyqnorm(Spec, Fnyq, verbose = TRUE, ...)

  \method{nyqnorm}{default} (Spec, Fnyq, verbose = TRUE,
    ...)

  \method{nyqnorm}{spec} (Spec, Fnyq, verbose = TRUE, ...)
}
\arguments{
  \item{Spec}{spectrum to normalize}

  \item{Fnyq}{nyquist frequency}

  \item{verbose}{logical; should messages be given?}

  \item{...}{(unused) additional parameters}
}
\value{
  original object, with values normalized accordingly
}
\description{
  Something
}
\author{
  A.J. Barbour <andy.barbour@gmail.com>
}
\seealso{
  \code{\link{psdcore}}
}
\keyword{normalization}
\keyword{prewhiten}
\keyword{spectrum-estimation}

