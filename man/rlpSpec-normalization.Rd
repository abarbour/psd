\name{rlpSpec-normalization}
\alias{normalization}
\alias{normalize}
\alias{normalize.default}
\alias{normalize.list}
\alias{normalize.spec}
\alias{rlpSpec-normalization}
\title{Normalization of the power spectral densities.}
\usage{
  normalize(Spec, Fsamp = 1, src = NULL, verbose = TRUE,
    ...)

  \method{normalize}{default} (Spec, Fsamp = 1, src = NULL,
    verbose = TRUE, ...)

  \method{normalize}{list} (Spec, Fsamp = 1, src = NULL,
    verbose = TRUE, ...)

  \method{normalize}{spec} (Spec, Fsamp = 1, src = NULL,
    verbose = TRUE, ...)
}
\arguments{
  \item{Spec}{spectrum to normalize}

  \item{Fsamp}{sampling frequency}

  \item{src}{character string; the source of the spectrum
  estimator}

  \item{verbose}{logical; should messages be given?}

  \item{...}{(unused) additional parameters}
}
\value{
  original object, with values normalized accordingly
}
\description{
  Normalization is vitally important in spectral analyses.
}
\section{Assumtions}{
  The normalizations performed here assume:
}
\author{
  A.J. Barbour <andy.barbour@gmail.com>
}
\seealso{
  \code{\link{psdcore}} \code{\link{spectral_properties}}
}
\keyword{normalization}
\keyword{prewhiten}
\keyword{spectrum-estimation}

