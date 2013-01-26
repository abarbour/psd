\name{pspectrum}
\alias{pspectrum}
\alias{pspectrum.default}
\title{Adaptive sine multitaper power spectral density estimation.}
\usage{
  pspectrum(x, x.frqsamp = 1, ntap_pilot = 5, niter = 3,
    verbose = TRUE, no.history = FALSE, plot = TRUE, ...)

  \method{pspectrum}{default} (x, x.frqsamp = 1,
    ntap_pilot = 5, niter = 3, verbose = TRUE,
    no.history = FALSE, plot = TRUE, ...)
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

  \item{no.history}{logical; Should the adaptive history
  \emph{not} be saved?}

  \item{plot}{logical; Should the results be plotted?}

  \item{...}{Optional parameters passed to
  \code{\link{riedsid}}}
}
\value{
  Object with class 'spec', invisibly.
}
\description{
  This is the primary function to be used in this package,
  and returns power spectral density estimates where the
  number of tapers at each frequency has been iteratively
  optimized (\code{niter} times).
}
\details{
  See the \strong{Adaptive estimation} section in the
  description of the \code{\link{rlpSpec-package}} for
  details regarding adaptive estimation.
}
\examples{
# \\dontrun{
NULL
# }
}
\author{
  A.J. Barbour <andy.barbour@gmail.com> adapted original by
  R.L. Parker.
}
\seealso{
  \code{\link{psdcore}}, \code{\link{riedsid}},
  \code{\link{rlpSpec-package}}
}
\keyword{numerical-derivative}
\keyword{riedel-sidorenko}
\keyword{spectrum-estimation}
\keyword{taper}
\keyword{taper-constraints}
\keyword{taper-weighting}

