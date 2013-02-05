\name{pspectrum}
\alias{pspectrum}
\alias{pspectrum.default}
\title{Adaptive sine multitaper power spectral density estimation.}
\usage{
  pspectrum(x, x.frqsamp = 1, ntap_pilot = 5, niter = 4,
    Nyquist.normalize = TRUE, verbose = TRUE,
    no.history = FALSE, plot = TRUE, ...)

  \method{pspectrum}{default} (x, x.frqsamp = 1,
    ntap_pilot = 5, niter = 4, Nyquist.normalize = TRUE,
    verbose = TRUE, no.history = FALSE, plot = TRUE, ...)
}
\arguments{
  \item{x}{vector; series to estimate PSD for.}

  \item{x.frqsamp}{scalar; the sampling rate (e.g. Hz) of
  the series \code{x}.}

  \item{ntap_pilot}{scalar; the number of sine tapers to
  use in the pilot spectrum estimation.}

  \item{niter}{scalar; the number of adaptive iterations to
  execute after the pilot spectrum.}

  \item{Nyquist.normalize}{logical; should the units be
  returned on Hz, rather than Nyquist?}

  \item{verbose}{logical; Should messages be given?}

  \item{no.history}{logical; Should the adaptive history
  \emph{not} be saved?}

  \item{plot}{logical; Should the results be plotted?}

  \item{...}{Optional parameters passed to
  \code{\link{riedsid}}}
}
\value{
  Object with class 'spec', invisibly. It also assigns the
  object to \code{"final_psd"} in the working environment.
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
# \dontrun{
##
## Adaptive multitaper PSD estimation
##
# load MAGSAT data
data(magsat)
#
#  use the defaults
pspectrum(X<-magsat$clean)
#
# turn off normalization
pspectrum(X, Nyquist.normalize=FALSE)
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
\keyword{tapers}
\keyword{tapers-constraints}
\keyword{tapers-weighting}

