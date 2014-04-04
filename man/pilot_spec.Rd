\name{pilot_spec}
\alias{pilot_spec}
\alias{pilot_spec.default}
\alias{pilot_spectrum}
\alias{spec.pilot}
\title{Calculate the pilot power spectral densities.}
\usage{
pilot_spec(x, x.frequency = 1, ntap = 7, remove.AR = 0, plot = FALSE,
  verbose = FALSE, ...)

\method{pilot_spec}{default}(x, x.frequency = 1, ntap = 7, remove.AR = 0,
  plot = FALSE, verbose = FALSE, ...)
}
\arguments{
  \item{x}{vector; the data series to find a pilot spectrum
  for}

  \item{x.frequency}{scalar; the sampling frequency (e.g.
  Hz) of the series}

  \item{ntap}{scalar; the number of tapers to apply during
  spectrum estimation}

  \item{remove.AR}{scalar; the max AR model to be removed
  from the data.}

  \item{plot}{logical; should a plot be created?}

  \item{verbose}{logical; should messages be given?}

  \item{...}{additional parameters passed to
  \code{\link{psdcore}}}
}
\value{
An object with class 'spec', invisibly.  It also assigns
the object to \code{"pilot_psd"} in the working
environment.
}
\description{
This PSD -- the pilot spectrum -- is used as the starting
point for the adaptive estimation routine.
}
\details{
A fixed number of tapers is applied across all frequencies
using \code{\link{psdcore}}, and subsequent
taper-refinements are based on the spectral derivatives of
this spectrum; hence, changes in the number of tapers can
affect how many adaptive stages may be needed (though there
are no formal convergence criteria to speak of).

The taper series of the returned spectrum is constrained
using \code{as.tapers(..., minspan=TRUE)}.

The default behaviour (\code{remove.AR <= 0}) is to remove
the standard linear model \eqn{[f(x) = \alpha x + \beta]}
from the data; however, the user can model the effect of an
autoregressive process by specifiying \code{remove.AR}.
}
\section{Removing an AR effect from the spectrum}{
  If \code{remove.AR > 0} the argument is used as
  \code{AR.max} in \code{\link{prewhiten}}, from which an
  AR-response spectrum is calculated using the best fitting
  model.

  If the value of \code{remove.AR} is too low the spectrum
  could become distorted, so use with care. \emph{Note,
  however, that the value of \code{remove.AR} will be
  restricted to within the range \eqn{[1,100]}.} If the AR
  order is much larger than this, it's unclear how
  \code{\link{prewhiten}} will perform and whether the AR
  model is appropriate.

  \emph{Note that this function does not produce a
  parametric spectrum estimation; rather, it will return
  the amplitude response of the best-fitting AR model as
  \code{spec.ar} would. \strong{Interpret these results
  with caution, as an AR response spectrum can be
  misleading.}}
}
\examples{
#RDEX#\dontrun{
require(psd)
##
## Pilot spectrum
##
data(magnet)
## simply calculate the pilot spectrum with a few tapers
plot(pilot_spec(xc <-  magnet$clean), log="dB", 
     main="Pilot PSDs for MAGNET and its AR-innovations (red)")
## remove the effect of an AR model works exceptionally
## well for the Project MAGNET data:
plot(pilot_spec(xc, remove.AR=10), log="dB", add=TRUE, col="red")
##
#RDEX#}
}
\author{
A.J. Barbour <andy.barbour@gmail.com>
}
\seealso{
\code{\link{psdcore}}, \code{\link{prewhiten}}

Documentation for \code{spec.ar}.
}

