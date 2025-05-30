% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/psd-package.R
\docType{data}
\name{wipp30}
\alias{wipp30}
\title{Water levels from borehole WIPP30}
\format{
A matrix with 13413 rows following 4 variables.

\describe{
\item{\code{time}}{Time (hours)}
\item{\code{wl}}{Water levels (psi)}
\item{\code{baro}}{Barometric pressure (psi)}
\item{\code{et}}{Earth tide gravity potential (nanometers/second^2)}
}
}
\description{
Observed water levels and barometric pressure from
well WIPP30 (WIPP: Waste Isolation Pilot Plant)
}
\details{
This is the dataset used in the multivariate PSD vignette
}
\examples{
data(wipp30)
summary(wipp30)
}
\references{
Toll, N.J., Rasmussen, T.C., (2007),
 Removal of Barometric Pressure Effects and Earth Tides from Observed Water Levels.
 \emph{Ground Water}, \strong{45}, 101–105,
 doi: 10.1111/j.1745-6584.2006.00254.x
}
\seealso{
\code{\link{pspectrum}}
}
\keyword{datasets}
