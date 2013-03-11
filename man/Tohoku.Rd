\docType{data}
\name{Tohoku}
\alias{Tohoku}
\title{Observations of teleseismic strains from the 2011 Tohoku earthquake.}
\format{A dataframe with 16000 observations on the following 15 variables.

\describe{
\item{\code{}}{}
\item{\code{}}{}
\item{\code{}}{}
\item{\code{}}{}
\item{\code{}}{}
\item{\code{}}{}
\item{\code{}}{}
\item{\code{}}{}
\item{\code{}}{}
\item{\code{}}{}
\item{\code{}}{}
\item{\code{}}{}
\item{\code{}}{}
\item{\code{}}{}
\item{\code{}}{}
}

and X attributes:
\describe{
\item{\code{}}{}
\item{\code{}}{}
\item{\code{}}{}
\item{\code{}}{}
\item{\code{}}{}
}}
\source{
  PBO High Frequency archive:

  \url{http://borehole.unavco.org/bsm/earthquakes/NeartheEastCoastofHonshuJapan_20110311}
}
\description{
  The \eqn{M_w 9} Tohoku earthquake happend on March 11,
  2011.  The seismic waves were recorded at stations across
  the globe, including by strainmeters in the Plate
  Boundary Observatory (PBO) borehole strainmeters.
}
\details{
  These data are for station B084, which is located
  approximately 8500 km away from the epicenter. Because
  this distance is large, the seismic waves didn't arrive
  at this station for more than 700 seconds after the
  origin time.  So there is a record of pre-seismic noise
  included, the timeseries extends 6784 seconds prior to
  the origin time, and 9215 seconds after.

  The data are classified with the \code{"epoch"} variable,
  which separates the series into pre-seismic and seismic
  data; this is defined relative to the predicted P-wave
  arrival time from a traveltime model.
}
\examples{
data(Tohoku)
str(Tohoku)
}
\references{
  USGS summary page:

  \url{http://earthquake.usgs.gov/earthquakes/eqinthenews/2011/usc0001xgp/}
}
\seealso{
  \code{\link{pspectrum}}, \code{\link{hfsnm}},
  \code{\link{magnet}}
}
\keyword{datasets}

