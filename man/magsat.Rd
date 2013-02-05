\docType{data}
\name{magsat}
\alias{magsat}
\title{A single line of MAGSAT horizontal field intensity.}
\format{A dataframe with 2048 observations on the following 4 variables.

\describe{
\item{\code{km}}{Relative along-track distance, in kilometers. The first observation is at zero kilometers.}
\item{\code{raw}}{Raw intensities, in nanoTesla.}
\item{\code{clean}}{Edited raw intensites, in nanoTesla}
\item{\code{mdiff}}{The difference between \code{clean} and \code{raw} intensities, in nanoTesla.}
}}
\source{
  Data subsetter:
  \url{http://ftpbrowser.gsfc.nasa.gov/magsat.html}
}
\description{
  The Magnetic Satellite (MAGSAT) mission provided a wealth
  of airborne-magnetometer data spanning the globe (Langel
  el al., 1982). This dataset represents a single track of
  horizontal field intensities (a very small subset of the
  full collection!).
}
\details{
  \subsection{Raw and Clean Sets}{ There are non-real data
  points in raw MAGSAT series; these are instrumental
  artefacts, and can severely affect power spectral density
  (PSD) estimates. A clean series has been included so that
  a comparison of PSDs may be made.

  Some command like \code{subset(magsat, abs(mdiff) > 0)}
  can be used to identify the rows where edits have been
  made. }
}
\examples{
data(magsat)
summary(magsat)
}
\references{
  Langel, R., G. Ousley, J. Berbert, J. Murphy, and M.
  Settle, (1982), The MAGSAT mission, \emph{Geophysical
  Research Letters}, \strong{9} (4), 243-245.
}
\keyword{datasets}
\keyword{magsat}

