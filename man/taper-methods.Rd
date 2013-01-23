\docType{methods}
\name{taper-methods}
\alias{lines}
\alias{lines.taper}
\alias{plot}
\alias{plot.taper}
\alias{points}
\alias{points.taper}
\alias{print}
\alias{print.summary.taper}
\alias{print.taper}
\alias{summary}
\alias{summary.taper}
\alias{taper}
\alias{taper-methods}
\title{Generic methods for objects with class 'taper'.}
\usage{
  \method{print}{taper} (x, ...)

  \method{summary}{taper} (object, ...)

  \method{print}{summary.taper} (x, ...)

  \method{lines}{taper} (x, lwd = 1.8, col = "red", ...)

  \method{points}{taper} (x, pch = "_", cex = 1, ...)

  \method{plot}{taper} (x,
    color.pal = c("Blues", "Spectral"), ylim = NULL, ...)
}
\arguments{
  \item{x}{taper object}

  \item{object}{taper object}

  \item{lwd}{line width (default is 1.8)}

  \item{col}{color of line (default is "red")}

  \item{pch}{point character (default is "_")}

  \item{cex}{point size (default is 1)}

  \item{color.pal}{color palette to use (choices are:
  "Blues","Spectral")}

  \item{ylim}{optional limits for y-axis}

  \item{...}{optional arguments}
}
\description{
  Generic methods for objects with class 'taper'.
}
\examples{
##
tap <- as.taper(c(1:49,50:0)+rnorm(1e2))
print(tap)
print(summary(tap))
plot(tap)
# no arithmetic methods
tap <- as.taper(tap/2)
lines(tap)
}
\author{
  A.J. Barbour <andy.barbour@gmail.com>
}
\seealso{
  \code{\link{as.taper}}, \code{\link{constrain_tapers}},
  \code{par}
}
\keyword{S3methods}
\keyword{methods}
\keyword{taper}

