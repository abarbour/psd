\docType{methods}
\name{tapers-methods}
\alias{as.data.frame.tapers}
\alias{data.frame.tapers}
\alias{lines.tapers}
\alias{plot.tapers}
\alias{points.tapers}
\alias{print.summary.tapers}
\alias{print.tapers}
\alias{summary.tapers}
\alias{tapers-methods}
\title{Generic methods for objects with class \code{'tapers'}.}
\usage{
  \method{as.data.frame}{tapers} (x, ...)

  \method{data.frame}{tapers} (x, ...)

  \method{print}{tapers} (x, ...)

  \method{summary}{tapers} (object, ...)

  \method{print}{summary.tapers} (x, ...)

  \method{lines}{tapers} (x, lwd = 1.8, col = "red", ...)

  \method{points}{tapers} (x, pch = "_", cex = 1, ...)

  \method{plot}{tapers} (x, xi = NULL,
    color.pal = c("Blues", "Spectral"), ylim = NULL,
    hv.lines = FALSE, ...)
}
\arguments{
  \item{x}{tapers object}

  \item{xi}{optional vector for indices of \code{x}}

  \item{object}{tapers object}

  \item{lwd}{line width (default is 1.8)}

  \item{col}{color of line (default is "red")}

  \item{pch}{point character (default is "_")}

  \item{cex}{point size (default is 1)}

  \item{color.pal}{color palette to use (choices are:
  "Blues","Spectral")}

  \item{ylim}{optional limits for y-axis}

  \item{hv.lines}{logical; should horizontal and vertival
  reference lines be plotted?}

  \item{...}{optional arguments}
}
\value{
  \code{plot} returns a list with names: \code{line.colors}
  (hex values)
}
\description{
  Generic methods for objects with class \code{'tapers'}.
}
\examples{
##
tap <- as.tapers(c(1:49,50:0)+rnorm(1e2))
print(tap)
print(summary(tap))
plot(tap)
# no arithmetic methods
tap <- as.tapers(tap/2)
lines(tap)
}
\author{
  A.J. Barbour <andy.barbour@gmail.com>
}
\seealso{
  \code{\link{as.tapers}}, \code{\link{constrain_tapers}},
  \code{par}
}
\keyword{S3methods}
\keyword{methods}
\keyword{tapers}

