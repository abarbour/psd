\name{taper-methods}
\alias{lines}
\alias{plot}
\alias{points}
\alias{print}
\alias{summary}
\alias{taper}
\alias{taper-methods}
\title{Generic methods for objects with class 'taper'.}
\usage{
  print.taper(x)

  summary.taper(x)

  print.summary.taper(x)

  lines.taper(x, lwd = 1.8, col = "red", ...)

  points.taper(x, pch = "_", cex = 1, ...)

  plot.taper(x, color.pal = c("Blues", "Spectral"),
    ylim = NULL, ...)
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
  Andrew Barbour <andy.barbour@gmail.com>
}
\seealso{
  \code{\link{as.taper}}, \code{\link{constrain_tapers}}
}
\keyword{S3methods}
\keyword{methods}
\keyword{taper}

