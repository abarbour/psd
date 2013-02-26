\docType{methods}
\name{spec-methods}
\alias{as.data.frame.spec}
\alias{data.frame.spec}
\alias{spec-methods}
\title{Generic methods for objects with class \code{'spec'}.}
\usage{
  \method{as.data.frame}{spec} (x, ...)

  \method{data.frame}{spec} (x, ...)
}
\arguments{
  \item{x}{spec object}

  \item{...}{optional arguments}
}
\description{
  Generic methods for objects with class \code{'spec'}.
}
\details{
  Objects with class \code{'spec'} are simply list objects.
  \code{as.data.frame} converts the list into a
  \code{'data.frame'} with individual columns for the
  frequency, PSD, and taper vectors; all other information
  will be retained as an attribute. \code{data.frame} is an
  alias.
}
\examples{
\dontrun{
require(psd)
##
## Objects with class 'spec'
##
set.seed(1234)
#
x <- spectrum(xn<-rnorm(10), plot=FALSE)
xdf <-as.data.frame(x)
str(xdf)
is.tapers(xdf$taper)
#
# tapers class is retained
#
x <- psdcore(xn)
xdf <- as.data.frame(x)
str(xdf)
is.tapers(xdf$taper)
}
}
\author{
  A.J. Barbour <andy.barbour@gmail.com>
}
\keyword{S3methods}
\keyword{methods}
\keyword{spec}

