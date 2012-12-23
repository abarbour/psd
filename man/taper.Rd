\name{taper}
\alias{taper}
\alias{taper-class}
\title{The 'taper' S4 class.}
\usage{
  taper(...)
}
\description{
  In this class the value of each position will be a
  non-zero, positive integer, since it represents the
  number of tapered sections to average.
}
\note{
  The prototypical S4 class has tapers==1, and length==1.
  Currently there are no \code{@slots}; this may change in
  the future.
}
\examples{
taper()
new("taper") # equivalent to taper()
print(ntap <- taper(1:10))
plot(ntap)
}
\author{
  Andrew Barbour
}

