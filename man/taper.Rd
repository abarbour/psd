\name{newTaper}
\alias{newTaper}
\alias{taper-class}
\title{The 'taper' S4 class.}
\usage{
  newTaper(...)
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
newTaper()
new("taper") # equivalent to newTaper()
print(ntap <- newTaper(1:10))
plot(ntap)
}
\author{
  Andrew Barbour <andy.barbour@gmail.com>
}

