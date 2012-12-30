\name{minspan}
\alias{minspan}
\alias{minspan.taper}
\title{minspan}
\usage{
  minspan(tapvec, nf = length(tapvec), ...)

  minspan.taper(tapvec, nf = length(tapvec), ...)
}
\arguments{
  \item{tapvec}{'taper' object; the number of tapers at
  each frequency}

  \item{nf}{scalar; number of positions or frequencies}

  \item{...}{(unused) optional argments}
}
\value{
  an object with class 'taper', with a constrained taper
  numbers
}
\description{
  Find the max span of taper object
}
\details{
  minimize the number of tapers: min(nf/2, 7*tapvec/5)
}
\note{
  The value of \code{nf} should represent the total number
  of tapers from the original spectrum; the default assumes
  the taper object represents that size.
}
\examples{
\dontrun{
ntap <- as.taper(1:10)
par(mfrow=c(2,1))
plot(ntap)
plot(minspan(ntap))
par(mfrow=c(1,1))
}
}
\author{
  Andrew Barbour <andy.barbour@gmail.com> ported original
  by R.L.Parker.
}
\seealso{
  \code{\link{constrain_tapers}}
}

