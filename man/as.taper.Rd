\name{as.taper}
\alias{as.taper}
\title{Coerce an object into a 'taper' object.}
\usage{
  as.taper(x, min_taper = 1, max_taper = max(x),
    setspan = FALSE)
}
\arguments{
  \item{x}{An object to set}

  \item{min_taper}{Set all values less than this to this.}

  \item{max_taper}{Set all values greater than this to
  this.}

  \item{setspan}{logical; should the taper object be passed
  through \code{\link{minspan}} before it is return?}
}
\description{
  In a tapered spectrum estimation algorithm, it is
  necessary to enforce rules on the number of tapers that
  may be applied.  For example, we cannot apply zero tapers
  (the result would be a raw periodogram) or one million
  tapers (that would be absurd).
}
\details{
  Formal requirements enforced by this function are:
  \itemize{ \item Non-zero. \item Integer values. \item
  Fewer than the half-length of the spectrum. }

  An object with S3 class 'taper' is created; this will
  have a minimum number of tapers in each position set by
  \code{min_taper}, and a maximum number of tapers in each
  position set by \code{max_taper}. If \code{minspan=TRUE},
  the bounded taper is fed through \code{\link{minspan}}
  which will restrict the maximum tapers to less than or
  equal to the half-length of the spectrum.

  Various classes can be coerced into a 'taper' object;
  those tested sofar include: scalar, vector, matrix,
  data.frame, and list.

  Multiple objects are concatenated into a single vector
  dimension.  For example, if the object is
  \code{list(x=c(1,2),y=c(3,4,5,0,1.1))} then the
  corresponding 'taper' object will be
  \code{1,2,3,4,5,1,1}, assuming \code{min_taper==1}.
}
\note{
  No support (yet) for \code{min_taper, max_taper} as
  vectors, although this could be desirable.
}
\examples{
is.taper(as.taper(1))
is.taper(as.taper(1:10))
is.taper(as.taper(matrix(1:10,ncol=1)))
as.taper(list(x=1:10,y=1:30)) # note dimensions
as.taper(x<-data.frame(x=1:10,y=10:19))
as.taper(x, min_taper=3, max_taper=10)
# class 'character' is in-coercible; raise error
as.taper(c("a","b"))
}
\author{
  Andrew Barbour <andy.barbour@gmail.com>
}
\seealso{
  \code{\link{is.taper}}
}
\keyword{S3methods}
\keyword{taper}

