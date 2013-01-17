\name{as.taper}
\alias{as.taper}
\title{Coerce an object into a 'taper' object.}
\usage{
  as.taper(x, min_taper = 1, max_taper = NULL,
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
  tapers (that would be absurd, and violate orthogonality
  conditions for any series less than two million terms
  long!).
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
  corresponding 'taper' objects for the following arguments
  are:

  \describe{ \item{\emph{defaults}}{\code{[1,2,3,4,5,1,1]}}
  \item{\code{setspan=TRUE}}{\code{[1,2,3,3,3,1,1]}}
  \item{\code{max_taper=5}}{\code{[1,2,3,4,5,1,1]}}
  \item{\code{max_taper=5,setspan=TRUE}}{\code{[1,2,3,3,3,1,1]}}
  }

  It should be clear that enabling \code{setspan} will only
  override \code{max_taper} should it be larger than the
  half-width of the series.
}
\note{
  No support (yet) for use of \code{min_taper,max_taper} as
  vectors, although this could be quite desirable.
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

