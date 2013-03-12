\name{as.tapers}
\alias{as.tapers}
\alias{tapers}
\title{Coerce an object into a \code{'tapers'} object.}
\usage{
  as.tapers(x, min_taper = 1, max_taper = NULL,
    setspan = FALSE)

  tapers(x, min_taper = 1, max_taper = NULL,
    setspan = FALSE)
}
\arguments{
  \item{x}{An object to set}

  \item{min_taper}{Set all values less than this to this.}

  \item{max_taper}{Set all values greater than this to
  this.}

  \item{setspan}{logical; should the tapers object be
  passed through \code{\link{minspan}} before it is
  return?}
}
\value{
  An object with class taper.
}
\description{
  In a tapered spectrum estimation algorithm, it is
  necessary to enforce rules on the number of tapers that
  may be applied.
}
\details{
  Formal requirements enforced by this function are:
  \itemize{ \item Non-zero. \item Integer values. \item
  Fewer than the half-length of the spectrum. } For
  example, we cannot apply zero tapers (the result would be
  a raw periodogram) or one million tapers (that would be
  absurd, and violate orthogonality conditions for any
  series less than two million terms long!).

  An object with S3 class \code{'tapers'} is created; this
  will have a minimum number of tapers in each position set
  by \code{min_taper}, and a maximum number of tapers in
  each position set by \code{max_taper}. If
  \code{minspan=TRUE}, the bounded taper is fed through
  \code{\link{minspan}} which will restrict the maximum
  tapers to less than or equal to the half-length of the
  spectrum.

  Various classes can be coerced into a \code{'tapers'}
  object; those tested sofar include: scalar, vector,
  matrix, data.frame, and list.

  Multiple objects are concatenated into a single vector
  dimension.

  Enabling \code{setspan} will only override
  \code{max_taper} should it be larger than the half-width
  of the series.
}
\note{
  No support (yet) for use of \code{min_taper,max_taper} as
  vectors, although this could be quite desirable.
}
\examples{
\dontrun{
require(psd)
##
## Objects with class 'tapers'
##
is.tapers(as.tapers(1))
is.tapers(as.tapers(1:10))
is.tapers(as.tapers(matrix(1:10,ncol=1)))
as.tapers(list(x=1:10,y=1:30)) # note dimensions
as.tapers(x<-data.frame(x=1:10,y=10:19))
as.tapers(x, min_taper=3, max_taper=10)
# class 'character' is in-coercible; raise error
try(as.tapers(c("a","b")), silent=TRUE)
}
}
\author{
  A.J. Barbour <andy.barbour@gmail.com>
}
\seealso{
  \code{\link{is.tapers}}
}
\keyword{S3methods}
\keyword{tapers}

