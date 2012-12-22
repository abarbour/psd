\name{as.taper}
\alias{as.taper}
\title{Coerces an object into a 'taper' object}
\usage{
  as.taper(x, min.taper = 1)
}
\arguments{
  \item{x}{An object to set}

  \item{min.taper}{Set all values less than this to this.}
}
\description{
  An object of class 'taper' is created, with a minimum
  number of tapers in each position set by
  \code{min.taper}.
}
\details{
  Various classes can be coerced into a 'taper' object;
  those tested sofar include: scalar, vector, matrix,
  data.frame, and list.  Multiple objects are concatenated
  into a single vector dimension.  For example, if the
  object is \code{list(x=c(1,2),y=c(3,4,5))} then the
  corresponding 'taper' will be \code{1,2,3,4,5}.
}
\examples{
is.taper(as.taper(1))
is.taper(as.taper(1:10))
is.taper(as.taper(matrix(1:10,ncol=1)))
is.taper(as.taper(list(x=1:10,y=1:30))) # note dimensions
is.taper(as.taper(data.frame(x=1:10,y=10:19)))
# class 'character' is in-coercible; raise error
as.taper(c("a","b"))
}
\seealso{
  \code{\link{is.taper}}, \code{\link{taper}}
}

