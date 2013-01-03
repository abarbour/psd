\name{is.rlpspec}
\alias{is.rlpspec}
\title{Reports whether x is an object with class 'rlpspec'}
\usage{
  is.rlpspec(x)
}
\arguments{
  \item{x}{An object to test}
}
\description{
  Reports whether x is an object with class 'rlpspec'
}
\examples{
is.rlpspec(1:10)
is.rlpspec(rlp <- newRlpspec())
is(1:10, "rlpspec")
is(rlp, "rlpspec")
}
\author{
  Andrew Barbour <andy.barbour@gmail.com>
}
\seealso{
  \code{\link{as.rlpspec}}
}

