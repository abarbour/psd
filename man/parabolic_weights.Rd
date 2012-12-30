\name{parabolic_weights}
\alias{parabolic_weights}
\alias{parabolic_weights.taper}
\alias{parabolic_weights_fast}
\title{parabolic_weights}
\usage{
  parabolic_weights(tapvec, tap.index = 1,
    vec.out = c("vertical", "horizontal"))

  parabolic_weights.taper(tapvec, tap.index = 1,
    vec.out = c("vertical", "horizontal"))

  parabolic_weights_fast(ntap)
}
\arguments{
  \item{tapvec}{'taper' object; the number of tapers at
  each frequency}

  \item{tap.index}{integer; the index of \code{tapvec} from
  which to find weights}

  \item{vec.out}{character; should the vector returned have
  long dimension as rows or columns}

  \item{ntap}{integer; the number of tapers to provide a
  weighting vector}
}
\value{
  an object with class 'matrix' with dimensions dependet on
  \code{vec.out}
}
\description{
  Calculate parabolic weighting factors
}
\author{
  Andrew Barbour <andy.barbour@gmail.com> ported original
  by R.L.Parker, and authored the optimized version.
}
\seealso{
  \code{\link{psdcore}}
}

