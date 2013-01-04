\name{parabolic_weights}
\alias{parabolic_weights}
\alias{parabolic_weights.taper}
\alias{parabolic_weights_fast}
\title{parabolic_weights}
\usage{
  parabolic_weights(tapvec, tap.index = 1L)

  parabolic_weights.taper(tapvec, tap.index = 1L)

  parabolic_weights_fast(ntap = 1L)
}
\arguments{
  \item{tapvec}{'taper' object; the number of tapers at
  each frequency}

  \item{tap.index}{integer; the index of \code{tapvec} from
  which to find weights}

  \item{ntap}{integer; the number of tapers to provide
  weightings for}
}
\value{
  a list with taper indices, and weighting parameters
}
\description{
  Calculate weighting factors for a series of tapers
}
\details{
  Weighting is calculated as follows:

  \deqn{n_T^2 - \frac{3 \cdot K_N^2}{2 \cdot n_T \cdot (n_T
  - 1/4) \cdot (n_T + 1)}}

  where \eqn{n_T} is the total number of tapers, and
  \eqn{K_N} is the integer sequence \eqn{[0,n_T-1]}
}
\author{
  Andrew Barbour <andy.barbour@gmail.com> ported original
  by R.L.Parker, and authored the optimized version.
}
\seealso{
  \code{\link{psdcore}}, \code{\link{riedsid}}
}

