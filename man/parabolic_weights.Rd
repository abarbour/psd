\name{parabolic_weights}
\alias{parabolic_weights}
\alias{parabolic_weights.tapers}
\alias{parabolic_weights_fast}
\alias{parabolic_weights_fast.default}
\alias{parabolic_weights_fast.tapers}
\title{Calculate parabolic weighting factors.}
\usage{
  parabolic_weights(tapvec, tap.index = 1L)

  \method{parabolic_weights}{tapers} (tapvec,
    tap.index = 1L)

  parabolic_weights_fast(ntap = 1L)

  \method{parabolic_weights_fast}{default} (ntap = 1L)
}
\arguments{
  \item{tapvec}{'tapers' object; the number of tapers at
  each frequency}

  \item{tap.index}{integer; the index of \code{tapvec} from
  which to find weights}

  \item{ntap}{integer; the number of tapers to provide
  weightings for.}
}
\value{
  A list with taper indices, and the weights \eqn{W_N}.
}
\description{
  The resampled spectrum involves summing weighted tapers;
  this produces the weighting factors.
}
\details{
  Weighting factors are calculated as follows: \deqn{W_N
  \equiv n_T^2 - \frac{3 K_N^2}{2 n_T (n_T - 1/4) (n_T +
  1)}} where \eqn{n_T} is the total number of tapers, and
  \eqn{K_N} is the integer sequence \eqn{[0,n_T-1]}
}
\author{
  A.J. Barbour <andy.barbour@gmail.com> adapted original by
  R.L.Parker, and authored the optimized version.
}
\seealso{
  \code{\link{psdcore}}, \code{\link{riedsid}}
}
\keyword{tapers}
\keyword{tapers-weighting}

