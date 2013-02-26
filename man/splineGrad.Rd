\name{splineGrad}
\alias{splineGrad}
\alias{splineGrad.default}
\title{Numerical derivatives of a series based on a weighted, smooth spline representation.}
\usage{
  splineGrad(dseq, dsig, plot.derivs = FALSE, ...)

  \method{splineGrad}{default} (dseq, dsig,
    plot.derivs = FALSE, ...)
}
\arguments{
  \item{dseq}{numeric; a vector of positions for
  \code{dsig}.}

  \item{dsig}{numeric; a vector of values (which will have
  a spline fit to them).}

  \item{plot.derivs}{logical; should the derivatives be
  plotted?}

  \item{...}{additional arguments passed to
  \code{smooth.spline}.}
}
\value{
  A matrix with columns representing \eqn{x, f(x), f'(x),
  f''(x)}.
}
\description{
  \code{splineGrad} computes the numerical derivatives of a
  spline representation of the input series;
  differentiation of spline curves is numerically
  efficient.

  With smoothing, the numerical instability for "noisy"
  data can be drastically reduced, since spline curves are
  inherently (at least) twice differentiable. See the
  \strong{Examples} for an illustration of this.
}
\examples{
\dontrun{
require(psd)
##
## Spline gradient
##
set.seed(1234)
x <- seq(0,5*pi,by=pi/64)
y <- cos(x) #**2
splineGrad(x, y, TRUE)
y <- y + rnorm(length(y), sd=.1)
# unfortunately, the presence of
# noise will affect numerical derivatives
splineGrad(x, y, TRUE)
# so change the smoothing used in smooth.spline
splineGrad(x, y, TRUE, spar=0.2)
splineGrad(x, y, TRUE, spar=0.6)
splineGrad(x, y, TRUE, spar=1.0)
##
}
}
\author{
  A.J. Barbour <andy.barbour@gmail.com>
}
\seealso{
  \code{smooth.spline}
}
\keyword{numerical-derivative}
\keyword{spline-gradient}
\keyword{utilities}

