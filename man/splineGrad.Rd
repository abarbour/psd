\docType{methods}
\name{splineGrad}
\alias{splineGrad}
\alias{spline_gradients}
\title{Numerical derivatives of a series based on a weighted, smooth spline representation.}
\usage{
  splineGrad(dseq, dsig, plot.derivs = FALSE, ...)

  splineGrad.default(dseq, dsig, plot.derivs = FALSE, ...)
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
  spline representation of the input series.

  Numerical instability can be reduced with this method,
  since spline curves are inherently (at least) twice
  differentiable.

  How does the first derivative of the first derivative of
  the original series compare to the second derivative of
  the original series? Apparently, there is no difference.
  Try something like this: \code{smspl.alt <-
  stats::smooth.spline(dseq, fsigderiv, ...); SPLFUN.alt <-
  stats::splinefun(smspl.alt$x, smspl.alt$y);
  fsigderiv2.alt <- SPLFUN.alt(dseq, deriv=1);
  print(all.equal(fsigderiv2,fsigderiv2.alt));}
}
\examples{
##
## Spline gradient
require(stats)
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
\seealso{
  \code{\link{smooth.spline}}
}

