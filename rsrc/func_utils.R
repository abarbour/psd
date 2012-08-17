###
### Various utility functions
###
##
#
.splineGrad.default <- function(dseq, dsig, plot.derivs=FALSE, ...){
  #
  # Use spline interpolation to help find an emprirical gradient
  # (reduces numerical instability)
  #
  # @dseq: the sequence (index) for @dsig, the signal
  # output is the same length as the input
  #
  require(stats, qraphics)
  #
  # create a weighted cubic spline func
  smspl <- stats::smooth.spline(dseq, dsig, ...)
  SPLFUN <- stats::splinefun(smspl$x, smspl$y)
  #   seq.rng <- range(dseq)
  #   from <- seq.rng[1]
  #   to <- seq.rng[2]
  #   n <- length(dseq)
  #
  # signal spline
  #fsig <<- SPLFUN(dseq)
  #   FD0 <- function(){graphics::curve(SPLFUN(x), from=from, to=to, n=n, add=NA)}
  #   fsig <<- FD0()
  #
  # first deriv
  #   FD1 <- function(){graphics::curve(SPLFUN(x,deriv=1), from=from, to=to, n=n, add=NA)}
  #   fsigderiv <<- FD1()
  fsigderiv <<- SPLFUN(dseq, deriv=1)
  #
  # second deriv
  #   FD2 <- function(){graphics::curve(SPLFUN(x,deriv=2), from=from, to=to, n=n, add=NA)}
  #   fsigderiv2 <<- FD2()
  fsigderiv2 <<- SPLFUN(dseq, deriv=1)
  ##
  toret <- data.frame(x=dseq, y=dsig, dydx=fsigderiv, d2yd2x=fsigderiv2)
  #
  if (plot.derivs){
    plot(dseq, dsig, cex=0.6, main="signal and weighted cubic spline derivs")
    lines(y ~ x, toret, col="red")
    lines(dydx ~ x, toret, col="green")
    lines(d2yd2x ~ x, toret, col="blue")
    legend("topleft", paste(c("weighted cubic spline fit","first deriv", "second deriv"), sep=''), col = 2:4, lty = 1)
  }
  return(invisible(toret))
}
##
