###
### Various utility functions
###
##
#' prewhiten a timeseries object
#' 
#' de-mean, de-trend (which also de-means), and (soon) fit an AR
#' model to the series
#'
#' @note the \code{AR.fit} option is not used (yet)
#'
#' @param tser  \code{ts} object
#' @param AR.fit boolean; FALSE
#' @param detrend  boolean; TRUE
#' @param demean  boolean; TRUE
#' @param plot  boolean; TRUE
#' @param verbose  boolean; TRUE
#'
#' @extends ts
prewhiten <- function(tser,
                      AR.max=0L,
                      detrend=TRUE,
                      demean=TRUE,
                      plot=TRUE,
                      verbose=TRUE) UseMethod("prewhiten")
prewhiten.ts <- function(tser,
                         AR.max=0L,
                         detrend=TRUE,
                         demean=TRUE,
                         plot=TRUE,
                         verbose=TRUE){
  # prelims
  stopifnot(is.ts(tser))
  require(zoo)
  # some other info... needed?
  sps <- frequency(tser)
  tstart <- start(tser)
  n.o <- length(tser)
  ttime <- sps*n.o
  if (AR.max > 0) {
    AR.max <- as.integer(max(1,AR.max))
    if (verbose) message("autoregressive model fit (returning innovations)")
    # solves Yule-Walker equations
    #http://svn.r-project.org/R/trunk/src/library/stats/R/ar.R
    arfit <<- stats::ar.yw(tser, 
                       aic=TRUE, 
                       order.max=AR.max, 
                       demean=demean)
    if (verbose) print(arfit)
    # ar returns a TS object
    tser.prew <- as.ts(zoo::na.locf(arfit$resid))
  }
  # data.frame with fit params
  if (AR.max < 1){
    fit.df <- data.frame(xr=seq.int(from=1, to=n.o, by=1), 
                         xc=rep.int(1, n.o), 
                         y=tser)
    if (detrend){
      if (verbose) message("detrending (and demeaning)")
      X <- as.matrix(stats::residuals( stats::lm(y ~ xr, fit.df)))
    } else if (demean) {
      if (verbose) message("demeaning")
      X <- as.matrix(stats::residuals( stats::lm(y ~ xc, fit.df)))
    } else {
      X <- tser
      if (verbose) warning("nothing was done to the timeseries object")
    }
    tser.prew <- stats::ts(X, frequency=sps, start=tstart)
  }
  if (plot) plot(ts.union(tser, tser.prew), 
                 yaxs="i", xaxs="i",
                 yax.flip=TRUE)
  return(invisible(tser.prew))
}


#' Reports whether x is a 'spec' object
#' @param x An object to test
#' @export
is.spec <- function(x) inherits(x, "spec")
#
.splineGrad.default <- function(dseq, dsig, plot.derivs=FALSE, ...){
  #
  # Use spline interpolation to help find an emprirical gradient
  # (reduces numerical instability)
  #
  # @dseq: the sequence (index) for @dsig, the signal
  # output is the same length as the input
  #
  require(stats, graphics)
  #
  # create a weighted cubic spline
  smspl <- stats::smooth.spline(dseq, dsig, ...)
  # and a function
  SPLFUN <- stats::splinefun(smspl$x, smspl$y)
  # ?splinefun:
  # splinefun returns a function with formal arguments x and deriv, 
  # the latter defaulting to zero. This function can be used to 
  # evaluate the interpolating cubic spline (deriv = 0), or its 
  # derivatives (deriv = 1, 2, 3) at the points x, where the spline 
  # function interpolates the data points originally specified. This 
  # is often more useful than spline.
  #
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
  fsigderiv2 <<- SPLFUN(dseq, deriv=2)
  # how does the first deriv of the first deriv compare to the second?
  #   smspl.alt <- stats::smooth.spline(dseq, fsigderiv, ...)
  #   SPLFUN.alt <- stats::splinefun(smspl.alt$x, smspl.alt$y)
  #   fsigderiv2.alt <<- SPLFUN.alt(dseq, deriv=1)
  #   print(all.equal(fsigderiv2,fsigderiv2.alt))
  # [1] TRUE
  #   plot(fsigderiv2,fsigderiv2.alt, asp=1)
  ##
  toret <- data.frame(x=dseq, y=dsig, 
                      dydx=fsigderiv, 
                      d2yd2x=fsigderiv2)
                      #d2yd2x.alt=fsigderiv2.alt)
  #
  if (plot.derivs){
    #     yl.u <- max(c(dsig,fsigderiv,fsigderiv2))#,fsigderiv2.alt))
    #     yl.l <- min(c(dsig,fsigderiv,fsigderiv2))#,fsigderiv2.alt))
    nr <- 3 # f, f', f''
    mar.multi <- c(2., 5.1, 2, 2.1)
    oma.multi <- c(6, 0, 5, 0)
    oldpar <- par(mar = mar.multi, oma = oma.multi, mfcol = c(nr, 1))
    on.exit(par(oldpar))
    par(las=1)
    plot(dseq, dsig, cex=0.6, pch=3,
         #ylim=1.1*c(yl.l,yl.u),
         xaxs="i", yaxs="i",
         xlab="x", ylab="f(x)",
         main="splineGrad: signal and weighted cubic-spline fit")
    lines(y ~ x, toret, col="grey", lwd=1)
    plot(dydx ~ x, toret, 
         xaxs="i", yaxs="i",
         main="first derivative",
         col="red", type="s", lwd=2.4)
    plot(d2yd2x ~ x, toret, 
         xaxs="i", yaxs="i",
         main="second derivative", xlab="x",
         col="blue", type="s", lwd=2.4, lty=3)
    #lines(d2yd2x.alt ~ x, toret, col="blue", type="s", lwd=2.4, lty=3)
    #     legend("topleft", 
    #            paste(c("weighted cubic spline fit","first deriv", "second deriv"), #, "deriv of first deriv"), 
    #                  sep=''), 
    #            col = c("grey","red","blue"), #,"blue"), 
    #            lty = c(rep(1,3)), #, 3),
    #            cex=0.7)
  }
  return(invisible(toret))
}
##
