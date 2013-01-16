#' Prewhiten a series.
#' 
#' Remove (optionally) mean, trend, and Auto Regressive (AR) model
#' from the original series.
#'
#' @section NA values:
#' \code{NA} values are allowed.  If present, the function
#' \code{\link[zoo]{na.locf}}, which stands for "Last Observation Carried Forward",
#' is used to impute them with real numbers.
#'
#' @name prewhiten
#' @import zoo
#' 
#' @param tser  vector; An object to prewhiten.
#' @param AR.max numeric; the maximum AR order to fit.
#' @param detrend  logical; Should a trend (and mean) be removed?
#' @param demean  logical; Should a mean value be removed?
#' @param plot  logical; Should the results be plotted?
#' @param verbose  logical; Should messages be printed?
#' @param x.fsamp sampling frequency (for non \code{ts} objects)
#' @param x.start start time of observations (for non \code{ts} objects)
#' @param ... variables passed to \code{\link{prewhiten.ts}} (for non \code{ts} objects)
#'
prewhiten <- function(tser, AR.max=0L, detrend=TRUE, demean=TRUE, plot=TRUE, verbose=TRUE) UseMethod("prewhiten")
#' @rdname prewhiten
#' @S3method prewhiten default
prewhiten.default <- function(tser, x.fsamp=1, x.start=c(1, 1), ...){
  Xts <- stats::ts(tser, frequency=x.fsamp, start=x.start)
  prewhiten(Xts, ...)
}
#' @rdname prewhiten
#' @S3method prewhiten ts
prewhiten.ts <- function(tser, AR.max=0L, detrend=TRUE, demean=TRUE, plot=TRUE, verbose=TRUE){
  # prelims
  stopifnot(stats::is.ts(tser))
  require(zoo)
  # some other info... needed?
  sps <- stats::frequency(tser)
  tstart <- stats::start(tser)
  n.o <- length(tser)
  ttime <- sps*n.o
  if (AR.max > 0) {
    AR.max <- as.integer(max(1,AR.max))
    if (verbose) message("autoregressive model fit (returning innovations)")
    # solves Yule-Walker equations
    #http://svn.r-project.org/R/trunk/src/library/stats/R/ar.R
    arfit <- stats::ar.yw(tser, aic=TRUE, order.max=AR.max, demean=demean)
    if (verbose) print(arfit)
    # ar returns a TS object
    tser.prew <- stats::as.ts(zoo::na.locf(arfit$resid))
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
  # optional plot
  if (plot) plot(stats::ts.union(tser, tser.prew), yaxs="i", xaxs="i", yax.flip=TRUE)
  return(invisible(tser.prew))
}
