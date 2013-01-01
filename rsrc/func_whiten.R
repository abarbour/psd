#' prewhiten a timeseries object
#' 
#' de-mean, de-trend (which also de-means), and (soon) fit an AR
#' model to the series
#'
#' @note the \code{AR.fit} option is not used (yet)
#'
#' @param tser  \code{ts} object
#' @param AR.fit logical; FALSE
#' @param detrend  logical; TRUE
#' @param demean  logical; TRUE
#' @param plot  logical; TRUE
#' @param verbose  logical; TRUE
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

prewhiten2 <- function(...) UseMethod(".whiten")
.whiten.default <- function(x) {
  ###
  # PORT of RLP's riedsid.m
  # abarbour
  # Dec 2011
  #
  # porting: Jan 3, 2012
  # testing: 
  #     Jan 3, 2012: 
  #        appears to work fine with 'filter' but not as well with 'open'
  #        noise<-rnorm(n,sd=0.2)
  #        signal<-.001*(1:n) + triang(n)+3*sin(pi*1:n/180+pi/4)+3*sin(2*pi*1:n/180+pi/4)
  ###
  #
  #  Applies a short convolution filter to time series x that
  #  flattens its spectrum (making it more nearly white).
  #  The resultant approximately white-noise time series makes
  #  nonstationarity and error spikes much more prominent.  
  #  The new series can be used for preliminary data screening.
  # 
  #  Method:  solves Yule-Walker equations using autocorrelation
  #  functions based on a sine-multitaper spectrum.
  #
  ##
  ## Args:	
  ##
  ## Returns:	
  ##
  ## TODO(abarbour):
  ##
  ## Prelims
  nord <- 4
  ntap <- 10
  nx <- length(x)
  if (nx < 10){ 
    print('Series too short to whiten') 
    return
  }
  ##
  #  Calculate psd with multitapers (fixed num tapers)
  # FIX function pointer:
  spec <- .psdcore.default(x, ntap, plotpsd=FALSE, as.spec=FALSE)$spec
  nf <- length(spec)        
  #  Fourier transform the psd to get the autocovariance function
  autoc <- Re(fft( c(spec[1]/2, spec[2:nf-1], matrix(0,nf-1,1)) ))
  autoc.n <- autoc / autoc[1]
  #  Form the Yule-Walker equations and solve them
  Toe <- toeplitz(autoc.n[1:nord])
  coeffs <- solve(Toe, autoc.n[2:(nord+1)])
  #  Convolve filter over original time series
  # 'filter' == 'same' (matlab), confirm [ ]
  whitex <- convolve(x-mean(x), c(-1,coeffs), type="filter")
  ##
  return(invisible(whitex))
}
# end whiten.default
# plot.qualcon or plot.whiten from show.whiten in suppfuncs [ ]
show.whiten <- function(wcrit, k1k2, x, x.w, limsc=2.5, pltlabs=TRUE) {
  ##  Display the series and its filtered version, interactively
  ##
  ## Args:  
  ##
  ## Returns:  
  ##
  ## TODO(abarbour):	
  ##     a) interative outlier identification or incorporate other package?
  ##     b) make this a plot method
  ##
  t1 <- k1k2[1]
  te <- k1k2[length(k1k2)]
  xplt <- x[k1k2]
  xplt <- xplt - median(xplt)
  xplt.w <- x.w[k1k2]
  x.wcriti <- which(abs(xplt.w) > wcrit) + t1 - 1
  x.wcrit <- x.w[x.wcriti]
  #
  plot(k1k2, xplt, type="s", col="gray", 
       ylim=c(-1*limsc*wcrit, limsc*wcrit),
       ylab="value",
       xlab="term",
       main="Outlier Inspector")
  
  abline(h=wcrit*c(1,-1), lty=3, col="black")
  abline(h=wcrit*c(2,-2), lty=3, col="red")
  if (pltlabs==TRUE){
    text(t1+25, -0.35*wcrit, "demeaned F(Xa:Xb)", pos=1, cex=0.8, col="dark gray")
    text(t1+17, -0.60*wcrit, "AR prew. F(X)", pos=1, cex=0.8, col="blue")
    text(t1+20, -1*wcrit, sprintf("threshold (%.3f)\ncrossings (+)",wcrit), pos=1, cex=0.8, col="black")
    text(t1+10, -2*wcrit, "2x thresh", pos=3, cex=0.8, col="red")
  }
  if (length(x.wcriti)>0){
    points(x.wcriti, x.wcrit, col="black", pch=3, cex=0.8)
  }
  lines(k1k2, xplt.w, type="s", col="blue")
  # display some info
  #   terms
  cat('Display of terms: ', t1, te, "\n")
  #   outliers
  if (length(x.wcriti>0)){
    cat("Threshold crossing indices:\n")
    toret <- t(rbind(x.wcriti,x.wcrit))
    print(toret)
    return(invisible(toret))
  } else {
    cat(sprintf("No critical outliers found for threshold of  %.03f\n",wcrit))
  }
} 
# end show.whiten
#
#