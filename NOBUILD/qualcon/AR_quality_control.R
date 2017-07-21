###
# PORT of RLP's qualcon.m/whiten.m pair
###
  
qualcon <-function(x, ...) UseMethod('qualcon')
qualcon.default <-function(x, crit.var=2, pm=100, sc=2.8, ...) {             
  ###
  #
  #  Searches for data errors indicated by the whitened series
  #  falling outside fixed limits.
  #  Plots +-100 points of x and whitened x around error section.
  #  badpoints <- number error regions detected
  #
  ###
  x.w <- ar_whiten(x)
  message("Critical variance multiplier: ", crit.var)
  wcrit <- crit.var * sqrt(var(x.w))
  par(ask=TRUE)
  
  lx <- length(x)
  toret <- show.whiten(wcrit, seq_len(lx), x, x.w, pltlabs=FALSE, ...)
  
  Icrit <- which(abs(x.w) > wcrit)
  
  ncrit <- length(Icrit)
  
  if (ncrit > 0){
	crits <- vector("list", ncrit)
    i2 <- 0
    badseg <- 0
    for (j  in  seq_len(ncrit)) {
      if (Icrit[j] >= i2){
      	badseg <- badseg + 1
        i1 <- Icrit[j] - pm
        i2 <- Icrit[j] + pm
        k1 <- max(1, i1)  
        k2 <- min(i2, lx)
        k1k2 <- seq.int(k1, k2)
        crits[[badseg]] <- show.whiten(wcrit, k1k2, x, x.w, limsc=sc, lbl=paste0("seg-",badseg))
      }
    }
  } else {
  	crits <- NA
  }
  toret <- Filter(Negate(is.null), crits)
  return(toret)
} 


ar_whiten <- function(x, ...) UseMethod("ar_whiten")
ar_whiten.default <- function(x, nord = 4, ntap = 10, ...) {
  ##
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
  library(psd)
  #
  nx <- length(x)
  if (nx < ntap){ 
    stop('Series too short (', nx, ' to whiten based on ', ntap, " tapers") 
  }
  ##
  #  Calculate psd with multitapers (fixed num tapers)
  p <- psdcore(x, ntapers=ntap, plot=FALSE) 
  spec <- p[['spec']]
  nf <- length(spec)        
  #  Fourier transform the psd to get the autocovariance function
  autoc <- Re(fft(c(spec[1]/2, spec[2:(nf-1)], matrix(0, (nf-1), 1))))
  autoc.n <- autoc / autoc[1]
  #  Form the Yule-Walker equations, and invert
  Toe <- toeplitz(autoc.n[1:nord])
  coeffs <- solve(Toe, autoc.n[2:(nord+1)])
  #  Convolve filter over original time series
  whitex <- convolve(x-mean(x), c(-1, coeffs), type="filter")
  ##
  return(invisible(whitex))
}

show.whiten <- function(wcrit, k1k2, x, x.w, limsc=1, lbl="-x-", pltlabs=TRUE) {
  ##
  #  Display the series and its filtered version, interactively
  ##
  ## TODO(abarbour):	
  ##     a) interative outlier identification or incorporate other package?
  ##     b) make this a plot method
  ##
  t1 <- k1k2[1]
  t1e <- k1k2[length(k1k2)]
  xplt <- x[k1k2]
  xplt <- xplt - median(xplt)
  xplt.w <- x.w[k1k2]
  x.wcriti <- which(abs(xplt.w) > wcrit) + t1 - 1
  x.wcrit <- x.w[x.wcriti]
  #
  plot(k1k2, xplt, col=NA,
  	xaxs='i', ylim=2 * limsc * wcrit * c(-1,1),
  	ylab="value", xlab="term", main="Outlier Inspector")
  abline(h=0)
  abline(h=wcrit*c(1,-1), lty=3, col="black")
  abline(h=wcrit*c(2,-2), lty=3, col="red")
  if (pltlabs==TRUE){
    text(t1, 0.20*c(1)*wcrit, "demeaned F(Xa:Xb)", pos=4, cex=0.8, col="grey40")
    text(t1, 0.40*c(1)*wcrit, "AR prew. F(X)", pos=4, cex=0.8, col="blue")
    text(t1, 0.9*c(-1,1)*wcrit, sprintf("threshold (%.3f) crossings (+)",wcrit), pos=4, cex=0.8, col="black")
    text(t1, 1.9*c(-1,1)*wcrit, "2x thresh", pos=4, cex=0.8, col="red")
  }
  lines(k1k2, xplt, type="l", col="grey")
  if (length(x.wcriti)>0){
    points(x.wcriti, x.wcrit, col="black", pch=3, cex=0.8)
  }
  lines(k1k2, xplt.w, type="l", col="blue")

  # display some info
  #   terms
  message('Display of terms: ', t1, " -- ", t1e)
  #   outliers
  toret <- if (length(x.wcriti>0)){
    message(length(x.wcriti), " threshold crossings")
    data.frame(Lbl=lbl, Ind=x.wcriti, Val=x.wcrit)
  } else {
    message(sprintf("No critical outliers found for threshold of  %.03f", wcrit))
    NA
  }
  return(toret)
} 


EXAMP <- function(...){
	try(dev.off())
	
	x <- rnorm(4e3)
	badinds <- c(500:520, 1500:1520)
	nbad <- length(badinds)
	x[badinds] <- runif(nbad, max=20) * rt(nbad, 2, 3)

	qualcon(x, ...)
}

run.examp <- FALSE
if (run.examp){
	set.seed(12444)
	print(EXAMP(sc=1.0, crit.var=2))
}

