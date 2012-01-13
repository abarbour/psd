psdcore <-function(x,  ntaper=1, ndecimate=1, plotpsd=TRUE, plotcolor="#000000") {
  ###
  # PORT of RLP's psdcore.m
  # abarbour
  # Dec 2011
  #
  # porting:  Jan 3, 2012 (partial)
  # testing:
  #           Jan 3, 2012: 
  #               have not checked decimation or interpolation interp1
  #               and persistent variables seem not to be an issue
  ###
  #
  #  Compute a spectral estimate of the power spectral density
  #  (psd) for the time series x using sine multitapers.
  #  Normalised to sampling interavl of 1.
  #  
  #  ntaper gives the number of tapers to be used at each frequency:
  #  if ntaper is a scalar, use same value at all freqs if a
  #  vector, use ntaper(j) sine tapers at frequency j. 
  #  If series length is n, psd is found at 1 + n/2 evenly spaced freqs
  #  if n is odd, x is truncted by 1.
  #  ndecimate: number of psds actually computed <- (1+n/2)/ndecimate
  #  these values are linearly interpolated into psd.
  #
  ##
  ## Args:	
  ##
  ## Returns:	
  ##
  ## TODO(abarbour):	
  ##
  require(signal, quietly=TRUE, warn.conflicts=FALSE)
  # for interp1
  #   require(clim.pact, quietly=T)
  #   for mod (just snaked and put in funcs.R)
  #
  # persistent [ ]
  #persistent fftz n nhalf varx
  #
  ###  When ntaper is a scalar, initialize
  lt <- length(ntaper)
  if (lt == 1){
    n.o <- length(x)
    n <<- n.o - n.o%%2    # Force series to be even in length (modulo division)
    #     n <- nextn(n.o)
    nhalf <<- n/2
    varx <<- var(x[1:n]) 
    ntap <- matrix(1, nhalf+1, 1)*ntaper  # Make a vector from scalar value
    #  Remove mean & pad with zeros
    z <- rbind( matrix(x[1:n],byrow=T) - mean(x[1:n]), matrix(0,n,1))
    #  Take double-length fft
    fftz <<- Re(fft(z))
  } else {
    ntap <- ntaper
  }
  ###  Select frequencies for PSD evaluation
  if  (lt > 1 && ndecimate > 1){
    # interp1 requires strict monotonicity (for solution stability)
    nsum <- cumsum(1/ntap)
    ns1<-nsum[1]
    tmp.x <- nhalf * (nsum - ns1) / (nsum[length(nsum)] - ns1)
    tmp.y <- seq(0, nhalf, by=1)
    tmp.xi <- seq(0, nhalf, by=ndecimate)
    # linear interplate (x,y) to (xi,yi)
    tmp.yi <- interp1(tmp.x, tmp.y, tmp.xi, method='linear', extrap=TRUE)
    f <- c(round(tmp.yi), nhalf)
    iuniq <- 1 + which(diff(f) > 0)
    f <- c(0, f[iuniq])   #  Remove repeat frequencies in the list
  } else {
    f <- seq(0, nhalf, by=1)
  }
  ###  Calculate the psd by averaging over tapered estimates
  nfreq <- length(f)
  psd <- matrix(0, nfreq, 1)
  ###  Loop over frequency
  for ( j in 1:nfreq ) {
     m <- f[j]
     tapers <- ntap[m+1]
     #  Sum over taper indexes weighting tapers parabolically
     if (tapers <= 0){
         k <- 1
       } else {
         k <- seq(1, tapers, by=1)
       }
     w <- rbind((tapers^2 - (k-1)^2) * (1.5/(tapers*(tapers-0.25)*(tapers+1))))
     # this is a distinction with order of operations and %% (2.14.0)
     # https://bugs.r-project.org/bugzilla3/show_bug.cgi?id=14771
     # (but correct in matlab's function call) 
     # So: enclose in parens or use mod in funcs2.R
     j1 <- mod((2*m + 2*n - k), 2*n)
     j2 <- mod((2*m + k), 2*n)
     f1 <- fftz[j1+1]
     f2 <- fftz[j2+1]
     psdv <- w %*% abs( f1 - f2 )^2
     psd[j,1] <- psdv
  }
  ##  Interpolate if necessary to uniform freq sampling
  if (length(ntaper) > 1 && ndecimate > 1){
    ## check [ ]
    tmp.x <- f
    tmp.y <- psd
    tmp.xi <- tmp.y
    tmp.yi <- interp1(tmp.x, tmp.y, tmp.xi, method='linear', extrap=TRUE)
    psd <- tmp.yi
  }
  ## Normalize by variance
  area <- (sum(psd) - psd[1]/2 - psd[length(psd)]/2)/nhalf  # 2*Trapezoid
  psd <- as.matrix((2*varx/area)*psd)
  psdtoplot <- 20*log10(psd[2:(nfreq-1)])
  frq <- seq(0, 0.5, length.out=nfreq)
  ftoplot <- frq[2:(nfreq-1)]
  ## Plot if desired
  if (plotpsd) {
    if (plotcolor=="#000000" || plotcolor==0){
      # initial plot (black)
      plot(ftoplot, psdtoplot, type="s",
           main="Adaptive Sine-multitaper PSD Estimation",
           xlab="Nyquist frequency",
           ylab="PSD, dB rel. Nyquist")
    } else {
      # so adaptive estimation may be visualized
      lines(ftoplot, psdtoplot, type="s", col=plotcolor)
    }
  }
  return(invisible(psd))
} 
# end psdcore
