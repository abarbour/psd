##
##  Default method for psdcore, which does the grunt work
##
.psdcore.default <-function(x,  
                           ntaper=1, ndecimate=1, 
                           plotpsd=TRUE, plotcolor="#000000",
                           xlims=c(0,0.5)) {
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
    # original series
    envAssign("len_orig", length(x))
    n.o <- envGet("len_orig")
    # Force series to be even in length (modulo division)
    envAssign("len_even", n.o - n.o%%2 )
    n.e <- envGet("len_even", psdenv)
    x_even <- x[1:n.e]
    envAssign("ser_orig", x)
    envAssign("ser_orig_even", x_even)
    # half length of even series
    envAssign("len_even_half", n.e/2)
    nhalf <- envGet("len_even_half")
    # variance of even series
    envAssign("ser_even_var", var(x_even))
    varx <- envGet("ser_even_var")
    # create uniform tapers
    ntap <- ones(nhalf+1)*ntaper
    ##  Remove mean & pad with zeros
    # convert to sweep [ ]
    z <- rbind( matrix(x_even, byrow=T) - mean(x_even), zeros(n.e) )
    ##  Take double-length fft
    envAssign("fft_even_demeaned_padded", Re(fft(z)))
    fftz <- envGet("fft_even_demeaned_padded")
  } else {
    ntap <- ntaper
    n.e <- envGet("len_even")
    nhalf <- envGet("len_even_half")
    varx <- envGet("ser_even_var")
    fftz <- envGet("fft_even_demeaned_padded")
  }
  ###  Select frequencies for PSD evaluation
  if  (lt > 1 && ndecimate > 1){
    # interp1 requires strict monotonicity (for solution stability)
    nsum <- cumsum(1/ntap)
    ns1 <- nsum[1]
    tmp.x <- nhalf * (nsum - ns1) / (nsum[length(nsum)] - ns1)
    tmp.y <- seq.int(0, nhalf, by=1)
    tmp.xi <- seq.int(0, nhalf, by=ndecimate)
    # linear interplate (x,y) to (xi,yi) where xi is decimated sequence
    tmp.yi <- signal::interp1(tmp.x, tmp.y, tmp.xi, method='linear', extrap=TRUE)
    f <- c(round(tmp.yi), nhalf)
    iuniq <- 1 + which(diff(f) > 0)
    f <- c(0, f[iuniq])   #  Remove repeat frequencies in the list
  } else {
    f <- seq.int(0, nhalf, by=1)
  }
  ###  Calculate the psd by averaging over tapered estimates
  nfreq <- length(f)
  #print(nfreq)
  psd <- zeros(nfreq)
  ###  Loop over frequency
  for ( j in 1:nfreq ) {
     m <- f[j]
     tapers <- ntap[m+1]
     #  Sum over taper indexes weighting tapers parabolically
     if (tapers <= 0){
         k <- 1
       } else {
         k <- seq.int(1, tapers, by=1)
       }
     w <- rbind((tapers^2 - (k-1)^2) * (1.5/(tapers*(tapers-0.25)*(tapers+1))))
     # this is a distinction with order of operations and %% (2.14.0)
     # https://bugs.r-project.org/bugzilla3/show_bug.cgi?id=14771
     # (but correct in matlab's function call) 
     # So: enclose in parens or use rlpSpec::mod.default
     j1 <- mod((2*m + 2*n.e - k), 2*n.e)
     j2 <- mod((2*m + k), 2*n.e)
     f1 <- fftz[j1+1]
     f2 <- fftz[j2+1]
     psdv <- w %*% abs( f1 - f2 )^2
     psd[j,1] <- psdv
  }
  ##  Interpolate if necessary to uniform freq sampling
  if (length(ntaper) > 1 && ndecimate > 1){
    ## check [ ]
    tmp.x <- f
    tmp.xi <- tmp.y
    tmp.y <- psd
    tmp.yi <- signal::interp1(tmp.x, tmp.y, tmp.xi, method='linear', extrap=TRUE)
    psd <- tmp.yi
  }
  ## Normalize by variance
  area <- (sum(psd) - psd[1]/2 - psd[length(psd)]/2)/nhalf  # 2*Trapezoid
  psd <- as.matrix((1*varx/area)*psd) #there was an apparently incorrect factor of 2 here
  psdtoplot <- 10*log10(psd[2:(nfreq-1)]) ## R uses 10* to scale to dB (why?)
  frq <- seq.int(0, 0.5, length.out=nfreq)
  ftoplot <- frq[2:(nfreq-1)]
  ## Plot if desired
  if (plotpsd) {
    if (plotcolor=="#000000" || plotcolor==0){
      # initial plot (black)
      plot(ftoplot, psdtoplot, type="s",
           main="Adaptive Sine-multitaper PSD Estimation",
           xlab="Nyquist frequency", xlim=xlims,
           ylab="PSD, dB rel. Nyquist")
    } else {
      # so adaptive estimation may be visualized
      lines(ftoplot, psdtoplot, type="s", col=plotcolor)
    }
  }
  return(invisible(psd))
} 
# end psdcore
###
