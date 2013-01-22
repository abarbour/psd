#[X] deprecate

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
    # original series and length
    envAssign("ser_orig", x)
    n.o <- envAssignGet("len_orig", length(x))
    # Force series to be even in length (modulo division)
    n.e <- envAssignGet("len_even", n.o - n.o%%2)
    # even-length series
    x_even <- envAssignGet("ser_orig_even", x[1:n.e])
    # half length of even series
    nhalf <- envAssignGet("len_even_half", n.e/2)
    # variance of even series
    varx <- envAssignGet("ser_even_var", var(x_even))
    # create uniform tapers
    ntap <- ones(nhalf+1)*ntaper
    ##  Remove mean & pad with zeros
    tmpx <- matrix(x_even, byrow=T)
    z <- rbind(sweep(tmpx, MARGIN=2, STATS=colMeans(tmpx), FUN="-", check.margin = FALSE), zeros(n.e))
    rm(tmpx)
    ##  Take double-length fft
    fftz <- envAssignGet("fft_even_demeaned_padded", Re(fft(z)))
  } else {
    ntap <- ntaper
    n.e <- envGet("len_even")
    nhalf <- envGet("len_even_half")
    varx <- envGet("ser_even_var")
    fftz <- envGet("fft_even_demeaned_padded")
  }
  if (!(is.taper(ntap))) ntap <- as.taper(ntap)
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
    #Sum over taper indexes weighting tapers parabolically
    m <- f[j]
    m2 <- m*2
    # parabolic weights, index m+1, column
    kW. <- parabolic_weights(ntap, tap.index=(m+1), vec.out="horizontal")
    k. <- seq.int(1,ntap(m+1),by=1)
    # there is a distinction with order of operations and %% (2.14.0)
    # https://bugs.r-project.org/bugzilla3/show_bug.cgi?id=14771
    # (but correct in matlab's function call) 
    # So: enclose in parens or use rlpSpec::mod.default
    n2.e <- 2*n.e
    j1 <- mod((m2 + n2.e - k.), n2.e)
    j2 <- mod((m2 + k.), n2.e)
    f1 <- fftz[j1+1]
    f2 <- fftz[j2+1]
    psdv <- kW. %*% abs( f1 - f2 )^2
    psd[j,1] <- psdv
  }
  ##  Interpolate if necessary to uniform freq sampling
  if (lt > 1 && ndecimate > 1){
    ## check [ ]
    tmp.x <- f
    tmp.xi <- tmp.y
    tmp.y <- psd
    tmp.yi <- signal::interp1(tmp.x, tmp.y, tmp.xi, method='linear', extrap=TRUE)
    psd <- tmp.yi
  }
  ## Normalize by variance
  area <- (sum(psd) - psd[1]/2 - psd[length(psd)]/2)/nhalf  # 2*Trapezoid
  psd <- as.matrix((2*varx/area)*psd) 
  #there was an apparently incorrect factor of 2 here <-- probably not
  # see notes/normalization.txt
  # dB_power = 10*log10(P1/P2)
  # 1 == 0 dB
  # 2 == 3 dB
  psdtoplot <- 10*log10(psd[2:(nfreq-1)]) ## R uses 10* to scale to power-dB
  frq <- seq.int(0, 0.5, length.out=nfreq)
  ftoplot <- frq[2:(nfreq-1)] # should fix this [ ]
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
