
##
##  Default method for psdcore, which does the grunt work
##
.dev_psdcore.default <-function(X,  
                           ntaper=1, ndecimate=1, 
                           plotpsd=TRUE, plotcolor="#000000",
                           xlims=c(0,0.5)) {
  stopifnot(is.vector(X))
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
    n.o <- envAssignGet("len_orig", length(X))
    # Force series to be even in length (modulo division)
    n.e <<- envAssignGet("len_even", n.o - n.o%%2 )
    x_even <<- X[1:n.e]
    envAssign("ser_orig", X)
    envAssign("ser_orig_even", x_even)
    # half length of even series
    nhalf <- envAssignGet("len_even_half", n.e/2)
    # variance of even series
    envAssign("ser_even_var", stats::var(x_even))
    varx <- envGet("ser_even_var")
    # create uniform tapers
    ntap <<- colvec(nrow=nhalf+1, val=ntaper)
    ##  Remove mean & pad with zeros
    # convert to sweep [ ]
    ##z <<- rbind( matrix(x_even, byrow=T) - mean(x_even), zeros(n.e) )
    z <<- matrix(c(x_even - base::mean(x_even), zeros(n.e)),ncol=1)
    ##  Take double-length fft
    # mvfft takes multicol matrix
    fftz <<- envAssignGet("fft_even_demeaned_padded", Re(stats::fft(z)))
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
    nsum <- base::cumsum(1/ntap)
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
    fj <<- f[j]
    tot.tapers <<- ntap[fj+1]
    if (tot.tapers <= 0){
      k. <<- 1
    } else {
      k. <<- seq.int(1, tot.tapers, by=1)
    }
    #  Sum over taper indexes weighting tapers parabolically
    W. <<- matrix(
      (tot.tapers*tot.tapers - (k.-1)*(k.-1))*3/2/tot.tapers/(tot.tapers-0.25)/(tot.tapers+1), 
      nrow=1)
    # this is a distinction with order of operations and %% (2.14.0)
    # https://bugs.r-project.org/bugzilla3/show_bug.cgi?id=14771
    # (but correct in matlab's function call) 
    # So: enclose in parens or use rlpSpec::mod.default
    m1. <<- 2*fj + 2*n.e - k.
    m2. <<- 2*n.e
    j1 <<- m1. %% m2.
    m1. <- 2*fj + k.
    j2 <<- m1. %% m2.
    f1 <<- fftz[j1+1]
    f2 <<- fftz[j2+1]
    af12. <<- abs( f1 - f2 )
    af122. <<- af12. * af12.
    psdv <<- W. %*% af122.
    psd[j,1] <- psdv
    
  }
  psd_F <<- psd
  ##  Interpolate if necessary to uniform freq sampling
  if (length(ntaper) > 1 && ndecimate > 1){
    ## check [ ]
    tmp.x <- f
    tmp.xi <- tmp.y
    tmp.y <- psd
    #psd_I <<- psd
    tmp.yi <- signal::interp1(tmp.x, tmp.y, tmp.xi, method="cubic") #, method='linear', extrap=TRUE)
    psd <- tmp.yi
    #psd_F <<- psd
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
      #plot(ftoplot, psdtoplot, type="s",
      plot(frq, 10*log10(psd), type="s",
           main="Adaptive Sine-multitaper PSD Estimation",
           xlab="Nyquist frequency", xlim=xlims,
           ylab="PSD, dB rel. Nyquist")
    } else {
      # so adaptive estimation may be visualized
      ##lines(ftoplot, psdtoplot, type="s", col=plotcolor)
    }
  }
  return(invisible(psd))
}
