## CHANGES WERE MADE TO FUNC_PSDCORE: reconcile changes here
##
##  Default method for psdcore, which does the grunt work
##
..dev_psdcore.default <-function(X.d, 
                                 X.frq=1, 
                                 ntaper=1, 
                                 ndecimate=1L,
                                 demean=TRUE, 
                                 detrend=TRUE,
                                 na.action = stats::na.fail,
                                 plotpsd=FALSE, 
                                 #plotcolor="#000000", xlims=c(0,0.5),
                                 as.spec=FALSE,
                                 ...
                                 ) {
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
  stopifnot(is.vector(X.d))
  series <- deparse(substitute(X.d))
  X.d <- na.action(stats::ts(X.d, frequency=X.frq))
  Nyq <- stats::frequency(X.d)/2
  #X.d <- as.matrix(X.d) # column mat
  ##
  ###  When ntaper is a scalar, initialize
  ##
  lt <- length(ntaper)
  if (lt == 1){
    #
    # 
    # original series
    n.o <- envAssignGet("len_orig", length(X.d))
    #
    X <- prewhiten(X.d, 
                   AR.fit=FALSE, 
                   detrend=detrend, 
                   demean=demean,
                   plot=FALSE,
                   verbose=FALSE)
    #
    # Force series to be even in length (modulo division)
    n.e <- envAssignGet("len_even", n.o - n.o%%2 )
    x_even <- as.matrix(X[1:n.e])
    envAssign("ser_orig", X)
    envAssign("ser_orig_even", x_even)
    # half length of even series
    nhalf <- envAssignGet("len_even_half", n.e/2)
    # variance of even series
    varx <- envAssignGet("ser_even_var", drop(stats::var(x_even)))
    # create uniform tapers
    ntap <- ntaper*ones(nhalf+1)
    ##  Remove mean & pad with zeros
    X.dem <- matrix(c(x_even, zeros(n.e)), ncol=1)
    ##  Take double-length fft
    # mvfft takes matrix (allos multicolumn)
    fftz <- envAssignGet("fft_even_demeaned_padded", stats::mvfft(X.dem))
  } else {
    #X <- X.d
    ntap <- ntaper
    n.e <- envGet("len_even")
    nhalf <- envGet("len_even_half")
    varx <- envGet("ser_even_var")
    fftz <- envGet("fft_even_demeaned_padded")
  }
  # if ntaper is a vector, this doesn't work [ ]
  ##
  if (!(is.taper(ntap))) ntap <- as.taper(ntap)
  ##
  ###  Select frequencies for PSD evaluation
  if  (lt > 1 && ndecimate > 1){
    stopifnot(!is.integer(ndecimate))
    message("decim stage 1")
    # interp1 requires strict monotonicity (for solution stability)
    nsum <- base::cumsum(1/ntap)
    ns1 <- nsum[1]
    tmp.x <- nhalf * (nsum - ns1) / (nsum[length(nsum)] - ns1)
    tmp.y <- seq.int(0, nhalf, by=1)
    tmp.xi <- seq.int(0, nhalf, by=ndecimate)
    # linear interplate (x,y) to (xi,yi) where xi is decimated sequence
    tmp.yi <- signal::interp1(tmp.x, tmp.y, tmp.xi, method='linear', extrap=TRUE)
    f <- c(round(tmp.yi), nhalf)
    #iuniq <- 1 + which(diff(f) > 0)
    f <- unique(c(0, f))   #  Remove repeat frequencies in the list
  } else {
    f <- seq.int(0, nhalf, by=1)
  }
  ##
  lt2 <- length(drop(ntaper))
  ##
  ###  Calculate the psd by averaging over tapered estimates
  nfreq <- length(f)
  ##
  ###  Loop over frequency
  ## change to lapply? []
  if (sum(ntaper) > 0) {
    psd <- zeros(nfreq)
    xfft <- Re(fftz)
    for ( j in 1:nfreq ) {
      #  Sum over taper indexes weighting tapers parabolically
      fj <- f[j]
      fj2 <- fj*2
      n2.e <- 2*n.e
      # parabolic weights, index m+1, column vec out
      KPW <<- parabolic_weights(ntap, tap.index=(fj+1), vec.out="horizontal")
      # this is a distinction with order of operations and %% (2.14.0)
      # https://bugs.r-project.org/bugzilla3/show_bug.cgi?id=14771
      # (but correct in matlab's function call) 
      # So: enclose in parens or use rlpSpec::mod.default
      m1. <- fj2 + n2.e - KPW$taper_seq
      j1 <- m1. %% n2.e
      m1. <- fj2 + KPW$taper_seq
      j2 <- m1. %% n2.e
      f1 <- xfft[j1+1]
      f2 <- xfft[j2+1]
      af12. <<- abs( f1 - f2 )
      af122. <- af12. * af12.
      #[ ] Error in KPW$taper_weights %*% af122. : 
      #requires numeric/complex matrix/vector arguments
      psdv <- KPW$taper_weights %*% as.colvec(af122.)
      psd[fj] <- drop(psdv)
      
    }
  } else {
    message("zero taper result == raw periodogram")
    xfft <- envGet("fft_even_demeaned_padded")
    ff <- xfft[1:nfreq]
    N0 <- envGet("len_orig")
    psd <- ff * Conj(ff) / N0
  }
  ##  Interpolate if necessary to uniform freq sampling
  if (lt > 1 && ndecimate > 1){
    ## check [ ]
    message("decim stage 2")
    tmp.x <- f
    tmp.xi <- tmp.y
    tmp.y <- psd
    tmp.yi <- signal::interp1(tmp.x, tmp.y, tmp.xi, method='linear', extrap=TRUE)
    psd <- tmp.yi
  }
  ##
  psd <- as.matrix(drop(Re(psd)))
  ## Normalize by variance
  #message("psd norm")
  trap.area <- (2*sum(psd) - psd[1] - psd[nrow(psd)])/2 # Trapezoidal rule
  bandwidth <- 1 / nhalf
  psd.norm <- drop(2 * trap.area / varx * bandwidth) #2*
  #there was an apparently incorrect factor of 2 here <-- probably not
  # see notes/normalization.tx
  # dB_power = 10*log10(P1/P2)
  # 1 == 0 dB
  # 2 == 3 dB
  #plot(psd,type="s")
  psd.n <- drop(psd / psd.norm)
  frq <- seq.int(0, 0.5, length.out=nfreq)
  # there seems to be an issue with f==0, so just extrapolate from the prev point
  psd.n <- as.matrix(
    exp(as.numeric(signal::interp1(frq[2:(nfreq-2)], 
                        log(psd.n[2:(nfreq-2)]), 
                        frq, 
                        method='linear', extrap=TRUE))))
  frq <- as.matrix(frq)
  
  pltpsd <- function(...){
    Xpg<-spec.pgram(X, log="no", pad=1, taper=0, detrend=F, demean=F, plot=F)
    opar <- par(no.readonly = TRUE)
    par(mar=rep(2,4), oma=rep(0,4))
    layout(matrix(c(1,2),ncol=1),c(1,2))
    plot(log(Xpg$freq), log10(Re(Xpg$spec)), col="red", type="l", main="spectra") #,ylim=5*c(-.5,1.5))
    lines(as.vector(log(frq)), as.vector(log10(psd.n)), type="l")
    legend("topright",c("spec.pgram","rlpSpec"),col=c("red","black"),lty=1)
    plot(X,type="l", main="spec series (de-trend/-mean if applied)")
    par(opar)
  }
  if (plotpsd) pltpsd(...)
  #plot(10*log10(psd),type="s")
  # 
  psd.out <- list(freq = frq, 
                  spec = psd.n, 
                  coh = NULL, 
                  phase = NULL, 
                  kernel = NA, 
                  df = NA, 
                  numfreq = nfreq,
                  bandwidth = bandwidth, 
                  n.used = envGet("len_even"), 
                  orig.n = envGet("len_orig"), 
                  series = series, 
                  snames = colnames(X), 
                  method = "Adaptive Sine Multitaper (rlpSpec)", 
                  taper = ntap, 
                  pad = TRUE, 
                  detrend = detrend, 
                  demean = demean)
  if (as.spec){class(psd.out) <- "spec"}
  return(invisible(psd.out))
}
