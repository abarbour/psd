##
##  Default method for psdcore, which does the grunt work
##
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
psdcore <- function(...) UseMethod(".psdcore")
#.devpsdcore <- function(x, ...) UseMethod("..dev_psdcore")
.psdcore.default <-function(X.d, 
                           X.frq=1, 
                           ntaper=as.taper(1), 
                           ndecimate=1L,
                           demean=TRUE, 
                           detrend=TRUE,
                           na.action = stats::na.fail,
                           first.last=TRUE,
                           Nyquist.normalize=TRUE,
                           plotpsd=FALSE,
                           as.spec=TRUE,
                           force_calc=FALSE,
                           ...
                           ) {
  #
  require(signal, quietly=TRUE, warn.conflicts=FALSE)
  # for interp1
  #   require(clim.pact, quietly=T)
  #   for mod (just snaked and put in funcs.R)
  #
  series <- deparse(substitute(X.d))
  X.d <- na.action(stats::ts(X.d, frequency=X.frq))
  Nyq <- stats::frequency(X.d)/2
  ##
  ###  When ntaper is a scalar, initialize
  ##
  lt <- length(ntaper)
  if ((lt == 1) | (force_calc)){
    #
    # 
    # original series
    n.o <- envAssignGet("len_orig", length(X.d))
    #
    X <- prewhiten(X.d, 
                   AR.max=0L, 
                   detrend=detrend, 
                   demean=demean,
                   plot=FALSE,
                   verbose=FALSE)
    #
    # Force series to be even in length (modulo division)
    n.e <- envAssignGet("len_even", n.o - n.o%%2 )
    X.even <- as.matrix(X[1:n.e])
    envAssign("ser_orig", X)
    envAssign("ser_orig_even", X.even)
    # half length of even series
    nhalf <- envAssignGet("len_even_half", n.e/2)
    # variance of even series
    varx <- envAssignGet("ser_even_var", drop(stats::var(X.even)))
    # create uniform tapers
    nt <- nhalf + 1
    if (lt < nt) {
      ntap <- ntaper*ones(nt) 
    } else {
      ntap <- ntaper[1:nt]
    }
    ##  Remove mean & pad with zeros
    X.dem <- matrix(c(X.even, zeros(n.e)), ncol=1)
    ##  Take double-length fft
    # mvfft takes matrix (allos multicolumn)
    fftz <- envAssignGet("fft_even_demeaned_padded", stats::mvfft(X.dem))
  } else {
    X <- X.d
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
  if (sum(ntaper) > 0) {
    psd <- zeros(nfreq)
    n2e <- 2*n.e
    Rfftz <- Re(fftz)
    PSDFUN <- function(fj, n2.e=n2e, ntaps=ntap, Xfft=Rfftz){
      # parabolic weights, index m+1, column vec out
      #print(c(fj,fj2,n2.e))
      NT <- ntaps[fj+1]
      #KPW <- parabolic_weights(ntaps, tap.index=(fj+1), vec.out="horizontal")
      KPW <- parabolic_weights_fast(NT)
      # this is a distinction with order of operations and %% (2.14.0)
      # https://bugs.r-project.org/bugzilla3/show_bug.cgi?id=14771
      # (but correct in matlab's function call) 
      # So: enclose in parens or use rlpSpec::mod.default
      fj2 <- 2*fj
      m1. <- fj2 + n2.e - KPW$taper_seq
      j1 <- m1. %% n2.e
      m1. <- fj2 + KPW$taper_seq
      j2 <- m1. %% n2.e
      f1 <- Xfft[j1+1]
      f2 <- Xfft[j2+1]
      af12. <- f1 - f2
      af122. <- af12. * af12.
      psdv <- KPW$taper_weights %*% matrix(af122., ncol=1)
      return(psdv)
    }
    #  easier to follow, but foreach solution is actually slower :(
    #require(foreach)
    #psd <- foreach::foreach(f.j=f[1:nfreq], .combine="c") %do% PSDFUN(fj=f.j)
    # vapply is much faster than even lapply
    psd <- vapply(X=f[1:nfreq], FUN=PSDFUN, FUN.VALUE=1.0)
  } else {
    message("zero taper result == raw periodogram")
    Xfft <- envGet("fft_even_demeaned_padded")
    ff <- Xfft[1:nfreq]
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
    # x y x_interp --> y_interp
    tmp.yi <- signal::interp1(tmp.x, tmp.y, tmp.xi, method='linear', extrap=TRUE)
    psd <- tmp.yi
  }
  ##
  stopifnot(!is.complex(psd))
  #psd <- as.rowvec(psd)
  ## Normalize by variance, 
  trap.area <- sum(psd) - psd[1]/2 - psd[length(psd)]/2 # Trapezoidal rule
  bandwidth <- 1 / nhalf
  timebp <- as.numeric(ntap/2)
  psd.n <- psd * (2 * varx / (trap.area * bandwidth))
  frq <- as.numeric(seq.int(0, 0.5, length.out=nfreq))
  #and (optionally) the Nyquist frequency so units will be in (units**2/Hz)
  if (Nyquist.normalize) {
    frq <- 2 * Nyq * frq
    psd.n <- Nyq * psd.n
  }
  # BUG: there seems to be an issue with f==0, & f[length(psd)]
  # so just extrapolate from the prev point
  if (first.last) psd.n <- (exp( signal::interp1(frq[2:(nfreq-1)], 
                                                          log(psd.n[2:(nfreq-1)]), 
                                                          frq, 
                                                          method='linear', 
                                                          extrap=TRUE) ))
  
  pltpsd <- function(...){
    Xpg <- spec.pgram(X, log="no", pad=1, taper=0, detrend=F, demean=F, plot=F)
    opar <- par(no.readonly = TRUE)
    par(mar=rep(2,4), oma=rep(0,4))
    layout(matrix(c(1,2),ncol=1),c(1,2))
    lpsd <- 10*log10(psd.n)
    lpgram <- 10*log10(Xpg$spec)
    r1 <- range(lpsd)
    r2 <- range(lpgram)
    plot(log10(Xpg$freq), lpgram, 
         col="red", type="l", main="spectra",
         ,ylim=c(min(r1,r2), max(r1,r2)))
    lines(log10(frq), lpsd, type="l")
    legend("bottomleft",c("spec.pgram","rlpSpec"),col=c("red","black"),lty=1)
    plot(X,type="l", main=sprintf("spec series ( dt:%s | dm:%s | f.l:%s )", demean, detrend, first.last))
    par(opar)
  }
  if (plotpsd) pltpsd(...)
  funcall<-paste(as.character(match.call()[]),collapse=" ") 
  psd.out <- list(freq = as.numeric(frq), 
                  spec = as.numeric(psd.n), 
                  coh = NULL, 
                  phase = NULL, 
                  kernel = NA, 
                  df = NA, 
                  numfreq = nfreq,
                  bandwidth = bandwidth, 
                  timebp=timebp,
                  n.used = envGet("len_even"), 
                  orig.n = envGet("len_orig"), 
                  series = series, 
                  snames = colnames(X), 
                  method = sprintf("Adaptive Sine Multitaper (rlpSpec)\n%s",funcall), 
                  taper = ntap, 
                  pad = TRUE, 
                  detrend = detrend, 
                  demean = demean)
  if (as.spec){class(psd.out) <- "spec"}
  return(invisible(psd.out))
}
