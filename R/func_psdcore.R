#' Multitaper power spectral density of a series
#'
#' Compute a power spectral denisty (PSD) estimate 
#' for the input series using sine multitapers.
#'  
#' The parameter \code{ntaper} specifies the number of sine tapers to be used 
#' at each frequency: if it's a scalar, the same number of tapers will be used
#' at every frequency; otherwise, use ntaper(j) sine tapers at frequency(j).
#'
#' The series length N is truncated, if necessary, so that 1+N/2 evenly spaced
#' frequencies are returned. 
#'
#' The parameter \code{ndecimate} specifies the number of psds actually 
#' computed, defined as \code{(1+n/2)/ndecimate}; other
#' values are found via linear interpolation.
#'
#' @note Decimation is not well tested as of this point (December 2012).
#'
#' @param X.d  the series to estimate a spectrum for 
#' @param X.frq  scalar; the sampling frequency (e.g. Hz)
#' @param ntaper  scalar, or vector; the number of tapers
#' @param ndecimate  scalar; decimation factor
#' @param demean  logical; should \code{X.d} be centered about the mean
#' @param detrend  logical; should \code{X.d} have a linear trend removed
#' @param na.action  function dealing with \code{NA} values
#' @param first.last  the extrapolates to give the zeroth and Nyquist frequency estimates
#' @param Nyquist.normalize  logical; should the units be returned in Hz, rather than Nyquist?
#' @param plotpsd  logical; should the estimate be shown compared to the \code{spec.pgram} estimate
#' @param as.spec  logical; should the object returned be of class 'spec'
#' @param force_calc  logical; force spectrum (used for development purposes)
#' @param ...  (unused); parameters passed to [ NULL ]
#'
#' @import signal
#' @name psdcore
#' @export
#' @author Andrew Barbour <andy.barbour@@gmail.com> ported original by R.L.Parker.
#' @seealso \code{\link{pspectrum}}, \code{\link{riedsid}}
#'
#' @examples
#' X.d <- rnorm(1e3)
#' plot(psdcore(X.d, ntaper=10), log="dB", ylim=10*c(-1,1))
#' psd.n <- psdcore(X.d, ntaper=10, Nyquist.normalize=FALSE)
#' lines(psd.n$freq, 10*log10(psd.n$spec), col="red") # note normalization
#' abline(h=c(0, 3), col=c("black","red"), lwd=2)
#'
#' # 10Hz sampling
#' plot(psdcore(X.d, X.frq=10, ntaper=10), log="dB", ylim=10*c(-0.3,1.7))
#' psd.n <- psdcore(X.d, X.frq=10, ntaper=10, Nyquist.normalize=FALSE)
#' lines(10*psd.n$freq, 10*log10(psd.n$spec), col="red") # note normalization
#' abline(h=c(10, 3), col=c("black","red"), lwd=2)
#' 
#' # if ntaper is a vector:
#' psdcore(X.d, ntaper=rep(10,length(X.d))
psdcore <- function(X.d, 
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
                    ...) UseMethod("psdcore")

#' @rdname psdcore
#' @S3method psdcore rlpspec
psdcore.rlpspec <- function(...){.NotYetImplemented()}

#' @rdname psdcore
#' @S3method psdcore default
psdcore.default <- function(X.d, 
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
    # get a set of all possible weights for the current taper-vector
    # then the function need only subset the master set
    # faster? YES
    KPWM <- parabolic_weights_fast(max(ntap))
    PSDFUN <- function(fj, n2.e=n2e, KPW=KPWM, ntaps=ntap, Xfft=Rfftz){
      # parabolic weights, index m+1, column vec out
      #print(c(fj,fj2,n2.e))
      ##NT <- ntaps[fj+1]
      #KPW <- parabolic_weights(ntaps, tap.index=(fj+1), vec.out="horizontal")
      ##KPW <- parabolic_weights_fast(NT)
      #
      # num tapers (for subsetting)
      NT <- ntaps[fj+1]
      # sequence
      Kseq <- KPW$taper_seq[1:NT]
      # weights
      Kwgt <- KPW$taper_weights[1:NT]
      #
      fj2 <- 2*fj
      m1. <- fj2 + n2.e - Kseq #KPW$taper_seq
      j1 <- m1. %% n2.e
      m1. <- fj2 + Kseq #KPW$taper_seq
      j2 <- m1. %% n2.e
      f1 <- Xfft[j1+1]
      f2 <- Xfft[j2+1]
      af12. <- f1 - f2
      af122. <- af12. * af12.
      #psdv <- KPW$taper_weights %*% matrix(af122., ncol=1)
      psdv <- Kwgt %*% matrix(af122., ncol=1)
      return(psdv)
    }
    # ** compiled code doesn't appear to help speed
    #     require(compiler)
    #     PSDFUNc <- cmpfun(PSDFUN)
    # ** foreach is easier to follow, but foreach solution is actually slower :(
    #     require(foreach)
    #     psd <- foreach::foreach(f.j=f[1:nfreq], .combine="c") %do% PSDFUN(fj=f.j)
    # ** vapply is much faster than even lapply
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
