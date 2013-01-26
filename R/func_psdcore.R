#' Multitaper power spectral density estimates of a series.
#'
#' Compute power spectral density (PSD) estimates
#' for the input series using sine multitapers.
#'  
#' @details
#' \subsection{Tapering}{
#' The parameter \code{ntaper} specifies the number of sine tapers to be used 
#' at each frequency: equal tapers at each frequency for a scalar; 
#' otherwise, use \code{ntaper(j)} sine tapers at \code{frequency(j)}.
#' }
#'
#' \subsection{Truncation}{
#' The series length \code{N} is truncated, if necessary, so that \code{1+N/2} evenly 
#' spaced frequencies are returned. 
#' }
#'
#' \subsection{Decimation}{
#' The parameter \code{ndecimate} determines the PSDs actually 
#' computed, defined as \code{(1+n/2)/ndecimate}; other
#' values are found via linear interpolation.
#' }
#'
#' \subsection{Sampling}{
#'  If \code{X.frq > 0} it's assumed the value represents \emph{frequency} (e.g. Hz).
#'  If \code{X.frq < 0} it's assumed the value represents \emph{interval} (e.g. seconds).
#' }
#'
#' @section Warning:
#' Decimation is not well tested as of this point.
#'
#' @param X.d  the series to estimate a spectrum for 
#' @param X.frq  scalar; the sampling information (see section Sampling)
#' @param ntaper  scalar, or vector; the number of tapers
#' @param ndecimate  scalar; decimation factor
#' @param demean  logical; should \code{X.d} be centered about the mean
#' @param detrend  logical; should \code{X.d} have a linear trend removed
#' @param na.action  function dealing with \code{NA} values
#' @param first.last  the extrapolates to give the zeroth and Nyquist frequency estimates
#' @param Nyquist.normalize  logical; should the units be returned in Hz, rather than Nyquist?
#' @param plotpsd  logical; should the estimate be shown compared to the \code{spec.pgram} estimate
#' @param as.spec  logical; should the object returned be of class 'spec'?
#' @param refresh  logical; ensure a free environment prior to execution
#' @param verbose logical; should messages be given?
#' @param ...  (unused) Optional parameters
#' @return An list object, invibibly.  If \code{as.spec=TRUE} then an object with class 'spec'.
#'
#' @name psdcore
#' @export
#' @keywords spectrum-estimation normalization prewhiten
#' @author A.J. Barbour <andy.barbour@@gmail.com> adapted original by R.L.Parker.
#' @seealso \code{\link{pspectrum}}, \code{\link{riedsid}}
#'
# @example x_examp/psdcore.R
psdcore <- function(X.d, X.frq=1, ntaper=as.taper(1), ndecimate=1L, demean=TRUE, detrend=TRUE, na.action = stats::na.fail, first.last=TRUE, Nyquist.normalize=TRUE, plotpsd=FALSE, as.spec=TRUE, refresh=FALSE, verbose=FALSE, ...) UseMethod("psdcore")
#' @rdname psdcore
#' @method psdcore default
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
                            refresh=FALSE,
                            verbose=FALSE,
                            ...
                           ) {
  #
  if (refresh) rlpSpec:::rlp_envClear(verbose=verbose)
  #
  series <- deparse(substitute(X.d))
  if (X.frq > 0){
    # value represents sampling frequency
    X.d <- na.action(stats::ts(X.d, frequency=X.frq))
    Nyq <- stats::frequency(X.d)/2
  } else if (X.frq < 0){
    # value is sampling interval
    X.frq <- abs(X.frq)
    X.d <- na.action(stats::ts(X.d, deltat=X.frq))
    Nyq <- 1/stats::deltat(X.d)/2
  } else {
    stop("bad sampling information")
  }
  #   X.d <- na.action(stats::ts(X.d, frequency=X.frq))
  #   Nyq <- stats::frequency(X.d)/2
  ##
  ###  When ntaper is a scalar, initialize
  ##
  # only one taper: usually means a first run
  lt <- length(ntaper)
  # onle one variable in the env (init): it hasn't been added to yet
  nenvar <- length(rlp_envStatus()$listing)
  if (lt == 1 | nenvar == 1 | refresh){
    # original series length
    n.o <- rlp_envAssignGet("len_orig", length(X.d))
    #
    X <- prewhiten(X.d, 
                   AR.max=0L, 
                   detrend=detrend, 
                   demean=demean,
                   plot=FALSE,
                   verbose=FALSE)
    #
    # Force series to be even in length (modulo division)
    n.e <- rlp_envAssignGet("len_even", n.o - n.o%%2 )
    X.even <- as.matrix(X[1:n.e])
    rlp_envAssign("ser_orig", X)
    rlp_envAssign("ser_orig_even", X.even)
    # half length of even series
    nhalf <- rlp_envAssignGet("len_even_half", n.e/2)
    # variance of even series
    varx <- rlp_envAssignGet("ser_even_var", drop(stats::var(X.even)))
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
    fftz <- rlp_envAssignGet("fft_even_demeaned_padded", stats::mvfft(X.dem))
  } else {
    X <- X.d
    ntap <- ntaper
    #stopifnot(length(X)==length(ntap))
    n.e <- rlp_envGet("len_even")
    nhalf <- rlp_envGet("len_even_half")
    varx <- rlp_envGet("ser_even_var")
    fftz <- rlp_envGet("fft_even_demeaned_padded")
  }
  # if ntaper is a vector, this doesn't work [ ]
  ##
  if (!(is.taper(ntap))) ntap <- as.taper(ntap, setspan=TRUE)
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
  lt2 <- length(drop(ntap))
  ##
  ###  Calculate the psd by averaging over tapered estimates
  nfreq <- length(f)
  ##
  if (sum(ntap) > 0) {
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
    #     PSDFUNc <- compiler::cmpfun(PSDFUN)
    # ** foreach is easier to follow, but foreach solution is actually slower :(
    #     psd <- foreach::foreach(f.j=f[1:nfreq], .combine="c") %do% PSDFUN(fj=f.j)
    # ** vapply is much faster than even lapply
    psd <- vapply(X=f[1:nfreq], FUN=PSDFUN, FUN.VALUE=1.0)
  } else {
    message("zero taper result == raw periodogram")
    Xfft <- rlp_envGet("fft_even_demeaned_padded")
    ff <- Xfft[1:nfreq]
    N0 <- rlp_envGet("len_orig")
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
  psd.n <- psd * (2 * varx / (trap.area * bandwidth))
  frq <- as.numeric(seq.int(0, 0.5, length.out=nfreq))
  #and (optionally) the Nyquist frequency so units will be in (units**2/Hz)
  if (Nyquist.normalize) {
    frq <- 2 * Nyq * frq
    psd.n <- Nyq * psd.n
  }
  ## timebp
  timebp <- as.numeric(ntap/2)
  ## bandwidth
  # http://biomet.oxfordjournals.org/content/82/1/201.full.pdf
  # half-width W = (K + 1)/{2(N + 1)}
  # effective bandwidth ~ 2 W (accurate for many spectral windows)
  mtap <- max(ntap)
  bandwidth <- bandwidth * (mtap + 1) 
  #
  # BUG: there seems to be an issue with f==0, & f[length(psd)]
  # so just extrapolate from the prev point
  indic <- 2:(nfreq-1)
  if (first.last) psd.n <- exp(signal::interp1(frq[indic], log(psd.n[indic]), frq, method='linear', extrap=TRUE))
  ##
  pltpsd <- function(Xser, frqs, psds, nyq, nNyq, detrend, demean, ...){
    Xser <- ts(Xser, frequency=1) # so we can normalize properly
    Xpg <- spec.pgram(Xser, log="no", pad=1, taper=0.2, detrend=detrend, demean=demean, plot=FALSE)
    if (nNyq) {
      Xpg$freq <- Xpg$freq * 2 * nyq
      Xpg$spec <- nyq * Xpg$spec
    }
    opar <- par(no.readonly = TRUE)
    par(mar=c(2, 3, 2.3, 1.2), oma=rep(2,4), las=1, tcl = -.3, mgp=c(2.2, 0.4, 0))
    #layout(matrix(c(1,2), ncol=1), c(1,2))
    layout(matrix(c(1,1,2,2,3,4), 3, 2, byrow = TRUE))
    ## Prelims:
    # rlp
    lfrq <- log10(frqs)
    rm(frqs)
    db_psd <- 10*log10(psds)
    rm(psds)
    # spec.pgram
    db_pgram <- 10*log10(Xpg$spec)
    lfrqp <- log10(Xpg$freq)
    rm(Xpg)
    r1 <- range(db_psd)
    r2 <- range(db_pgram)
    ## Spectra, in decibels
    plot(lfrqp, db_pgram, col="red", type="l", 
         main="Naive and Multitaper PSD",
         xaxs="i", xlab="",
         ylab="dB, units^2 * delta",
         ylim=c(min(r1,r2), max(r1,r2)))
    mtext("frequency, log10 1 / delta", side=1, line=1.5)
    #mtext(, cex=0.6)
    abline(h=3.01*c(-1,0,1), v=log10(nyq), col="dark gray", lwd=c(0.8,2,0.8), lty=c(4,3,4))
    lines(lfrq, db_psd, type="l")
    legend("bottomleft",c("20% cosine","rlpSpec"), col=c("red","black"), lty=1, lwd=2, cex=0.9)
    ## tapers
    plot(ntap, lfrq)
    ## original series
    plot(Xser, type="l", ylab="units", xlab="", xaxs="i", main="Modified series")
    mtext("index", side=1, line=1.5)
    mtext(sprintf("( dt:%s | dm:%s | f.l:%s )", demean, detrend, first.last), cex=0.6)
    ## autocorrelation
    acf(Xser, main="")
    mtext("lag", side=1, line=1.5)
    ## reset params
    par(opar)
  }
  ## Plot it
  if (plotpsd) pltpsd(Xser=X, frqs=frq, psds=psd.n, nyq=Nyq, nNyq=Nyquist.normalize, 
                      detrend=detrend, demean=demean, ...)
  ##
  funcall<-paste(as.character(match.call()[]),collapse=" ") 
  psd.out <- list(freq = as.numeric(frq), 
                  spec = as.numeric(psd.n), 
                  coh = NULL, 
                  phase = NULL, 
                  kernel = NULL, 
                  df = 2*mtap, #mtap-1, # must be a scalar for plot.spec to give conf ints:
                  # Percival and Walden eqn (370b)
                  numfreq = nfreq,
                  bandwidth = bandwidth, 
                  n.used = rlp_envGet("len_even"), 
                  orig.n = rlp_envGet("len_orig"), 
                  series = series, 
                  snames = colnames(X), 
                  method = sprintf("Adaptive Sine Multitaper (rlpSpec)\n%s",funcall), 
                  taper = ntap, 
                  pad = 1, 
                  detrend = detrend, 
                  demean = demean,
                  timebp=timebp,
                  nyquist.frequency=Nyq,
                  nyquist.normalized=Nyquist.normalize
                  )
  if (as.spec){class(psd.out) <- "spec"}
  return(invisible(psd.out))
}
