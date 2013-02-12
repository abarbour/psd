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
#'  If \code{X.frq} is NULL, the value is assumed to be 1, unless \code{X.d} is a 'ts' object.
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
#' @param preproc  logical; should \code{X.d} have a linear trend removed?
#' @param na.action  the function to deal with \code{NA} values
#' @param first.last  the extrapolates to give the zeroth and Nyquist frequency estimates
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
#' @example inst/Examples/rdex_psdcore.R
psdcore <- function(X.d, X.frq=NULL, ntaper=as.tapers(1), ndecimate=1L, preproc=TRUE, na.action = stats::na.fail, first.last=TRUE, plotpsd=FALSE, as.spec=TRUE, refresh=FALSE, verbose=FALSE, ...) UseMethod("psdcore")
#' @rdname psdcore
#' @method psdcore default
#' @S3method psdcore default
psdcore.default <- function(X.d, 
                            X.frq=NULL, 
                            ntaper=as.tapers(1), 
                            ndecimate=1L,
                            preproc=TRUE,
                            na.action = stats::na.fail,
                            first.last=TRUE,
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
  if (is.null(X.frq)){
    # make some assumptions about the sampling rate
    X.frq <- 1
    if (is.ts(X.d)){
      X.frq <- stats::frequency(X.d)
    }
    if (verbose) message(sprintf("Sampling frequency assumed to be  %f", X.frq))
  } else if (X.frq > 0){
    # value represents sampling frequency
    X.d <- na.action(stats::ts(X.d, frequency=X.frq))
  } else if (X.frq < 0){
    # value is sampling interval
    X.d <- na.action(stats::ts(X.d, deltat=abs(X.frq)))
  } else {
    stop("bad sampling information")
  }
  # sampling and nyquist
  X.frq <- stats::frequency(X.d)
  Nyq <- X.frq/2
  ##
  ###  When ntaper is a scalar, initialize
  ##
  # only one taper: usually means a first run
  lt <- length(ntaper)
  # onle one variable in the env (init): it hasn't been added to yet
  nenvar <- length(rlpSpec:::rlp_envStatus()$listing)
  if (lt == 1 | nenvar == 1 | refresh){
    # original series length
    n.o <- rlpSpec:::rlp_envAssignGet("len_orig", length(X.d))
    #
    X <- X.d
    if (preproc) X <- prewhiten(X, AR.max=0L, detrend=TRUE, plot=FALSE, verbose=verbose)$prew_lm
    #
    # Force series to be even in length (modulo division)
    # nextn(factors=2) ?
    n.e <- rlpSpec:::rlp_envAssignGet("len_even", n.o - n.o %% 2 )
    X.even <- as.matrix(X[seq_len(n.e)]) #1:n.e])
    rlpSpec:::rlp_envAssign("ser_orig", X)
    rlpSpec:::rlp_envAssign("ser_orig_even", X.even)
    # half length of even series
    nhalf <- rlpSpec:::rlp_envAssignGet("len_even_half", n.e/2)
    # variance of even series
    varx <- rlpSpec:::rlp_envAssignGet("ser_even_var", drop(stats::var(X.even)))
    # create uniform tapers
    nt <- nhalf + 1
    if (lt < nt) {
      ntap <- ntaper * ones(nt) 
    } else {
      ntap <- ntaper[seq_len(nt)] #1:nt]
    }
    ## zero pad and take double-length fft
    # fftw is faster (becomes apparent for long series)
    fftz <- fftw::FFT(as.numeric(c(X.even, zeros(n.e))))
    fftz <- rlpSpec:::rlp_envAssignGet("fft_even_demeaned_padded", fftz)
  } else {
    if (verbose){warning("Working environment *not* refreshed. Results may be bogus.")}
    X <- X.d
    ntap <- ntaper
    #stopifnot(length(X)==length(ntap))
    n.e <- rlpSpec:::rlp_envGet("len_even")
    nhalf <- rlpSpec:::rlp_envGet("len_even_half")
    varx <- rlpSpec:::rlp_envGet("ser_even_var")
    fftz <- rlpSpec:::rlp_envGet("fft_even_demeaned_padded")
  }
  #
  # if ntaper is a vector, this doesn't work [ ] ?
  ##
  # if the user wants a raw periodogram: by all means
  DOAS <- FALSE
  if (lt == 1){
    if (ntaper > 0) DOAS <- TRUE
  } else {
    if (!(is.tapers(ntap))) DOAS <- TRUE
  }
  ##TMP <<- ntap
  if (DOAS) ntap <- as.tapers(ntap, setspan=TRUE)
  
  ## interpolation
  ###  Select frequencies for PSD evaluation
  if  (lt > 1 && ndecimate > 1){
    stopifnot(!is.integer(ndecimate))
    if (verbose) message("decim stage 1")
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
    f <- base::seq.int(0, nhalf, by=1)
  }
  ##
  ###  Calculate the psd by averaging over tapered estimates
  nfreq <- length(f)
  NF <- seq_len(nfreq) # faster or slower than 1:nfreq?
  ##
  if (sum(ntap) > 0) {
    psd <- 0 * NF
    n2e <- n.e * 2
    # get a set of all possible weights for the current taper-vector
    # then the function need only subset the master set
    # faster? YES
    KPWM <- parabolic_weights_fast(max(ntap))
    PSDFUN <- function(fj, n2.e=n2e, KPW=KPWM, ntaps=ntap, Xfft=fftz){
      # number tapers (for subsetting)
      NT <- base::seq_len(ntaps[fj+1])
      # sequence
      Kseq <- KPW$taper_seq[NT]
      # weights
      #Kwgt <- KPW$taper_weights[NT]
      # Resampling weighted spectral values:
      fj2 <- fj * 2
      #m1. <- fj2 + n2.e - Kseq
      j1 <- (fj2 + n2.e - Kseq) %% n2.e
      #m1. <- fj2 + Kseq
      j2 <- (fj2 + Kseq) %% n2.e
      #f1 <- Xfft[j1+1]
      #f2 <- Xfft[j2+1]
      af12. <- Xfft[j1+1] - Xfft[j2+1]
      af122. <- af12. * af12. # will be complex, so use abs:
      psdv <- KPW$taper_weights[NT] %*% base::matrix(base::abs(af122.), ncol=1)
      return(psdv)
    }
    # ** compiled code doesn't appear to help speed
    #     PSDFUNc <- compiler::cmpfun(PSDFUN)
    # ** foreach is easier to follow, but foreach solution is actually slower :(
    #     psd <- foreach::foreach(f.j=f[1:nfreq], .combine="c") %do% PSDFUN(fj=f.j)
    # ** vapply is much faster than even lapply
    psd <- vapply(X=f[NF], FUN=PSDFUN, FUN.VALUE=double(1))
  } else {
    if (verbose) message("zero taper result == raw periodogram")
    Xfft <- rlpSpec:::rlp_envGet("fft_even_demeaned_padded")
    ff <- Xfft[NF]
    N0. <- rlpSpec:::rlp_envGet("len_orig")
    psd <- base::abs(ff * base::Conj(ff)) / N0.
  }
  ##  Interpolate if necessary to uniform freq sampling
  if (lt > 1 && ndecimate > 1){
    ## check [ ]
    if (verbose) message("decim stage 2")
    tmp.x <- f
    tmp.xi <- tmp.y
    tmp.y <- psd
    # x y x_interp --> y_interp
    tmp.yi <- signal::interp1(tmp.x, tmp.y, tmp.xi, method='linear', extrap=TRUE)
    psd <- tmp.yi
  }
  ##
  # should not be complex at this point!!
  stopifnot(!is.complex(psd))
  #psd <- as.rowvec(psd)
  ## Normalize by variance, 
  trap.area <- base::sum(psd) - psd[1]/2 - psd[length(psd)]/2 # Trapezoidal rule
  bandwidth <- 1 / nhalf
  ## normalize to two sided spectrum
  psd.n <- psd * (2 * varx / (trap.area * bandwidth))
  ## Nyquist frequencies
  frq <- as.numeric(base::seq.int(0, Nyq, length.out=nfreq)) # was just 0.5
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
  indic <- base::seq_len(nfreq - 2) + 1 
  #2:(nfreq-1)
  if (first.last) psd.n <- base::exp(signal::interp1(frq[indic], base::log(psd.n[indic]), frq, method='linear', extrap=TRUE))
  ##
  pltpsd <- function(Xser, frqs, psds, taps, nyq, detrend, demean, ...){
    fsamp <- frequency(Xser)
    stopifnot(fsamp==X.frq)
    Xpg <- spec.pgram(Xser, log="no", pad=1, taper=0.2, detrend=detrend, demean=detrend, plot=FALSE)
    # frequencies are appropriate,
    # but spectrum is normed for double-sided whereas rlpSpec single-sided; hence,
    # factor of 2
    Xpg <- normalize(Xpg, fsamp, "spectrum", verbose=FALSE)
    ##
    opar <- par(no.readonly = TRUE)
    par(mar=c(2, 3, 2.3, 1.2), oma=rep(2,4), las=1, tcl = -.3, mgp=c(2.2, 0.4, 0))
    #layout(matrix(c(1,2), ncol=1), c(1,2))
    layout(matrix(c(1,1,2,2,3,4), 3, 2, byrow = TRUE))
    ## Prelims:
    # rlp
    lfrq <- log10(frqs); rm(frqs)
    db_psd <- dB(psds/fsamp); rm(psds) # quick normalization
    # spec.pgram
    lfrqp <- log10(Xpg$freq)
    db_pgram <- dB(Xpg$spec); rm(Xpg)
    # plotting
    r1 <- range(db_psd, na.rm=TRUE, finite=TRUE)
    r2 <- range(db_pgram, na.rm=TRUE, finite=TRUE)
    ylims <- round(c(min(r1, r2), max(r1, r2)), 1)
    r1 <- range(lfrqp, na.rm=TRUE, finite=TRUE)
    r2 <- range(lfrq, na.rm=TRUE, finite=TRUE)
    xlims <- round(c(min(r1, r2), max(r1, r2)), 1) + .1*c(-1,1)
    ## Spectra, in decibels
    plot(lfrqp, db_pgram, col="red", type="l", 
         main="Naive and Multitaper PSD",
         xaxs="i", xlab="", xlim=xlims,
         ylab="dB, units^2 * delta", ylim=ylims)
    mtext("log10 frequency", side=1, line=1.6)
    #mtext(, cex=0.6)
    abline(h=3.01*c(-1,0,1), v=log10(nyq), col="dark gray", lwd=c(0.8,2,0.8), lty=c(4,3,4))
    lines(lfrq, db_psd, type="l")
    legend("bottomleft",c("spec.pgram (20% cosine taper)",sprintf("psdcore (max %i tapers)",max(taps))), col=c("red","black"), lty=1, lwd=2, cex=0.9)
    ## tapers
    if (is.tapers(taps)){
      plot(taps, lfrq, xlim=xlims, xaxs="i")
    } else {
      plot(lfrq, taps, type="h", xlim=xlims, xaxs="i")
    }
    ## original series
    plot(Xser, type="l", ylab="units", xlab="", xaxs="i", main="Modified series")
    mtext("index", side=1, line=1.5)
    mtext(sprintf("( dt+dm: %s | f.l: %s )", preproc, first.last), cex=0.4)
    ## autocorrelation
    acf(Xser, main="")
    mtext("lag", side=1, line=1.5)
    ## reset params
    par(opar)
  }
  ## Plot it
  if (plotpsd) pltpsd(Xser=X, frqs=frq, psds=psd.n, taps=ntap, nyq=Nyq, detrend=preproc, demean=preproc, ...)
  ##
  funcall <- sprintf("psdcore (dem.+detr. %s f.l. %s refr. %s)", preproc, first.last, refresh) 
  ## paste(as.character(match.call()[]),collapse=" ") 
  psd.out <- list(freq = as.numeric(frq), 
                  spec = as.numeric(psd.n), 
                  coh = NULL, 
                  phase = NULL, 
                  kernel = NULL, 
                  # must be a scalar for plot.spec to give conf ints:
                  df = 2 * mtap, # 2 DOF per taper, Percival and Walden eqn (370b)
                  numfreq = nfreq,
                  bandwidth = bandwidth, 
                  n.used = rlpSpec:::rlp_envGet("len_even"), 
                  orig.n = rlpSpec:::rlp_envGet("len_orig"), 
                  series = series, 
                  snames = colnames(X), 
                  method = sprintf("Sine multitaper\n%s",funcall), 
                  taper = ntap, 
                  pad = 1, # always!
                  detrend = preproc, 
                  demean = preproc,
                  timebp=timebp,
                  nyquist.frequency=Nyq
                  )
  if (as.spec){class(psd.out) <- "spec"}
  return(invisible(psd.out))
}
