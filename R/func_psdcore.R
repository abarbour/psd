#' Multitaper power spectral density estimates of a series.
#'
#' Compute power spectral density (PSD) estimates
#' for the input series using sine multitapers.
#'  
#' @details
#' \subsection{Tapering}{
#' The parameter \code{ntaper} specifies the number of sine tapers to be used 
#' at each frequency: equal tapers at each frequency for a scalar; 
#' otherwise, use \code{ntaper[j]} sine tapers at \code{frequency[j]}.
#' }
#'
#' \subsection{Truncation}{
#' The series, with length \eqn{N}, is necessarily truncated so that \eqn{1+N/2} evenly 
#' spaced frequencies are returned.  This truncation makes the series length ``highly composite",
#' which the discrete Fourier transform (DFT) is most efficient.
#' The vignette "fftw" (accessed with \code{vignette("fftw",package="psd")}) shows
#' how the performance of a DFT can be affected by series length.
#' }
#'
#' \subsection{Decimation}{
#' The parameter \code{ndecimate} determines the number of PSD estimates actually 
#' computed.  This number is defined as a fraction of the truncated length, \eqn{(1+N/2)/n_d}.
#' Linear interpolation is used.
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
#' The \code{first.last} parameter is a workaround for potential bug (under investigation), 
#' which causes
#' the power at the zero and Nyquist frequencies to have anomalously low values.  
#' This argument enables using
#' linear \emph{extrapolation} to correct these values.
#' \strong{The feature will be deprecated if the supposed bug is both identified and fixed.}
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
#' @return An list object, invisibly.  If \code{as.spec=TRUE} then an object with class \code{spec};
#' otherwise the list object will have information similar to a \code{spec} object, but with 
#' a few additional fields.
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
  if (refresh) psd::psd_envRefresh(verbose=verbose)
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
  nenvar <- length(psd::psd_envStatus()$listing)
  if (lt == 1 | nenvar == 1 | refresh){
    # original series length
    n.o <- psd::psd_envAssignGet("len_orig", length(X.d))
    #
    X <- X.d
    if (preproc) X <- prewhiten(X, AR.max=0L, detrend=TRUE, plot=FALSE, verbose=verbose)$prew_lm
    #
    # Force series to be even in length (modulo division)
    # nextn(factors=2) ?
    n.e <- psd::psd_envAssignGet("len_even", n.o - n.o %% 2 )
    X.even <- as.matrix(X[seq_len(n.e)])
    psd::psd_envAssign("ser_orig", X)
    psd::psd_envAssign("ser_orig_even", X.even)
    # half length of even series
    nhalf <- psd::psd_envAssignGet("len_even_half", n.e/2)
    # variance of even series
    varx <- psd::psd_envAssignGet("ser_even_var", drop(stats::var(X.even)))
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
    fftz <- psd::psd_envAssignGet("fft_even_demeaned_padded", fftz)
  } else {
    if (verbose){warning("Working environment *not* refreshed. Results may be bogus.")}
    X <- X.d
    ntap <- ntaper
    #stopifnot(length(X)==length(ntap))
    n.e <- psd::psd_envGet("len_even")
    nhalf <- psd::psd_envGet("len_even_half")
    varx <- psd::psd_envGet("ser_even_var")
    fftz <- psd::psd_envGet("fft_even_demeaned_padded")
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
  ###  Calculate the PSD by averaging over tapered estimates
  nfreq <- length(f)
  NF <- seq_len(nfreq) # faster or slower than 1:nfreq?
  ##
  if (sum(ntap) > 0) {
    PSD <- 0 * NF
    n2e <- n.e * 2
    fjs <- f[NF]
    fjs2 <- fjs * 2
    # get a set of all possible weights for the current taper-vector
    # then the function need only subset the master set
    # faster? YES
    # BUT this is wrong, because it normalizes by the max taper!!
    #KPWM <- parabolic_weights_fast(max(ntap))
    
    PSDFUN <- function(fj, fi=fj+1, fj2=fjs2[fi], n2.e=n2e, ntaps=ntap, Xfft=fftz){
      # number tapers (for subsetting)
      NT <- ntaps[fi]
      #NT <- base::seq_len(NT)
      # sequence
      KPW <- parabolic_weights_fast(NT)
      Kseq <- KPW$taper_seq
      # weights
      Kwgt <- KPW$taper_weights
      #rm(KPW)
      # Resampling weighted spectral values:
      j1 <- (fj2 + n2.e - Kseq) %% n2.e
      j2 <- (fj2 + Kseq) %% n2.e
      # select spectra
      af12. <- abs(Xfft[j1+1] - Xfft[j2+1])
      af12. <- matrix(af12. * af12., ncol=1)
      # reweighted spectra
      PSDv <- Kwgt %*% af12.
      return(PSDv)
    }
    # ** compiled code doesn't appear to help speed
    #     PSDFUNc <- compiler::cmpfun(PSDFUN)
    # ** foreach is easier to follow, but foreach solution is actually slower :(
    #     PSD <- foreach::foreach(f.j=f[1:nfreq], .combine="c") %do% PSDFUN(fj=f.j)
    # ** vapply is much faster than even lapply
    PSD <- vapply(X=fjs, FUN=PSDFUN, FUN.VALUE=double(1))
  } else {
    if (verbose) message("zero taper result == raw periodogram")
    Xfft <- psd::psd_envGet("fft_even_demeaned_padded")
    ff <- Xfft[NF]
    N0. <- psd::psd_envGet("len_orig")
    PSD <- base::abs(ff * base::Conj(ff)) / N0.
  }
  ##  Interpolate if necessary to uniform freq sampling
  if (lt > 1 && ndecimate > 1){
    ## check [ ]
    if (verbose) message("decim stage 2")
    tmp.x <- f
    tmp.xi <- tmp.y
    tmp.y <- PSD
    # x y x_interp --> y_interp
    tmp.yi <- signal::interp1(tmp.x, tmp.y, tmp.xi, method='linear', extrap=TRUE)
    PSD <- tmp.yi
  }
  ##
  # should not be complex at this point!!
  stopifnot(!is.complex(PSD))
  #PSD <- as.rowvec(PSD)
  ## Normalize by variance, 
  trap.area <- base::sum(PSD) - PSD[1]/2 - PSD[length(PSD)]/2 # Trapezoidal rule
  bandwidth <- 1 / nhalf
  ## normalize to one sided spectrum
  PSD.n <- PSD * (2 * varx / (trap.area * bandwidth))
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
  if (first.last){ 
    # BUG: there may be an issue with f==0, & f[length(PSD)]
    # so just extrapolate from the prev point
    indic <- base::seq_len(nfreq - 2) + 1 
    PSD.n <- base::exp(signal::interp1(frq[indic], base::log(PSD.n[indic]), frq, method='linear', extrap=TRUE))
    if (verbose) message("first.last=TRUE: Zero and Nyquist frequencies were extrapolated")
  }
  ##
  funcall <- sprintf("psdcore (dem.+detr. %s f.l. %s refr. %s)", preproc, first.last, refresh) 
  ##
  pltpsd <- function(Xser, frqs, PSDS, taps, nyq, detrend, demean, ...){
    fsamp <- frequency(Xser)
    stopifnot(fsamp==X.frq)
    Xpg <- spec.pgram(Xser, log="no", pad=1, taper=0.2, detrend=detrend, demean=detrend, plot=FALSE)
    # frequencies are appropriate,
    # but spectrum is normed for double-sided whereas psd single-sided; hence,
    # factor of 2
    Xpg <- normalize(Xpg, fsamp, "spectrum", verbose=FALSE)
    ##
    opar <- par(no.readonly = TRUE)
    par(mar=c(2, 3, 2.3, 1.2), oma=rep(2,4), las=1, tcl = -.3, mgp=c(2.2, 0.4, 0))
    #layout(matrix(c(1,2), ncol=1), c(1,2))
    layout(matrix(c(1,1,2,2,3,4), 3, 2, byrow = TRUE))
    ## Prelims:
    # psd
    lfrq <- log10(frqs); rm(frqs)
    db_PSD <- dB(PSDS/fsamp); rm(PSDS) # quick normalization
    # spec.pgram
    lfrqp <- log10(Xpg$freq)
    db_pgram <- dB(Xpg$spec); rm(Xpg)
    # plotting
    r1 <- range(db_PSD, na.rm=TRUE, finite=TRUE)
    r2 <- range(db_pgram, na.rm=TRUE, finite=TRUE)
    ylims <- round(c(min(r1, r2), max(r1, r2)), 1)
    r1 <- range(lfrqp, na.rm=TRUE, finite=TRUE)
    r2 <- range(lfrq, na.rm=TRUE, finite=TRUE)
    xlims <- round(c(min(r1, r2), max(r1, r2)), 1) + .1*c(-1,1)
    ## Spectra, in decibels
    plot(lfrqp, db_pgram, col="red", type="l", 
         main=funcall,
         xaxs="i", xlab="", xlim=xlims,
         ylab="dB, units^2 * delta", ylim=ylims)
    mtext("log10 frequency", side=1, line=1.6)
    #mtext(, cex=0.6)
    abline(h=3.01*c(-1,0,1), v=log10(nyq), col="dark gray", lwd=c(0.8,2,0.8), lty=c(4,3,4))
    lines(lfrq, db_PSD, type="l")
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
  if (plotpsd) pltpsd(Xser=X, frqs=frq, PSDS=PSD.n, taps=ntap, nyq=Nyq, detrend=preproc, demean=preproc, ...)
  ##
 ## paste(as.character(match.call()[]),collapse=" ") 
  PSD.out <- list(freq = as.numeric(frq), 
                  spec = as.numeric(PSD.n), 
                  coh = NULL, 
                  phase = NULL, 
                  kernel = NULL, 
                  # must be a scalar for plot.spec to give conf ints:
                  df = 2 * mtap, # 2 DOF per taper, Percival and Walden eqn (370b)
                  numfreq = nfreq,
                  bandwidth = bandwidth, 
                  n.used = psd::psd_envGet("len_even"), 
                  orig.n = psd::psd_envGet("len_orig"), 
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
  if (as.spec){class(PSD.out) <- c("spec","amt")}
  return(invisible(PSD.out))
}
