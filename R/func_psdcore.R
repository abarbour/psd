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
#' @param ndecimate  now ignored
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
#' @seealso \code{\link{pspectrum}}, \code{\link{riedsid}}, \code{\link{parabolic_weights}}
#'
#' @example inst/Examples/rdex_psdcore.R
psdcore <- function(X.d, ...) UseMethod("psdcore")

#' @rdname psdcore
#' @export
#psdcore.ts <- function(X.d, ...) .NotYetImplemented()

#' @rdname psdcore
#' @export
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
  if (ndecimate != 1){
    # force a warning if the user wanted to use decimation, which is no longer supported
    warning('Support for decimation has been removed.')
  }
  
  # get a clean environment
  if (refresh) psd_envRefresh(verbose=verbose)
  
  # named series (?)
  series <- deparse(substitute(X.d))
  
  if (is.null(X.frq)){
    # make an assumption about the sampling rate
    X.frq <- ifelse(is.ts(X.d), stats::frequency(X.d), 1)
    if (verbose) message(sprintf("Sampling frequency assumed to be  %f", X.frq))
  } 
  
  X.d <- if (X.frq > 0){
    # value represents sampling frequency
    na.action(stats::ts(X.d, frequency=X.frq))
  } else if (X.frq < 0){
    # value is sampling interval
    na.action(stats::ts(X.d, deltat=abs(X.frq)))
  } else {
    stop("bad sampling information")
  }
  
  # Refresh sampling rate, and get Nyquist frequency, tapers, and status
  X.frq <- stats::frequency(X.d)
  Nyq <- X.frq/2
  len_tapseq <- length(ntaper)
  #  only one variable in the env (init) means it hasn't been added to yet
  nenvar <- length(psd_envStatus()[['listing']])
  
  ops <- getOption("psd.ops")
  stopifnot(!is.null(ops))
  
  evars <- ops[['names']]
  
  #"fft_even_demeaned_padded"
  
  # initialize fft and other things, since this usually means a first run
  fftz <- if ( len_tapseq == 1 | nenvar == 1 | refresh ){
    
    # original series length
    n.o <- psd_envAssignGet(evars[['n.orig']], length(X.d))
    
    X <- psd_envAssignGet(evars[['series.orig']], {
      if (preproc){
        # option for fast-detrend only [ ]
        prewhiten(X.d, AR.max=0L, detrend=TRUE, plot=FALSE, verbose=verbose)$prew_lm
      } else {
        X.d
      }
    })
    
    # Force series to be even in length (modulo division)
    #
    n.e <- psd_envAssignGet(evars[['n.even']], modulo_floor_rcpp(n.o))
    even_seq <- seq_len(n.e)
    X.even <- psd_envAssignGet(evars[['series.even']], as.matrix(X[even_seq]))
    
    # half length of even series
    nhalf <- psd_envAssignGet(evars[['n.even.half']], n.e/2)
    
    # variance of even series
    varx <- psd_envAssignGet(evars[['var.even']], drop(stats::var(X.even)))
    
    # create uniform tapers
    kseq <- psd_envAssignGet(evars[['last.taper']], {
      if (len_tapseq == 1){
        message("tap c A ", len_tapseq, " ", nhalf)
        print(ntaper)
        rep.int(ntaper, nhalf+1)
      } else {
		message("tap c B ")
        tmptap <- ntaper[nhalf+1] # if length < nhalf + 1 the remnants will be NA
        tmptap[is.na(tmptap)] <- ops[['tapmin']]
      }
    })
    
    ## zero pad and take double-length fft (fftw is faster for very long series)
    padded <- as.numeric(c(X.even, zeros(n.e)))
    padded.fft <- psd_envAssignGet(evars[['fft']], fftw::FFT(padded))
    
    # Fix first value of fft -- always basically zero -- this will get rid 
	# of the bug, and prevent needing first-last extrapolation, I think!
	if (first.last){
		n.fft <- length(padded.fft)
		x.extrap <- c(2:(n.fft - 1))
		xi.extrap <- seq_len(n.fft)
		if (verbose) message("zero-frequency fft was extrapolated")
		padded.fft <- psd_envAssignGet(evars[['fft.extrap']], {
			signal::interp1(x = x.extrap, y=padded.fft[x.extrap], xi = xi.extrap, method='linear', extrap = TRUE)
			})
	}
    padded.fft
    
  } else {
    
    if (verbose) warning("Working environment *not* refreshed. Results may be bogus.")
    
    X <- X.d
    kseq <- ntaper
    #stopifnot(length(X)==length(kseq))
    n.e <- psd_envGet(evars[['n.even']])
    nhalf <- psd_envGet(evars[['n.even.half']])
    varx <- psd_envGet(evars[['var.even']])
    
    psd_envGet(evars[['fft']])
  }
  
  # TODO: if ntaper is a vector, this doesn't work [ ] ?
  
  # Switch: multitaper if TRUE, periodogram if FALSE
  DOMT <- if (len_tapseq == 1){
    ifelse(ntaper > 0, TRUE, FALSE)
  } else {
    ifelse(all(kseq > 0), TRUE, FALSE)
  }
  
  ###  Select frequencies for PSD evaluation
  message("f creation")
  f <- base::seq.int(0, nhalf, by=1)
  nfreq <- length(f)
  
  ###  Calculate the PSD by averaging over tapered estimates  
  PSD <- psd_envAssignGet(evars[["last.psdcore"]], {
    
    if (DOMT){
      
      if (verbose) message("multitaper psd")
      
      ## resample fft with taper sequence and quadratic weighting
      kseq <- as.integer(as.tapers(kseq, setspan=TRUE))
      #    this is where the majority of the work goes on:
      reff <- try(resample_fft_rcpp(fftz, kseq, verbose=verbose))
	  #    check status:
      if (inherits(reff,'try-error')){
      	stop("Could not resample fft... inspect with psd_envGet(",evars[['fft']],"), etc.")
      } else {
        reff[['psd']]
      }
      
    } else {
      
      if (verbose) message("raw periodogram")
      Xfft <- psd_envGet(evars[['fft']])
      ff <- Xfft[seq_len(nfreq)]
      N0. <- psd_envGet(evars[['n.orig']])
      # 
      # if the user wants it, then by all means
      # ... but force a warning on them
      warning("careful with these zero-taper results!")
      base::abs(ff * base::Conj(ff)) / N0.
      
    }
  })
  
  # should not be complex at this point!
  stopifnot(!is.complex(PSD))
  
  npsd <- length(PSD)
  nonfin <- is.infinite(PSD)
  
  if (any(nonfin)) PSD <- replace(PSD, nonfin, NA)
  
  ## Nyquist frequencies
  frq <- as.numeric(base::seq.int(0, Nyq, length.out=npsd))
  
  # First-last extrapolation
  # in spec.pgram: pgram[1, i, j] <- 0.5*(pgram[2, i, j] + pgram[N, i, j])
#   if (first.last){ 
#     #
#     # should this go before area calculation?? Yes.
#     #
#     # BUG (still appears -- Feb 2015): there may be an issue with f==0, & f[length(PSD)]
#     # so just extrapolate from the prev point
#     #
#     pre.extrap <- PSD[c(1:4, (npsd-2):npsd)]
#     indic <- seq_len(nfreq - 2) + 1
#     PSD <- base::exp(signal::interp1(frq[indic], base::log(PSD[indic]), frq, method='linear', extrap=TRUE))
#     post.extrap <- PSD[c(1:4, (npsd-2):npsd)]
#     #
#     print(rbind(pre.extrap, post.extrap, extrap.diff=pre.extrap - post.extrap))
#     #
#     if (verbose) message("psd estimates at zero and Nyquist were extrapolated")
#   }
  
  ## Normalize by variance using trapezoidal rule and convert to one-sided spectrum
  trap.area <- base::sum(PSD, na.rm=TRUE) - mean(PSD[c(1,npsd)], na.rm=TRUE)
  
  print(dB(varx))
  print(dB((trap.area / nhalf)))
  print(dB(varx / (trap.area / nhalf)))
  
  PSD <- 2 * PSD * varx / (trap.area / nhalf)
  
  funcall <- sprintf("psdcore (dem.+detr. %s f.l. %s refr. %s)", preproc, first.last, refresh) 
  
  ## Plot result locally
  # TODO: cleanup
  # TODO: class and method (?)
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
  
  ## Plot results
  if (plotpsd) pltpsd(Xser=X, frqs=frq, PSDS=PSD, taps=kseq, nyq=Nyq, detrend=preproc, demean=preproc, ...)
  
  ## Return results
  mtap <- max(kseq, na.rm=TRUE)
  PSD.out <- list(freq = as.numeric(frq), 
                  spec = as.numeric(PSD), 
                  coh = NULL, 
                  phase = NULL, 
                  kernel = NULL, 
                  # must be a scalar for plot.spec to give conf ints:
                  df = 2 * mtap, # 2 DOF per taper, Percival and Walden eqn (370b)
                  numfreq = nfreq,    
                  ## bandwidth
                  # http://biomet.oxfordjournals.org/content/82/1/201.full.pdf
                  # half-width W = (K + 1)/{2(N + 1)}
                  # effective bandwidth ~ 2 W (accurate for many spectral windows)
                  bandwidth = (mtap + 1) / nhalf, 
                  n.used = psd_envGet(evars[['n.even']]), 
                  orig.n = psd_envGet(evars[['n.orig']]), 
                  series = series, 
                  snames = colnames(X), 
                  method = sprintf("Sine multitaper\n%s",funcall), 
                  taper = kseq, 
                  pad = TRUE, # always!
                  detrend = preproc, # always true?
                  demean = preproc,
                  timebp = as.numeric(kseq/2),
                  nyquist.frequency = Nyq
  )
  if (as.spec){class(PSD.out) <- c("spec","amt")}
  return(invisible(PSD.out))
}
