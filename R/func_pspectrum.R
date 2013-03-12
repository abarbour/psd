#' Adaptive sine multitaper power spectral density estimation.
#' 
#' This is the primary function to be used in this package, and returns
#' power spectral density estimates where the number of tapers at each
#' frequency has been iteratively optimized (\code{niter} times).
#'
#' See the \strong{Adaptive estimation} section in the description of
#' the \code{\link{psd-package}} for details regarding adaptive estimation.
#'
#' @name pspectrum
#' @export
#' @author A.J. Barbour <andy.barbour@@gmail.com> adapted original by R.L. Parker.
#' @seealso \code{\link{psdcore}}, \code{\link{pilot_spec}}, \code{\link{riedsid}}, \code{\link{prewhiten}}
#' @keywords spectrum-estimation riedel-sidorenko tapers tapers-constraints tapers-weighting numerical-derivative
#' 
#' @param x vector; series to estimate PSD for.
#' @param x.frqsamp scalar; the sampling rate (e.g. Hz) of the series \code{x}.
#' @param ntap_pilot scalar; the number of sine tapers to use in the pilot spectrum estimation.
#' @param niter scalar; the number of adaptive iterations to execute after the pilot spectrum.
#' @param AR logical; should the effects of an AR model be removed from the pilot spectrum?
#' @param Nyquist.normalize  logical; should the units be returned on Hz, rather than Nyquist?
#' @param verbose logical; Should messages be given?
#' @param no.history logical; Should the adaptive history \emph{not} be saved?
#' @param plot logical; Should the results be plotted?
#' @param ... Optional parameters passed to \code{\link{riedsid}}
#' @return Object with class 'spec', invisibly. It also assigns the object to
#' \code{"final_psd"} in the working environment.
#'
#' @example inst/Examples/rdex_pspectrum.R
pspectrum <- function(x, x.frqsamp=1, ntap_pilot=7, niter=5, AR=FALSE, Nyquist.normalize=TRUE, verbose=TRUE, no.history=FALSE, plot=FALSE, ...) UseMethod("pspectrum")
#' @rdname pspectrum
#' @method pspectrum default
#' @S3method pspectrum default
pspectrum.default <- function(x, x.frqsamp=1, ntap_pilot=7, niter=5, AR=FALSE, Nyquist.normalize=TRUE, verbose=TRUE, no.history=FALSE, plot=FALSE, ...){
  stopifnot(length(x)>1)
  #
  adapt_message <- function(stage, dvar=NULL){
    stopifnot(stage>=0)
    if (stage==0){
      stage <- paste(stage,"est. (pilot)")
    } else {
      stage <- paste(stage, sprintf("est. (Ave. S.V.R. %.01f dB)", dB(dvar)))
    }
    message(sprintf("Stage  %s ", stage))
  }
  #
  niter <- abs(niter)
  plotpsd_ <- FALSE
  for (stage in 0:niter){
    if (stage==0){
      if (verbose) adapt_message(stage)
      # --- setup the environment ---
      psd:::psd_envRefresh(verbose=verbose)
      # --- pilot spec ---
      # ** normalization is here:
      if (niter==0){
        if (plot) plotpsd_ <- TRUE
      }
      ##
      ordAR <- ifelse(AR, 100, 0)
      pilot_spec(x, x.frequency=x.frqsamp, ntap=ntap_pilot, 
                 remove.AR=ordAR, verbose=verbose, plot=plotpsd_)
      # ensure it's in the environment
      Pspec <- psd:::psd_envGet("pilot_psd")
      xo <- psd:::psd_envAssignGet("original_series", x)
      # starting spec variance
      dvar.o <- varddiff(Pspec$spec)
      # --- history ---
      save_hist <- ifelse(niter < 10, TRUE, FALSE)
      if (no.history) save_hist <- FALSE
      if (save_hist){
        psd:::new_adapt_history(niter)
        psd:::update_adapt_history(0, Pspec$taper, Pspec$spec, Pspec$freq)
      }
      #xo <- 0 # to prevent passing orig data back/forth
    } else {
      # enforce no verbosity
      rverb <- ifelse(stage > 0, FALSE, TRUE)
      ## calculate optimal tapers
      kopt <- riedsid(Pspec, verbose=rverb, ...)
      stopifnot(exists('kopt'))
      ## reapply to spectrum
      if (stage==niter & plot){
        plotpsd_ <- TRUE
        #xo <- x
        #rm(x)
      }
      # preproc done in pilot_spec
      Pspec <- psdcore(X.d=xo, X.frq=x.frqsamp, ntaper=kopt, 
                       preproc=FALSE, plotpsd=plotpsd_, verbose=FALSE)
      if (verbose) if (verbose) adapt_message(stage, vardiff(Pspec$spec, double.diff=TRUE)/dvar.o)
      ## update history
      if (save_hist) psd:::update_adapt_history(stage, Pspec$taper, Pspec$spec)
    }
  }
  if (Nyquist.normalize) Pspec <- normalize(Pspec, x.frqsamp, "psd", verbose=verbose)
  return(invisible(psd:::psd_envAssignGet("final_psd", Pspec)))
}

#' Calculate the pilot power spectral densities.
#'
#' This PSD -- the pilot spectrum -- is used as the starting point
#' for the adaptive estimation routine.
#'
#' A fixed number
#' of tapers is applied across all frequencies using \code{\link{psdcore}}, and
#' subsequent taper-refinements are based on the spectral derivatives
#' of this spectrum; hence, changes in the number of tapers can affect
#' how many adaptive stages may be needed (though there are no formal convergence
#' criteria to speak of).
#'
#' The taper series of the returned spectrum is constrained using
#' \code{as.tapers(..., minspan=TRUE)}.
#'
#' The default behaviour (\code{remove.AR <= 0}) is to remove the standard linear 
#' model \eqn{[f(x) = \alpha x + \beta]} from the data; however,
#' the user can remove the effect of an autoregressive process by specifiying
#' \code{remove.AR}.
#'
#' @section Removing an AR effect from the spectrum:
#' If \code{remove.AR > 0} the argument is used as \code{AR.max} in 
#' \code{\link{prewhiten}}, from which an AR spectrum is calculated using
#' the best fitting model; this is removed from the spectrum calculated
#' with the full data.
#'
#' If the value of \code{remove.AR} is too low the spectrum 
#' could become distorted,
#' so use with care.
#' \emph{Note, however, that the 
#' value of \code{remove.AR} will be restricted to within the 
#' range \eqn{[1,100]}.}
#' If the AR order is much larger than this, it's unclear how \code{\link{prewhiten}}
#' will perform.
#'
#' @name pilot_spec
#' @aliases pilot_spectrum spec.pilot
#' @export
#' @author A.J. Barbour <andy.barbour@@gmail.com>
#' @seealso \code{\link{psdcore}}, \code{\link{prewhiten}}
#'
#' @param x  vector; the data series to find a pilot spectrum for
#' @param x.frequency  scalar; the sampling frequency (e.g. Hz) of the series
#' @param ntap  scalar; the number of tapers to apply during spectrum estimation
#' @param remove.AR  scalar; the max AR model to be removed from the data.
#' @param plot  logical; should a plot be created?
#' @param verbose  logical; should messages be given?
#' @param ...  additional parameters passed to \code{\link{psdcore}}
#' @return An object with class 'spec', invisibly.  It also assigns the object to
#' \code{"pilot_psd"} in the working environment.
#'
#' @example inst/Examples/rdex_pilotspec.R
pilot_spec <- function(x, x.frequency=1, ntap=7, remove.AR=0, plot=FALSE, verbose=FALSE, ...) UseMethod("pilot_spec")
#' @rdname pilot_spec
#' @method pilot_spec default
#' @S3method pilot_spec default
pilot_spec.default <- function(x, x.frequency=1, ntap=7, remove.AR=0, plot=FALSE, verbose=FALSE, ...){
  stopifnot(length(ntap)==1)
  stopifnot(length(remove.AR)==1)
  if (is.ts(x)) x.frequency <- stats::frequency(x)
  # setup a universal calculator
  PSDFUN <- function(X.., Xf.., Xk.., AR=FALSE){
    toret <- psdcore(X.., Xf.., Xk.., 
                     preproc=FALSE, 
                     first.last=!AR, 
                     as.spec=TRUE, 
                     refresh=TRUE, 
                     verbose=FALSE)
    return(toret)
  }
  # preprocess
  REMAR <- FALSE
  if (remove.AR > 0){
    REMAR <- TRUE
    # restrict to within [1,100]
    remove.AR <- max(1, min(100, abs(remove.AR)))
  }
  xprew <- prewhiten(x, AR.max=remove.AR, detrend=TRUE, 
                     impute=TRUE, plot=FALSE, verbose=verbose)
  if (REMAR){
    # AR fit
    ordAR <- xprew$ardfit$order
    if (ordAR==0){
      warning("AR(0) was the highest model found!\n\t\tConsider fitting a linear model instead ( remove.AR = 0 ).")
    } else {
      if (verbose) message(sprintf("removed AR(%s) effects from the spectrum", ordAR))
    }
    xar <- xprew$prew_ar
    # PSD of the AR fit
    Pspec_ar <- PSDFUN(xar, x.frequency, ntap, AR=TRUE)
    arvar <- var(Pspec_ar$spec)
    Pspec_ar$spec <- Pspec_ar$spec / (mARs <- mean(Pspec_ar$spec))
  }
  #
  # Initial spectrum:
  Pspec <- PSDFUN(xprew$prew_lm, x.frequency, ntap, AR=FALSE)
  num_frq <- length(Pspec$freq)
  num_tap <- length(Pspec$taper)
  stopifnot(num_tap <= num_frq)
  Ptap <- Pspec$taper
  # generate a series, if necessary
  if (num_tap < num_frq) Ptap <- rep.int(Ptap[1], num_frq)
  # return tapers object
  Pspec$taper <- as.tapers(Ptap, setspan=TRUE)
  ##
  ## remove the spectrum of the AR process
  if (REMAR){
    stopifnot(exists("Pspec_ar"))
    if (verbose) message(sprintf("Removing AR(%s) effects from spectrum", ordAR))
    Ospec <- Pspec
    Pspec$spec <- Pspec$spec / Pspec_ar$spec
    # reup the spectrum
    psd:::psd_envAssignGet("AR_psd", Pspec_ar)
  }
  if (plot){
    if (REMAR){
      par(las=1)
      plot(Ospec, log="dB", col="red", main="Pilot spectrum estimation")
      mtext(sprintf("(with AR(%s) correction)", ordAR), line=0.4)
      Pspec_ar$spec <- Pspec_ar$spec * mARs
      plot(Pspec_ar, log="dB", col="blue", add=TRUE)
      plot(Pspec, log="dB", add=TRUE, lwd=2)
      legend("bottomleft", 
             c("original PSD",
               sprintf("AR-innovations PSD\n(mean %.01f +- %.01f dB)", dB(mARs), dB(sqrt(arvar))/4),
               "AR-corrected PSD"), 
             lwd=2, col=c("red","blue","black"))
    } else {
      plot(Pspec, log="dB", main="Pilot spectrum estimation")
    }
  }
  return(invisible(psd:::psd_envAssignGet("pilot_psd", Pspec)))
}
