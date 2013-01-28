#' Adaptive sine multitaper power spectral density estimation.
#' 
#' This is the primary function to be used in this package, and returns
#' power spectral density estimates where the number of tapers at each
#' frequency has been iteratively optimized (\code{niter} times).
#'
#' See the \strong{Adaptive estimation} section in the description of
#' the \code{\link{rlpSpec-package}} for details regarding adaptive estimation.
#'
#' @name pspectrum
#' @export
#' @author A.J. Barbour <andy.barbour@@gmail.com> adapted original by R.L. Parker.
#' @seealso \code{\link{psdcore}}, \code{\link{riedsid}}, \code{\link{rlpSpec-package}}
#' @keywords spectrum-estimation riedel-sidorenko tapers tapers-constraints tapers-weighting numerical-derivative
#' 
#' @param x vector; series to estimate PSD for.
#' @param x.frqsamp scalar; the sampling rate (e.g. Hz) of the series \code{x}.
#' @param ntap_pilot scalar; the number of sine tapers to use in the pilot spectrum estimation.
#' @param niter scalar; the number of adaptive iterations to execute after the pilot spectrum.
#' @param verbose logical; Should messages be given?
#' @param no.history logical; Should the adaptive history \emph{not} be saved?
#' @param plot logical; Should the results be plotted?
#' @param ... Optional parameters passed to \code{\link{riedsid}}
#' @return Object with class 'spec', invisibly.
#' @example x_examp/pspec.R
pspectrum <- function(x, x.frqsamp=1, ntap_pilot=5, niter=3, verbose=TRUE, no.history=FALSE, plot=TRUE, ...) UseMethod("pspectrum")
#' @rdname pspectrum
#' @method pspectrum default
#' @S3method pspectrum default
pspectrum.default <- function(x, x.frqsamp=1, ntap_pilot=5, niter=3, verbose=TRUE, no.history=FALSE, plot=TRUE, ...){
  stopifnot(length(x)>1)
  xo <- rlpSpec:::rlp_envAssignGet("original_series", x)
  #
  adapt_message <- function(stage, dvar=NULL){
    stopifnot(stage>=0)
    if (stage==0){
      stage <- paste(stage,"est. (pilot)")
    } else {
      stage <- paste(stage, sprintf("est. (A.V.R. %.01f dB)", dB(dvar)))
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
      rlp_initEnv(refresh=TRUE, verbose=verbose)
      # --- pilot spec ---
      # normalization is there
      Pspec <- pilot_spec(x=xo, x.frequency=x.frqsamp, ntap=ntap_pilot)
      dvar.o <- vardiff(Pspec$spec)
      # --- history ---
      save_hist <- ifelse(niter < 10, TRUE, FALSE)
      if (no.history) save_hist <- FALSE
      if (save_hist){
        new_adapt_history(niter)
        update_adapt_history(0, Pspec$taper, Pspec$spec, Pspec$freq)
      }
      xo <- 0 # to prevent passing orig data back/forth
    }
    rverb <- ifelse(stage>0, FALSE, TRUE)
    ## calculate optimal tapers
    kopt <- riedsid(Pspec, verbose=rverb, ...)
    stopifnot(exists('kopt'))
    ## reapply to spectrum
    if (stage==niter & plot){
      plotpsd_ <- TRUE
      xo <- x
      rm(x)
    }
    ##print(plotpsd_)
    Pspec <- psdcore(X.d=xo, X.frq=x.frqsamp, ntaper=kopt, plotpsd=plotpsd_, verbose=FALSE)
    if (verbose) if (verbose) adapt_message(stage, vardiff(Pspec$spec)/dvar.o)
    ## update history
    if (save_hist) update_adapt_history(stage, Pspec$taper, Pspec$spec)
  }
  return(invisible(Pspec))
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
#' A mean value and linear trend are removed from the series prior to
#' estimation, and no decimation is performed.
#'
#' The taper series of the returned spectrum is constrained using
#' \code{as.tapers(..., minspan=TRUE)}.
#'
#' @name pilot_spec
#' @aliases pilot_spectrum spec.pilot
#' @export
#' @author A.J. Barbour <andy.barbour@@gmail.com>
#' @seealso \code{\link{psdcore}}, \code{\link{as.tapers}}, \code{\link{minspan}}
#'
#' @param x vetor; the data series to find a pilot spectrum for
#' @param x.frequency scalar; the sampling frequency (e.g. Hz) of the series
#' @param ntap scalar; the number of tapers to apply during spectrum estimation
#' @param ... additional parameters passed to \code{\link{psdcore}}
#' @return An object with class 'spec', invisibly.  It also assigns it to
#' \code{"final_psd"} in the working environment.
pilot_spec <- function(x, x.frequency=1, ntap=5, ...) UseMethod("pilot_spec")
#' @rdname pilot_spec
#' @method pilot_spec default
#' @S3method pilot_spec default
pilot_spec.default <- function(x, x.frequency=1, ntap=5, ...){
  stopifnot(length(ntap)==1)
  if (is.ts(x)) x.frequency <- stats::frequency(x)
  # initial spectrum:
  Pspec <- psdcore(X.d=x, X.frq=x.frequency, ntaper=ntap, 
                   ndecimate=1L, 
                   demean=TRUE, detrend=TRUE, 
                   first.last=TRUE,
                   Nyquist.normalize=TRUE,
                   as.spec=TRUE, 
                   refresh=TRUE, ...)
  num_frq <- length(Pspec$freq)
  num_tap <- length(Pspec$taper)
  stopifnot(num_tap <= num_frq)
  Ptap <- Pspec$taper
  # generate a series, if necessary
  if (num_tap < num_frq) Ptap <- rep.int(Ptap[1], num_frq)
  # return tapers object
  Pspec$taper <- as.tapers(Ptap, setspan=TRUE)
  #
  return(invisible(rlpSpec:::rlp_envAssignGet("final_psd", Pspec)))
}
