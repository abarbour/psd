#' Adaptive sine multitaper power spectral density estimation.
#' 
#'
#'
#' @name pspectrum
#' @export
#' @author Andrew Barbour <andy.barbour@@gmail.com> based on R.L. Parker's algorithm.
#' @seealso \code{\link{psdcore}}, \code{\link{riedsid}}
#' 
#' @param x vector; series to estimate PSD for.
#' @param x.frqsamp scalar; the sampling rate (e.g. Hz) of the series \code{x}.
#' @param ntap_pilot scalar; the number of sine tapers to use in the pilot spectrum estimation.
#' @param niter scalar; the number of adaptive iterations to execute after the pilot spectrum.
#' @param verbose logical; Should messages be given?
#' @param ... optional arguments (unused)
#' @return Object with class 'spec'.
#' 
pspectrum <- function(x, x.frqsamp=1, ntap_pilot=10, niter=4, verbose=TRUE, ...) UseMethod("pspectrum")
#' @rdname pspectrum
#' @S3method pspectrum default
pspectrum.default <- function(x, x.frqsamp=1, ntap_pilot=10, niter=2, verbose=TRUE, ...){
  stopifnot(length(x)>1)
  #
  adapt_message <- function(stage){
    stopifnot(stage>=0)
    if (stage==0){stage <- paste(stage,"(pilot)")}
    message(sprintf("Stage  %s  estimation", stage))
  }
  #
  for (stage in 0:abs(niter)){
    if (verbose) adapt_message(stage)
    if (stage==0){
      # --- setup the environment ---
      rlp_initEnv(refresh=TRUE,verbose=verbose)
      # --- pilot spec ---
      Pspec <- pilot_spec(x=x, x.frequency=x.frqsamp, ntap=ntap_pilot)
      # --- history ---
      save_hist <- ifelse(niter < 10, TRUE, FALSE)
      if (save_hist){
        new_adapt_history(niter)
        update_adapt_history(0, Pspec$taper, Pspec$spec, Pspec$freq)
      }
      x <- 0 # to prevent passing orig data back/forth
    }
    ## calculate optimal tapers
    kopt <- riedsid(Pspec$spec, Pspec$taper)
    stopifnot(exists('kopt'))
    ## reapply to spectrum
    Pspec <- psdcore(X.d=x, ntaper=kopt)
    ## update history
    if (save_hist) update_adapt_history(stage, Pspec$taper, Pspec$spec)
  }
  return(Pspec)
}

#' Calculate the pilot power spectral densities.
#'
#' This PSD -- the pilot spectrum -- is used as the figurative starting point
#' for the adaptive estimation routine.
#' A fixed number
#' of tapers is applied across all frequencies using \code{\link{psdcore}}, and
#' subsequent taper-refinements are based on the spectral derivatives
#' of this spectrum; hence, changes in the number of tapers can affect
#' how many adaptive stages may be needed (though there are no formal convergence
#' criteria to speak of).
#' 
#' @note A mean value and linear trend are removed from the series prior to
#' estimation, and no decimation is performed  
#' The taper series of the returned spectrum is constrained using
#' \code{as.taper(..., minspan=TRUE)}.
#'
#' @name pilot_spec
#' @aliases pilot_spectrum spec.pilot
#' @export
#' @author Andrew Barbour <andy.barbour@@gmail.com>
#' @seealso \code{\link{psdcore}}
#'
#' @param x vetor; the data series to find a pilot spectrum for
#' @param x.frequency scalar; the sampling frequency (e.g. Hz) of the series
#' @param ntap scalar; the number of tapers to apply during spectum estimation
#' @param ... additional parameters
#' @return An object with class 'spec'
#'
pilot_spec <- function(...) UseMethod("pilot_spec")
#' @rdname pilot_spec
#' @S3method pilot_spec default
pilot_spec.default <- function(x, x.frequency=1, ntap=10, ...){
  stopifnot(length(ntap)==1)
  if (is.ts(x)) x.frequency <- stats::frequency(x)
  # initial spectrum:
  Pspec <- psdcore(X.d=x, X.frq=x.frequency, ntaper=ntap, 
                   ndecimate=1L, demean=TRUE, detrend=TRUE, as.spec=TRUE) #, ...)
  num_frq <- length(Pspec$freq)
  num_tap <- length(Pspec$taper)
  stopifnot(num_tap <= num_frq)
  Ptap <- Pspec$taper
  # generate a series, if necessary
  if (num_tap < num_frq) Ptap <- rep.int(Ptap[1], num_frq)
  # return taper object
  Pspec$taper <- as.taper(Ptap, setspan=TRUE)
  #
  return(Pspec)
}