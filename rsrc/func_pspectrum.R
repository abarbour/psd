#' Calculate the pilot spectrum.
#'
#' The pilot spectrum calculated here is found using a fixed number
#' of tapers across all frequencies using \code{\link{psdcore}}.  Any
#' adaptive taper-refinements are based on the spectral derivatives
#' of this spectrum; hence, changes in the number of tapers can affect
#' how many adaptive stages may be needed.
#' 
#' @note A mean value and linear are removed from the series prior to spectrum
#' estimation, and no decimation is performed on the taper optimization.  
#' The taper series of the returned spectrum is constrained using
#' \code{\link{as.taper(..., minspan=TRUE)}}.
#'
#' @name pilot_spec
#' @aliases pilot_spectrum spec.pilot
#' @export
#' @author Andrew Barbour <andy.barbour@@gmail.com> based on R.L. Parker's algorithm.
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

#' Initialize a nested-list object to store the taper-optimization iterations.
#'
#' The top names are
#' \itemize{
#' \item stg_kopt.  Sequential taper vectors.
#' \item stg_psd.  Sequential power spectral density vectors.
#' \item freq. The frequencies for each set of \code{stg_kopt} and \code{stg_psd}.
#' }
#'
#' @note The object is stored as \code{'histlist'} in the \code{.rlpSpecEnv} environment. 
#' At any point the list may be accessed with \code{\link{rlp_envGet("histlist")}}.
#'
#' @name new_adapt_history
#' @export
#'
#' @param adapt_stages scalar; The number of adaptive iterations to save (excluding pilot spectrum).
#' @return the empty list
#' @seealso \code{\link{rlp_envGet}}
new_adapt_history <- function(adapt_stages){
  stopifnot(length(adapt_stages)==1)
  histlist <- vector("list", 3) # freq, list-tap, list-psd
  names(histlist) <- c("freq", "stg_kopt", "stg_psd")
  num_pos <- 1 + adapt_stages # pilot + adapts
  histlist[[2]] <- histlist[[3]] <- vector("list", adapt_stages+1)
  rlp_envAssignGet("histlist", histlist)
}

#' Update the adaptive estimation history list.
#'
#' @name update_adapt_history
#' @export
#'
#' @param stage x
#' @param ntap x
#' @param psd x
#' @param freq x
#' @return
#'
update_adapt_history <- function(stage, ntap, psd, freq=NULL){
  histlist <- rlp_envGet("histlist")
  # stage==0 <--> index==1
  stg_ind <- stage+1
  nulfrq <- is.null(freq)
  if (!nulfrq) histlist$freq <- freq
  histlist$stg_kopt[[stg_ind]] <- ntap
  histlist$stg_psd[[stg_ind]] <- psd
  if (is.null(histlist$freq) & stage>0) warning("freqs absent despite non-pilot stage update")
  rlp_envAssignGet("histlist",histlist)
}

#' Message printing utility for adapt stages
#' @param stage scalar; the current stage.
#' @return NULL
adapt_message <- function(stage){
  stopifnot(stage>=0)
  if (stage==0){stage <- paste(stage,"(pilot)")}
  message(sprintf("Stage  %s  estimation", stage))
}

#' Adaptive sine multitaper power spectral density estimation.
#' 
#' @name pspectrum
#' @export
#' @seealso \code{\link{psdcore}}
#' 
#' @param x x
#' @param x.frqsamp x
#' @param ntap_pilot x
#' @param niter x
#' @param verbose x
#' @param ... x
#' @return Object with class 'spec'.
#' 
pspectrum <- function(x, x.frqsamp, ntap_pilot=10, niter=2, verbose=TRUE, ...) UseMethod("pspectrum")
#' @rdname pspectrum
#' @S3method pspectrum default
pspectrum.default <- function(x, x.frqsamp, ntap_pilot=10, niter=2, verbose=TRUE, ...){
  for (stage in 0:niter){
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