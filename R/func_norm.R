#' @title Normalization of the power spectral densities.
#'
#' @description Normalization is vitally important in spectral
#' analyses.
#' @section Assumtions:
#' The normalizations performed here assume:
#'
#' @name rlpSpec-normalization
#' @rdname rlpSpec-normalization
#' @aliases normalization
#' @keywords spectrum-estimation normalization prewhiten
#' @author A.J. Barbour <andy.barbour@@gmail.com>
#'
#' @param Spec spectrum to normalize
#' @param Fsamp sampling frequency
#' @param src character string; the source of the spectrum estimator
#' @param verbose logical; should messages be given?
#' @param ... (unused) additional parameters
#' @return original object, with values normalized accordingly
#'
#' @seealso \code{\link{psdcore}} \code{\link{spectral_properties}}
NULL
 
#' @rdname rlpSpec-normalization
#' @aliases normalize
#' @export
normalize <- function(Spec, Fsamp=1, src=NULL, verbose=TRUE, ...) UseMethod("normalize")
#' @rdname rlpSpec-normalization
#' @aliases normalize.default
#' @method normalize default
#' @S3method normalize default
normalize.default <- function(Spec, Fsamp=1, src=NULL, verbose=TRUE, ...){
  .NotYetImplemented()
}
#' @rdname rlpSpec-normalization
#' @aliases normalize.list
#' @method normalize list
#' @S3method normalize list
normalize.list <- function(Spec, Fsamp=1, src=NULL, verbose=TRUE, ...){
  stopifnot(exists("freq", where=Spec) & exists("spec", where=Spec))
  class(Spec) <- "spec"
  Spec <- normalize(Spec, Fsamp, src, verbose, ...)
  class(Spec) <- "list"
  return(Spec)
}
#' @rdname rlpSpec-normalization
#' @aliases normalize.spec
#' @method normalize spec
#' @S3method normalize spec
normalize.spec <- function(Spec, Fsamp=1, src=NULL, verbose=TRUE, ...){
  stopifnot(is.spec(Spec))
  #
  if (Fsamp > 0){
    # value represents sampling frequency
    Fsamp <- Fsamp
  } else if (Fsamp < 0){
    # value is sampling interval
    Fsamp <- abs(1/Fsamp)
  } else {
    stop("bad sampling information")
  }
  #
  # assume its from spectrum
  rlp <- switch(src <- match.arg(src, c("spectrum","rlpspec")), 
                rlpspec=TRUE, spectrum=FALSE)
  if (rlp){
    ptyp <- "single"
    # spectrum is from rlpspec, and is single-sided
    Spec$spec <- Spec$spec / Fsamp
  } else {
    ptyp <- "double"
    # spectrum is from spectrum or others, double sided
    Spec$spec <- Spec$spec * 2
  }
  if (verbose) message(sprintf("Normalized  %s-sided PSD  (%s)  to sampling-freq.  %s", ptyp, src, Fsamp))
  return(invisible(Spec))
}
#