#' @title Normalization of the power spectral densities.
#'
#' @description Something
#'
#' @name rlpSpec-normalization
#' @rdname rlpSpec-normalization
#' @aliases normalization
#' @keywords spectrum-estimation normalization prewhiten
#' @author A.J. Barbour <andy.barbour@@gmail.com>
#'
#' @param Spec spectrum to normalize
#' @param Fnyq nyquist frequency
#' @param verbose logical; should messages be given?
#' @param ... (unused) additional parameters
#' @return original object, with values normalized accordingly
#'
#' @seealso \code{\link{psdcore}}
NULL
 
#' @rdname rlpSpec-normalization
#' @aliases psdnorm
#' @export
psdnorm <- function(Spec, Fnyq, verbose=TRUE, ...) UseMethod("nyqnorm")
#' @rdname rlpSpec-normalization
#' @aliases psdnorm.default
#' @method psdnorm default
#' @S3method psdnorm default
psdnorm.default <- function(Spec, Fnyq, verbose=TRUE, ...){
  NULL
}
#' @rdname rlpSpec-normalization
#' @aliases psdnorm.spec
#' @method psdnorm spec
#' @S3method psdnorm spec
psdnorm.spec <- function(Spec, Fnyq, verbose=TRUE, ...){
  stopifnot(is.spec(Spec))
  src <- "spec.pgram"
  is.rlp <- exists("nyquist.normalized",where=Spec)
  if (is.rlp) src <- "rlpSpec"
  return(Spec)
}

##
## Nyquist normalization
##

#' @rdname rlpSpec-normalization
#' @aliases nyqnorm
#' @export
nyqnorm <- function(Spec, Fnyq, verbose=TRUE, ...) UseMethod("nyqnorm")
#' @rdname rlpSpec-normalization
#' @aliases nyqnorm.default
#' @method nyqnorm default
#' @S3method nyqnorm default
nyqnorm.default <- function(Spec, Fnyq, verbose=TRUE, ...){
  NULL
}
#' @rdname rlpSpec-normalization
#' @aliases nyqnorm.spec
#' @method nyqnorm spec
#' @S3method nyqnorm spec
nyqnorm.spec <- function(Spec, Fnyq, verbose=TRUE, ...){
  stopifnot(is.spec(Spec))
  src <- "spec.pgram"
  is.rlp <- exists("nyquist.normalized",where=Spec)
  if (is.rlp) src <- "rlpSpec"
  ##
  if (is.rlp & !(Spec$nyquist.normalized)){
    # is from psdcore, and hasn't been normalized
    
  } else if (!is.rlp) {
    # it's from spectrum
    Spec$freq <- Spec$freq * 2 * Fnyq
    Spec$spec <- Spec$spec * Fnyq
  } else {
    if (verbose) message("Spectrum is already normalized.")
  }
  if (verbose) message(sprintf("Normalized  %s PSD  for Nyquist  %s",src, Fnyq))
  return(Spec)
}
