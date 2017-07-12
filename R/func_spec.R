#' Generic methods for objects with class \code{'spec'}
#'
#' @details Objects with class \code{'spec'} are simply lists with spectral estimates and parameters 
#' \code{as.data.frame} converts the list into a \code{'data.frame'} with individual
#' columns for the frequency, PSD, and taper vectors; 
#' all other information will be retained as a list in the attributes.
#'
#' @name spec-methods
#' @author A.J. Barbour
#' @rdname spec-methods
#' @docType methods
#' 
#' @param x a \code{'spec'} object
#' @param y optional coordinate vector for the y-axis
#' @param type character; the type of plot
#' @param ... optional arguments
#' 
#' @example inst/Examples/rdex_spec.R
NULL

#' @rdname spec-methods
#' @aliases lines.spec
#' @export
lines.spec <- function(x, y=NULL, type = 'l', ...){
  plot(x, add=TRUE, type=type, ...)
}

#' @rdname spec-methods
#' @export
spec_details <- function(x, ...){
  x[-which(names(x) %in% c('freq','spec','taper'))]
}

#' @rdname spec-methods
#' @aliases as.data.frame.spec
#' @export
as.data.frame.spec <- function(x, ...){
  # get the meat-n-potatoes
  xdf <- data.frame(freq=x[['freq']], spec=x[['spec']], taper = x[['taper']])
  # keep everything but the mean-n-potatoes as an attribute
  # trouble is attributes are not always retained
  attr(xdf, 'spec.details') <- spec_details(x)
  return(xdf)
}

#' @rdname spec-methods
#' @export
as.matrix.spec <- function(x, ...){
  as.matrix(as.data.frame.spec(x), ...)
}

#' @rdname spec-methods
#' @export
as.list.spec <- function(x, ...){
  return(unclass(x))
}
