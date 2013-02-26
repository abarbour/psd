##
## Functions for objects with class spec
##
#' @title Generic methods for objects with class \code{'spec'}.
#'
#' @details Objects with class \code{'spec'} are simply list objects. 
#' \code{as.data.frame} converts the list into
#' a \code{'data.frame'} with individual
#' columns for the frequency, PSD, and taper vectors; 
#' all other information will be retained as an attribute.
#' \code{data.frame} is an alias.
#'
#' @keywords methods S3methods spec
#' @name spec-methods
#' @author A.J. Barbour <andy.barbour@@gmail.com>
#' @rdname spec-methods
#' @docType methods
#' @param x spec object
#' @param ... optional arguments
#' @example inst/Examples/rdex_spec.R
NULL

#' @rdname spec-methods
#' @name as.data.frame.spec
#' @method as.data.frame spec
#' @S3method as.data.frame spec
as.data.frame.spec <- function(x, ...){
  # [1]  "freq"      "spec"      "coh"       "phase"     "kernel"    
  #      "df"        "bandwidth" "n.used"    "orig.n"   
  # [10] "series"    "snames"    "method"    "taper"     "pad"       
  #      "detrend"   "demean" 
  xdf <- as.data.frame(x[1:2])
  attr(xdf,"coh") <- x$coh
  attr(xdf,"phase") <- x$phase
  attr(xdf,"kernel") <- x$kernel
  attr(xdf,"df") <- x$df
  attr(xdf,"bandwidth") <- x$bandwidth
  attr(xdf,"n.used") <- x$n.used
  attr(xdf,"orig.n") <- x$orig.n
  attr(xdf,"series") <- x$series
  attr(xdf,"snames") <- x$snames
  attr(xdf,"method") <- x$method
  xdf$taper <- x$taper
  attr(xdf,"pad") <- x$pad
  attr(xdf,"detrend") <- x$detrend
  attr(xdf,"demean") <- x$demean
  return(xdf)
}
#' @rdname spec-methods
#' @name data.frame.spec
#' @method data.frame spec
#' @S3method data.frame spec
data.frame.spec <- as.data.frame.spec
