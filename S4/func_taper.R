#' The 'taper' S4 class.
#'
#' In this class the value of each position will be a non-zero, positive
#' integer, since it represents the number of tapered sections to average.
#'
# \section{Slots}{
#   \describe{
#     \item{\code{tapers}:}{Object of class \code{"integer"}, containing 
#     a vector of tapers.}
#   }
# }
#'
#' @note  The prototypical S4 class has tapers==1, and length==1.
#' Currently there are no \code{@@slots}; this may change in the future.
#'
#' @name newTaper
#' @rdname taper
#' @aliases taper-class
#' @exportClass taper
#' @author Andrew Barbour <andy.barbour@@gmail.com>
#' @examples
#' newTaper()
#' new("taper") # equivalent to newTaper()
#' print(ntap <- newTaper(1:10))
#' plot(ntap)
newTaper <- setClass("taper",
                     # if slots, add 'taper="integer",...
                     representation=representation("integer"),
                     prototype = 1L)

#' S4 from S3
#' @rdname taper-methods
#' @name print
#' @export
#' @docType methods
setMethod("print", signature("taper"), print.taper)