#' This package provides tools to
#' 
#' The main function to be used 
#' is \code{\link{some_func}}
#' 
#' There are also two helper functions included: 
# \describe{
# \item{\code{\link{some_other_func}}}{ to do something.}
# }
#'
#' @section Scientific background:
#'
#' A bunch of stuff, and inline equation \eqn{r}, and a newline equation:
#' \deqn{
#' \frac{\partial^2 s}{\partial r^2}= 0
#' }
#' and some more.
#' 
#' And more.
#'
#' @docType package
#' @name rlpSpec-package
#' @aliases rlpSpec, rlpspec-package, spec.rlp
#' @title Adaptively estimate power spectral densities of an optimally tapered series.
#' 
#' @author Andrew Barbour <andy.barbour@@gmail.com> 
#' 
#' @import matrixStats
#' @import RColorBrewer
#' @import signal
#' @useDynLib rlpSpec
#' 
#' @references Riedel, K. S., & Sidorenko, A. (1995), 
#' Minimum bias multiple taper spectral estimation,
#' \emph{Signal Processing, IEEE Transactions on}, \strong{43}(1), 188--195.
#
#' @references Riedel, K. S. (1996),
#' Adaptive smoothing of the log-spectrum with multiple tapering,
#' \emph{Signal Processing, IEEE Transactions on}, \strong{44}(7), 1794--1800.
#'
#' @references Walden, A. T. WALDEN, and  E. J. McCoy, and D. B. Percival (1995),
#' The effective bandwidth of a multitaper spectral estimator,
#' \emph{Biometrika}, \strong{82}(1), 201--214,
#' \url{http://biomet.oxfordjournals.org/content/82/1/201}
#'
# problem with %7E or # at end? (BOTH)
#' @references \url{http://igppweb.ucsd.edu/\%7Eparker/Software/\#PSD}
#'
# @seealso \code{\link{pspectrum}}, \code{\link{psdcore}}
#'  
NULL
