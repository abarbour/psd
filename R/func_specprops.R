#' Calculate properties of multitaper power spectral density estimates
#'
#' @description
#' Various spectral properties may be computed from the vector of tapers, and
#' if necessary the sampling frequency.
#' 
#' @details
#' Parameter Details:
#' \subsection{Uncertainty}{
#' See \code{\link{spec_confint}} for details.
#' }
#'
#' \subsection{Resolution}{
#' The frequency resolution depends on the number of tapers (\eqn{K}), and
#' is found from 
#' \deqn{\frac{K \cdot f_N}{N_f}} 
#' where \eqn{f_N} is the Nyquist
#' frequency and \eqn{N_f} is the 
#' number of frequencies estimated.
#' }
#'
#' \subsection{Degrees of Freedom}{
#' There are two degrees of freedom for each taper \eqn{K}:
#' \deqn{\nu = 2 K}
#' }
#'
#' \subsection{Bandwidth}{
#' The bandwidth of a multitaper estimate depends on the number of
#' tapers.
#' Following Walden et al (1995) the effective bandwidth is \eqn{\approx 2W}
#' where
#' \deqn{W = \frac{K + 1}{2 N}} 
#(N+1)}}
#' and \eqn{N} is the number of terms in the series, which makes \eqn{N \cdot W} the
#' approximate time-bandwidth product.
#' }
#' 
#' @author A.J. Barbour
#' @name spectral_properties
#' @export
#' @seealso \code{\link{spec_confint}}, \code{\link{psd-package}}
#'
#' @param x object to calculate spectral properties for; or a vector of number of tapers
#' @param f.samp numeric; the sampling frequency (e.g. Hz) of the series the tapers are for
#' @param n.freq integer; the number of frequencies of the original spectrum 
#' (if \code{NULL} the length of the tapers object is assumed to be the number)
#' @param p numeric; the coverage probability, bound within \eqn{[0,1)}
#' @param db.ci logical; should the uncertainty confidence intervals be returned as decibels?
#' @param ... additional arguments
#' 
#' @return A list with the following properties (and names):
#' \itemize{
#' \item{\code{taper}: the number of tapers}
#' \item{\code{stderr.chi .upper, .lower, .median}: results returned from \code{\link{spec_confint}}}
#' \item{\code{resolution}: effective spectral resolution}
#' \item{\code{dof}: degrees of freedom; will be slightly inaccurate for single-taper periodograms}
#' \item{\code{bw}: effective bandwidth of the spectrum}
#' }
#'
#' @example inst/Examples/rdex_spectralproperties.R
spectral_properties <- function(x, ...) UseMethod("spectral_properties")

#' @rdname spectral_properties
#' @export
spectral_properties.spec <- function(x, ...){
  frq <- x[['freq']]
  n.freq <- length(frq)
  f.samp <- 2 * frq[n.freq]
  tapvec <- x[['taper']]
  # correct for case this is not an adaptive spectrum (then tapers is in % length)
  # this will still lead to too low and uncertainty, since d.o.f. is 1.79 instead of 2
  if (!is.amt(x)) tapvec <- ceiling(tapvec) 
  spectral_properties(tapvec, f.samp, n.freq, ...)
}

#' @rdname spectral_properties
#' @export
spectral_properties.tapers <- function(x, ...){
  spectral_properties(as.vector(x), ...)
}

#' @rdname spectral_properties
#' @export
spectral_properties.default <- function(x, f.samp=1, n.freq=NULL, p=0.95, db.ci=FALSE, ...){
  K <- as.vector(x)
  if (is.null(n.freq)) n.freq <- length(K)
  Nyquist <- f.samp/2
    ## Deg Freedom: PW93 Ch7 343
  DOF <- 2 * K
  ## Half-width
  W <- (K+1)/(2*n.freq)
  ## Effective bandwidth ~ 2 W (accurate for many spectral windows)
  BW <- 2 * W
  ## Resolution
  Resolu <- 2 * BW
  ## Uncertainty CI
  StdErrCI <- .spec_confint(DOF, p, as.db=db.ci)
  ##
  return(data.frame(taper=K, stderr.chi=StdErrCI, resolution=Resolu, dof=DOF, bw=BW))
}

#' Confidence intervals for multitaper power spectral density estimates 
#'
#' @details
#' The errors are estimated 
#' from the number of degrees of freedom \eqn{\nu} by evaluating
#' the \eqn{\chi_{p,\nu}^{2}(\nu,\nu)} distribution for an optional 
#' coverage probability \eqn{p} (defaulting to \eqn{p=0.95}).  
#' Additionally, the
#' \eqn{p=0.5} values and an approximation from \eqn{1/\sqrt{\nu - 1}}
#' are returned.
#'
#' A more 
#' sophisticated (and complicated) approach would be to
#' estimate via jack-knifing (Prieto et al 2007), but this is not yet
#' made available.
#'
#' Additive uncertainties \eqn{\delta S} are returned, such that 
#' the spectrum with confidence interval is \eqn{S \pm \delta S}.
#'
#' @author A.J. Barbour; some code modified from the \code{spec.ci} function inside \code{stats::plot.spec}
#' @name spec_confint
#' @export
#' @seealso \code{\link{spectral_properties}}, \code{\link{psd-package}}, \code{stats::plot.spec}, \code{\link{dB}}
#' 
#' @param x object to calculate spectral properties
#' @param dof numeric; the degrees of freedom \eqn{\nu}
#' @param p numeric; the coverage probability \eqn{p}, bound within \eqn{[0,1)}
#' @param as.db logical; should the values be returned as decibels?
#' @param ... additional arguments
#' 
#' @return A \code{data.frame} with the following properties (and names):
#' \itemize{
#' \item{\code{lower}: Based on upper tail probabilities (\eqn{p})}
#' \item{\code{upper}: Based on lower tail probabilities (\eqn{1-p})}
#' \item{\code{median}: Based on lower tail probabilities (\eqn{p=0.5})}
#' \item{\code{approx}: Approximation based on \eqn{1/\sqrt(\nu - 1)}.}
#' }
#' @example inst/Examples/rdex_confint.R
spec_confint <- function(x, ...) UseMethod("spec_confint")

#' @rdname spec_confint
#' @export
spec_confint.spec <- function(x, ...){
  res <- if (is.amt(x)){
    # 'tapers' object
    taps <- x[['taper']]
    if (!is.tapers(taps)) taps <- as.tapers(taps)
    taps
  } else {
    # degrees of freedom
    x[['df']]
  }
  spec_confint(res)
}

#' @rdname spec_confint
#' @export
spec_confint.tapers <- function(x, ...){
  # two degrees of freedom per taper 
  dof <- 2 * as.vector(x)
  .spec_confint(dof, ...)
}

#' @rdname spec_confint
#' @export
spec_confint.default <- function(x, ...){
    # assumes x is the number of degrees of freedom
    dof <- as.vector(x)
    .spec_confint(dof, ...)
}

#' @rdname spec_confint
#' @export
.spec_confint <- function(dof, p = 0.95, as.db=FALSE, ...) {
  # Mostly from spec.ci, lifted from plot.spec
  if (p < 0 || p >= 1) stop("coverage probability out of range [0,1)")
  ptail <- (1 - p)
  # qchisq gives distribution function
  # if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]
  upper.qp <- 1 - ptail * pchisq(dof, dof, lower.tail = FALSE)
  med.qp <-  0.5 * pchisq(dof, dof, lower.tail = TRUE)
  lower.qp <- ptail * pchisq(dof, dof, lower.tail = TRUE)
  ndof <- length(dof)
  # qchisq gives quantile function
  # spec.ci calculates with 1/(spec/dof)
  ci.ul <- 1 / (qchisq(c(upper.qp, lower.qp, med.qp), dof)/dof)
  # dS = S*ci.ul/dof
  ## heuristically tuned to approximate the median distribution of Chi^2 uncertainties
  approx <- 1 + 1/sqrt(dof-1)
  ci <- data.frame(lower=ci.ul[1:ndof], 
                   upper=ci.ul[(ndof+1):(2*ndof)], 
                   median=ci.ul[(2*ndof+1):(3*ndof)],
                   approx=approx)
  if (as.db) ci <- dB(ci)
  return(ci)
}
