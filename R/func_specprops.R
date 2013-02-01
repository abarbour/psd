#' Calculate spectral properties such as standard error and resolution.
#'
#' Various spectral properties may be computed from the vector of tapers, and
#' if necessary the sampling frequency.
#'
#' @section Parameter Details:
#' \subsection{Uncertainty}{
#' The errors are estimated in the simplest way, 
#' from the number of degrees of freedom; a more 
#' sophisticated (and complicated) approach is to
#' estimate via jack-knifing (Prieto et al 2007)
#' which is not yet available.
#'
#' Here the standard error \eqn{\delta S} is returned so \eqn{\delta S \cdot S} 
#' represents spectral uncertainty.
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
#' There are two degrees of freedom for each taper.
#' }
#'
#' \subsection{Bandwidth}{
#' The bandwidth of a multitaper estimate depends on the number of
#' tapers.
#' Following Walden et al (1995) the effective bandwidth is \eqn{\approx 2W}
#' where
#' \deqn{W = \frac{K + 1}{2N}} 
#(N+1)}}
#' and \eqn{N} is the number of terms in the series, which makes \eqn{N \cdot W} the
#' approximate time-bandwidth product.
#' }
#' 
#' @author A.J. Barbour <andy.barbour@@gmail.com>
#' @name spectral_properties
#' @export
#' @keywords properties tapers resolution uncertainty degrees-of-freedom bandwidth
#' @seealso \code{\link{spec_confint}}, \code{\link{rlpSpec-package}}
#'
#' @param tapvec object with class \code{tapers} or \code{spec}
#' @param f.samp scalar; the sampling frequency (e.g. Hz) of the series the tapers are for
#' @param n.freq scalar; the number of frequencies of the original spectrum (if \code{NULL} the length of the tapers object is assumed to be the number)
#' @param p numeric; the coverage probability, bound within \eqn{[0,1)}
#' @param  db.ci logical; should the uncertainty confidence intervals be returned as decibels?
#' @param ... additional arguments (unused)
#' @return A list with the following properties (and names):
#' \itemize{
#' \item{\code{taper}: The original taper vector.}
#' \item{\code{stderr}: The standard error of the spectrum.}
#' \item{\code{stderr.chi, .upper, .lower, .median}:A data frame with the results from \code{\link{spec_confint}}.}
#' \item{\code{resolution}: The effective spectral resolution.}
#' \item{\code{dof}: The number of degrees of freedom.}
#' \item{\code{bw}: The effective bandwidth of the spectrum.}
#' }
spectral_properties <- function(tapvec, f.samp=1, n.freq=NULL, p=0.95, db.ci=FALSE, ...) UseMethod("spectral_properties")
#' @rdname spectral_properties
#' @aliases spectral_properties.spec
#' @method spectral_properties spec
#' @S3method spectral_properties spec
spectral_properties.spec <- function(tapvec, ...){
  stopifnot(is.spec(Pspec <- tapvec))
  n.freq <- length(Pspec$freq)
  f.samp <- 2*Pspec$freq[n.freq]
  tapvec <- Pspec$taper
  if (!is.tapers(tapvec)) tapvec <- as.tapers(tapvec)
  spectral_properties(tapvec, f.samp, n.freq, ...)
}
#' @rdname spectral_properties
#' @aliases spectral_properties.tapers
#' @method spectral_properties tapers
#' @S3method spectral_properties tapers
spectral_properties.tapers <- function(tapvec, f.samp=1, n.freq=NULL, p=0.95, db.ci=FALSE, ...){
  stopifnot(is.tapers(tapvec))
  K <- unclass(tapvec)
  Nyquist <- f.samp/2
  if (is.null(n.freq)) n.freq <- length(tapvec)
  #Var <- 10 / K / 12
  ## Deg Freedom: PW93 Ch7 343
  Dof <- 2 * K
  #
  ##BW <- K / n.freq
  ##Resolu <- BW * Nyquist
  #
  ## Bandwidth
  # Walden et al
  # half-width W = (K + 1)/{2(N + 1)}
  # effective bandwidth ~ 2 W (accurate for many spectral windows)
  W <- (K+1)/(2*n.freq)
  BW <- 2 * W
  ## Resolution
  Resolu <- 2 * BW
  ## Uncertainty CI -- how does it compare to StdErr
  StdErrCI <- spec_confint(Dof, p, as.db=db.ci)
  ##
  return(data.frame(taper=K, stderr.chi=StdErrCI, resolution=Resolu, dof=Dof, bw=BW))
}

#' Calculate confidence intervals for a given number of degrees of freedom.
#'
#' In a multitaper scheme, the degrees of freedom is two per taper.
#'
#' @author A.J. Barbour <andy.barbour@@gmail.com> modified from the 
#' \code{spec.ci} function inside
#' \code{stats::plot.spec}.
#' @name spec_confint
#' @export
#' @seealso \code{\link{spectral_properties}}, \code{\link{rlpSpec-package}}, \code{plot.spec}
#' @param dof numeric; the degrees of freedom
#' @param p numeric; the coverage probability, bound within \eqn{[0,1)}
#' @param as.db logical; should the values be returned as decibels?
#' @return A \code{data.frame} with the following properties (and names):
#' \itemize{
#' \item{\code{upper}: The original taper vector.}
#' \item{\code{lower}: The standard error of the spectrum.}
#' \item{\code{median}: The standard error of the spectrum from \code{\link{spec_confint}}.}
#' }
#' @keywords properties tapers uncertainty degrees-of-freedom
spec_confint <- function(dof, p = 0.95, as.db=FALSE) UseMethod("spec_confint")
#' @rdname spec_confint
#' @aliases spec_confint.spec
#' @method spec_confint spec
#' @S3method spec_confint spec
spec_confint.spec <- function(dof, p = 0.95, as.db=FALSE){
  stopifnot(is.spec(dof))
  dof <- dof$df
  spec_confint(dof, p)
}
#' @rdname spec_confint
#' @aliases spec_confint.tapers
#' @method spec_confint tapers
#' @S3method spec_confint tapers
spec_confint.tapers <- function(dof, p = 0.95, as.db=FALSE){
  stopifnot(is.tapers(dof))
  # two degrees of freedom per taper 
  dof <- 2 * unclass(dof)
  spec_confint(dof, p)
}
#' @rdname spec_confint
#' @aliases spec_confint.default
#' @method spec_confint default
#' @S3method spec_confint default
spec_confint.default <- function(dof, p = 0.95, as.db=FALSE) {
  # Mostly from spec.ci, lifted from plot.spec
  if (p < 0 || p >= 1) stop("coverage probability out of range [0,1)")
  ptail <- (1 - p)
  # qchisq gives distribution function
  # if TRUE (default), probabilities are P[X ≤ x], otherwise, P[X > x]
  upper.qp <- 1 - ptail * pchisq(dof, dof, lower.tail = FALSE)
  med.qp <-  0.5 * pchisq(dof, dof, lower.tail = TRUE)
  lower.qp <- ptail * pchisq(dof, dof, lower.tail = TRUE)
  ndof <- length(dof)
  # qchisq gives quantile function
  # spec.ci calculates with 1/(spec/dof)
  ci.ul <- 1 / (qchisq(c(upper.qp, lower.qp, med.qp), dof)/dof)
  # dS = S*ci.ul/dof
  ## heuristically tuned to approximate the median distribution of Chi^2 uncertainties
  approx <- 1+1/sqrt(dof-1)
  ci <- data.frame(upper=ci.ul[1:ndof], 
                   lower=ci.ul[(ndof+1):(2*ndof)], 
                   median=ci.ul[(2*ndof+1):(3*ndof)],
                   approx=approx)
  if (as.db) ci <- dB(ci)
  return(ci)
}