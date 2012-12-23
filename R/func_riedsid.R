#' Riedel & Sidorenko taper optimization
#' 
#' Estimates optimal number of tapers at each frequency of
#' given spec, based on Riedel-Sidorenko MSE recipe, and other
#' tweaks by RLP.
#' 
#' @title riedsid
#' @export
#' 
#' @param spec vector; the spectral values used to optimize taper numbers
#' @param ntaper scalar or vector; number of tapers to apply optimization
#' @param tapseq vector; representing positions or frequencies (same length as spec)
#' @param c.method string; constraint method to use
#' @param ... optional argments passed to \code{constrain_tapers}
#' 
#' @seealso \code{\link{psdcore}}, \code{\link{constrain_tapers}}
# @example
# riedsid(rnorm(10), 10)
# riedsid(rnorm(10), c(1:5,5:1))
riedsid <- function(spec, ntaper, tapseq=NULL, c.method=NULL, ...) UseMethod("riedsid")

#' @rdname riedsid
#' @S3method riedsid default
#' @return \code{NULL}
riedsid.default <- function(spec, ntaper, tapseq=NULL, c.method=NULL, ...) {
  ## spectral values
  spec <- as.vector(spec)
  # num freqs
  nf <- envAssignGet("num_freqs", length(spec))
  # prelims
  eps <- .Machine$double.eps # was: 1e-78  #  A small number to protect against zeros
  Ones <- ones(nf) # was rowvec, now col 
  Zeros <- zeros(nf) # was rowvec, now col
  # vectorize initial estimate
  if (length(ntaper)==1){
    ntap <- ntaper*Ones
  } else {
    ntap <- ntaper
  }
  if (!is.taper(ntap)) ntap <- as.taper(ntap)
  #plot(ntap)
  # find the minimum by column for 1/2 nf, 7/5 ntap
  # rowMins produces a rowvec of rowwise minimums; convert to colvec
  nspan <- minspan(ntap, nf)
  # The spectral gradients should be in log-space, so
  # create a log spec, and pad to handle begnning and end values
  nadd <- 1 + max(nspan)
  Y <- c(spec[nadd:2], spec, spec[(nf-1):(nf-nadd)])
  Y[Y <= 0] <- eps
  lY <- log(Y)
  dY <- d2Y <- Zeros
  #
  if (is.null(tapseq) | (length(tapseq) != length(spec))){
    #kseq <- seq.int(from=1, to=nf, by=1)
    xfreq <- frequency(spec) #1 # frequency(x) <==> sps, Hz
    Nspec <- floor(nf/2)
    kseq <- seq.int(from = xfreq/nf, by = xfreq/nf, length.out = Nspec)
  } else {
    kseq <- tapseq # sort?
  }
  #
  # Smooth spectral derivatives
  #
  for (  j  in  1:nf ) {
    j1 <- j - nspan[j] + nadd - 1
    j2 <- j + nspan[j] + nadd - 1
    #  Over an interval proportional to taper length, fit a least
    #  squares quadratic to Y to estimate smoothed 1st, 2nd derivs
    jr <- j1:j2           # rowvec
    u <- jr - (j1 + j2)/2 # rowvec 
    u2 <- u*u             # rowvec
    L <- j2-j1+1          # constant
    L2 <- L*L             # constant
    LL2 <- L*L2           # constant
    LL2L <- LL2 - L       # constant
    uzero <- (L2 - 1)/12  # constant
    # first deriv
    dY[j] <- u %*% lY[jr] * 12 / LL2L
    # second deriv (?)
    d2Y[j] <- (u2 - uzero) %*% lY[jr] * 360 / LL2L / (L2-4)
  }
  #
  #  R <- spec"/spec <- Y" + (Y')^2  2nd form preferred for consistent smoothing
  #
  #  Riedel-Sidorenko recipe (eq 13): 
  #       kopt <- (12*abs(spec ./ d2spec)).^0.4 
  #  but parabolic weighting in psdcore requires: 
  #               (480)^0.2*abs(spec ./ d2spec).^0.4
  #  Original form:  kopt <- 3.428*abs(spec ./ d2spec).^0.4
  #
  # the optimal number of tapers (in an MSE sense):
  kopt <- as.taper( 3.437544 / abs(eps + d2Y + dY*dY) ^ 0.4 )
  kopt.bound <- constrain_tapers(kopt, kseq, c.method, ...)
  ##
  return(kopt.bound)
} 
# end riedsid.default

#' @description Constrain the number of tapers applied
#' 
#' @details Refines the number of tapers; the method which it does this
#' may be chosen by the user.
#' 
#' @title constrain_tapers
#' @aliases constrain_tapers
#' @export
#' @seealso \code{\link{riedsid}}, \code{\link{ctap_simple}}
#' @param tapvec vector; the number of tapers at each frequency
#' @param tapseq vector; positions or frequencies -- necessary for smoother methods
#' @param constraint.method method to use for constraints on tapers numbers
#' @param min.tapers integer; the minimum number of tapers
#' @param verbose boolean; should warnings and messages be given?
#' @param ... optional argments passed to the constraint.method
constrain_tapers <- function(tapvec, 
                             tapseq=NULL,
                             constraint.method=c("simple.slope",
                                                 "markov.chain",
                                                 "loess.smooth",
                                                 "friedman.smooth",
                                                 "none"),
                             min.tapers=1,
                             verbose=TRUE, 
                             ...){
  stopifnot(is.taper(tapvec))
  stopifnot(min.tapers>=0)
  # choose the appropriate method to apply taper constraints
  cmeth <- match.arg(constraint.method) 
  if (cmeth=="none"){
    if (verbose) warning("no taper optimization constraints applied")
    tapvec.adj <- tapvec
  } else {
    if (verbose) message(sprintf("Constraining tapers with ** %s ** method",cmeth))
    CTAPFUN <- switch(cmeth,
                      "simple.slope"=ctap_simple,
                      "markov.chain"=ctap_markov,
                      "loess.smooth"=ctap_loess,
                      "friedman.smooth"=ctap_friedman)
    tapvec.adj <- CTAPFUN(tapvec, tapseq, ...)
  }
  # MAX/MIN bounds
  # set the maximim tapers: Never average over more than the length of the spectrum!
  maxtap <- round(length(tapvec)/2)
  # ensure the minimum is below
  mintap <- min(min.tapers, maxtap)
  # sort for posterity
  tapbounds <- sort(c(mintap, maxtap))
  mintap <- tapbounds[1]
  maxtap <- tapbounds[2]
  tapvec.adj[tapvec.adj < mintap] <- mintap
  tapvec.adj[tapvec.adj > maxtap] <- maxtap
  return(tapvec.adj)
}
