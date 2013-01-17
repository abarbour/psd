#' Riedel & Sidorenko taper optimization
#' 
#' Estimates optimal number of tapers at each frequency of
#' given psd, based on Riedel-Sidorenko MSE recipe, and other
#' tweaks by RLP.
#' 
#' @title riedsid
#' @export
#' @keywords taper taper-constraints riedel-sidorenko
#' @author Andrew Barbour <andy.barbour@@gmail.com> ported original by R.L.Parker.
#' 
#' @param psd vector; the spectral values used to optimize taper numbers
#' @param pspec object with class 'spec'
#' @param ntaper scalar or vector; number of tapers to apply optimization
#' @param tapseq vector; representing positions or frequencies (same length as psd)
#' @param c.method string; constraint method to use
#' @param ... optional argments passed to \code{constrain_tapers}
#' 
#' @seealso \code{\link{psdcore}}, \code{\link{constrain_tapers}}
#' @examples
#' riedsid(rnorm(10), 10)
#' riedsid(rnorm(10), c(1:5,5:1))
riedsid <- function(psd, ntaper, tapseq=NULL, c.method=NULL, ...) UseMethod("riedsid")

#' @rdname riedsid
#' @S3method riedsid spec
riedsid.spec <- function(pspec, ...){
  psd <- pspec$spec
  ntap <- pspec$taper
  riedsid(psd, ntap, ...)
  #.NotYetImplemented()
}

#' @rdname riedsid
#' @S3method riedsid default
riedsid.default <- function(psd, ntaper, tapseq=NULL, c.method=NULL, ...) {
  ## spectral values
  psd <- as.vector(psd)
  # num freqs
  nf <- rlp_envAssignGet("num_freqs", length(psd))
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
  
  # Set the number of tapers to within the range: 1/2 nf, 7/5 ntap
  # rowMins produces a rowvec of rowwise minimums; convert to colvec
  nspan <- minspan(ntap, nf)
  # The spectral gradients should be in log-space, so
  # create a log spec, and pad to handle begnning and end values
  nadd <- 1 + max(nspan)
  Y <- c(psd[nadd:2], psd, psd[(nf-1):(nf-nadd)])
  Y[Y <= 0] <- eps
  lY <- log(Y)
  dY <- d2Y <- Zeros
  #
  if (is.null(tapseq) | (length(tapseq) != length(psd))){
    #kseq <- seq.int(from=1, to=nf, by=1)
    xfreq <- frequency(psd) #1 # frequency(x) <==> sps, Hz
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
  #  R <- psd"/psd <- Y" + (Y')^2  2nd form preferred for consistent smoothing
  #
  #  Riedel-Sidorenko recipe (eq 13): 
  #       kopt <- (12*abs(psd ./ d2psd)).^0.4 
  #  but parabolic weighting in psdcore requires: 
  #               (480)^0.2*abs(psd ./ d2psd).^0.4
  #  Original form:  kopt <- 3.428*abs(psd ./ d2psd).^0.4
  #
  # the optimal number of tapers (in an MSE sense):
  kopt <- as.taper( 3.437544 / abs(eps + d2Y + dY*dY) ^ 0.4 )
  kopt.bound <- constrain_tapers(kopt, kseq, c.method, ...)
  ##
  return(kopt.bound)
} 
# end riedsid.default
