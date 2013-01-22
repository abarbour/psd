#' Constrained, optimal tapers using the Riedel & Sidorenko--Parker method.
#' 
#' Estimates the
#' optimal number of tapers at each frequency of
#' given PSD, using a modified Riedel-Sidorenko MSE recipe (RS-RLP).
#' 
#' The optimization is as follows. First, weighted derivatives of the 
#' input PSD are computed.
#' Using those derivates the optimal number of tapers is found through the 
#' RS-RLP formulation.
#' Constraints are then placed on the practicable number of tapers.
#'
#' Set \code{constrained=FALSE} to turn off taper constraints.
#'
#' @export
#' @keywords taper taper-constraints riedel-sidorenko
#' @author A.J. Barbour <andy.barbour@@gmail.com> adapted original by R.L. Parker.
#' 
#' @param psd vector or class 'spec'; the spectral values used to optimize taper numbers
#' @param ntaper scalar or vector; number of tapers to apply optimization
#' @param tapseq vector; representing positions or frequencies (same length as psd)
#' @param constrained logical; should the taper constraints be applied to the optimum tapers?
#' @param c.method string; constraint method to use if \code{constrained=TRUE}
#' @param ... optional argments passed to \code{\link{constrain_tapers}}
#' @return Object with class 'taper'.
#' 
#' @seealso \code{\link{constrain_tapers}}, \code{\link{psdcore}}
#' @example x_examp/ried.R
riedsid <- function(psd, ntaper, tapseq=NULL, constrained=TRUE, c.method=NULL, ...) UseMethod("riedsid")

#' @rdname riedsid
#' @S3method riedsid spec
riedsid.spec <- function(psd, ntaper=psd$taper, tapseq=NULL, constrained=TRUE, c.method=NULL, ...){
  stopifnot(is.spec(psd))
  psd <- psd$spec
  riedsid(psd, ntaper, ...)
  #.NotYetImplemented()
}

#' @rdname riedsid
#' @S3method riedsid default
riedsid.default <- function(psd, ntaper, tapseq=NULL, constrained=TRUE, c.method=NULL, ...) {
  ## spectral values
  psd <- as.vector(psd)
  # num freqs
  nf <- rlp_envAssignGet("num_freqs", length(psd))
  # prelims
  eps <- .Machine$double.eps**2 # was: 1e-78  #  A small number to protect against zeros
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
  DFUN <- function(j, 
                   j1=j-nspan[j]+nadd-1, 
                   j2=j+nspan[j]+nadd-1, 
                   jr=j1:j2, 
                   logY=lY[jr], 
                   dEps=eps){
    u <- jr - (j1 + j2)/2 # rowvec 
    u2 <- u*u             # rowvec
    L <- j2-j1+1          # constant
    L2 <- L*L             # constant
    LL2 <- L*L2           # constant
    LL2L <- LL2 - L       # constant
    uzero <- (L2 - 1)/12  # constant
    # first deriv
    dY <- u %*% logY * 12 / LL2L
    # second deriv
    d2Y <- (u2 - uzero) %*% logY * 360 / LL2L / (L2-4)
    return(c(fdY2=dY*dY, fd2Y=d2Y, fdEps=dEps))
  }
  DX <- 1:nf
  RStmp <- rlp_envAssignGet("spectral_derivatives", vapply(X=DX, FUN=DFUN, FUN.VALUE=c(1,1,1)))
  #[1,] fdY2
  #[2,] fd2Y
  #[3,] fdEps
  kopt <- as.taper( 3.437544 / abs(colSums(RStmp)) ** 0.4 )
  rm(RStmp)
  
  #
  #   for (  j  in  1:nf ) {
  #     j1 <- j - nspan[j] + nadd - 1
  #     j2 <- j + nspan[j] + nadd - 1
  #     #  Over an interval proportional to taper length, fit a least
  #     #  squares quadratic to Y to estimate smoothed 1st, 2nd derivs
  #     jr <- j1:j2           # rowvec
  #     u <- jr - (j1 + j2)/2 # rowvec 
  #     u2 <- u*u             # rowvec
  #     L <- j2-j1+1          # constant
  #     L2 <- L*L             # constant
  #     LL2 <- L*L2           # constant
  #     LL2L <- LL2 - L       # constant
  #     uzero <- (L2 - 1)/12  # constant
  #     # first deriv
  #     dY[j] <- u %*% lY[jr] * 12 / LL2L
  #     # second deriv (?)
  #     d2Y[j] <- (u2 - uzero) %*% lY[jr] * 360 / LL2L / (L2-4)
  #   }
  #   #
  #   #  R <- psd"/psd <- Y" + (Y')^2  2nd form preferred for consistent smoothing
  #   #
  #   #  Riedel-Sidorenko recipe (eq 13): 
  #   #       kopt <- (12*abs(psd ./ d2psd)).^0.4 
  #   #  but parabolic weighting in psdcore requires: 
  #   #               (480)^0.2*abs(psd ./ d2psd).^0.4
  #   #  Original form:  kopt <- 3.428*abs(psd ./ d2psd).^0.4
  #   #
  #   # the optimal number of tapers (in an MSE sense):
  #   kopt_old <- as.taper( 3.437544 / abs(eps + d2Y + dY*dY) ^ 0.4 )
  #   #
  #print(all.equal(kopt_old,kopt)) # TRUE!
  ##
  ## Constrain tapers
  if (constrained) kopt <- constrain_tapers(kopt, kseq, c.method, ...)
  ##
  return(kopt)
} 
# end riedsid.default
