#' Constrained, optimal tapers using the Riedel & Sidorenko--Parker method.
#' 
#' Estimates the
#' optimal number of tapers at each frequency of
#' given PSD, using a modified Riedel-Sidorenko MSE recipe (RS-RLP).
#' 
#' @details
#' The optimization is as follows. First, weighted derivatives of the 
#' input PSD are computed.
#' Using those derivates the optimal number of tapers is found through the 
#' RS-RLP formulation.
#' Constraints are then placed on the practicable number of tapers.
#'
#' \subsection{Taper constraints}{
#' The parameter \code{c.method} provides an option to change the method
#' of taper constraints.  A description of each may be found in 
#' the documentation for \code{\link{constrain_tapers}}.
#'
#' Once can use \code{constrained=FALSE} to turn off all taper constraints; this
#' could lead to strange behavior though.
#' }
#'
#' \subsection{Spectral derivatives}{
#' The parameter \code{Deriv.method} determines which method is used
#' to estimate derivatives.
#' \itemize{
#' \item{\code{"local_qls"}}{ (\strong{default}) uses quadratic weighting and
#' local least-squares estimation; 
#' then, \code{Local.loss} can alter slightly the weighting to make the derivatives more
#' or less succeptible to changes in spectral values. Can be slower than \code{"spg"}.}
#' \item{\code{"spg"}}{ uses \code{\link{splineGrad}}; then, additional arguments
#' may be passed to
#' control the smoothness of the derivatives
#' (e.g \code{spar} in \code{smooth.spline}).}
#' }
#' }
#'
#' @section Warning:
#' The \code{"spg"} can become numerically unstable, and it's not clean when it will
#' be preferred to the \code{"local_qls"} method other than for efficiency's sake.
#'
#' @export
#' @keywords tapers tapers-constraints riedel-sidorenko
#' @author A.J. Barbour <andy.barbour@@gmail.com> adapted original by R.L. Parker.
#' 
#' @param PSD vector or class 'spec'; the spectral values used to optimize taper numbers
#' @param ntaper scalar or vector; number of tapers to apply optimization
#' @param tapseq vector; representing positions or frequencies (same length as \code{PSD})
#' @param Deriv.method character string; choice of gradient estimation method 
#' @param Local.loss string; sets how sensitive the spectral derivatives are
#' @param constrained logical; should the taper constraints be applied to the optimum tapers?
#' @param c.method string; constraint method to use if \code{constrained=TRUE}
#' @param verbose logical; should messages be printed?
#' @param ... optional argments passed to \code{\link{constrain_tapers}}
#' @return Object with class 'tapers'.
#' 
#' @seealso \code{\link{constrain_tapers}}, \code{\link{psdcore}}, \code{smooth.spline}
#' @example inst/Examples/rdex_riedsid.R
riedsid <- function(PSD, ntaper, 
                    tapseq=NULL, 
                    Deriv.method=c("local_qls","spg"),
                    Local.loss=c("Optim","Less","More"),
                    constrained=TRUE, c.method=NULL,
                    verbose=TRUE, ...) UseMethod("riedsid")

#' @rdname riedsid
#' @aliases riedsid.spec
#' @method riedsid spec
#' @S3method riedsid spec
riedsid.spec <- function(PSD, ...){
  stopifnot(is.spec(PSD))
  Pspec <- PSD$spec
  Ntap <- PSD$taper
  Tapseq <- PSD$freq
  riedsid(PSD=Pspec, ntaper=Ntap, tapseq=Tapseq, ...)
  #.NotYetImplemented()
}

#' @rdname riedsid
#' @method riedsid default
#' @S3method riedsid default
riedsid.default <- function(PSD, ntaper, 
                            tapseq=NULL, 
                            Deriv.method=c("local_qls","spg"),
                            Local.loss=c("Optim","Less","More"),
                            constrained=TRUE, c.method=NULL,
                            verbose=TRUE, ...) {
  ## spectral values
  PSD <- as.vector(PSD)
  # num freqs
  nf <- psd:::psd_envAssignGet("num_freqs", length(PSD))
  # prelims
  eps <- 1e-78 
  # .Machine$double.eps #  A small number to protect against zeros
  # vectorize initial estimate
  Zeros <- zeros(nf)
  nt <- length(ntaper)
  if (nt==1){
    ntap <- ntaper + Zeros
  } else {
    ntap <- ntaper
  }
  if (!is.tapers(ntap)) ntap <- as.tapers(ntap)
  
  # Set the number of tapers to within the range: 1/2 nf, 7/5 ntap
  # rowMins produces a rowvec of rowwise minimums; convert to colvec
  nspan <- minspan(ntap, nf)
  # The spectral gradients should be in log-space, so
  # create a log spec, and pad to handle begnning and end values
  nadd <- 1 + max(nspan)
  Y <- c(PSD[nadd:2], PSD, PSD[(nf-1):(nf-nadd)])
  Y[Y <= 0] <- eps
  lY <- log10(Y) # log in matlab is log_10
  dY <- d2Y <- Zeros
  #
  if (is.null(tapseq) | (length(tapseq) != nf)){
    kseq <- seq.int(from=0, to=1/2, length.out=length(PSD))
  } else {
    kseq <- tapseq
  }
  #
  # Smooth spectral derivatives
  #
  lsmeth <- switch(match.arg(Deriv.method), local_qls=TRUE, spg=FALSE)
  stopifnot(exists("lsmeth"))
  if (lsmeth){
    dsens <- switch(match.arg(Local.loss), Optim=12, More=6, Less=24) # 12 is optim
    DFUN <- function(j, 
                     j1=j-nspan[j]+nadd-1, 
                     j2=j+nspan[j]+nadd-1, 
                     jr=j1:j2, 
                     logY=lY[jr], 
                     dEps=eps,
                     CC=dsens){
      u <- jr - (j1 + j2)/2 # rowvec 
      u2 <- u*u             # rowvec
      L <- j2-j1+1          # constant
      L2 <- L*L             # constant
      LL2 <- L*L2           # constant
      LL2L <- LL2 - L       # constant
      #CC <- 12 # (orig)
      uzero <- (L2 - 1)/CC  # constant
      # first deriv
      dY <- u %*% logY * CC / LL2L
      # second deriv
      d2Y <- (u2 - uzero) %*% logY * 360 / LL2L / (L2-4)
      return(c(fdY2=dY*dY, fd2Y=d2Y, fdEps=dEps))
    }
    DX <- seq_len(nf) #1:nf
    RSS <- vapply(X=DX, FUN=DFUN, FUN.VALUE=c(1,1,1))
    attr(RSS, which="lsderiv") <- lsmeth
    RSS <- psd:::psd_envAssignGet("spectral_derivatives.ls", RSS)
    RSS <- abs(colSums(RSS))
    # sums:
    #[ ,1] fdY2
    #[ ,2] fd2Y
    #[ ,3] fdEps
    msg <- "local quadratic regression"
  } else {
    RSS <- splineGrad(dseq=log10(0.5+kseq), #seq.int(0,.5,length.out=length(PSD))), 
                      dsig=log10(PSD),
                      plot.derivs=FALSE, ...) #, spar=1)
    attr(RSS, which="lsderiv") <- lsmeth
    RSS <- psd:::psd_envAssignGet("spectral_derivatives", RSS) 
    #returns log
    RSS[,2:4] <- 10**RSS[,2:4]
    RSS <- abs(eps + RSS[,4] + RSS[,3]**2)
    msg <- "weighted cubic spline"
  }
  if (verbose) message(sprintf("Using spectral derivatives from  %s", msg))
  ##
  # (480)^0.2 == 3.437544
  KC <- 3.437544
  ##
  #(480)^0.2*abs(PSD/d2psd)^0.4
  # Original form:  kopt = 3.428*abs(PSD ./ d2psd).^0.4;
  # kopt = round( 3.428 ./ abs(eps + d2Y + dY.^2).^0.4 );
  ##
  stopifnot(exists("RSS"))
  kopt <- as.tapers( KC / (RSS ** 0.4) ) #/
  rm(RSS)
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
  #   #  R <- PSD"/PSD <- Y" + (Y')^2  2nd form preferred for consistent smoothing
  #   #
  #   #  Riedel-Sidorenko recipe (eq 13): 
  #   #       kopt <- (12*abs(PSD ./ d2psd)).^0.4 
  #   #  but parabolic weighting in psdcore requires: 
  #   #               (480)^0.2*abs(PSD ./ d2psd).^0.4
  #   #  Original form:  kopt <- 3.428*abs(PSD ./ d2psd).^0.4
  #   #
  #   # the optimal number of tapers (in an MSE sense):
  #   kopt_old <- as.tapers( 3.437544 / abs(eps + d2Y + dY*dY) ^ 0.4 ) 
  #   #
  #print(all.equal(kopt_old,kopt)) # TRUE!
  ##
  ## Constrain tapers
  stopifnot(diff(length(kopt), length(kseq))==0)
  if (constrained) kopt <- constrain_tapers(tapvec=kopt, 
                                            tapseq=kseq, 
                                            constraint.method=c.method, 
                                            verbose=verbose)
  ##
  return(kopt)
} 
# end riedsid.default
