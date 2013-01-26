#' Coerce an object into a 'taper' object.
#'
#' In a tapered spectrum estimation algorithm, it is
#' necessary to enforce rules on the number of tapers
#' that may be applied.
#'
#' Formal requirements enforced by this function are:
#' \itemize{
#' \item Non-zero.
#' \item Integer values.
#' \item Fewer than the half-length of the spectrum.
#' }
#' For example, we cannot apply
#' zero tapers (the result would be a raw periodogram)
#' or one million tapers (that would be absurd, and
#' violate orthogonality
#' conditions for any series less than two million terms long!).
#' 
#' An object with S3 class 'taper' is created;
#' this will have
#' a minimum number of tapers in each position
#' set by \code{min_taper}, and
#' a maximum number of tapers in each position
#' set by \code{max_taper}.
#' If \code{minspan=TRUE}, the bounded taper is fed through \code{\link{minspan}}
#' which will restrict the maximum tapers to less than or equal to
#' the half-length of the spectrum.
#'
#' Various classes can be coerced into a 'taper' object; those
#' tested sofar include: scalar, vector, matrix, data.frame, 
#' and list.  
#'
#' Multiple objects are concatenated into a single
#' vector dimension.  
#' 
#' Enabling \code{setspan} will only override
#' \code{max_taper} should it be larger than the half-width of the series.
#'
# @section Example of For example, if the object is 
# \code{list(x=c(1,2),y=c(3,4,5,0,1.1))} then the corresponding 'taper'
# objects for the following arguments are:
#
# \describe{
# \item{\emph{defaults}}{\code{[1,2,3,4,5,1,1]}}
# \item{\code{setspan=TRUE}}{\code{[1,2,3,3,3,1,1]}}
# \item{\code{max_taper=5}}{\code{[1,2,3,4,5,1,1]}}
# \item{\code{max_taper=5,setspan=TRUE}}{\code{[1,2,3,3,3,1,1]}}
# }
#'
#' @note No support (yet) for use of \code{min_taper,max_taper} as vectors, although
#' this could be quite desirable.
#'
#' @keywords taper S3methods
#' @param x An object to set
#' @param min_taper Set all values less than this to this.
#' @param max_taper Set all values greater than this to this.
#' @param setspan logical; should the taper object be passed through \code{\link{minspan}} before it is return?
#' @export
#' @author A.J. Barbour <andy.barbour@@gmail.com>
#' @seealso \code{\link{is.taper}}
#' @example x_examp/taper.R
as.taper <- function(x, min_taper=1, max_taper=NULL, setspan=FALSE){
  # taper should be non-zero integer, since it represents the
  # number of tapered sections to average; hence, floor.
  # pmin/pmax.int are fast versions of
  x <- as.vector(unlist(x))
  if (is.null(max_taper)) max_taper <- max(x)
  #print(summary(unclass(x)))
  stopifnot(min_taper*max_taper >= 1 & max_taper >= min_taper & !(is.character(x)))
  x <- as.integer(pmin.int(max_taper, pmax.int(min_taper, floor(x))))
  #x[x < min_taper] <- min_taper
  #   > as.integer(as.matrix(data.frame(x=1:10,y=10:19)))
  #   [1]  1  2  3  4  5  6  7  8  9 10 10 11 12 13 14 15 16 17 18 19
  class(x) <- "taper"
  if (setspan) x <- minspan(x)
  return(x)
}

###
###  Generic methods
###

#' @title Generic methods for objects with class 'taper'.
#' @keywords methods S3methods taper
#' @name taper-methods
#' @author A.J. Barbour <andy.barbour@@gmail.com>
#' @aliases taper
#' @rdname taper-methods
#' @docType methods
#' @import RColorBrewer
#'
#' @seealso \code{\link{as.taper}}, \code{\link{constrain_tapers}}, \code{par}
#' @param x taper object
#' @param xi optional vector for indices of \code{x}
#' @param object taper object
#' @param lwd line width (default is 1.8)
#' @param col color of line (default is "red")
#' @param pch point character (default is "_")
#' @param cex point size (default is 1)
#' @param color.pal color palette to use (choices are: "Blues","Spectral")
#' @param ylim optional limits for y-axis
#' @param ... optional arguments
#' @examples
#' ##
#' tap <- as.taper(c(1:49,50:0)+rnorm(1e2))
#' print(tap)
#' print(summary(tap))
#' plot(tap)
#' # no arithmetic methods
#' tap <- as.taper(tap/2)
#' lines(tap)
NULL

#' @rdname taper-methods
#' @name print
#' @aliases print.taper
#' @method print taper
#' @S3method print taper
print.taper <- function(x, ...){
  stopifnot(is.taper(x))
  xh <- paste(as.character(head(x)), collapse=" ")
  xt <- paste(as.character(tail(x)), collapse=" ")
  cat(sprintf("taper object:\n\thead:  %s\n\t\t...\n\ttail:  %s\n",xh,xt))
}

#' @rdname taper-methods
#' @name summary
#' @aliases summary.taper
#' @method summary taper
#' @S3method summary taper
summary.taper <- function(object, ...){
  stopifnot(is.taper(object))
  toret <- summary.default(object)
  class(toret) <- "summary.taper"
  return(toret)
}

#' @rdname taper-methods
#' @name print
#' @aliases print.summary.taper
#' @method print summary.taper
#' @S3method print summary.taper
print.summary.taper <- function(x, ...){
  cat("summary of tapers:\n")
  print(summary(x))
}

#' @rdname taper-methods
#' @name lines
#' @aliases lines.taper
#' @method lines taper
#' @S3method lines taper
lines.taper <- function(x, lwd=1.8, col="red", ...){
  stopifnot(is.taper(x))
  nt <- length(x)
  xi <- 1:nt
  #mx <- max(x)
  graphics::lines(xi, x, lwd=lwd, col=col, ...)
}

#' @rdname taper-methods
#' @name points
#' @aliases points.taper
#' @method points taper
#' @S3method points taper
points.taper <- function(x, pch="_", cex=1, ...){
  stopifnot(is.taper(x))
  nt <- length(x)
  xi <- 1:nt
  #mx <- max(x)
  graphics::points(xi, x, pch=pch, cex=cex, ...)
}

#' @rdname taper-methods
#' @name plot
#' @aliases plot.taper
#' @method plot taper
#' @S3method plot taper
plot.taper <- function(x, xi=NULL, color.pal=c("Blues","Spectral"), ylim=NULL, ...){
  stopifnot(is.taper(x))
  #if isS4(x) x <- x@tapers
  nt <- length(x)
  if (is.null(xi)){
    xi <- 1:nt
  }
  stopifnot(length(xi)==nt)
  mx <- max(x)
  pal <- match.arg(color.pal)
  npal <- switch(pal, RdYlBu=11, Spectral=11, Blues=9)
  pal.col <- RColorBrewer::brewer.pal(npal, pal)
  cols <- grDevices::colorRampPalette(pal.col)(mx)
  if (is.null(ylim)) ylim <- 1.15*c(0.5, mx)
  graphics::plot.default(xi, x,
               ylab="number of tapers",
               xlab="taper index",
               ylim=ylim, yaxs="i", 
               #xlim=c(-1, nt+2), 
                         xaxs="i",
               lwd=1.8,
               type="h",
               col=cols[x],
               ...)
  graphics::lines.default(xi,x,col="black",lwd=0.7)
  # plot log2 multiples as horiz lines
  hl <- 2**(1:round(log2(mx)))
  graphics::abline(h=hl,lty=1,lwd=0.6,col="black")
  vl <- c(1, nt)
  graphics::abline(v=vl,lty=3,lwd=2,col="blue")
}

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
#' @param tapvec object with class taper
#' @param f.samp scalar; the sampling frequency (e.g. Hz) of the series the tapers are for
#' @param n.freq scalar; the number of frequencies of the original spectrum (if \code{NULL} the length of the taper object is assumed to be the number)
#' @param ... additional arguments (unused)
#' @return A list with the following properties (and names):
#' \itemize{
#' \item{\code{taper}: The original taper vector.}
#' \item{\code{stderr}: The standard error of the spectrum.}
#' \item{\code{resolution}: The effective spectral resolution.}
#' \item{\code{dof}: The number of degrees of freedom.}
#' \item{\code{bw}: The effective bandwidth of the spectrum.}
#' }
#' @export
#' @keywords properties taper resolution uncertainty degrees-of-freedom bandwidth
spectral_properties <- function(tapvec, f.samp=1, n.freq=NULL, ...) UseMethod("spectral_properties")
#' @rdname spectral_properties
#' @aliases spectral_properties.spec
#' @method spectral_properties spec
#' @S3method spectral_properties spec
spectral_properties.spec <- function(tapvec, f.samp=1, n.freq=NULL, ...) .NotYetImplemented()
#' @rdname spectral_properties
#' @aliases spectral_properties.taper
#' @method spectral_properties taper
#' @S3method spectral_properties taper
spectral_properties.taper <- function(tapvec, f.samp=1, n.freq=NULL, ...){
  stopifnot(is.taper(tapvec))
  K <- unclass(tapvec)
  Nyquist <- f.samp/2
  if (is.null(n.freq)) n.freq <- length(tapvec)
  ## Resolution
  Resolu <- K * Nyquist / n.freq
  ## Uncertainty
  StdErr <- 1 / sqrt(K / 1.2)
  #Var <- 10 / K / 12
  ## Deg Freedom
  Dof <- 2 * K
  ## Bandwidth
  # Walden et al
  # half-width W = (K + 1)/{2(N + 1)}
  # effective bandwidth ~ 2 W (accurate for many spectral windows)
  BW <- K/n.freq
  ##
  return(data.frame(taper=K, stderr=StdErr, resolution=Resolu, dof=Dof, bw=BW))
}

###
###  Weighting methods
###

#' Calculate parabolic weighting factors.
#'
#' The resampled spectrum involves summing weighted tapers; this produces
#' the weighting factors.
#'
#' Weighting factors are calculated as follows:
#' \deqn{W_N \equiv n_T^2 - \frac{3  K_N^2}{2 n_T (n_T - 1/4) (n_T + 1)}}
#' where \eqn{n_T} is the total number of tapers, and 
#' \eqn{K_N} is the integer sequence \eqn{[0,n_T-1]} 
#'
#' @export
#' @keywords taper taper-weighting
#' @author A.J. Barbour <andy.barbour@@gmail.com> adapted original by R.L.Parker,
#' and authored the optimized version.
#' @seealso \code{\link{psdcore}}, \code{\link{riedsid}}
#'
#' @param tapvec 'taper' object; the number of tapers at each frequency
#' @param tap.index integer; the index of \code{tapvec} from which to find weights
#' @param ntap integer; the number of tapers to provide weightings for.
#' @return A list with taper indices, and the weights \eqn{W_N}.
parabolic_weights <- function(tapvec, tap.index=1L) UseMethod("parabolic_weights")
#' @rdname parabolic_weights
#' @aliases parabolic_weights.taper
#' @method parabolic_weights taper
#' @S3method parabolic_weights taper
parabolic_weights.taper <- function(tapvec, tap.index=1L){
  stopifnot(is.taper(tapvec) | ((tap.index > 0L) & (tap.index <= length(tapvec))))
  kWeights <- parabolic_weights_fast(tapvec[as.integer(tap.index)])
  return(kWeights)
}
#' @rdname parabolic_weights
#' @aliases parabolic_weights_fast
#' @export
#' @keywords taper taper-weighting
parabolic_weights_fast <- function(ntap=1L) UseMethod("parabolic_weights_fast")
#' @rdname parabolic_weights
#' @aliases parabolic_weights_fast.taper
#' @method parabolic_weights_fast default
#' @S3method parabolic_weights_fast default
parabolic_weights_fast.default <- function(ntap=1L){
  ntap <- max(1, as.integer(ntap[1]))
  kseq <- seq.int(from=0, to=ntap-1, by=1)
  lk <- length(kseq)
  K2 <- kseq * kseq # vector
  NT2 <- ntap * ntap # scalar
  NT3 <- NT2 * ntap # scalar
  #w = (tapers^2 - (k-1).^2)*(1.5/(tapers*(tapers-0.25)*(tapers+1)));
  TW <- matrix((NT2 - K2) * 3/(2*NT3 + NT2*3/2 - ntap/2), ncol=lk)
  return(list(taper_seq=kseq+1, taper_weights=TW ))
}

###
###  Constraint methods
###

#' @title Taper constraint methods.
#'
#' @description
#' In the Riedel-Sidorenko recipe, the number of optimal tapers
#' at each frequency is strongly dependent on the first and
#' second derivatives of the spectrum. It is crucial to enforce
#' constraints on the number of actual tapers applied; this is
#' because the derivatives of "noisy" series can be bogus.
#'
#' @keywords taper taper-constraints riedel-sidorenko
#' @rdname taper-constraints
#' @name taper-constraints
NULL

#' @description \code{\link{minspan}} sets the maximum span a taper object
#' may have, which is necessary because it would be nonsense to
#' have more tapers than the length of the series. 
#' 
#' @details \code{\link{minspan}} bounds the number of tapers between
#' the minimum of the half-length of the series, and 7/5 times
#' the tapers.  In code this would look something like: 
#' \code{min(length(tapvec)/2, 7*tapvec/5)}
#'
#' @rdname taper-constraints
#' @title minspan
#' @export
#' @keywords taper taper-constraints
#' @seealso \code{\link{splineGrad}}, \code{\link{riedsid}}
#'
#' @author A.J. Barbour <andy.barbour@@gmail.com> and R.L.Parker. 
#' AJB adapted some of RLP's original code,
#' and wrote the main function in \code{\link{ctap_simple}} for dynamic loading C-code.
#' The main function used by \code{\link{ctap_markov}} is from Morhac (2008).
#'
minspan <- function(tapvec, ...) UseMethod("minspan")
#' @rdname taper-constraints
#' @aliases minspan.taper
#' @method minspan taper
#' @S3method minspan taper
minspan.taper <- function(tapvec, ...){
  stopifnot(is.taper(tapvec))
  nspan <- as.taper(7*tapvec/5, max_taper=length(tapvec)/2)
  return(nspan)
}

#' @description \code{\link{constrain_tapers}} refines the number of tapers; 
#' the method by which it does this is chosen with the \code{constraint.method}
#' parameter. See \strong{Constraint methods} section for descriptions of each method.
#' Below is a summary of the function associated with each \code{constraint.method}:
#' \itemize{
#'   \item \code{'simple.slope'} uses \code{\link{ctap_simple}}
#'   \item \code{'markov.chain'} uses \code{\link{ctap_markov}}
#'   \item \code{'loess.smooth'} uses \code{\link{ctap_loess}}
#'   \item \code{'friedman.smooth'} uses \code{\link{ctap_friedman}}
#'   \item \code{'none'} returns unbounded tapers.
#' }
#' 
#' @section Details of Constraint Methods:
#' \subsection{via first differencing (default)}{
#'  \code{\link{ctap_simple}} is the default, and preferred constraint method.
#' The algortihm uses first-differencing to modify the number
#' of tapers in the previous position.  Effectively, the constraint
#' is based on a causal, 1st-order Finite Impulse-response Filter (FIR) 
#' which makes the method
#' sensitive to rapid changes in the number of tapers; naturally,
#' smoother spectra tend to produce less fluctuation in taper numbers, 
#' which makes this well suited for adaptive processing. 
#'
#' This produces, generally, the most
#' stable results, meaning repeatedly running the constraint will not change values
#' other than on the first execution; the same cannot be said for the other
#' methods.
#'
#' In pure-R this algorithm can be very slow; however, here we have
#' included it as dynamically loaded c-code so it it reasonably fast.
#' }
#'
#' \subsection{via Markov Chain}{
#' \code{\link{ctap_markov}} uses a Markov Chain, based on the theory of 
#' quantum-well probability chains, which are
#' commonly used in gamma-ray spectroscopy.
#'
#' The main function behind this method is from Morhac (2008): \code{SpectrumSmoothMarkov}.
#' This calculates the probability that the number of tapers would have
#' changed (from it's previous value); it is very fast.  Details of the theory 
#' behind this algorithm may be found in Morhac (2008) and Silagadze (1996).
#' }
#' 
#' \subsection{via LOESS smoothing}{
#' \code{\link{ctap_loess}} uses \code{loess} to smooth the taper vector; is
#' can be very slow thanks to quadratic scaling.
#' }
#'
#' \subsection{via Friedman super-smoothing}{
#' \code{\link{ctap_friedman}} uses \code{supsmu}, the Friedman super-smoother.
#' }
#'
#' @section Warning:
#' \code{\link{ctap_markov}} can produce "unstable"
#' results in the sense that for
#' successive application on taper vectors, even modest sized serially-correlated
#' peaks tends to sharpen; hence, this method should be used with care, unless the 
#' intention is to specifically
#' enhance peaks.  The \code{'chain.width'} parameter controls the broadness
#' of the a priori distribution.  As a rule of thumb: the smaller the parameter is, 
#' the shorter the tails become.
#'
#' \code{\link{ctap_markov}} results tend to be strongly dependent on
#' the tuning parameters given to \code{loess} (for obvious reasons); hence, 
#' some effort should be given to understand
#' their effect, and/or re-tuning them if needed.
#'
#' \code{\link{ctap_friedman}} results are generally poor in my opinion; 
#' hence, the method may be removed in future releases.
#'
#' @rdname taper-constraints
#' @title constrain_tapers
#' @aliases constrain_tapers taper-constraints
#' @export
#' @import Peaks
#' @keywords taper taper-constraints
#' @param tapvec 'taper' object; the number of tapers at each frequency
#' @param tapseq vector; positions or frequencies -- necessary for smoother methods
#' @param constraint.method  character; method to use for constraints on tapers numbers
## @param min_tapers integer; the minimum number of tapers
#' @param verbose logical; should warnings and messages be given?
#' @param maxslope integer; constrain based on this maximum first difference
#' @param chain.width  scalar; the width the MC should consider for the change probability
#' @param normalize logical; should the refined tapers be normalized?
#' @param loess.span  scalar; the span used in \code{loess}
#' @param loess.degree  scalar; the polynomial degree
#' @param smoo.span  scalar; fraction of the observations in the span of the running lines smoother
#' @param smoo.bass  scalar; controls the smoothness of the fitted curve 
#' @param ... optional arguments (unused)
#' @return An object with class 'taper'.
#'
#' @references Morhac, M. (2008), Peaks: Peaks, \emph{R package}, \strong{version 0.2}
#' @references Silagadze, Z.K. (1996), A new algorithm for automatic photopeak searches,
#' \emph{Nucl. Instrum. Meth. A}, \strong{376} 451,
#' \url{http://arxiv.org/abs/hep-ex/9506013}
#'
#' @examples
#' \dontrun{
#' ## compare all the methods:
#' demo("ctap")
#' }
constrain_tapers <- function(tapvec, tapseq=NULL,
                             constraint.method=c("simple.slope",
                                                 "markov.chain",
                                                 "loess.smooth",
                                                 "friedman.smooth",
                                                 "none"),
                             verbose=TRUE, ...) UseMethod("constrain_tapers")
#' @rdname taper-constraints
#' @aliases constrain_tapers.taper
#' @method constrain_tapers taper
#' @S3method constrain_tapers taper
constrain_tapers.taper <- function(tapvec, tapseq=NULL,
                                   constraint.method=c("simple.slope",
                                                       "markov.chain",
                                                       "loess.smooth",
                                                       "friedman.smooth",
                                                       "none"),
                                   verbose=TRUE, ...){
  stopifnot(is.taper(tapvec))
  # choose the appropriate method to apply taper constraints
  cmeth <- match.arg(constraint.method) 
  if (cmeth=="none"){
    if (verbose) warning("no taper optimization constraints applied")
    tapvec.adj <- tapvec
  } else {
    if (verbose) message(sprintf("Constraining tapers with  ...  %s  ...  method",cmeth))
    CTAPFUN <- switch(cmeth,
                      "simple.slope"=ctap_simple,
                      "markov.chain"=ctap_markov,
                      "loess.smooth"=ctap_loess,
                      "friedman.smooth"=ctap_friedman)
    tapvec.adj <- CTAPFUN(tapvec, tapseq, ...)
  }
  # MAX/MIN bounds
  # set the maximim tapers: Never average over more than the length of the spectrum!
  #   maxtap <- round(length(tapvec)/2)
  #   # ensure the minimum is below
  #   mintap <- min(1, maxtap)
  #   # sort for posterity
  #   tapbounds <- sort(c(mintap, maxtap))
  #   mintap <- tapbounds[1]
  #   maxtap <- tapbounds[2]
  #   tapvec.adj[tapvec.adj < mintap] <- mintap
  #   tapvec.adj[tapvec.adj > maxtap] <- maxtap
  tapvec.adj <- as.taper(tapvec.adj, min_taper=1, max_taper=round(length(tapvec.adj)/2))
  return(tapvec.adj)
}

#' @rdname taper-constraints
#' @aliases constrain_taper_simple_slope ctap_simple
#' @export
#' @keywords taper taper-constraints
ctap_simple <- function(tapvec, tapseq=NA, maxslope=1, ...) UseMethod("ctap_simple")
#' @rdname taper-constraints
#' @aliases ctap_simple.taper
#' @method ctap_simple taper
#' @S3method ctap_simple taper
ctap_simple.taper <- function(tapvec, tapseq=NA, maxslope=1, ...){
  stopifnot(is.taper(tapvec))
  # tapseq not needed
  # 'taper' object gives integer values, but code requires real
  tapvec <- as.numeric(tapvec)
  maxslope <- as.numeric(maxslope)
  # c-code used for speed up of forward+backward operations
  tapvec.adj <- as.taper(.Call("rlp_constrain_tapers", tapvec, maxslope, PACKAGE = "rlpSpec"))
  return(tapvec.adj)
}

#' @rdname taper-constraints
#' @aliases constrain_taper_markov_chain ctap_markov
#' @export
#' @keywords taper taper-constraints
ctap_markov <- function(tapvec, tapseq=NA, chain.width=round(5*length(tapvec)), ...) { UseMethod("ctap_markov") }
#' @rdname taper-constraints
#' @aliases ctap_markov.taper
#' @method ctap_markov taper
#' @S3method ctap_markov taper
ctap_markov.taper <- function(tapvec, tapseq=NA, chain.width=round(5*length(tapvec)), normalize=TRUE, ...){
  stopifnot(is.taper(tapvec))
  # tapseq not needed
  # bound the chain.width
  MC.win <- max(1, round(chain.width))
  tapvec.adj <- as.taper(Peaks::SpectrumSmoothMarkov(as.numeric(tapvec), MC.win))
  if (normalize) tapvec.adj <- as.taper(tapvec.adj*max(tapvec)/max(tapvec.adj))
  return(tapvec.adj)
}

#' @rdname taper-constraints
#' @aliases constrain_taper_loess_smooth ctap_loess
#' @export
#' @keywords taper taper-constraints
ctap_loess <- function(tapvec, tapseq=NULL, loess.span=.3, loess.degree=1, verbose=TRUE, ...){ UseMethod("ctap_loess") }
#' @rdname taper-constraints
#' @aliases ctap_loess.taper
#' @method ctap_loess taper
#' @S3method ctap_loess taper
ctap_loess.taper <- function(tapvec, tapseq=NULL, loess.span=.3, loess.degree=1, verbose=TRUE, ...){
  stopifnot(is.taper(tapvec))
  # having an x-sequence is absolutely critical to obtaining useful results
  if (is.null(tapseq)){
    tapseq <- 1:length(tapvec)
    if (verbose) warning("Generated a position sequence; results may be bogus.")
  }
  lt <- length(tapvec)
  if (verbose & lt > 1e4) warning("Loess-method has quadratic memory scaling (1e3 pt -> 10 Mb)...")
  trc <- ifelse(lt >= 1e3, "approximate", "exact")
  loe <- stats::loess(y ~ x, 
                      data.frame(x=tapseq, y=as.numeric(tapvec)), 
                      span=loess.span, degree=loess.degree,
                      control = loess.control(trace.hat = trc))
  tapvec.adj <- as.taper(stats::predict(loe))
  return(tapvec.adj)
}

#' @rdname taper-constraints
#' @aliases constrain_taper_friedman_smooth ctap_friedman
#' @aliases ctap_friedman.taper
#' @export
#' @keywords taper taper-constraints
ctap_friedman <- function(tapvec, tapseq=NULL, smoo.span=.3, smoo.bass=2, verbose=TRUE, ...){ UseMethod("ctap_friedman") }
#' @rdname taper-constraints
#' @method ctap_friedman taper
#' @S3method ctap_friedman taper
ctap_friedman.taper <- function(tapvec, tapseq=NULL, smoo.span=.3, smoo.bass=2, verbose=TRUE, ...){
  stopifnot(is.taper(tapvec))
  # having an x-sequence is absolutely critical to obtaining useful results
  if (is.null(tapseq)){
    tapseq <- 1:length(tapvec)
    if (verbose) warning("Generated a position sequence; results may be bogus.")
  }
  # Friedman's super smoother.
  # Meh.
  tapvec.adj <- as.taper(stats::supsmu(tapseq, as.numeric(tapvec), span=smoo.span, bass=smoo.bass)$y)
  return(tapvec.adj)
}
###
