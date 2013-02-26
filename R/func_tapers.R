#' Coerce an object into a \code{'tapers'} object.
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
#' An object with S3 class \code{'tapers'} is created;
#' this will have
#' a minimum number of tapers in each position
#' set by \code{min_taper}, and
#' a maximum number of tapers in each position
#' set by \code{max_taper}.
#' If \code{minspan=TRUE}, the bounded taper is fed through \code{\link{minspan}}
#' which will restrict the maximum tapers to less than or equal to
#' the half-length of the spectrum.
#'
#' Various classes can be coerced into a \code{'tapers'} object; those
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
# \code{list(x=c(1,2),y=c(3,4,5,0,1.1))} then the corresponding \code{'tapers'}
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
#' @keywords tapers S3methods
#' @param x An object to set
#' @param min_taper Set all values less than this to this.
#' @param max_taper Set all values greater than this to this.
#' @param setspan logical; should the tapers object be passed through \code{\link{minspan}} before it is return?
#' @export
#' @return An object with class taper.
#' @author A.J. Barbour <andy.barbour@@gmail.com>
#' @seealso \code{\link{is.tapers}}
#' @example inst/Examples/rdex_tapers.R
as.tapers <- function(x, min_taper=1, max_taper=NULL, setspan=FALSE){
  # taper should be non-zero integer, since it represents the
  # number of tapered sections to average; hence, floor.
  # pmin/pmax.int are fast versions of
  x <- as.vector(unlist(x))
  stopifnot(!(is.character(x)))
  if (is.null(max_taper)) max_taper <- ceiling(max(x))
  stopifnot(min_taper*max_taper >= 1 & max_taper >= min_taper)
  #
  x <- as.integer(pmin.int(max_taper, pmax.int(min_taper, floor(x))))
  #
  class(x) <- "tapers"
  
  attr(x, "n_taper_limits") <- c(min_taper, max_taper)
  attr(x, "taper_positions") <- NA
  #
  if (setspan) x <- minspan(x)
  #
  attr(x, "span_was_set") <- setspan
  attr(x, "n_taper_limits_orig") <- c(min_taper, max_taper)
  #
  return(x)
}
#' @rdname as.tapers
#' @name tapers
#' @export
tapers <- as.tapers

###
###  Generic methods
###

#' @title Generic methods for objects with class \code{'tapers'}.
#' @keywords methods S3methods tapers
#' @name tapers-methods
#' @author A.J. Barbour <andy.barbour@@gmail.com>
#' @rdname tapers-methods
#' @docType methods
#' @import RColorBrewer
#'
#' @seealso \code{\link{as.tapers}}, \code{\link{constrain_tapers}}, \code{par}
#' @param x tapers object
#' @param xi optional vector for indices of \code{x}
#' @param object tapers object
#' @param lwd line width (default is 1.8)
#' @param col color of line (default is "red")
#' @param pch point character (default is "_")
#' @param cex point size (default is 1)
#' @param color.pal color palette to use (choices are: "Blues","Spectral")
#' @param ylim optional limits for y-axis
#' @param hv.lines logical; should horizontal and vertival reference lines be plotted?
#' @param ... optional arguments
#' @return \code{plot} returns a list with names: \code{line.colors} (hex values)
#' @examples
#' ##
#' tap <- as.tapers(c(1:49,50:0)+rnorm(1e2))
#' print(tap)
#' print(summary(tap))
#' plot(tap)
#' # no arithmetic methods
#' tap <- as.tapers(tap/2)
#' lines(tap)
NULL

#' @rdname tapers-methods
#' @name as.data.frame.tapers
#' @method as.data.frame tapers
#' @S3method as.data.frame tapers
as.data.frame.tapers <- function(x, ...){
  df <- as.data.frame.numeric(x)
  names(df) <- "n.tapers"
  return(df)
}
#' @rdname tapers-methods
#' @name data.frame.tapers
#' @method data.frame tapers
#' @S3method data.frame tapers
data.frame.tapers <- as.data.frame.tapers

#' @rdname tapers-methods
#' @name print
#' @aliases print.tapers
#' @method print tapers
#' @S3method print tapers
print.tapers <- function(x, ...){
  stopifnot(is.tapers(x))
  xh <- paste(as.character(head(x)), collapse=" ")
  xt <- paste(as.character(tail(x)), collapse=" ")
  cat(sprintf("'tapers' object: num. tapers applied by index\n\thead:  %s\n\t\t...\n\ttail:  %s\n",xh,xt))
}

#' @rdname tapers-methods
#' @name summary
#' @aliases summary.tapers
#' @method summary tapers
#' @S3method summary tapers
summary.tapers <- function(object, ...){
  stopifnot(is.tapers(object))
  toret <- summary.default(object)
  class(toret) <- "summary.tapers"
  return(toret)
}

#' @rdname tapers-methods
#' @name print
#' @aliases print.summary.tapers
#' @method print summary.tapers
#' @S3method print summary.tapers
print.summary.tapers <- function(x, ...){
  cat("summary of tapers:\n")
  print(summary(x))
}

#' @rdname tapers-methods
#' @name lines
#' @aliases lines.tapers
#' @method lines tapers
#' @S3method lines tapers
lines.tapers <- function(x, lwd=1.8, col="red", ...){
  stopifnot(is.tapers(x))
  nt <- length(x)
  xi <- 1:nt
  #mx <- max(x)
  graphics::lines(xi, x, lwd=lwd, col=col, ...)
}

#' @rdname tapers-methods
#' @name points
#' @aliases points.tapers
#' @method points tapers
#' @S3method points tapers
points.tapers <- function(x, pch="_", cex=1, ...){
  stopifnot(is.tapers(x))
  nt <- length(x)
  xi <- 1:nt
  #mx <- max(x)
  graphics::points(xi, x, pch=pch, cex=cex, ...)
}

#' @rdname tapers-methods
#' @name plot
#' @aliases plot.tapers
#' @method plot tapers
#' @S3method plot tapers
plot.tapers <- function(x, xi=NULL, color.pal=c("Blues","Spectral"), ylim=NULL, hv.lines=FALSE, ...){
  stopifnot(is.tapers(x))
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
  if (hv.lines){
    # plot log2 multiples as horiz lines
    hl <- 2**(1:round(log2(mx)))
    graphics::abline(h=hl,lty=1,lwd=0.6,col="black")
    vl <- c(1, nt)
    graphics::abline(v=vl,lty=3,lwd=2,col="blue")
  }
  return(invisible(list(line.colors=cols)))
}

###
###  Weighting methods
###

#' Calculate parabolic weighting factors.
#'
#' The resampled spectrum involves summing weighted tapers; this produces
#' the weighting factors.
#'
#' If one has a \code{tapers} object, specify the \code{taper.index} to
#' produce a sequence of weights up to the value at that index; the user
#' is likely to never need to use this function.
#'
#' Weighting factors are calculated as follows:
#' \deqn{W_N \equiv n_T^2 - \frac{3  K_N^2}{2 n_T (n_T - 1/4) (n_T + 1)}}
#' where \eqn{n_T} is the total number of tapers, and 
#' \eqn{K_N} is the integer sequence \eqn{[0,n_T-1]} 
#'
#' @export
#' @keywords tapers tapers-weighting
#' @author A.J. Barbour <andy.barbour@@gmail.com> adapted original by R.L.Parker,
#' and authored the optimized version.
#' @seealso \code{\link{psdcore}}, \code{\link{riedsid}}
#'
#' @param tapvec \code{'tapers'} object; the number of tapers at each frequency
#' @param tap.index integer; the index of \code{tapvec} from which to produce a sequence of weights for
#' @param ntap integer; the number of tapers to provide weightings for.
#' @return A list with taper indices, and the weights \eqn{W_N}.
#'
#' @example inst/Examples/rdex_parabolicweights.R
parabolic_weights <- function(tapvec, tap.index=1L) UseMethod("parabolic_weights")
#' @rdname parabolic_weights
#' @aliases parabolic_weights.tapers
#' @method parabolic_weights tapers
#' @S3method parabolic_weights tapers
parabolic_weights.tapers <- function(tapvec, tap.index=1L){
  stopifnot(is.tapers(tapvec) | ((tap.index > 0L) & (tap.index <= length(tapvec))))
  kWeights <- parabolic_weights_fast(tapvec[as.integer(tap.index)])
  return(kWeights)
}
#' @rdname parabolic_weights
#' @aliases parabolic_weights_fast
#' @export
#' @keywords tapers tapers-weighting
parabolic_weights_fast <- function(ntap=1L) UseMethod("parabolic_weights_fast")
#' @rdname parabolic_weights
#' @aliases parabolic_weights_fast.tapers
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
#' @keywords tapers tapers-constraints riedel-sidorenko
#' @rdname tapers-constraints
#' @name tapers-constraints
NULL

#' @description \code{\link{minspan}} sets the maximum span a tapers object
#' may have, which is necessary because it would be nonsense to
#' have more tapers than the length of the series. 
#' 
#' @details \code{\link{minspan}} bounds the number of tapers to within
#' the minimum of
#' either the maximum number of tapers found in the object, 
#' or the half-length of the series.
#'
### the minimum of the half-length of the series, and 7/5 times
### the tapers.  In code this would look something like: 
### \code{min(length(tapvec)/2, 7*tapvec/5)}
#'
#' @rdname tapers-constraints
#' @title minspan
#' @export
#' @keywords tapers tapers-constraints
#' @seealso \code{\link{splineGrad}}, \code{\link{riedsid}}
#'
#' @author A.J. Barbour <andy.barbour@@gmail.com> and R.L.Parker. 
#' AJB adapted some of RLP's original code,
#' and wrote the main function in \code{\link{ctap_simple}} for dynamic loading C-code.
#' The main function used by \code{\link{ctap_markov}} is from Morhac (2008).
#'
minspan <- function(tapvec, ...) UseMethod("minspan")
#' @rdname tapers-constraints
#' @aliases minspan.tapers
#' @method minspan tapers
#' @S3method minspan tapers
minspan.tapers <- function(tapvec, ...){
  stopifnot(is.tapers(tapvec))
  ##tapvec <- 7*tapvec/5
  maxtap <- min(max(tapvec), length(tapvec)/2)
  nspan <- as.tapers(tapvec, min_taper=1, max_taper=maxtap, setspan=FALSE)
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
#' \code{\link{ctap_loess}} results tend to be strongly dependent on
#' the tuning parameters given to \code{loess} (for obvious reasons); hence, 
#' some effort should be given to understand
#' their effect, and/or re-tuning them if needed.
#'
#' \code{\link{ctap_friedman}} results are generally poor in my opinion; 
#' hence, the method may be removed in future releases.
#'
#' @rdname tapers-constraints
#' @title constrain_tapers
#' @aliases constrain_tapers tapers-constraints
#' @export
#' @import Peaks
#' @keywords tapers tapers-constraints
#' @param tapvec \code{'tapers'} object; the number of tapers at each frequency
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
#' @return An object with class \code{'tapers'}.
#'
#' @references Morhac, M. (2008), Peaks: Peaks, \emph{R package}, \strong{version 0.2}
#' @references Silagadze, Z.K. (1996), A new algorithm for automatic photopeak searches,
#' \emph{Nucl. Instrum. Meth. A}, \strong{376} 451,
#' \url{http://arxiv.org/abs/hep-ex/9506013}
#'
#' @example inst/Examples/rdex_constraintapers.R
constrain_tapers <- function(tapvec, tapseq=NULL,
                             constraint.method=c("simple.slope",
                                                 "markov.chain",
                                                 "loess.smooth",
                                                 "friedman.smooth",
                                                 "none"),
                             verbose=TRUE, ...) UseMethod("constrain_tapers")
#' @rdname tapers-constraints
#' @aliases constrain_tapers.tapers
#' @method constrain_tapers tapers
#' @S3method constrain_tapers tapers
constrain_tapers.tapers <- function(tapvec, tapseq=NULL,
                                   constraint.method=c("simple.slope",
                                                       "markov.chain",
                                                       "loess.smooth",
                                                       "friedman.smooth",
                                                       "none"),
                                   verbose=TRUE, ...){
  stopifnot(is.tapers(tapvec))
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
  tapvec.adj <- as.tapers(tapvec.adj, min_taper=1, max_taper=round(length(tapvec.adj)/2))
  return(tapvec.adj)
}

#' @rdname tapers-constraints
#' @aliases constrain_taper_simple_slope ctap_simple
#' @export
#' @keywords tapers tapers-constraints
ctap_simple <- function(tapvec, tapseq=NA, maxslope=1, ...) UseMethod("ctap_simple")
#' @rdname tapers-constraints
#' @aliases ctap_simple.tapers
#' @method ctap_simple tapers
#' @S3method ctap_simple tapers
ctap_simple.tapers <- function(tapvec, tapseq=NA, maxslope=1, ...){
  stopifnot(is.tapers(tapvec))
  # tapseq not needed
  # \code{'tapers'} object gives integer values, but code requires real
  tapvec <- as.numeric(tapvec)
  maxslope <- as.numeric(maxslope)
  # c-code used for speed up of forward+backward operations
  tapvec.adj <- as.tapers(.Call("rlp_constrain_tapers", tapvec, maxslope, PACKAGE="psd"))
  return(tapvec.adj)
}

#' @rdname tapers-constraints
#' @aliases constrain_taper_markov_chain ctap_markov
#' @export
#' @keywords tapers tapers-constraints
ctap_markov <- function(tapvec, tapseq=NA, chain.width=round(5*length(tapvec)), ...) { UseMethod("ctap_markov") }
#' @rdname tapers-constraints
#' @aliases ctap_markov.tapers
#' @method ctap_markov tapers
#' @S3method ctap_markov tapers
ctap_markov.tapers <- function(tapvec, tapseq=NA, chain.width=round(5*length(tapvec)), normalize=TRUE, ...){
  stopifnot(is.tapers(tapvec))
  # tapseq not needed
  # bound the chain.width
  MC.win <- max(1, round(chain.width))
  tapvec.adj <- as.tapers(Peaks::SpectrumSmoothMarkov(as.numeric(tapvec), MC.win))
  if (normalize) tapvec.adj <- as.tapers(tapvec.adj*max(tapvec)/max(tapvec.adj))
  return(tapvec.adj)
}

#' @rdname tapers-constraints
#' @aliases constrain_taper_loess_smooth ctap_loess
#' @export
#' @keywords tapers tapers-constraints
ctap_loess <- function(tapvec, tapseq=NULL, loess.span=.3, loess.degree=1, verbose=TRUE, ...){ UseMethod("ctap_loess") }
#' @rdname tapers-constraints
#' @aliases ctap_loess.tapers
#' @method ctap_loess tapers
#' @S3method ctap_loess tapers
ctap_loess.tapers <- function(tapvec, tapseq=NULL, loess.span=.3, loess.degree=1, verbose=TRUE, ...){
  stopifnot(is.tapers(tapvec))
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
  tapvec.adj <- as.tapers(stats::predict(loe))
  return(tapvec.adj)
}

#' @rdname tapers-constraints
#' @aliases constrain_taper_friedman_smooth ctap_friedman
#' @aliases ctap_friedman.tapers
#' @export
#' @keywords tapers tapers-constraints
ctap_friedman <- function(tapvec, tapseq=NULL, smoo.span=.3, smoo.bass=2, verbose=TRUE, ...){ UseMethod("ctap_friedman") }
#' @rdname tapers-constraints
#' @method ctap_friedman tapers
#' @S3method ctap_friedman tapers
ctap_friedman.tapers <- function(tapvec, tapseq=NULL, smoo.span=.3, smoo.bass=2, verbose=TRUE, ...){
  stopifnot(is.tapers(tapvec))
  # having an x-sequence is absolutely critical to obtaining useful results
  if (is.null(tapseq)){
    tapseq <- 1:length(tapvec)
    if (verbose) warning("Generated a position sequence; results may be bogus.")
  }
  # Friedman's super smoother.
  # Meh.
  tapvec.adj <- as.tapers(stats::supsmu(tapseq, as.numeric(tapvec), span=smoo.span, bass=smoo.bass)$y)
  return(tapvec.adj)
}
###
