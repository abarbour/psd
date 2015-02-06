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
  x <- as.integer(pmin.int(max_taper, pmax.int(min_taper, round(x))))
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
#' @export
as.data.frame.tapers <- function(x, ...){
  df <- as.data.frame.numeric(x)
  names(df) <- "n.tapers"
  return(df)
}
#' @rdname tapers-methods
#' @export
data.frame.tapers <- as.data.frame.tapers

#' @rdname tapers-methods
#' @aliases print.tapers
#' @export
print.tapers <- function(x, ...){
  stopifnot(is.tapers(x))
  xh <- paste(as.character(head(x)), collapse=" ")
  xt <- paste(as.character(tail(x)), collapse=" ")
  cat(sprintf("'tapers' object: num. tapers applied by index\n\thead:  %s\n\t\t...\n\ttail:  %s\n",xh,xt))
}

#' @rdname tapers-methods
#' @aliases summary.tapers
#' @export
summary.tapers <- function(object, ...){
  stopifnot(is.tapers(object))
  toret <- summary.default(object)
  class(toret) <- "summary.tapers"
  return(toret)
}

#' @rdname tapers-methods
#' @aliases print.summary.tapers
#' @export
print.summary.tapers <- function(x, ...){
  cat("summary of tapers:\n")
  print(summary(x))
}

#' @rdname tapers-methods
#' @aliases lines.tapers
#' @export
lines.tapers <- function(x, lwd=1.8, col="red", ...){
  stopifnot(is.tapers(x))
  nt <- length(x)
  xi <- 1:nt
  #mx <- max(x)
  graphics::lines(xi, x, lwd=lwd, col=col, ...)
}

#' @rdname tapers-methods
#' @aliases points.tapers
#' @export
points.tapers <- function(x, pch="_", cex=1, ...){
  stopifnot(is.tapers(x))
  nt <- length(x)
  xi <- 1:nt
  #mx <- max(x)
  graphics::points(xi, x, pch=pch, cex=cex, ...)
}

#' @rdname tapers-methods
#' @aliases plot.tapers
#' @export
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
#' Weighting factors, \eqn{W}{w}, are calculated as follows:
#' \deqn{
#'  W \equiv \frac{6 (n^2 - K^2)}{n (4 * n - 1) (n + 1)}
#' }{
#'  w = 6 (k^2 - (K-1)^2) / (k (4 * k - 1) (k + 1))
#' }
#' where \eqn{n}{k} is the total number of tapers, and 
#' \eqn{K}{K} is the integer sequence \eqn{[0,n-1]}{[0,K-1]} 
#'
#' @export
#' @author A.J. Barbour <andy.barbour@@gmail.com> adapted original by R.L.Parker,
#' and authored the optimized version.
#' @seealso \code{\link{psdcore}}, \code{\link{riedsid}}
#'
#' @param tapvec \code{'tapers'} object; the number of tapers at each frequency
#' @param tap.index integer; the index of \code{tapvec} from which to produce a sequence of weights for
#' @param ntap integer; the number of tapers to provide weightings for.
#' @return A list with taper indices, and the weights \eqn{W_N}{w}.
#'
#' @example inst/Examples/rdex_parabolicweights.R
parabolic_weights <- function(tapvec, tap.index=1L) UseMethod("parabolic_weights")
#' @rdname parabolic_weights
#' @aliases parabolic_weights.tapers
#' @export
parabolic_weights.tapers <- function(tapvec, tap.index=1L){
  stopifnot(is.tapers(tapvec) | ((tap.index > 0L) & (tap.index <= length(tapvec))))
  kWeights <- parabolic_weights_fast(tapvec[as.integer(tap.index)])
  return(kWeights)
}
#' @rdname parabolic_weights
#' @export
parabolic_weights_fast <- function(ntap=1L) UseMethod("parabolic_weights_fast")

#' @rdname parabolic_weights
#' @export
parabolic_weights_fast.default <- function(ntap=1L){
  # Must be long-integer, otherwise overflow for ntap > 1e3
  K <- as.double(ntap)
  kseq <- seq_len(ntap) - 1
  #
  ksq <- kseq * kseq # vector
  K2 <- K * K   	   # scalar
  K3 <- K2 * K 		   # scalar
  #
  # orig: w = (tapers^2 - (k-1).^2) * (1.5/(tapers*(tapers-0.25)*(tapers+1)));
  # or:       (tapers^2 - (k-1).^2) * 3/(2*K3 + K2*3/2 - K/2)
  # or:
  w <- (K2 - ksq) * 6 / ( 4 * K3  +  3 * K2  -  K )
  #
  toret <- list(taper_seq = matrix(kseq + 1, ncol=ntap),
                taper_weights = matrix(w, ncol=ntap)
  )
  return(toret)
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
#'
minspan <- function(tapvec, ...) UseMethod("minspan")
#' @rdname tapers-constraints
#' @aliases minspan.tapers
#' @export
minspan.tapers <- function(tapvec, ...){
  stopifnot(is.tapers(tapvec))
  tapvec <- 7*tapvec/5
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
#'   \item \code{'loess.smooth'} uses \code{\link{ctap_loess}}
#'   \item \code{'none'} returns unbounded tapers.
#' }
#' 
#' @section Details of Constraint Methods:
#' 
#' \subsection{via first differencing (the default)}{
#' \code{\link{ctap_simple}} is the preferred constraint method.
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
#' \subsection{via LOESS smoothing}{
#' \code{\link{ctap_loess}} uses \code{loess} to smooth the taper vector; is
#' can be very slow thanks to quadratic scaling.
#' }
#'
#' @section Warning:
#'
#' \code{\link{ctap_loess}} results tend to be strongly dependent on
#' the tuning parameters given to \code{loess} (for obvious reasons); hence, 
#' some effort should be given to understand
#' their effect, and/or re-tuning them if needed.
#'
#' @rdname tapers-constraints
#' @title constrain_tapers
#' @aliases constrain_tapers tapers-constraints
#' @export
#' @keywords tapers tapers-constraints
#' @param tapvec \code{'tapers'} object; the number of tapers at each frequency
#' @param tapseq vector; positions or frequencies -- necessary for smoother methods
#' @param constraint.method  character; method to use for constraints on tapers numbers
#' @param verbose logical; should warnings and messages be given?
#' @param maxslope integer; constrain based on this maximum first difference
#' @param loess.span  scalar; the span used in \code{loess}
#' @param loess.degree  scalar; the polynomial degree
#' @param ... optional arguments (unused)
#' @return An object with class \code{'tapers'}.
#'
#' @example inst/Examples/rdex_constraintapers.R
constrain_tapers <- function(tapvec, tapseq=NULL,
                             constraint.method=c("simple.slope",
                                                 "loess.smooth",
                                                 "none"),
                             verbose=TRUE, ...) UseMethod("constrain_tapers")
#' @rdname tapers-constraints
#' @aliases constrain_tapers.tapers
#' @export
constrain_tapers.tapers <- function(tapvec, tapseq=NULL,
                                   constraint.method=c("simple.slope",
                                                       "loess.smooth",
                                                       "none"),
                                   verbose=TRUE, ...){
  stopifnot(is.tapers(tapvec))
  # choose the appropriate method to apply taper constraints
  cmeth <- match.arg(constraint.method) 
  tapvec.adj <- if (cmeth=="none"){
    if (verbose) warning("no taper optimization constraints applied")
    tapvec
  } else {
    if (verbose) message(sprintf("Constraining tapers with  ...  %s  ...  method",cmeth))
    CTAPFUN <- switch(cmeth,
                      "simple.slope"=ctap_simple,
                      "loess.smooth"=ctap_loess)
    CTAPFUN(tapvec, tapseq, ...)
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
#' @export
ctap_simple <- function(tapvec, tapseq=NA, maxslope=1, ...) UseMethod("ctap_simple")
#' @rdname tapers-constraints
#' @aliases ctap_simple.tapers
#' @export
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
#' @export
ctap_loess <- function(tapvec, tapseq=NULL, loess.span=.3, loess.degree=1, verbose=TRUE, ...){ UseMethod("ctap_loess") }
#' @rdname tapers-constraints
#' @export
ctap_loess.tapers <- function(tapvec, tapseq=NULL, loess.span=.3, loess.degree=1, verbose=TRUE, ...){
  stopifnot(is.tapers(tapvec))
  # having an x-sequence is absolutely critical to obtaining useful results
  if (is.null(tapseq)){
    tapseq <- seq_along(tapvec)
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
ctap_markov <- function() UseMethod("ctap_markov")
#' @rdname tapers-constraints
ctap_markov.tapers <- function() .Defunct("ctap_simple", package="psd")

#' @rdname tapers-constraints
ctap_friedman <- function(){ UseMethod("ctap_friedman") }
ctap_friedman.tapers <- function()  .Defunct("ctap_simple", package="psd")
