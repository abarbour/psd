#' The 'taper' S4 class.
#'
#' In this class the value of each position will be a non-zero, positive
#' integer, since it represents the number of tapered sections to average.
#'
# \section{Slots}{
#   \describe{
#     \item{\code{tapers}:}{Object of class \code{"integer"}, containing 
#     a vector of tapers.}
#   }
# }
#'
#' @note  The prototypical S4 class has tapers==1, and length==1.
#' Currently there are no \code{@@slots}; this may change in the future.
#'
#' @name taper
#' @rdname taper
#' @aliases taper-class
#' @exportClass taper
#' @author Andrew Barbour <andy.barbour@@gmail.com>
#' @examples
#' taper()
#' new("taper") # equivalent to taper()
#' print(ntap <- taper(1:10))
#' plot(ntap)
taper <- setClass("taper",
                  # if slots, add 'taper="integer",...
                  representation=representation("integer"),
                  prototype = 1L)
         
###
###  Check/conversion methods
###

#' Reports whether x is a 'taper' object
#' @param x An object to test
#' @export
#' @author Andrew Barbour <andy.barbour@@gmail.com>
#' @seealso \code{\link{as.taper}}, \code{\link{taper}}
#' @examples
#' is.taper(1:10)
#' is.taper(as.taper(1:10))
is.taper <- function(x) inherits(x, "taper")

#' Coerces an object into a 'taper' object
#'
#' An object of class 'taper' is created, with
#' a minimum number of tapers in each position
#' set by \code{min.taper}.
#'
#' Various classes can be coerced into a 'taper' object; those
#' tested sofar include: scalar, vector, matrix, data.frame, 
#' and list.  
#'
#' Multiple objects are concatenated into a single
#' vector dimension.  For example, if the object is 
#' \code{list(x=c(1,2),y=c(3,4,5,0,1.1))} then the corresponding 'taper'
#' object
#' will be \code{1,2,3,4,5,1,1}, assuming \code{min.taper==1}.
#'
#' @param x An object to set
#' @param min.taper Set all values less than this to this.
#' @export
#' @author Andrew Barbour <andy.barbour@@gmail.com>
#' @seealso \code{\link{is.taper}}, \code{\link{taper}}
#' @examples
#' is.taper(as.taper(1))
#' is.taper(as.taper(1:10))
#' is.taper(as.taper(matrix(1:10,ncol=1)))
#' is.taper(as.taper(list(x=1:10,y=1:30))) # note dimensions
#' is.taper(as.taper(data.frame(x=1:10,y=10:19)))
#' # class 'character' is in-coercible; raise error
#' as.taper(c("a","b"))
as.taper <- function(x, min.taper=1){
  # taper should be non-zero integer, since it represents the
  # number of tapered sections to average; hence, ceiling
  # then integerize
  stopifnot(!(is.character(x)) | (min.taper>0))
  x <- unlist(x)
  x[x < min.taper] <- min.taper
  #   > as.integer(as.matrix(data.frame(x=1:10,y=10:19)))
  #   [1]  1  2  3  4  5  6  7  8  9 10 10 11 12 13 14 15 16 17 18 19
  x <- as.integer(ceiling(x))
  class(x) <- "taper"
  return(x)
}

###
###  Generic methods
###

#' @title Generic methods for 'taper' objects
#' @keywords methods generic
#' @name taper-methods
#' @author Andrew Barbour <andy.barbour@@gmail.com>
#' @aliases taper
#' @rdname taper-methods
#' @seealso \code{\link{is.taper}}, \code{\link{as.taper}}
NULL

#' @rdname taper-methods
#' @name print
#' @S3method print taper
print.taper <- function(x){
  stopifnot(is.taper(x))
  xh <- paste(as.character(head(x)), collapse=" ")
  xt <- paste(as.character(tail(x)), collapse=" ")
  cat(sprintf("taper object:\n\thead:  %s\n\t\t...\n\ttail:  %s\n",xh,xt))
}
#S4 from S3
#' @rdname taper-methods
#' @name print
#' @export
#' @docType methods
setMethod("print", signature("taper"), print.taper)

#' @rdname taper-methods
#' @name summary
#' @S3method summary taper
summary.taper <- function(x){
  stopifnot(is.taper(x))
  toret <- summary.default(x)
  class(toret) <- "summary.taper"
  return(toret)
}

#' @rdname taper-methods
#' @name print
# @aliases print.summary.taper
#' @S3method print summary.taper
print.summary.taper <- function(x){
  cat("summary of tapers:\n")
  print(summary(x))
}

#' @rdname taper-methods
#' @name plot
#' @S3method plot taper
plot.taper <- function(x, color.pal=c("Blues","Spectral"), ...){
  stopifnot(is.taper(x))
  #if isS4(x) x <- x@tapers
  require(graphics, grDevices, RColorBrewer)
  nt <- length(x)
  xi <- 1:nt
  mx <- max(x)
  pal <- match.arg(color.pal)
  npal <- switch(pal, RdYlBu=11, Spectral=11, Blues=9)
  pal.col <- RColorBrewer::brewer.pal(npal, pal)
  cols <- grDevices::colorRampPalette(pal.col)(mx)
  graphics::plot.default(xi, x,
               main="",
               ylab="number of tapers",
               xlab="taper index",
               ylim=1.15*c(0.5, mx), yaxs="i", 
               xlim=c(-1, nt+2), xaxs="i",
               lwd=1.8,
               type="h",
               col=cols[x],
               ...)
  lines(xi,x,col="black",lwd=0.9)
  # plot log2 multiples as horiz lines
  hl <- 2**(1:round(log2(mx)))
  graphics::abline(h=hl,lty=1,lwd=0.6,col="black")
  vl <- c(1, nt)
  graphics::abline(v=vl,lty=3,lwd=2,col="blue")
}

#' Calculate weighting factors for a series of tapers
#'
#' Weighting is calculated as follows:
#'
#' \deqn{n_T^2 - \frac{3 \cdot K_N^2}{2 \cdot n_T \cdot (n_T - 1/4) \cdot (n_T + 1)}}
#'
#' where \eqn{n_T} is the total number of tapers, and 
#' \eqn{K_N} is the integer sequence \eqn{[0,n_T-1]} 
#'
# w = (tapers^2 - (k-1).^2)*(1.5/(tapers*(tapers-0.25)*(tapers+1)));
#'
#' @title parabolic_weights
#' @export
#' @author Andrew Barbour <andy.barbour@@gmail.com> ported original by R.L.Parker,
#' and authored the optimized version.
#' @seealso \code{\link{psdcore}}, \code{\link{riedsid}}
#'
#' @param tapvec 'taper' object; the number of tapers at each frequency
#' @param tap.index integer; the index of \code{tapvec} from which to find weights
#' @return a list with taper indices, and weighting parameters
parabolic_weights <- function(tapvec, tap.index=1L) UseMethod("parabolic_weights")

#' @rdname parabolic_weights
#' @S3method parabolic_weights taper
parabolic_weights.taper <- function(tapvec, tap.index=1L){
  stopifnot(is.taper(tapvec) | ((tap.index > 0L) & (tap.index <= length(tapvec))))
  kWeights <- parabolic_weights_fast(tapvec[as.integer(tap.index)])
  return(kWeights)
}

#' @param ntap integer; the number of tapers to provide weightings for
#' @rdname parabolic_weights
#' @export
parabolic_weights_fast <- function(ntap=1L){
  ntap <- max(1, as.integer(ntap))
  kseq <- seq.int(from=0, to=ntap-1, by=1)
  lk <- length(kseq)
  K2 <- kseq * kseq # vector
  NT2 <- ntap * ntap # scalar
  NT3 <- NT2 * ntap # scalar
  TW <- matrix((NT2 - K2) * 3/(2*NT3 + NT2*3/2 - ntap/2), ncol=lk)
  return(list(taper_seq=kseq+1, taper_weights=TW))
}

###
###  Constraint methods
###

#' @description Find the max span of taper object
#'
#' @details minimize the number of tapers: min(nf/2, 7*tapvec/5)
#'
#' @note The value of \code{nf} should represent the total number of
#' tapers from the original spectrum; the default assumes the taper
#' object represents that size.
#'
#' @title minspan
#' @export
#' @seealso \code{\link{constrain_tapers}}
#' @author Andrew Barbour <andy.barbour@@gmail.com> ported original by R.L.Parker.
#'
#' @param tapvec 'taper' object; the number of tapers at each frequency
#' @param nf scalar; number of positions or frequencies
#' @param ... (unused) optional argments
#' @return an object with class 'taper', with a constrained taper numbers
#' @examples
#' \dontrun{
#' ntap <- as.taper(1:10)
#' par(mfrow=c(2,1))
#' plot(ntap)
#' plot(minspan(ntap))
#' par(mfrow=c(1,1))
#' }
minspan <- function(tapvec, 
                    nf=length(tapvec), 
                    ...) UseMethod("minspan")

#' @rdname minspan
#' @S3method minspan taper
minspan.taper <- function(tapvec, nf=length(tapvec), ...){
  stopifnot(is.taper(tapvec))
  require(matrixStats)
  Ones <- ones(nf)
  nspan <- as.taper(matrixStats::rowMins(cbind(Ones*nf/2, 7*as.matrix(tapvec)/5)))
  return(nspan)
}

#' @description Apply constraints on the number of tapers
#' 
#' @details Refines the number of tapers; the method by which it does this
#' may be chosen by the user.
#' 
#' @title constrain_tapers
#' @aliases constrain_tapers
#' @export
#' @author Andrew Barbour <andy.barbour@@gmail.com> ported original by R.L.Parker to R and C.
#' @seealso \code{\link{ctap_simple}}, \code{\link{ctap_markov}}, \code{\link{riedsid}}
#'
#' @param tapvec vector; the number of tapers at each frequency
#' @param tapseq vector; positions or frequencies -- necessary for smoother methods
#' @param constraint.method method to use for constraints on tapers numbers
#' @param min.tapers integer; the minimum number of tapers
#' @param verbose logical; should warnings and messages be given?
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
                             ...) UseMethod("constrain_tapers")

#' @rdname constrain_tapers
#' @S3method constrain_tapers taper
constrain_tapers.taper <- function(tapvec, 
                                   tapseq=NULL,
                                   constraint.method=c("simple.slope",
                                                       "markov.chain",
                                                       "loess.smooth",
                                                       "friedman.smooth",
                                                       "none"),
                                   min.tapers=1,
                                   verbose=TRUE, 
                                   ...){
  stopifnot(is.taper(tapvec) | (min.tapers>=0))
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


#' @description Constrain tapers with first differencing.
#' 
#' @details This is the default, and preferred constraint method.
#' The algortihm uses first-differencing to modify the number
#' of tapers in the previous position.  Effectively, the constraint
#' is based on a causal, 1st-order Finite Impulse-response Filter (FIR) 
#' which makes the method
#' sensitive to rapid changes in the number of tapers; naturally,
#' smoother spectra tend to produce less fluctuation in taper numbers, 
#' which makes this well suited for adaptive processing. 
#'
#' In pure-R this algorithm can be very slow; however, here we have
#' included it as dynamically loaded c-code so it it reasonably fast.
#' 
#' @note The results obtained by \strong{\code{'simple.slope'}} are generally the most
#' stable, meaning repeatedly running the constraint will not change values
#' other than on the first execution; the same cannot be said for the other
#' methods.
#'
#' @title ctap_simple
#' @aliases constrain_taper_simple_slope
#' @export
#' @author Andrew Barbour <andy.barbour@@gmail.com> ported original by R.L.Parker to R and C.
#'
#' @seealso \code{\link{constrain_tapers}}
#' @param tapvec 'taper' object; the number of tapers at each frequency
#' @param tapseq (unused) vector; positions or frequencies -- necessary for smoother methods
#' @param maxslope (\code{'simple.slope'}) integer; constrain based on this maximum first difference
#' @param ... (unused) optional argments
#' @return an object with class 'taper', with a constrained taper numbers
#' @examples
#' \dontrun{
#' ntap <- as.taper(1:10)
#' ctap_simple(ntap)
#' }
ctap_simple <- function(tapvec,
                        tapseq=NA, 
                        maxslope=1,
                        ...) UseMethod("ctap_simple")

#' @rdname ctap_simple
#' @S3method ctap_simple taper
ctap_simple.taper <- function(tapvec,
                              tapseq=NA, 
                              maxslope=1,
                              ...){
  stopifnot(is.taper(tapvec))
  # tapseq not needed
  # 'taper' object gives integer values, but code requires real
  tapvec <- as.numeric(tapvec)
  maxslope <- as.numeric(maxslope)
  # c-code used for speed up of forward+backward operations
  # until it's packaged, need to dynamic load:
  owd <- getwd()
  src <- "/Users/abarbour/nute.processing/development/rlpSpec/src"
  #src <- "/Users/abarbour/kook.processing/R/dev/packages/rlpSpec/src"
  setwd(src)
  system("rm riedsid.*o")
  system("R CMD SHLIB riedsid.c")
  dyn.load("riedsid.so")
  setwd(owd)
  tapvec.adj <- as.taper(.Call("rlp_constrain_tapers", tapvec, maxslope))
  ##as.matrix(.Call("rlp_constrain_tapers", tapvec, maxslope, PACKAGE = "rlpSpec"))
  return(tapvec.adj)
}

#' @description Constrain tapers with a Markov Chain, based on quantum-well probability.
#' 
#' @details This method uses a Markov Chain, based on the theory of 
#' quantum-well probability, commonly used in gamma-ray spectroscopy.
#'
#' The main function behind this method is \link[Peaks]{SpectrumSmoothMarkov} 
#' which calculates the probability that the number of tapers would have
#' changed (from it's previous value); it is very fast.  Details of the theory 
#' behind this algorithm may be found in Morhac (2008) and Silagadze (1996).
#'
#' @note This algorithm can produce ``unstable"
#' results in the sense that for
#' successive application on taper vectors, even modest sized serially-correlated
#' peaks tends to sharpen; hence, this method should be used with care, unless the 
#' intention is to specifically
#' enhance peaks.  The \code{'chain.width'} parameter controls the broadness
#' of the a priori distribution.  As a rule of thumb: the smaller the parameter is, 
#' the shorter the tails become.
#'
#' @author Andrew Barbour <andy.barbour@@gmail.com> adapted the algorithm from Morhac (2008)
#' for use here.
#'
#' @references Morhac, M. (2008), Peaks: Peaks, \emph{R package}, \strong{version 0.2}
#' @references Silagadze, Z.K. (1996), A new algorithm for automatic photopeak searches,
#' \emph{Nucl. Instrum. Meth. A}, \strong{376} 451.
#' @references \url{http://arxiv.org/abs/hep-ex/9506013}
#'
#' @title ctap_markov
#' @aliases constrain_taper_markov_chain
#' @export
#' @seealso \code{\link{constrain_tapers}}, \code{\link[Peaks]{SpectrumSmoothMarkov}}
#' @param tapvec 'taper' object; the number of tapers at each frequency
#' @param tapseq (unused) vector; positions or frequencies -- necessary for smoother methods
#' @param chain.width  scalar; the width the MC should consider for the change probability
#' @param ... (unused) optional argments
#' @return an object with class 'taper', with a constrained taper numbers
#' @examples
#' \dontrun{
#' ntap <- as.taper(1:10)
#' ctap_markov(ntap)
#' }
ctap_markov <- function(tapvec, 
                        tapseq=NA, 
                        chain.width=round(5*length(tapvec)),
                        ...) { UseMethod("ctap_markov") }

#' @rdname ctap_markov
#' @S3method ctap_markov taper
ctap_markov.taper <- function(tapvec, 
                              tapseq=NA, 
                              chain.width=round(5*length(tapvec)),
                              ...){
  stopifnot(is.taper(tapvec))
  require("Peaks")
  # tapseq not needed
  # bound the chain.width
  MC.win <- max(1, round(chain.width))
  tapvec.adj <- as.taper(Peaks::SpectrumSmoothMarkov(as.numeric(tapvec), MC.win))
  return(tapvec.adj)
}

#' @description Constrain tapers with the Loess smoother.
#' 
#' @details This method uses \code{stats::loess} to smooth the taper vector; is
#' can be very slow thanks to quadratic scaling.
#'
#' @note The results obtained can be strongly dependent on
#' the \code{loess} tuning parameters; hence, some effort should be given to understand
#' their effect, and/or re-tuning the default parameters.
#'
#' @title ctap_loess
#' @aliases constrain_taper_loess_smooth
#' @export
#' @seealso \code{\link{constrain_tapers}}, \code{\link{loess}}
#' @author Andrew Barbour <andy.barbour@@gmail.com>
#'
#' @param tapvec 'taper' object; the number of tapers at each frequency
#' @param tapseq vector; positions or frequencies -- necessary for smoother methods
#' @param loess.span  scalar; the span used in \code{loess}
#' @param loess.degree  scalar; the polynomial degree
#' @param verbose logical; should warnings and messages be given?
#' @param ... optional argments; unused
#' @return an object with class 'taper', with a constrained taper numbers
#' @examples
#' \dontrun{
#' ntap <- as.taper(1:10)
#' ctap_loess(ntap)
#' }
ctap_loess <- function(tapvec, 
                       tapseq=NULL, 
                       loess.span=.3, 
                       loess.degree=1,
                       verbose=TRUE,
                       ...){ UseMethod("ctap_loess") }

#' @rdname ctap_loess
#' @S3method ctap_loess taper
ctap_loess.taper <- function(tapvec,
                             tapseq=NULL, 
                             loess.span=.3, 
                             loess.degree=1,
                             verbose=TRUE,
                             ...){
  stopifnot(is.taper(tapvec))
  # having an x-sequence is absolutely critical to obtaining useful results
  if (is.null(tapseq)){
    tapseq <- 1:length(tapvec)
    if (verbose) warning("Generated a position sequence; results may be bogus.")
  }
  require(stats)
  if (verbose) warning("Loess-method has quadratic memory scaling (1e3 pt -> 10 Mb)...")
  trc <- ifelse(length(tapvec)>=1000, "approximate", "exact")
  loe <- stats::loess(y ~ x, 
                      data.frame(x=tapseq, y=as.numeric(tapvec)), 
                      span=loess.span, degree=loess.degree,
                      control = loess.control(trace.hat = trc))
  tapvec.adj <- as.taper(stats::predict(loe))
  return(tapvec.adj)
}

#' @description Constrain tapers with the Friedman 'super-smoother'.
#' 
#' @details This method uses the base Fiedman super-smoother.
#' 
#' @note The results obtained by \strong{\code{'friedman.smooth'}} are generally poor; 
#' hence, the method may be removed in future releases.
#' 
#' @title ctap_friedman
#' @aliases constrain_taper_friedman_smooth
#' @export
#' @seealso \code{\link{constrain_tapers}}, \code{\link{supsmu}}
#' @author Andrew Barbour <andy.barbour@@gmail.com>
#'
#' @param tapvec 'taper' object; the number of tapers at each frequency
#' @param tapseq vector; positions or frequencies -- necessary for smoother methods
#' @param smoo.span  scalar; fraction of the observations in the span of the running lines smoother
#' @param smoo.bass  scalar; controls the smoothness of the fitted curve 
#' @param verbose logical; should warnings and messages be given?
#' @param ... optional argments; unused
#' @return an object with class 'taper', with a constrained taper numbers
#' @examples
#' \dontrun{
#' ntap <- as.taper(1:10)
#' ctap_friedman(ntap)
#' }
ctap_friedman <- function(tapvec, 
                          tapseq=NULL, 
                          smoo.span=.3, 
                          smoo.bass=2, 
                          verbose=TRUE,
                          ...){ UseMethod("ctap_friedman") }

#' @rdname ctap_friedman
#' @S3method ctap_friedman taper
ctap_friedman.taper <- function(tapvec, 
                                tapseq=NULL, 
                                smoo.span=.3, 
                                smoo.bass=2, 
                                verbose=TRUE,
                                ...){
  stopifnot(is.taper(tapvec))
  # having an x-sequence is absolutely critical to obtaining useful results
  if (is.null(tapseq)){
    tapseq <- 1:length(tapvec)
    if (verbose) warning("Generated a position sequence; results may be bogus.")
  }
  require(stats)
  # Friedman's super smoother.
  # Meh.
  tapvec.adj <- as.taper(stats::supsmu(tapseq, as.numeric(tapvec), span=smoo.span, bass=smoo.bass)$y)
  return(tapvec.adj)
}
###
