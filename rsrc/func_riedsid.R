#' Riedel & Sidorenko taper optimization
#' 
#' Estimates optimal number of tapers at each frequency of
#' given spec, based on Riedel-Sidorenko MSE recipe, and other
#' tweaks by RLP.
#' 
#' @title riedsid
#' @export
#' 
#' @seealso \code{\link{psdcore}}, \code{\link{constrain_tapers}}
#' @example
#' riedsid(rnorm(10), 10)
#' riedsid(rnorm(10), c(1:5,5:1))
riedsid <- function(...) UseMethod(".riedsid")

#' @return \code{NULL}
#' @rdname riedsid
#' @docType methods
#' @method riedsid default
#' @S3method riedsid default
.riedsid.default <- function(spec, ntaper, c.method=NULL) {
  ###
  # PORT of RLP's riedsid.m
  # abarbour
  # Dec 2011
  #
  # porting:  Jan 3, 2012
  # testing:
  ###
  #
  #  Estimates optimal number of tapers at each frequency of
  #  given spec, based on Riedel-Sidorenko MSE recipe and other
  #  tweaks due to RLP.  
  #  ntaper is the vector of taper lengths in the previous iteration.
  #
  ##
  ## Args:	
  ##
  ## Returns:	
  ##
  ## TODO(abarbour):
  ##

  ##message("\t\ttaper optimization")
  ##
  eps <- .Machine$double.eps # was: 1e-78  #  A small number to protect against zeros
  # spectral values
  spec <- as.vector(spec)
  # num freqs
  nf <- envAssignGet("num_freqs", length(spec))
  # prelims
  Ones <- ones(nf) # was rowvec, now col 
  Zeros <- zeros(nf) # was rowvec, now col
  # vectorize initial estimate
  if (length(ntaper)==1){
    ntap <- ntaper*Ones
  } else {
    ntap <- ntaper
  }
  # find the minimum by column for 1/2 nf, 7/5 ntap
  # ones is colvec, ntap is colvec
  require(matrixStats)
  # rowMins produces a rowvec of rowwise minimums; convert to colvec
  nspan <- matrix(matrixStats::rowMins(cbind(nf*Ones/2, 7*ntap/5)), ncol=1)
  
  #  Create log spec, and pad to handle begnning and end values
  nadd <- 1 + max(nspan)
  # [spec(nadd:-1:2); spec; spec(nf-1:-1:nf-nadd)]
  Y <- log(eps + c(spec[nadd:2], spec, spec[(nf-1):(nf-nadd)]))
  #  R <- spec"/spec <- Y" + (Y')^2  2nd form preferred for consistent smoothing
  dY <- d2Y <- Zeros # was row: t(Zeros)
  nadd1 <- nadd-1
  for (  j  in  1:nf ) {
    j1 <- j-nspan[j]+nadd1 
    j2 <- j+nspan[j]+nadd1 
    #  Over an interval proportional to taper length, fit a least
    #  squares quadratic to Y to estimate smoothed 1st, 2nd derivs
    jr <- j1:j2
    u <- jr - (j1 + j2)/2 # rowvec 
    u2 <- u*u             # rowvec
    L <- j2-j1+1          # constant
    L2 <- L*L             # constant
    LL2 <- L*L2           # constant
    LL2L <- LL2 - L       # constant
    uzero <- (L2 - 1)/12  # constant
    # first deriv
    dY[j] <- u %*% Y[jr] * 12 / LL2L
    # second deriv (?)
    d2Y[j] <- (u2 - uzero) %*% Y[jr] * 360/LL2L/(L2-4)
  }
  #
  #  Riedel-Sidorenko recipe (eq 13): kopt <- (12*abs(spec ./ d2spec)).^0.4 but
  #  parabolic weighting in psdcore requires: (480)^0.2*abs(spec./d2spec)^0.4
  #
  #  Original form:  kopt <- 3.428*abs(spec ./ d2spec).^0.4
  kopt <- round( 3.428 / abs(eps + d2Y + dY*dY) ^ 0.4 )
  #
  kseq <- seq.int(1, nf, by=1)
  print(c.method)
  kopt.bound <- constrain_tapers(kopt, kseq, c.method)
  ##
  return(invisible(kopt.bound))
} 
# end riedsid.default

#' Constrain the number of tapers applied
#' @export
#' @seealso \code{\link{riedsid}}
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
  stopifnot(min.tapers>=0)
  # choose the appropriate method to apply taper constraints
  cmeth <- match.arg(constraint.method) 
  if (cmeth=="none"){
    if (verbose) warning("no taper optimization constraints applied")
    tapvec.adj <- tapvec
  } else {
    if (verbose) message(sprintf("Constraining tapers with ** %s ** method",cmeth))
    ctapmeth <- paste(".ctap", switch(cmeth, 
                                      simple.slope="simple",
                                      "markov.chain"="markov",
                                      "loess.smooth"="loess",
                                      "friedman.smooth"="friedman"), sep="_")
    # execute alias method-function
    TAPFUN <- function(tapvec, tapseq, ...) UseMethod(ctapmeth)
    tapvec.adj <- TAPFUN(tapvec, tapseq, ...)
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


###
###  Constraint methods
###


#' Constrain tapers with slopes from first differencing
#' @return \code{NULL}
#' @rdname constrain_tapers
#' @docType methods
#' @method constrain_tapers default
#' @S3method constrain_tapers default
#' @param tapvec vector; the number of tapers for each frequency
#' @param tapseq (unused)
#' @param maxslope integer; constrain based on this maximum first difference
.ctap_simple.default <- function(tapvec, 
                                 tapseq=NULL, 
                                 maxslope=1L){
  # tapseq not needed
  tapvec <- as.numeric(tapvec)
  maxslope <- as.numeric(maxslope) #max(1L, ceiling(maxslope)))
  # c-code used for speed up of forward+backward operations
  # until it's packaged, need to dynamic load:
  owd <- getwd()
  src <- "/Users/abarbour/kook.processing/R/dev/packages/rlpSpec/src"
  setwd(src)
  system("rm riedsid.*o")
  system("R CMD SHLIB riedsid.c")
  dyn.load("riedsid.so")
  setwd(owd)
  tapvec.adj <- as.matrix(.Call("rlp_constrain_tapers", tapvec, maxslope))
  #as.matrix(.Call("rlp_constrain_tapers", tapvec, maxslope, PACKAGE = "rlpSpec"))
  return(tapvec.adj)
}

#' Constrain tapers with a Markov Chain, based on quantum-well probability
#' @return \code{NULL}
#' @rdname constrain_tapers
#' @docType methods
#' @method constrain_tapers default
#' @S3method constrain_tapers default
#' @param tapvec vector; the number of tapers for each frequency
#' @param (unused)
#' @param chain.width the width the MC should consider for the change probability
.ctap_markov.default <- function(tapvec, 
                                 tapseq=NULL, 
                                 chain.width=3*length(tapvec)){
  #stopifnot(("Peaks" %in% rownames(installed.packages())))
  require("Peaks")
  # tapseq not needed
  MC.win <- max(1, round(chain.width))
  tapvec.adj <- as.matrix(Peaks::SpectrumSmoothMarkov(tapvec, MC.win))
  return(tapvec.adj)
}

#' Constrain tapers with the Loess smoother (slow)
#' @return \code{NULL}
#' @rdname constrain_tapers
#' @docType methods
#' @method constrain_tapers default
#' @S3method constrain_tapers default
#' @param tapvec vector; the number of tapers for each frequency
#' @param tapseq vector; positions (or frequencies) necessary for smoother methods
#' @param loess.span span used in loess smoothing
#' @param loess.degree the polynomial degree
#' @seealso \code{\link{loess}}
.ctap_loess.default <- function(tapvec, 
                                tapseq, 
                                loess.span=.3, 
                                loess.degree=1){
  # having an x-sequence is absolutely critical to obtaining useful results
  if (is.null(tapseq)){
    tapseq <- 1:length(tapvec)
    warning("Generated a position sequence; results may be bogus.")
  }
  require(stats)
  warning("Loess-method has quadratic memory scaling (1e3 pt -> 10 Mb)...")
  trc <- ifelse(length(tapvec)>=1000, "approximate", "exact")
  loe <- stats::loess(y ~ x, 
                      data.frame(x=tapseq, y=tapvec), 
                      span=loess.span, degree=loess.degree,
                      control = loess.control(trace.hat = trc))
  tapvec.adj <- as.matrix(stats::predict(loe))
  return(tapvec.adj)
}

#' Constrain tapers with the Friedman 'super-smoother'
#' @return \code{NULL}
#' @rdname constrain_tapers
#' @docType methods
#' @method constrain_tapers default
#' @S3method constrain_tapers default
#' @param tapvec vector; the number of tapers for each frequency
#' @param tapseq vector; positions (or frequencies) necessary for smoother methods
#' @param smoo.span fraction of the observations in the span of the running lines smoother
#' @param smoo.bass controls the smoothness of the fitted curve 
#' @seealso \code{\link{supsmu}}
.ctap_friedman.default <- function(tapvec, 
                                   tapseq, 
                                   smoo.span=.3, 
                                   smoo.bass=2){
  # having an x-sequence is absolutely critical to obtaining useful results
  if (is.null(tapseq)){
    tapseq <- 1:length(tapvec)
    warning("Generated a position sequence; results may be bogus.")
  }
  require(stats)
  # Friedman's super smoother.
  # Meh.
  tapvec.adj <- as.matrix(stats::supsmu(tapseq, tapvec, span=smoo.span, bass=smoo.bass)$y)
  return(tapvec.adj)
}
###