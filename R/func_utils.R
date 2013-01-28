#' @title Various utility functions.
#'
#' @description Various utility functions are used:
#'
#' @keywords methods S3methods utilities
#' @author A.J. Barbour <andy.barbour@@gmail.com>
#' @rdname rlpSpec-utilities
#' @name rlpSpec-utilities
NULL

#' @description \code{vardiff} reeturns the variance of the first difference of the series.
#' @rdname rlpSpec-utilities
#' @export
#' @keywords utilities first-difference variance
#' @seealso \code{\link{rlpSpec-package}}
#' @param Xd object to difference
#' @return numeric
vardiff <- function(Xd){stats::var(diff(Xd))}

#' @description \code{dB} returns an object converted to decibels.
#' @details Decibels are defined as \deqn{10 \log{}_{10} x}.
#' @rdname rlpSpec-utilities
#' @param Rat numeric; A ratio to convert to decibels (\code{dB}).
#' @param invert logical; assumes \code{Rat} is already in decibels, so return ratio
#' @param pos.only logical; if \code{invert=FALSE}, sets negative or zero values to NA

#' @return numeric
#' @export
#' @aliases decibels db
#' @keywords utilities normalization decibel
#' @seealso \code{\link{rlpSpec-package}}
#' @examples
#' dB(1) # signal is equal <--> zero dB
dB <- function(Rat, invert=FALSE, pos.only=TRUE){
  CC <- 10
  if (invert) {
    toret <- 10**(Rat/CC)
  } else {
    if (pos.only) Rat[Rat <= 0] <- NA
    toret <- CC*log10(Rat)
  }
  return(toret)
}

#' @description \code{char2envir} converts a character string of an environment 
#' name to an evaluated name; whereas, \code{envir2char} converts an environment 
#' name to a character string.
#' 
#' @note \code{char2envir} ensures the \code{envchar} object is a character, 
#' so that something
#' is not unintentionally evaluated; \code{envir2char} simply deparses the
#' object name.
#' 
#' @rdname rlpSpec-utilities
#' @export
#' @keywords utilities environment
#' @param envchar An object with class 'character'.
#' @param envir An object of class 'environment'.
#' @return \code{char2envir} returns the result of evaluating the object: an 
#' environment object; \code{envir2char} returns the result of deparsing the 
#' environment name: a character string.
#' @seealso \code{\link{rlpSpec-package}}
#' @examples
#' ##
#' ## Evaluate character strings
#' print(.rlpenv)
#' char2envir(.rlpenv)
#' char2envir("some nonexistent environment") # error
char2envir <- function(envchar){
  stopifnot(is.character(envchar))
  eval(as.name(envchar))
}
#' @rdname rlpSpec-utilities
#' @export
#' @keywords utilities environment
#' @examples
#' # and environment objects:
#' print(.GlobalEnv)
#' envir2char(.GlobalEnv)
#' envir2char(.rlpSpecEnv)
#' char2envir(some_nonexistent_environment) # error
#' ##
#' ##
envir2char <- function(envir){
  stopifnot(is.environment(envir))
  deparse(substitute(envir))
}

#' @description \code{vector_reshape} reshapes a vector into another vector.
#' @rdname rlpSpec-utilities
#' @name vector_reshape
#' @param x  An object to reshape (\code{vector_reshape}).
#' @param vec.shape  choice between horizontally-long or vertically-long vector.
#' @return \code{vector_reshape} returns a "reshaped" vector, meaning it has
#' had it's dimensions changes so that it has either one row 
#' (if \code{vec.shape=="horizontal"}), or one column (\code{"vertical"}).
#' @export
#' @keywords utilities vector-manipulation matrix-manipulation
vector_reshape <- function(x, vec.shape=c("horizontal","vertical")) UseMethod("vector_reshape")
#' @rdname rlpSpec-utilities
#' @method vector_reshape default
#' @S3method vector_reshape default
vector_reshape.default <- function(x, vec.shape=c("horizontal","vertical")){
  x <- as.vector(x)
  vec.shape <- match.arg(vec.shape)
  nrow <- switch(vec.shape, "horizontal"=1, "vertical"=length(x))
  return(matrix(x, nrow=nrow))
}

#' @description \code{colvec} returns the object as a vertically long vector; whereas
#' \code{rowvec} returns the object as a horizontally long vector.
#' @details \code{colvec, rowvec} are simple wrapper functions to \code{vector_reshape}.
#' @rdname rlpSpec-utilities
#' @export
#' @keywords utilities vector-manipulation matrix-manipulation
colvec <- function(x) vector_reshape(x, "vertical")

#' @rdname rlpSpec-utilities
#' @aliases as.rowvec
#' @export
#' @keywords utilities vector-manipulation matrix-manipulation
rowvec <- function(x) vector_reshape(x, "horizontal")

#' @description \code{is.spec} reports whether an object has class S3 class 'spec', as
#' would one returned by, for example, \code{spectrum}.
#' @rdname rlpSpec-utilities
#' @param Obj  An object to test for class inheritance.
#' @return \code{is.spec} and \code{is.tapers} both return
#' logicals about whether or not the object does have class 'spec' or 'tapers', 
#' respectively
#' @export
#' @keywords utilities inherits is
#' @examples
#' ## Check for spec object:
#' # quick power spectral density
#' x <- rnorm(1e2, sd=10)
#' psd <- spectrum(x, plot=FALSE)
#' # return is class 'spec'
#' is.spec(psd) # TRUE
#' #
#' # but the underlying structure is just a list
#' psd <- unclass(psd)
#' is.spec(psd) # FALSE
#' ##
is.spec <- function(Obj) inherits(Obj, "spec")

#' @description \code{is.tapers} reports whether an object has S3 class 'tapers', as
#' would one returned by, for example, \code{\link{as.tapers}}.
#' @export
#' @keywords utilities inherits is
#' @rdname rlpSpec-utilities
#' @seealso \code{\link{as.tapers}}
#' @examples
#' ## check for tapers
#' is.tapers(x)
#' is.tapers(as.tapers(x))
#' ##
is.tapers <- function(Obj) inherits(Obj, "tapers")

#' Numerical derivatives of a series based on a weighted, smooth spline representation.
#' 
#' @description \code{splineGrad} computes the numerical derivatives of a spline 
#' representation of the input series; differentiation of spline curves is 
#' numerically efficient.
#' 
#' With smoothing, the numerical instability for "noisy" data can be drastically
#' reduced, 
#' since spline curves are
#' inherently (at least) twice differentiable. See the \strong{Examples} for
#' an illustration of this.
#' 
#' How does the first derivative of the first derivative of the
#' original series compare to the second derivative of the original series? 
#' Apparently, there is no difference.  
#'
# Try something like this:
# \code{smspl.alt <- stats::smooth.spline(dseq, fsigderiv, ...);
#  SPLFUN.alt <- stats::splinefun(smspl.alt$x, smspl.alt$y);
#  fsigderiv2.alt <- SPLFUN.alt(dseq, deriv=1);
#  print(all.equal(fsigderiv2,fsigderiv2.alt));}
#' 
#' @author A.J. Barbour <andy.barbour@@gmail.com>
#' @name splineGrad
#' @param dseq  numeric; a vector of positions for \code{dsig}.
#' @param dsig  numeric; a vector of values (which will have a spline fit to them).
#' @param plot.derivs  logical; should the derivatives be plotted?
#' @param ... additional arguments passed to \code{smooth.spline}.
#' @return A matrix with columns representing \eqn{x, f(x), f'(x), f''(x)}.
#' @export
#' @keywords utilities spline-gradient numerical-derivative
#' @seealso \code{smooth.spline}
#' @example x_examp/splinegrad.R
splineGrad <- function(dseq, dsig, plot.derivs=FALSE, ...) UseMethod("splineGrad")
#' @rdname splineGrad
#' @method splineGrad default
#' @S3method splineGrad default
splineGrad.default <- function(dseq, dsig, plot.derivs=FALSE, ...){
  #
  # Use spline interpolation to help find an emprirical gradient
  # (reduces numerical instability)
  #
  # @dseq: the sequence (index) for @dsig, the signal
  # output is the same length as the input
  #
  # create a weighted cubic spline
  smspl <- stats::smooth.spline(dseq, dsig, ...)
  # and a function
  SPLFUN <- stats::splinefun(smspl$x, smspl$y)
  # ?splinefun:
  # splinefun returns a function with formal arguments x and deriv, 
  # the latter defaulting to zero. This function can be used to 
  # evaluate the interpolating cubic spline (deriv = 0), or its 
  # derivatives (deriv = 1, 2, 3) at the points x, where the spline 
  # function interpolates the data points originally specified. This 
  # is often more useful than spline.
  #
  #   seq.rng <- range(dseq)
  #   from <- seq.rng[1]
  #   to <- seq.rng[2]
  #   n <- length(dseq)
  #
  # signal spline
  #fsig <<- SPLFUN(dseq)
  #   FD0 <- function(){graphics::curve(SPLFUN(x), from=from, to=to, n=n, add=NA)}
  #   fsig <<- FD0()
  #
  # first deriv
  #   FD1 <- function(){graphics::curve(SPLFUN(x,deriv=1), from=from, to=to, n=n, add=NA)}
  #   fsigderiv <<- FD1()
  fsigderiv <- SPLFUN(dseq, deriv=1)
  #
  # second deriv
  #   FD2 <- function(){graphics::curve(SPLFUN(x,deriv=2), from=from, to=to, n=n, add=NA)}
  #   fsigderiv2 <<- FD2()
  fsigderiv2 <- SPLFUN(dseq, deriv=2)
  # how does the first deriv of the first deriv compare to the second?
  #   smspl.alt <- stats::smooth.spline(dseq, fsigderiv, ...)
  #   SPLFUN.alt <- stats::splinefun(smspl.alt$x, smspl.alt$y)
  #   fsigderiv2.alt <<- SPLFUN.alt(dseq, deriv=1)
  #   print(all.equal(fsigderiv2,fsigderiv2.alt))
  # [1] TRUE
  #   plot(fsigderiv2,fsigderiv2.alt, asp=1)
  ##
  toret <- data.frame(x=dseq, y=dsig, 
                      dydx=fsigderiv, 
                      d2yd2x=fsigderiv2)
                      #d2yd2x.alt=fsigderiv2.alt)
  #
  if (plot.derivs){
    #     yl.u <- max(c(dsig,fsigderiv,fsigderiv2))#,fsigderiv2.alt))
    #     yl.l <- min(c(dsig,fsigderiv,fsigderiv2))#,fsigderiv2.alt))
    nr <- 3 # f, f', f''
    mar.multi <- c(2., 5.1, 2, 2.1)
    oma.multi <- c(6, 0, 5, 0)
    oldpar <- par(mar = mar.multi, oma = oma.multi, mfcol = c(nr, 1))
    on.exit(par(oldpar))
    par(las=1)
    plot(dseq, dsig, cex=0.6, pch=3,
         #ylim=1.1*c(yl.l,yl.u),
         xaxs="i", yaxs="i",
         xlab="x", ylab="f(x)",
         main=sprintf("splineGrad: signal and weighted cubic-spline fit (spar = %g)",smspl$spar))
    lines(y ~ x, toret, col="dark grey", lwd=1.5)
    plot(dydx ~ x, toret, 
         xaxs="i", yaxs="i",
         main="first derivative",
         col="red", type="s", lwd=2.4)
    plot(d2yd2x ~ x, toret, 
         xaxs="i", yaxs="i",
         main="second derivative", xlab="x",
         col="blue", type="s", lwd=2.4, lty=3)
    #lines(d2yd2x.alt ~ x, toret, col="blue", type="s", lwd=2.4, lty=3)
    #     legend("topleft", 
    #            paste(c("weighted cubic spline fit","first deriv", "second deriv"), #, "deriv of first deriv"), 
    #                  sep=''), 
    #            col = c("grey","red","blue"), #,"blue"), 
    #            lty = c(rep(1,3)), #, 3),
    #            cex=0.7)
  }
  return(invisible(toret))
}
##

#' @description \code{na_mat} populates a matrix of specified dimensions 
#' with \code{NA} values.
#' @rdname rlpSpec-utilities
#' @param nrow integer; the number of rows to create.
#' @param ncol integer; the number of columns to create (default 1).
#' @return \code{na_mat} returns a matrix of dimensions \code{(nrow,ncol)} with
#' \code{NA} values, the representation of which is set by \code{NA_real_}
#' @export
#' @keywords utilities vector-creation matrix-creation
#' @examples
#' ## matrix and vector creation:
#' # NA matrix
#' nd <- 5
#' na_mat(nd)
#' na_mat(nd,nd-1)
na_mat <- function(nrow, ncol=1) UseMethod("na_mat")
#' @rdname rlpSpec-utilities
#' @method na_mat default
#' @S3method na_mat default
na_mat.default <- function(nrow, ncol=1){matrix(NA_real_, nrow, ncol)}
#' @description \code{zeros} populate a column-wise matrix with zeros; whereas,
#' \code{ones} populates a column-wise matrix with ones
#' @rdname rlpSpec-utilities
#' @export 
#' @keywords utilities vector-creation matrix-creation
# params described by na_mat
#' @return For \code{zeros} or \code{ones} respectively, a matrix vector 
#' with \code{nrow} zeros or ones.
#' @examples
#' # zeros
#' zeros(nd)
#' zeroes(nd)
zeros <- function(nrow) UseMethod("zeros")
#' @rdname rlpSpec-utilities
#' @method zeros default
#' @S3method zeros default
zeros.default <- function(nrow){stopifnot(!is.null(nrow)); matrix(rep.int(0, nrow), nrow=nrow)}

#' @rdname rlpSpec-utilities
#' @export
#' @keywords utilities vector-creation matrix-creation
#' @examples
#' # and ones
#' ones(nd)
#' ##
#' ##
ones <- function(nrow) UseMethod("ones")
#' @rdname rlpSpec-utilities
#' @method ones default
#' @S3method ones default
ones.default <- function(nrow){stopifnot(!is.null(nrow)); matrix(rep.int(1, nrow), nrow=nrow)}

#' @description \code{mod} finds the modulo division of X and Y.
#' 
#' @details Modulo division has higher order-of-operations ranking than other
#' arithmetic operations; hence, \code{x + 1 \%\% y} is equivalent to
#' \code{x + (1 \%\% y)} which can produce confusing results. \code{mod}
#' is simply a series of \code{trunc} commands which
#' reduces the chance for unintentionally erroneous results.
#' 
#' @note The performance of \code{mod} has not been tested against the 
#' \code{\%\%} arithmetic method -- it may or may not be slower for large
#' numeric vectors.
#' 
#' @references For \code{\link{mod}}: see Peter Dalgaard's explanation of 
#' the non-bug (#14771) I raised (instead I should've asked it on R-help): 
#' \url{https://bugs.r-project.org/bugzilla3/show_bug.cgi?id=14771\#c2}
#' 
#' @param X numeric; the "numerator" of the modulo division
#' @param Y numeric; the "denominator" of the modulo division
#' @return \code{mod} returns the result of a modulo division, which is 
#' equivalent to \code{(X) \%\% (Y)}.
#' @export
#' @keywords utilities modulo-division arithmetic-operations
#' @rdname rlpSpec-utilities
#' @aliases modulo
#' @example x_examp/mod.R
mod <- function(X, Y) UseMethod("mod")
#' @rdname rlpSpec-utilities
#' @method mod default
#' @S3method mod default
mod.default <- function(X, Y){
  stopifnot(is.numeric(c(X, Y)))
  ## modulo division
  X1 <- trunc( trunc(X/Y) * Y)
  Z <- trunc(X) - X1
  return(Z)
}
