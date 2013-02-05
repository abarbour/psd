###
### Various utility functions
###
##
#' convert a character string of an environment name to an evaluated name
#' @rdname rlpspec-utilities
#' @param envchar object with class 'character'
#' @return The result of evaluating the object: an environment object
#' @examples
#' print(.rlpenv)
#' char2envir(.rlpenv) #
#' char2envir("some nonexistent environment") # error
char2envir <- function(envchar){
  stopifnot(is.character(envchar))
  eval(as.name(envchar))
}
#' convert an environment name to a character string
#' @rdname rlpspec-utilities
#' @param envir object of class 'environment'
#' @return The result of evaluating the object: an environment object
envir2char <- function(envir){
  stopifnot(is.environment(envir))
  deparse(substitute(envir))
}


vector_shape <- function(x, ...) UseMethod(".reshape_vector")
.reshape_vector.default <- function(x, vec.shape=c("horizontal","vertical")){
  x <- as.vector(x)
  vec.shape <- match.arg(vec.shape)
  nrow <- switch(vec.shape, "horizontal"=1, "vertical"=length(x))
  return(matrix(x, nrow=nrow))
}
as.colvec <- function(x) vector_shape(x, "vertical")
as.rowvec <- function(x) vector_shape(x, "horizontal")


#' Reports whether x is an object of class 'spec'
#' @param x An object to test
#' @export
is.spec <- function(x) inherits(x, "spec")

splineGrad <- function(x, ...) UseMethod(".splineGrad")
.splineGrad.default <- function(dseq, dsig, plot.derivs=FALSE, ...){
  #
  # Use spline interpolation to help find an emprirical gradient
  # (reduces numerical instability)
  #
  # @dseq: the sequence (index) for @dsig, the signal
  # output is the same length as the input
  #
  require(stats, graphics)
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
  fsigderiv <<- SPLFUN(dseq, deriv=1)
  #
  # second deriv
  #   FD2 <- function(){graphics::curve(SPLFUN(x,deriv=2), from=from, to=to, n=n, add=NA)}
  #   fsigderiv2 <<- FD2()
  fsigderiv2 <<- SPLFUN(dseq, deriv=2)
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
         main="splineGrad: signal and weighted cubic-spline fit")
    lines(y ~ x, toret, col="grey", lwd=1)
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
##
na_mat <- function(x, ...) UseMethod(".na_mat")
.na_mat.default <- function(nrow, ncol=1){matrix(NA, nrow, ncol)}
##
zeros <- function(x, ...) UseMethod(".zeros")
.zeros.default <- function(nrow){matrix(rep.int(0, nrow),nrow=nrow)}
##
ones <- function(x, ...) UseMethod(".ones")
.ones.default <- function(nrow){matrix(rep.int(1, nrow),nrow=nrow)}
##
mod <- function(x, ...) UseMethod(".mod")
.mod.default <- function(x,y){
  ## modulo division
  ## R %/% requires strict consideration of order of operations whereas this is
  ## function internal and thus less prone to error, but perhaps less efficient?
  ##
  ## Args:  x	val 1
  ##		    y	val 2
  ##
  ## Returns:	modulo division of x vs y
  ##
  x1 <- trunc(trunc(x/y)*y)
  z <- trunc(x)-x1
  z
}
# end mod
###