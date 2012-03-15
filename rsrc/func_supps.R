###
###  Default methods for supplementary functions 
###    (to make porting from Matlab to R tractable and supplement)
###
# TODO(abarbour):
##
initEnv.default <- function(refresh=FALSE, ...){
  # initialize the psd environment
  env <- "psdenv"
  if(!exists(env) | refresh){
    psdenv <<- new.env(parent=baseenv(), ...)
    cat(sprintf("\t>>>> ** %s ** environment initialized\n",env))
  } else if (!refresh) {
    cat(sprintf("\t>>>> The ** %s ** environment is already initialized.\n",env))
  }
}
envGet.default <- function(variable, value, envir=psdenv){
  ## return contents on envir::variable
  get(variable, envir=envir)
}
envAssign.default <- function(variable, value, envir=psdenv){
  ## set contents of envir::variable to value
  assign(variable, value, envir=envir)
}
##
nas.default <- function(nrow, ncol=1){matrix(NA, nrow, ncol)}
##
zeros.default <- function(nrow, ncol=1){matrix(0, nrow, ncol)}
##
ones.default <- function(nrow, ncol=1){matrix(1, nrow, ncol)}
##
mod.default <- function(x,y){
  ## modulo division
  ## R %/% requires strict consideration of order of operations whereas this is
  ## function internal and thus less prone to error, but perhaps less efficient?
  ##
  ## Args:	x	val 1
  ##		y	val 2
  ##
  ## Returns:	modulo division of x vs y
  ##
  x1 <- trunc(trunc(x/y)*y)
  z <- trunc(x)-x1
  z
}
# end mod
###