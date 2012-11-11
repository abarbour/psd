###
###  Default methods for supplementary functions 
###    (to make porting from Matlab to R tractable and supplement)
###
# TODO(abarbour):
##
.rlp_initEnv.default <- function(refresh=FALSE, verbose=TRUE, ...){
  # initialize the psd environment
  env <- "psdenv"
  if(!exists(env) | refresh){
    psdenv <<- new.env(parent=baseenv(), ...)
    if (verbose) message(sprintf("\t>>>> ** %s ** environment initialized",env))
  } else if (!refresh) {
    if (verbose) message(sprintf("\t>>>> The ** %s ** environment is already initialized",env))
  }
}
.rlp_envList.default <- function(envir=psdenv){
  ## return listing of envir::variable
  ls(envir=envir)
}
.rlp_envGet.default <- function(variable, envir=psdenv){
  ## return contents on envir::variable
  get(variable, envir=envir)
}
.rlp_envAssign.default <- function(variable, value, envir=psdenv){
  ## set contents of envir::variable to value
  assign(variable, value, envir=envir)
}
.rlp_envAssignGet.default <- function(variable, value, envir=psdenv){
  ## set contents of envir::variable to value
  envAssign(variable, value, envir=envir)
  envGet(variable, envir=envir)
}
##
.nas.default <- function(nrow, ncol=1){matrix(NA, nrow, ncol)}
##
.colvec.default <- function(nrow, val){matrix(val, nrow=nrow, ncol=1)}
.rowvec.default <- function(ncol, val){matrix(val, nrow=1, ncol=ncol)}
##
.zeros.default <- function(nrow){colvec(nrow, 0)}
##
.ones.default <- function(nrow){colvec(nrow, 1)}
##
.mod.default <- function(x,y){
  ## modulo division
  ## R %/% requires strict consideration of order of operations whereas this is
  ## function internal and thus less prone to error, but perhaps less efficient?
  ##
  ## Args:	x	val 1
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