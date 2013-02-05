###
###  Default methods for supplementary functions 
###    (to make porting from Matlab to R tractable and supplement)
###
# TODO(abarbour):
##


envClear <- function(...) rlp_initEnv(refresh=TRUE, ...)

envAssign <- function(x, ...) UseMethod(".rlp_envAssign")
envGet <- function(x, ...) UseMethod(".rlp_envGet")
envAssignGet <- function(x, ...) UseMethod(".rlp_envAssignGet")

rlp_envStatus <- function(...) UseMethod(".rlp_envStatus")
.rlp_envStatus.default <- function(envir=.rlpenv){
  #rlp_initEnv(refresh=FALSE, verbose=FALSE)
  return(list(env.name=envir, 
              obvious.exists=exists(envir), 
              listing=rlp_envList(envir),
              init.stamp=Sys.time()))
  
}

rlp_initEnv <- function(x, ...) UseMethod(".rlp_initEnv")
.rlp_initEnv.default <- function(envir=.rlpenv, refresh=FALSE, verbose=TRUE, ...){
  # initialize the psd calculation environment
  if( !exists(envir) | refresh ){
    assign(x=envir, 
           value=new.env(parent=baseenv(), ...), 
           envir=.GlobalEnv)
    if (verbose) {
      msg <- "initialized"
      if (refresh){ msg <- "refreshed"}
      message(sprintf("\tenvironment  ** %s **  %s", envir, msg))
    }
  } else if (!refresh) {
    if (verbose) message(sprintf("\t** %s ** is already initialized: try 'refresh=TRUE' to clear", envir))
  }
  rlp_envStatus(envir)
}

rlp_envList <- function(x, ...) UseMethod(".rlp_envList")
.rlp_envList.default <- function(envir=.rlpenv){
  ## return listing of envir::variable
  print(envir)
  if (!is.environment(envir)) envir <- as.name(envir)
  print(envir)
  ls(envir=envir)
}

.rlp_envGet.default <- function(variable, envir=.rlpenv){
  ## return contents on envir::variable
  if (!is.environment(envir)) envir <- as.symbol(envir)
  get(variable, envir=envir)
}
.rlp_envAssign.default <- function(variable, value, envir=.rlpenv){
  ## set contents of envir::variable to value
  if (!is.environment(envir)) envir <- as.symbol(envir)
  assign(variable, value, envir=envir)
}
.rlp_envAssignGet.default <- function(variable, value, envir=.rlpenv){
  ## set contents of envir::variable to value
  envAssign(variable, value, envir=envir)
  envGet(variable, envir=envir)
}
##
.reshape_vector.default <- function(x, vec.shape=c("horizontal","vertical")){
  x <- as.vector(x)
  vec.shape <- match.arg(vec.shape)
  nrow <- switch(vec.shape, "horizontal"=1, "vertical"=length(x))
  return(matrix(x, nrow=nrow))
}
##
.nas.default <- function(nrow, ncol=1){matrix(NA, nrow, ncol)}
##
.zeros.default <- function(nrow){matrix(rep.int(0, nrow),nrow=nrow)}
##
.ones.default <- function(nrow){matrix(rep.int(1, nrow),nrow=nrow)}
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
