###
###  Default methods for supplementary functions 
###    (to make porting from Matlab to R tractable and supplement)
###
# TODO(abarbour):
##
.rlp_envStatus.default <- function(){
  initEnv(refresh=FALSE, verbose=FALSE)
}
.rlp_initEnv.default <- function(env="psdenv", refresh=FALSE, verbose=TRUE, ...){
  # initialize the psd environment
  if( !exists(env) | refresh ){
    assign(x=env, value=new.env(parent=baseenv(), ...), envir=.GlobalEnv)
    if (verbose) {
      msg <- "initialized"
      if (refresh){ msg <- "refreshed"}
      message(sprintf("\tenvironment  ** %s **  %s", env, msg))
    }
  } else if (!refresh) {
    if (verbose) message(sprintf("\t** %s ** is already initialized: try 'refresh=TRUE' to clear", env))
  }
  return(list(env.name=env, 
              obvious.exists=exists(env), 
              listing=envList(),
              init.stamp=Sys.time()))
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
