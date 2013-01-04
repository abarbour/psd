###
###  Default methods for supplementary functions 
###    (to make porting from Matlab to R tractable and supplement)
###
# TODO(abarbour):
##

rlp_initEnv <- function(x, ...) UseMethod(".rlp_initEnv")
.rlp_initEnv.default <- function(envir=.rlpenv, refresh=FALSE, verbose=TRUE, ...){
  # initialize the psd calculation environment
  if (!exists(envir)){ new_env <- TRUE } else {new_env <- FALSE}
  if( new_env | refresh ){
    assign(x=envir, 
           value=new.env(parent=baseenv(), ...), 
           envir=.GlobalEnv)
    msg <- "initialized"
    if (verbose) {
      if (refresh & !new_env){ msg <- "refreshed"}
      message(sprintf("\tenvironment  ** %s **  %s", envir, msg))
    }
  } else if (!refresh) {
    if (verbose) message(sprintf("\t** %s ** is already initialized: try 'refresh=TRUE' to clear", envir))
  }
  rlp_envStatus(envir)
}

rlp_envClear <- function(...) rlp_initEnv(refresh=TRUE, ...)

rlp_envStatus <- function(...) UseMethod(".rlp_envStatus")
.rlp_envStatus.default <- function(envir=.rlpenv){
  #rlp_initEnv(refresh=FALSE, verbose=FALSE)
  if (is.environment(envir)) envir <- substitute(deparse(envir))
  return(list(env.name=envir, 
              obvious.exists=exists(envir), 
              listing=rlp_envList(envir),
              init.stamp=Sys.time()))
  
}

rlp_envList <- function(x, ...) UseMethod(".rlp_envList")
.rlp_envList.default <- function(envir=.rlpenv){
  ## return listing of envir::variable
  if (is.character(envir)) envir <- char2envir(envir)
  ls(envir=envir)
}


rlp_envGet <- function(x, ...) UseMethod(".rlp_envGet")
.rlp_envGet.default <- function(variable, envir=.rlpenv){
  ## return contents on envir::variable
  if (is.character(envir)) envir <- char2envir(envir)
  get(variable, envir=envir)
}

rlp_envAssign <- function(x, ...) UseMethod(".rlp_envAssign")
.rlp_envAssign.default <- function(variable, value, envir=.rlpenv){
  ## set contents of envir::variable to value
  if (is.character(envir)) envir <- char2envir(envir)
  assign(variable, value, envir=envir)
}

rlp_envAssignGet <- function(x, ...) UseMethod(".rlp_envAssignGet")
.rlp_envAssignGet.default <- function(variable, value, envir=.rlpenv){
  ## set contents of envir::variable to value
  envAssign(variable, value, envir=envir)
  envGet(variable, envir=envir)
}

