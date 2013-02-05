#' @title Various environment manipulation functions.
#' 
#' @description The computation of \emph{adaptive} power spectral density estimates 
#' requires bookkeeping and non-destructive manipulation of variables.  
#' The functions here are mainly convenience wrappers
#' designed to maintain variable separation from the 
#' \code{.GlobalEnv} environment so that no innocent variable is destroyed in
#' the process of iteratively computing spectra.
#' \strong{The user should generally not be using the \emph{setters}; hence, only the
#' \emph{getters} exist in the namespace.}
#' 
#' @rdname rlpSpec-environment
#' @name rlpSpec-environment
#' @docType methods
#' @param envir  character string for the environment used.  See section \strong{Defaults}
#' @section Defaults and Initialization:
#' By default, these functions all use \code{envir=.rlpenv} to set the enviroment-name
#' string; \code{.rlpenv} is set when attaching the package. 
#' 
#' If the environment has not yet been initialized (it should never need to be prior
#' to running \code{pspectrum}) \code{rlp_initEnv} should be used.
#' If a fresh environment is desired, and the environment already exists, 
#' \code{rlp_envClear} (which is
#' really just an alias for \code{rlp_initEnv(refresh=TRUE)}) can be used.
#'
#' One could set \code{.rlpenv} 
#' to another string, if a different environment is desired.
#'
#' @section Assigning and Retieving:
#' \code{rlp_envAssign} and \code{rlp_envGet} perform the assignments and retrieval
#' of objects in the environment.  A convenience function, \code{rlp_envAssignGet},
#' is included so that both assignment and retrieval may be performed at the same
#' time.  This ensures the assignment has succeeded, and the returned value is
#' not from the \code{.GlobalEnv} or any other environment.
#'
#' @seealso \code{\link{rlpSpec-utilities}}, \code{\link{char2envir}}, \code{\link{pspectrum}}
#' @example inst/Examples/rdex_rlpenv.R
NULL

# # @description \code{envParamSet} sets up pointers for the \code{.rlpenv} environment,
# # in the \code{.GlobalEnv}.
# # @rdname rlpSpec-environment
# # @name rlp_envParamSet
# # @export
# rlp_envParamSet <- function() assign(".rlpenv", ".rlpSpecEnv", envir=.GlobalEnv)

#' @description \code{rlp_initEnv} initializes the \code{.rlpenv} environment with
#' an option to clear the contents (if the environment already exists).
#' @note \code{rlp_initEnv} will not re-initialize the enviroment, unless told to 
#' do so with \code{refresh=TRUE}.
#' @rdname rlpSpec-environment
#' @name rlp_initEnv
#' @param refresh logical; should the contents of the environment be trashed?
#' @param verbose logical; should messages be given?
#' @return \code{rlp_initEnv} returns (invisibly) the result of \code{rlp_envStatus()}.
#' @seealso \code{new.env}, \code{baseenv}
rlp_initEnv <- function(envir=.rlpenv, refresh=FALSE, verbose=TRUE, ...) {
  # initialize the psd calculation environment
  if (exists(envir)){ new_env <- FALSE } else {new_env <- TRUE}
  msg <- "no action"
  if( new_env | refresh ){
    assign(x=envir, 
           value=new.env(parent=globalenv(), ...),  #baseenv()
           envir=.GlobalEnv)
    msg <- "initialized"
    if (refresh & !new_env){ msg <- "refreshed"}
    if (verbose) {
      message(sprintf("\tenvironment  ** %s **  %s", envir, msg))
    }
  } else if (!refresh) {
    if (verbose) message(sprintf("\t** %s ** is already initialized: try 'refresh=TRUE' to clear", envir))
  }
  rlp_envAssign("init", sprintf("%s at %s", msg, Sys.time()))
  return(invisible(rlp_envStatus(envir)))
}

#' @description \code{rlp_envClear} clears the contents of the environment.
#' @note \code{rlp_envClear} does \emph{not} remove the environment--simply the assignments within it.
#' @rdname rlpSpec-environment
#' @name rlp_envClear
rlp_envClear <- function(...) rlp_initEnv(refresh=TRUE, ...)

#' @description \code{rlp_envStatus} returns a list of some information regarding
#' the status of the environment.
#' @rdname rlpSpec-environment
#' @name rlp_envStatus
#' @export
rlp_envStatus <- function(envir=.rlpenv){
  #rlp_initEnv(refresh=FALSE, verbose=FALSE)
  if (is.environment(envir)) envir <- substitute(deparse(envir))
  #is_init <- ifelse(exists(envir))
  return(list(env_name=envir, 
              obviously_exists=exists(envir), 
              listing=rlp_envList(envir),
              env_init=rlp_envGet("init"),
              env_status_stamp=Sys.time() ))
  
}

#' @description \code{rlp_envList} returns a listing of the assignments.
#' @rdname rlpSpec-environment
#' @name rlp_envList
#' @export
rlp_envList <- function(envir=.rlpenv){
  ## return listing of envir::variable
  stopifnot(is.character(envir))
  ENV <- char2envir(envir)
  ls(envir=ENV, all.names=TRUE)
}

#' @description \code{rlp_envGet} returns a the value of \code{variable}.
#' @rdname rlpSpec-environment
#' @name rlp_envGet
#' @param variable character; the name of the variable to get or assign
#' @return the object represented by \code{variable} in the \code{rlpSpec} environment.
#' @export
rlp_envGet <- function(variable, envir=.rlpenv){
  ## return contents on envir::variable
  stopifnot(is.character(variable) & is.character(envir))
  ENV <- char2envir(envir)
  if (!exists(variable, envir=ENV)){
    return(NULL)
  } else {
    return(get(variable, envir=ENV))
  }
}

#' @description \code{rlp_envAssign} assigns \code{value} to \code{variable}, but does not return it.
#' @rdname rlpSpec-environment
#' @name rlp_envAssign
#' @param value character; the name of the variable to assign
rlp_envAssign <- function(variable, value, envir=.rlpenv){
  ## set contents of envir::variable to value
  stopifnot(is.character(variable) & is.character(envir))
  ENV <- char2envir(envir)
  assign(variable, value, envir=ENV)
}

#' @description \code{rlp_envAssignGet} both assigns and returns a value.
#' @rdname rlpSpec-environment
#' @name rlp_envAssignGet
# placing these at the end for orderliness.
#' @param ... For \code{rlp_envClear}: arguments passed to \code{rlp_initEnv}. For \code{rlp_initEnv}: arguments passed to \code{new.env}
rlp_envAssignGet <- function(variable, value, envir=.rlpenv){
  ## set contents of envir::variable to value
  rlp_envAssign(variable, value, envir=envir)
  rlp_envGet(variable, envir=envir)
}

#' @description \code{new_adapt_history} initializes a nested-list object to store the 
#' data from each iteration.
#' @rdname rlpSpec-environment
#'
#' @section Adaptive History:
#' The list object for historical adapt-data may be accessed with \code{\link{get_adapt_history}}.
#' The top names of the returned list are
#' \describe{
#' \item{\code{stg_kopt}}{Sequential taper vectors.}
#' \item{\code{stg_psd}}{Sequential power spectral density vectors.}
#' \item{\code{freq}}{The frequencies for each set of \code{stg_kopt} and \code{stg_psd}.}
#' }
#' @param adapt_stages scalar; The number of adaptive iterations to save (excluding pilot spectrum).
new_adapt_history <- function(adapt_stages){
  stopifnot(length(adapt_stages)==1)
  histlist <- vector("list", 3) # freq, list-tap, list-psd
  names(histlist) <- c("freq", "stg_kopt", "stg_psd")
  num_pos <- 1 + adapt_stages # pilot + adapts
  histlist[[2]] <- histlist[[3]] <- vector("list", adapt_stages+1)
  rlp_envAssignGet("histlist", histlist)
}

#' @export
#' @rdname rlpSpec-environment
get_adapt_history <- function() rlp_envGet("histlist")

#' @description \code{update_adapt_history} Updates the adaptive estimation history list.
#' @rdname rlpSpec-environment
#' @param stage scalar; the current stage of the adaptive estimation procedure
#' @param ntap vector; the tapers
#' @param psd vector; the power spectral densities
#' @param freq vector; the frequencies
update_adapt_history <- function(stage, ntap, psd, freq=NULL){
  histlist <- get_adapt_history()
  # stage==0 <--> index==1
  stg_ind <- stage+1
  nulfrq <- is.null(freq)
  if (!nulfrq) histlist$freq <- freq
  histlist$stg_kopt[[stg_ind]] <- ntap
  histlist$stg_psd[[stg_ind]] <- psd
  if (is.null(histlist$freq) & stage>0) warning("freqs absent despite non-pilot stage update")
  rlp_envAssignGet("histlist",histlist)
}
