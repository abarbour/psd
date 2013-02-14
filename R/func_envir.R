#' @title Various environment manipulation functions.
#' 
#' @description The computation of \emph{adaptive} power spectral density estimates 
#' requires bookkeeping and non-destructive manipulation of variables.  
#' The functions here are mainly convenience wrappers
#' designed to maintain variable separation from the 
#' \code{.GlobalEnv} environment so that no innocent variable is destroyed in
#' the process of iteratively computing spectra.
#' \strong{The user should generally not be using the \emph{setters} even though
#' all functions exist in the namespace.}
#' 
#' @rdname rlpSpec-environment
#' @name rlpSpec-environment
#' @docType methods
#' @section Defaults and Initialization:
#' By default, these functions all use the global pointer and name set by \code{rlp_envGlobals}.
#' One can use \code{get_rlp_env_pointer()} and \code{get_rlp_env_name()} to access them if
#' needed.
#'
#' If the environment has not yet been initialized (it should never need to be prior
#' to running \code{pspectrum} unless it was destroyed) \code{rlp_initEnv} should be used.
#' If a fresh environment is desired, and the environment already exists, 
#' \code{rlp_envClear} (which is
#' really just an alias for \code{rlp_initEnv(refresh=TRUE)}) can be used.
#'
#' @section Assigning and Retieving:
#' \code{rlp_envAssign} and \code{rlp_envGet} perform the assignments and retrieval
#' of objects in the environment.  A convenience function, \code{rlp_envAssignGet},
#' is included so that both assignment and retrieval may be performed at the same
#' time.  This ensures the assignment has succeeded, and the returned value is
#' not from the \code{.GlobalEnv} or any other environment.
#'
#' @section Getters and Setters:
#' The functions here can be classified whether the get, or set variables in the
#' environment (noted if in the global environment); some do both.  
#' Others make no modifications to the environment.
#' 
#' \subsection{Getter}{
#' \itemize{
#' \item{\code{get_adapt_history}}{}
#' \item{\code{get_rlp_env_name}}{global}
#' \item{\code{get_rlp_env_pointer}}{global}
#' \item{\code{rlp_envGet}}{}
#' \item{\code{rlp_envList}}{}
#' \item{\code{rlp_envStatus}}{both}
#' }
#' }
#'
#' \subsection{Setter}{
#' \itemize{
#' \item{\code{new_adapt_history}}{} #S
#' \item{\code{rlp_envAssign}}{}
#' \item{\code{rlp_envGlobals}}{global}
#' }
#' }
#' 
#' \subsection{Getter and Setter}{
#' \itemize{
#' \item{\code{rlp_envAssignGet}}{}
#' \item{\code{rlp_envClear}}{both}
#' \item{\code{rlp_initEnv}}{both}
#' \item{\code{update_adapt_history}}{}
#' }
#' }
#' 
#' @seealso \code{\link{rlpSpec-utilities}}, \code{\link{char2envir}}, \code{\link{pspectrum}}
#' @example inst/Examples/rdex_rlpenv.R
NULL

#' @description \code{rlp_envGlobals} sets up the environment's pointer, and name
#' strings in the \code{.GlobalEnv}.
#' @rdname rlpSpec-environment
#' @name rlp_envGlobals
#' @param envpoint  character; the pointer to the environment
#' @param envname  character; the name of the environment
#' @export
rlp_envGlobals <- function(envpoint=".rlpenv", envname=".rlpSpecEnv") assign(envpoint, envname, envir=.GlobalEnv)
#' @description \code{get_rlp_env_pointer} is a convenience wrapper to get the environment pointer.
#' @rdname rlpSpec-environment
#' @name get_rlp_env_pointer
#' @export
get_rlp_env_pointer <- function() as.name(formals(rlp_envGlobals)$envpoint)
#' @description \code{get_rlp_env_name} is a convenience wrapper to get the environment name.
#' @rdname rlpSpec-environment
#' @name get_rlp_env_name
#' @export
get_rlp_env_name <- function() eval(get_rlp_env_pointer())

#' @description \code{rlp_initEnv} initializes the environment with
#' an option to clear the contents (if the environment already exists).
#' @note \code{rlp_initEnv} will not re-initialize the enviroment, unless told to 
#' do so with \code{refresh=TRUE}.
#' @rdname rlpSpec-environment
#' @name rlp_initEnv
#' @param refresh logical; should the contents of the environment be trashed?
#' @param verbose logical; should messages be given?
#' @return \code{rlp_initEnv} returns (invisibly) the result of \code{rlp_envStatus()}.
#' @seealso \code{new.env}, \code{baseenv}
#' @export
rlp_initEnv <- function(refresh=FALSE, verbose=TRUE, ...) {
  # env params
  new_env <- FALSE
  envir <- get_rlp_env_pointer()
  envname <- get_rlp_env_name()
  #
  if (!exists(envname)){
    rlpSpec:::rlp_envGlobals()
    new_env <- TRUE
  }
  # initialize the psd calculation environment
  msg <- "no action"
  if( new_env | refresh ){
    assign(x=envname, 
           value=new.env(parent=globalenv(), ...),  #baseenv()
           envir=.GlobalEnv)
    msg <- "initialized"
    if (refresh & !new_env){ msg <- "refreshed"}
    if (verbose) {
      message(sprintf("\tenvironment  ** %s **  %s", envname, msg))
    }
  } else if (!refresh) {
    if (verbose) message(sprintf("\t** %s ** is already initialized: try 'refresh=TRUE' to clear", envname))
  }
  rlpSpec:::rlp_envAssign("init", sprintf("%s at %s", msg, Sys.time()))
  return(invisible(rlpSpec:::rlp_envStatus()))
}

#' @description \code{rlp_envClear} clears the contents of the environment.
#' @note \code{rlp_envClear} does \emph{not} remove the environment--simply the assignments within it.
#' @rdname rlpSpec-environment
#' @name rlp_envClear
rlp_envClear <- function(...) rlpSpec:::rlp_initEnv(refresh=TRUE, ...)

#' @description \code{rlp_envStatus} returns a list of some information regarding
#' the status of the environment.
#' @rdname rlpSpec-environment
#' @name rlp_envStatus
#' @export
rlp_envStatus <- function(){
  envir <- get_rlp_env_pointer()
  envname <- get_rlp_env_name()
  return(list(env_name=envname, 
              env_pointer=envir,
              env_obviously_exists=exists(envname), 
              env_is_env=is.environment(char2envir(envname)), 
              listing=rlpSpec:::rlp_envList(),
              env_init=rlpSpec:::rlp_envGet("init"),
              env_status_stamp=Sys.time() ))
  
}

#' @description \code{rlp_envList} returns a listing of the assignments.
#' @rdname rlpSpec-environment
#' @name rlp_envList
#' @export
rlp_envList <- function(){
  ## return listing of envir::variable
  ENV <- char2envir(get_rlp_env_name())
  ls(envir=ENV, all.names=TRUE)
}

#' @description \code{rlp_envGet} returns a the value of \code{variable}.
#' @rdname rlpSpec-environment
#' @name rlp_envGet
#' @param variable character; the name of the variable to get or assign
#' @return the object represented by \code{variable} in the \code{rlpSpec} environment.
#' @export
rlp_envGet <- function(variable){
  ## return contents on envir::variable
  ENV <- char2envir(get_rlp_env_name())
  if (!exists(variable, envir=ENV)){
    warning(sprintf("Variable  '%s'  not found! Use rlp_envList()", variable))
    return(NULL)
  } else {
    return(get(variable, envir=ENV))
  }
}

#' @description \code{rlp_envAssign} assigns \code{value} to \code{variable}, but does not return it.
#' @rdname rlpSpec-environment
#' @name rlp_envAssign
#' @param value character; the name of the variable to assign
#' @export
rlp_envAssign <- function(variable, value){
  ## set contents of envir::variable to value
  ENV <- char2envir(get_rlp_env_name())
  assign(variable, value, envir=ENV)
}

#' @description \code{rlp_envAssignGet} both assigns and returns a value.
#' @rdname rlpSpec-environment
#' @name rlp_envAssignGet
# placing these at the end for orderliness.
#' @param ... For \code{rlp_envClear}: arguments passed to \code{rlp_initEnv}. For \code{rlp_initEnv}: arguments passed to \code{new.env}
#' @export
rlp_envAssignGet <- function(variable, value){
  ## set contents of envir::variable to value
  rlpSpec:::rlp_envAssign(variable, value)
  rlpSpec:::rlp_envGet(variable)
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
#' @export
new_adapt_history <- function(adapt_stages){
  stopifnot(length(adapt_stages)==1)
  histlist <- vector("list", 3) # freq, list-tap, list-psd
  names(histlist) <- c("freq", "stg_kopt", "stg_psd")
  num_pos <- 1 + adapt_stages # pilot + adapts
  histlist[[2]] <- histlist[[3]] <- vector("list", adapt_stages+1)
  rlpSpec:::rlp_envAssignGet("histlist", histlist)
}
#' @export
#' @rdname rlpSpec-environment
get_adapt_history <- function() rlpSpec:::rlp_envGet("histlist")

#' @description \code{update_adapt_history} Updates the adaptive estimation history list.
#' @rdname rlpSpec-environment
#' @param stage scalar; the current stage of the adaptive estimation procedure
#' @param ntap vector; the tapers
#' @param psd vector; the power spectral densities
#' @param freq vector; the frequencies
#' @export
update_adapt_history <- function(stage, ntap, psd, freq=NULL){
  histlist <- get_adapt_history()
  # stage==0 <--> index==1
  stg_ind <- stage+1
  nulfrq <- is.null(freq)
  if (!nulfrq) histlist$freq <- freq
  histlist$stg_kopt[[stg_ind]] <- ntap
  histlist$stg_psd[[stg_ind]] <- psd
  if (is.null(histlist$freq) & stage>0) warning("freqs absent despite non-pilot stage update")
  rlpSpec:::rlp_envAssignGet("histlist",histlist)
}
