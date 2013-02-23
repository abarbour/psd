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
#' @rdname psd-environment
#' @name psd-environment
#' @docType methods
#' @section Defaults and Initialization:
#' By default, these functions all use the global pointer and name set by \code{psd_envGlobals()}.
#' One can use \code{get_psd_env_pointer()} and \code{get_psd_env_name()} to access them if
#' needed.
#'
#' If the environment has not yet been initialized (it should never need to be prior
#' to running \code{pspectrum} unless it was destroyed) \code{psd_initEnv} should be used.
#' If a fresh environment is desired, and the environment already exists, 
#' \code{psd_envClear} (which is
#' really just an alias for \code{psd_initEnv(refresh=TRUE)}) can be used.
#'
#' @section Assigning and Retieving:
#' \code{psd_envAssign} and \code{psd_envGet} perform the assignments and retrieval
#' of objects in the environment.  A convenience function, \code{psd_envAssignGet},
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
#' \item{\code{get_psd_env_name}}{ (global)}
#' \item{\code{get_psd_env_pointer}}{ (global)}
#' \item{\code{psd_envGet}}{}
#' \item{\code{psd_envList}}{}
#' \item{\code{psd_envStatus}}{both}
#' }
#' }
#'
#' \subsection{Setter}{
#' \itemize{
#' \item{\code{new_adapt_history}}{}
#' \item{\code{psd_envAssign}}{}
#' \item{\code{psd_envGlobals}}{ (global)}
#' }
#' }
#' 
#' \subsection{Getter and Setter}{
#' \itemize{
#' \item{\code{psd_envAssignGet}}{}
#' \item{\code{psd_envClear}}{ (both)}
#' \item{\code{psd_initEnv}}{ (both)}
#' \item{\code{update_adapt_history}}{}
#' }
#' }
#' 
#' @seealso \code{\link{psd-utilities}}, \code{\link{char2envir}}, \code{\link{pspectrum}}
#' @example inst/Examples/rdex_psdenv.R
NULL

#' @description \code{psd_envGlobals} sets up the environment's pointer, and name
#' strings in the \code{.GlobalEnv} environment.  See the Usage section for
#' the names of the default pointer and environment.
#' @rdname psd-environment
#' @name psd_envGlobals
#' @param envpoint  character; the pointer to the environment
#' @param envname  character; the name of the environment
#' @export
psd_envGlobals <- function(envpoint=".psdenv", envname=".PsdSpecEnv") assign(envpoint, envname, envir=.GlobalEnv)
#' @description \code{get_psd_env_pointer} is a convenience wrapper to get the environment pointer.
#' @rdname psd-environment
#' @name get_psd_env_pointer
#' @export
get_psd_env_pointer <- function() as.name(formals(psd_envGlobals)$envpoint)
#' @description \code{get_psd_env_name} is a convenience wrapper to get the environment name.
#' @rdname psd-environment
#' @name get_psd_env_name
#' @export
get_psd_env_name <- function() eval(get_psd_env_pointer())

#' @description \code{psd_initEnv} initializes the environment with
#' an option to clear the contents (if the environment already exists).
#' @note \code{psd_initEnv} will not re-initialize the enviroment, unless told to 
#' do so with \code{refresh=TRUE}.
#' @rdname psd-environment
#' @name psd_initEnv
#' @param refresh logical; should the contents of the environment be trashed?
#' @param verbose logical; should messages be given?
#' @return \code{psd_initEnv} returns (invisibly) the result of \code{psd_envStatus()}.
#' @seealso \code{new.env}, \code{baseenv}
#' @export
psd_initEnv <- function(refresh=FALSE, verbose=TRUE, ...) {
  # env params
  new_env <- FALSE
  envir <- get_psd_env_pointer()
  envname <- get_psd_env_name()
  #
  if (!exists(envname)){
    psd:::psd_envGlobals()
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
  psd:::psd_envAssign("init", sprintf("%s at %s", msg, Sys.time()))
  return(invisible(psd:::psd_envStatus()))
}

#' @description \code{psd_envClear} clears the contents of the environment.
#' @note \code{psd_envClear} does \emph{not} remove the environment--simply the assignments within it.
#' @rdname psd-environment
#' @name psd_envClear
psd_envClear <- function(...) psd:::psd_initEnv(refresh=TRUE, ...)

#' @description \code{psd_envStatus} returns a list of some information regarding
#' the status of the environment.
#' @rdname psd-environment
#' @name psd_envStatus
#' @export
psd_envStatus <- function(){
  envir <- get_psd_env_pointer()
  envname <- get_psd_env_name()
  return(list(env_name=envname, 
              env_pointer=envir,
              env_obviously_exists=exists(envname), 
              env_is_env=is.environment(char2envir(envname)), 
              listing=psd:::psd_envList(),
              env_init=psd:::psd_envGet("init"),
              env_status_stamp=Sys.time() ))
  
}

#' @description \code{psd_envList} returns a listing of the assignments.
#' @rdname psd-environment
#' @name psd_envList
#' @export
psd_envList <- function(){
  ## return listing of envir::variable
  ENV <- char2envir(get_psd_env_name())
  ls(envir=ENV, all.names=TRUE)
}

#' @description \code{psd_envGet} returns a the value of \code{variable}.
#' @rdname psd-environment
#' @name psd_envGet
#' @param variable character; the name of the variable to get or assign
#' @return the object represented by \code{variable} in the \code{psd} environment.
#' @export
psd_envGet <- function(variable){
  ## return contents on envir::variable
  ENV <- char2envir(get_psd_env_name())
  if (!exists(variable, envir=ENV)){
    warning(sprintf("Variable  '%s'  not found! Use psd_envList()", variable))
    return(NULL)
  } else {
    return(get(variable, envir=ENV))
  }
}

#' @description \code{psd_envAssign} assigns \code{value} to \code{variable}, but does not return it.
#' @rdname psd-environment
#' @name psd_envAssign
#' @param value character; the name of the variable to assign
#' @export
psd_envAssign <- function(variable, value){
  ## set contents of envir::variable to value
  ENV <- char2envir(get_psd_env_name())
  assign(variable, value, envir=ENV)
}

#' @description \code{psd_envAssignGet} both assigns and returns a value.
#' @rdname psd-environment
#' @name psd_envAssignGet
# placing these at the end for orderliness.
#' @param ... For \code{psd_envClear}: arguments passed to \code{psd_initEnv}. For \code{psd_initEnv}: arguments passed to \code{new.env}
#' @export
psd_envAssignGet <- function(variable, value){
  ## set contents of envir::variable to value
  psd:::psd_envAssign(variable, value)
  psd:::psd_envGet(variable)
}

#' @description \code{new_adapt_history} initializes a nested-list object to store the 
#' data from each iteration.
#' @rdname psd-environment
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
  psd:::psd_envAssignGet("histlist", histlist)
}
#' @export
#' @rdname psd-environment
get_adapt_history <- function() psd:::psd_envGet("histlist")

#' @description \code{update_adapt_history} Updates the adaptive estimation history list.
#' @rdname psd-environment
#' @param stage scalar; the current stage of the adaptive estimation procedure
#' @param ntap vector; the tapers
#' @param PSD vector; the power spectral densities
#' @param freq vector; the frequencies
#' @export
update_adapt_history <- function(stage, ntap, PSD, freq=NULL){
  histlist <- get_adapt_history()
  # stage==0 <--> index==1
  stg_ind <- stage+1
  nulfrq <- is.null(freq)
  if (!nulfrq) histlist$freq <- freq
  histlist$stg_kopt[[stg_ind]] <- ntap
  histlist$stg_psd[[stg_ind]] <- PSD
  if (is.null(histlist$freq) & stage>0) warning("freqs absent despite non-pilot stage update")
  psd:::psd_envAssignGet("histlist",histlist)
}
