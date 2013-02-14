\docType{methods}
\name{rlpSpec-environment}
\alias{get_adapt_history}
\alias{get_rlp_env_name}
\alias{get_rlp_env_pointer}
\alias{new_adapt_history}
\alias{rlpSpec-environment}
\alias{rlp_envAssign}
\alias{rlp_envAssignGet}
\alias{rlp_envClear}
\alias{rlp_envGet}
\alias{rlp_envGlobals}
\alias{rlp_envList}
\alias{rlp_envStatus}
\alias{rlp_initEnv}
\alias{update_adapt_history}
\title{Various environment manipulation functions.}
\usage{
  rlp_envGlobals(envpoint = ".rlpenv",
    envname = ".rlpSpecEnv")

  get_rlp_env_pointer()

  get_rlp_env_name()

  rlp_initEnv(refresh = FALSE, verbose = TRUE, ...)

  rlp_envClear(...)

  rlp_envStatus()

  rlp_envList()

  rlp_envGet(variable)

  rlp_envAssign(variable, value)

  rlp_envAssignGet(variable, value)

  new_adapt_history(adapt_stages)

  get_adapt_history()

  update_adapt_history(stage, ntap, psd, freq = NULL)
}
\arguments{
  \item{envpoint}{character; the pointer to the
  environment}

  \item{envname}{character; the name of the environment}

  \item{refresh}{logical; should the contents of the
  environment be trashed?}

  \item{verbose}{logical; should messages be given?}

  \item{variable}{character; the name of the variable to
  get or assign}

  \item{value}{character; the name of the variable to
  assign}

  \item{...}{For \code{rlp_envClear}: arguments passed to
  \code{rlp_initEnv}. For \code{rlp_initEnv}: arguments
  passed to \code{new.env}}

  \item{adapt_stages}{scalar; The number of adaptive
  iterations to save (excluding pilot spectrum).}

  \item{stage}{scalar; the current stage of the adaptive
  estimation procedure}

  \item{ntap}{vector; the tapers}

  \item{psd}{vector; the power spectral densities}

  \item{freq}{vector; the frequencies}
}
\value{
  \code{rlp_initEnv} returns (invisibly) the result of
  \code{rlp_envStatus()}.

  the object represented by \code{variable} in the
  \code{rlpSpec} environment.
}
\description{
  The computation of \emph{adaptive} power spectral density
  estimates requires bookkeeping and non-destructive
  manipulation of variables. The functions here are mainly
  convenience wrappers designed to maintain variable
  separation from the \code{.GlobalEnv} environment so that
  no innocent variable is destroyed in the process of
  iteratively computing spectra. \strong{The user should
  generally not be using the \emph{setters} even though all
  functions exist in the namespace.}

  \code{rlp_envGlobals} sets up the environment's pointer,
  and name strings in the \code{.GlobalEnv}.

  \code{get_rlp_env_pointer} is a convenience wrapper to
  get the environment pointer.

  \code{get_rlp_env_name} is a convenience wrapper to get
  the environment name.

  \code{rlp_initEnv} initializes the environment with an
  option to clear the contents (if the environment already
  exists).

  \code{rlp_envClear} clears the contents of the
  environment.

  \code{rlp_envStatus} returns a list of some information
  regarding the status of the environment.

  \code{rlp_envList} returns a listing of the assignments.

  \code{rlp_envGet} returns a the value of \code{variable}.

  \code{rlp_envAssign} assigns \code{value} to
  \code{variable}, but does not return it.

  \code{rlp_envAssignGet} both assigns and returns a value.

  \code{new_adapt_history} initializes a nested-list object
  to store the data from each iteration.

  \code{update_adapt_history} Updates the adaptive
  estimation history list.
}
\note{
  \code{rlp_initEnv} will not re-initialize the enviroment,
  unless told to do so with \code{refresh=TRUE}.

  \code{rlp_envClear} does \emph{not} remove the
  environment--simply the assignments within it.
}
\section{Defaults and Initialization}{
  By default, these functions all use the global pointer
  and name set by \code{rlp_envGlobals}. One can use
  \code{get_rlp_env_pointer()} and
  \code{get_rlp_env_name()} to access them if needed.

  If the environment has not yet been initialized (it
  should never need to be prior to running \code{pspectrum}
  unless it was destroyed) \code{rlp_initEnv} should be
  used. If a fresh environment is desired, and the
  environment already exists, \code{rlp_envClear} (which is
  really just an alias for
  \code{rlp_initEnv(refresh=TRUE)}) can be used.
}

\section{Assigning and Retieving}{
  \code{rlp_envAssign} and \code{rlp_envGet} perform the
  assignments and retrieval of objects in the environment.
  A convenience function, \code{rlp_envAssignGet}, is
  included so that both assignment and retrieval may be
  performed at the same time.  This ensures the assignment
  has succeeded, and the returned value is not from the
  \code{.GlobalEnv} or any other environment.
}

\section{Getters and Setters}{
  The functions here can be classified whether the get, or
  set variables in the environment (noted if in the global
  environment); some do both. Others make no modifications
  to the environment.

  \subsection{Getter}{ \itemize{
  \item{\code{get_adapt_history}}{}
  \item{\code{get_rlp_env_name}}{global}
  \item{\code{get_rlp_env_pointer}}{global}
  \item{\code{rlp_envGet}}{} \item{\code{rlp_envList}}{}
  \item{\code{rlp_envStatus}}{both} } }

  \subsection{Setter}{ \itemize{
  \item{\code{new_adapt_history}}{} #S
  \item{\code{rlp_envAssign}}{}
  \item{\code{rlp_envGlobals}}{global} } }

  \subsection{Getter and Setter}{ \itemize{
  \item{\code{rlp_envAssignGet}}{}
  \item{\code{rlp_envClear}}{both}
  \item{\code{rlp_initEnv}}{both}
  \item{\code{update_adapt_history}}{} } }
}

\section{Adaptive History}{
  The list object for historical adapt-data may be accessed
  with \code{\link{get_adapt_history}}. The top names of
  the returned list are \describe{
  \item{\code{stg_kopt}}{Sequential taper vectors.}
  \item{\code{stg_psd}}{Sequential power spectral density
  vectors.} \item{\code{freq}}{The frequencies for each set
  of \code{stg_kopt} and \code{stg_psd}.} }
}
\examples{
\dontrun{
##
## rlpSpec working environment
##
# Get some status information about the .rlpenv environment
rlp_envStatus()
#
# Get a list of all variables in .rlpenv
rlp_envList()
#
# Pull the variable "init" into .GlobalEnv
print(x <- rlp_envGet("init"))
#
# pull the adaptive history into .GlobalEnv
get_adapt_history()
}
}
\seealso{
  \code{\link{rlpSpec-utilities}},
  \code{\link{char2envir}}, \code{\link{pspectrum}}

  \code{new.env}, \code{baseenv}
}

