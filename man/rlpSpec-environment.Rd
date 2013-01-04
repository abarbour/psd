\docType{methods}
\name{rlpSpec-environment}
\alias{rlpSpec-environment}
\alias{rlp_envAssign}
\alias{rlp_envAssignGet}
\alias{rlp_envGet}
\alias{rlp_envList}
\alias{rlp_envStatus}
\alias{rlp_initEnv}
\title{Various environment manipulation functions.}
\usage{
  rlp_initEnv(envir = .rlpenv, refresh = FALSE,
    verbose = TRUE, ...)

  rlp_initEnv.default(envir = .rlpenv, refresh = FALSE,
    verbose = TRUE, ...)

  rlp_envClear(...)

  rlp_envStatus(envir = .rlpenv)

  rlp_envStatus.default(envir = .rlpenv)

  rlp_envList(envir = .rlpenv)

  rlp_envList.default(envir = .rlpenv)

  rlp_envGet(variable, envir = .rlpenv)

  rlp_envGet.default(variable, envir = .rlpenv)

  rlp_envAssign(variable, value, envir = .rlpenv)

  rlp_envAssign.default(variable, value, envir = .rlpenv)

  rlp_envAssignGet(variable, value, envir = .rlpenv)

  rlp_envAssignGet.default(variable, value,
    envir = .rlpenv)
}
\arguments{
  \item{envir}{character string for the environment used.
  See section \strong{Defaults}}

  \item{refresh}{logical; should the contents of the
  environment be trashed?}

  \item{verbose}{logical; should messages be given?}

  \item{variable}{character; the name of the variable to
  get or assign}

  \item{value}{character; the name of the variable to
  assign}

  \item{...}{For \code{rlp_envClear}: arguments passed to
  \code{rlp_initEnv}}

  \item{...}{For \code{rlp_initEnv}: arguments passed to
  \code{new.env}}
}
\value{
  \code{rlp_initEnv} returns (invisibly) the result of
  \code{rlp_envStatus}.
}
\description{
  The computation of \emph{adaptive} power spectral density
  estimates requires bookkeeping and non-destructive
  manipulation of variables. The functions here are mainly
  convenience wrappers designed to maintain variable
  separation from the \code{.GlobalEnv} environment so that
  no innocent variable is destroyed in the process of
  iteratively computing spectra.

  \code{rlp_initEnv} initializes the \code{.rlpenv}
  environment with an option to clear the contents (if the
  environment already exists).

  \code{rlp_envClear} clears the contents of the
  environment.

  \code{rlp_envStatus} returns a list of some information
  regarding the status of the environment.

  \code{rlp_envList} returns a listing of the assignments.

  \code{rlp_envGet} returns a the value of \code{variable}.

  \code{rlp_envAssign} assigns \code{value} to
  \code{variable}, but does not return it.

  \code{rlp_envAssignGet} both assigns and returns a value.
}
\note{
  \code{rlp_initEnv} will not re-initialize the enviroment,
  unless told to do so with \code{refresh=TRUE}.

  \code{rlp_envClear} does \emph{not} remove the
  environment--simply the assignments within it.
}
\section{Defaults and Initialization}{
  By default, these functions all use \code{envir=.rlpenv}
  to set the enviroment-name string; \code{.rlpenv} is set
  when attaching the package.

  If the environment has not yet been initialized (it
  should never need to be prior to running
  \code{pspectrum}) the command \code{rlp_initEnv()} should
  be used. If a fresh environment is desired, and the
  environment already exists, the command
  \code{rlp_envClear()} (which is really just an alias for
  \code{rlp_initEnv(refresh=TRUE)}) can be used.

  One could set \code{.rlpenv} to another string, if a
  different environment is desired.
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
\seealso{
  \code{\link{rlpSpec-utilities}},
  \code{\link{char2envir}}, \code{\link{pspectrum}}

  \code{\link{new.env}}, \code{\link{baseenv}},
  \code{\link{rlp_envStatus}}

  \code{\link{rlp_initEnv}}
}

