### Define default methods for functions
pspectrum <- function(x, ...) UseMethod(".pspectrum")
psdcore <- function(x, ...) UseMethod(".psdcore")
.devpsdcore <- function(x, ...) UseMethod("..dev_psdcore")
riedsid <- function(x, ...) UseMethod(".riedsid")
whiten <- function(x, ...) UseMethod(".whiten")
qualcon <- function(x, ...) UseMethod(".qualcon")
boundslope <- function(x, ...) UseMethod(".boundslope")
# in suppfuncs
initEnv <- function(x, ...) UseMethod(".rlp_initEnv")
envList <- function(x, ...) UseMethod(".rlp_envList")
envAssign <- function(x, ...) UseMethod(".rlp_envAssign")
envGet <- function(x, ...) UseMethod(".rlp_envGet")
envAssignGet <- function(x, ...) UseMethod(".rlp_envAssignGet")
nas <- function(x, ...) UseMethod(".nas")
mod <- function(x, ...) UseMethod(".mod")
colvec <- function(x, ...) UseMethod(".colvec")
rowvec <- function(x, ...) UseMethod(".rowvec")
ones <- function(x, ...) UseMethod(".ones")
zeros <- function(x, ...) UseMethod(".zeros")
#
splineGrad <- function(x, ...) UseMethod(".splineGrad")
###
