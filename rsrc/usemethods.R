### Define default methods for functions
pspectrum <- function(x, ...) UseMethod(".pspectrum")
psdcore <- function(x, ...) UseMethod(".psdcore")
.devpsdcore <- function(x, ...) UseMethod("..dev_psdcore")

whiten <- function(x, ...) UseMethod(".whiten")
qualcon <- function(x, ...) UseMethod(".qualcon")
boundslope <- function(x, ...) UseMethod(".boundslope")
# in suppfuncs
initEnv <- function(x, ...) UseMethod(".rlp_initEnv")
envStatus <- function(...) UseMethod(".rlp_envStatus")
envClear <- function(...) initEnv(refresh=TRUE, ...)
envList <- function(x, ...) UseMethod(".rlp_envList")
envAssign <- function(x, ...) UseMethod(".rlp_envAssign")
envGet <- function(x, ...) UseMethod(".rlp_envGet")
envAssignGet <- function(x, ...) UseMethod(".rlp_envAssignGet")
nas <- function(x, ...) UseMethod(".nas")
mod <- function(x, ...) UseMethod(".mod")
reshape_vector <- function(x, ...) UseMethod(".reshape_vector")
as.colvec <- function(x) reshape_vector(x, "vertical")
as.rowvec <- function(x) reshape_vector(x, "horizontal")
ones <- function(x, ...) UseMethod(".ones")
zeros <- function(x, ...) UseMethod(".zeros")
#
splineGrad <- function(x, ...) UseMethod(".splineGrad")
###
