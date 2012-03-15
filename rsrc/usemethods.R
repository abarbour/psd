### Define default methods for functions
pspectrum <- function(x, ...) UseMethod("pspectrum")
psdcore <- function(x, ...) UseMethod("psdcore")
riedsid <- function(x, ...) UseMethod("riedsid")
whiten <- function(x, ...) UseMethod("whiten")
qualcon <- function(x, ...) UseMethod("qualcon")
boundslope <- function(x, ...) UseMethod("boundslope")
# in suppfuncs
initEnv <- function(x, ...) UseMethod("initEnv")
envAssign <- function(x, ...) UseMethod("envAssign")
envGet <- function(x, ...) UseMethod("envGet")
nas <- function(x, ...) UseMethod("nas")
mod <- function(x, ...) UseMethod("mod")
ones <- function(x, ...) UseMethod("ones")
zeros <- function(x, ...) UseMethod("zeros")
###