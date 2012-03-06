### Define default methods for functions
pspectrum <- function(x, ...) UseMethod("pspectrum")
psdcore <- function(x, ...) UseMethod("psdcore")
riedsid <- function(x, ...) UseMethod("riedsid")
whiten <- function(x, ...) UseMethod("whiten")
qualcon <- function(x, ...) UseMethod("qualcon")
boundslope <- function(x, ...) UseMethod("boundslope")
# in suppfuncs
mod <- function(x, ...) UseMethod("mod")
###