.First.lib <- function(lib, pkg) {
  library.dynam("rlpSpec", pkg, lib)
}

.onLoad <- function(...) {
  NULL
}

# executed after .onLoad is executed
.onAttach <- function(...) {
  ## environment string equivalent to .rlpenv <- ".rlpSpecEnv"
  assign(".rlpenv", ".rlpSpecEnv", envir=.GlobalEnv)
  ## is this really needed .onAttach?
  ## rlp_initEnv(envir=.rlpenv, refresh=TRUE, verbose=TRUE)
  
  packageStartupMessage(
    sprintf("Loaded rlpSpec (%s) -- Adaptive multitaper spectrum estimation.",
            utils:::packageVersion("rlpSpec")))
}

.onUnload <- function(libpath)
{
  library.dynam.unload("rlpSpec", libpath)
}