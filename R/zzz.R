# ##
# Anything needed for the functioning of the namespace should 
# be handled at load/unload times by the .onLoad and .onUnload 
# hooks. For example, DLLs can be loaded (unless done by a 
#                                         useDynLib directive in 
#                                         the ‘NAMESPACE’ file) 
# and initialized in .onLoad and unloaded in .onUnload. 
# Use .onAttach only for actions that are needed only when the 
# package becomes visible to the user (for example a start-up message) 
# or need to be run after the package environment has been created.
# ##
## .First.lib is a defunct(ion)
# .First.lib <- function(lib, pkg) {
#   library.dynam("rlpSpec", pkg, lib)
# }

##
.onLoad <- function(...) {
  ## DLL
  #library.dynam("rlpSpec", pkg, lib)
  # useDynLib(rlpSpec) in NAMESPACE though... any conflicts?
  ## env pointer
  assign(".rlpenv", ".rlpSpecEnv", envir=.GlobalEnv)
}
.onUnload <- function(libpath)
{
  if (exists(.rlpenv, envir=.GlobalEnv)) rm(.rlpenv, envir=.GlobalEnv)
  library.dynam.unload("rlpSpec", libpath)
}

##
# executed after .onLoad is executed, once the namespace is visible to user
.onAttach <- function(...) {
  ##
  rlpSpec:::rlp_initEnv(envir=.rlpenv, refresh=TRUE, verbose=FALSE)
  .rlpSpecEnv$init <- paste(.rlpSpecEnv$init,"(upon attach)")
  ##
  packageStartupMessage(
    sprintf("Loaded rlpSpec (%s) -- Adaptive multitaper spectrum estimation.",
            utils:::packageVersion("rlpSpec")))
}
.Last.lib <- function(...){
  NULL
}
