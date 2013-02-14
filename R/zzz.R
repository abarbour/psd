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
##
.onLoad <- function(...) {
  ## DLL
  #library.dynam("rlpSpec", pkg, lib)
  # useDynLib(rlpSpec) in NAMESPACE though... any conflicts?
}

.onUnload <- function(libpath)
{
  library.dynam.unload("rlpSpec", libpath)
}
##
# executed after .onLoad is executed, once the namespace is visible to user
.onAttach <- function(...) {
  ##
  ## env pointer
  rlpSpec:::rlp_envGlobals()
  ## initialize it
  rlpSpec:::rlp_initEnv(refresh=TRUE, verbose=FALSE)
  ## add some info
  rlpSpec:::rlp_envAssign("init", paste(rlpSpec:::rlp_envGet("init"), "(upon attach)"))
  ##
  packageStartupMessage(
    sprintf("Loaded rlpSpec (%s) -- Adaptive multitaper spectrum estimation.",
            utils:::packageVersion("rlpSpec")))
}
.Last.lib <- function(...){
  NULL
}
