.First.lib <- function(lib, pkg) {
  library.dynam("rlpSpec", pkg, lib)
}
.onLoad <- function(...) {
packageStartupMessage(
sprintf(
"Loaded rlpSpec (%s) -- Adaptive multitaper spectrum estimation.",
utils::packageVersion("rlpSpec")
)
)
}
