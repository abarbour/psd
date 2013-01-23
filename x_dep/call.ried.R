rm(list=ls())
system("rm ctap_simple.*o")
system("R CMD SHLIB ctap_simple.c")
dyn.load("ctap_simple.so")
source("../R/func_utils.R")
sim.ntaps <- function(){
	slp1 <- seq(1,5,by=1) - 1
	slp2 <- 2*slp1 - 2
	slp5 <- 5*slp1 - 5
	fwd <- as.numeric(
	c(slp1, slp2+5, slp5+11)
	)
	return(c(fwd, rev(fwd)))
}
ntaps.orig <- sim.ntaps()
#slopes <- as.numeric(.splineGrad.default(1:length(ntaps), ntaps)$dydx)
maxslope <- round(0.8)
(ntaps.adj <- .Call("rlp_constrain_tapers", ntaps.orig, maxslope))
all.equal(length(ntaps.orig),length(ntaps.adj))
#
plt <- function(nt=sim.ntaps(), nta=ntaps.adj){
  message("plotting.")
  plot(nt)
  lines(nta, col="red", lty=3)
}
plt()
