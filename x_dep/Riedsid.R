ConstrainTapers <- function(ntap, maxslope=1){
  # calc slopes
  slopes <- splineGrad(ntap, deriv=1)
  #
  # "flatten" based on maxslope
  p <- .Call("r_constrain_tapers",
  	  as.vector(slopes), 
	  as.integer(maxslope), 
	  PACKAGE="rlpSpec")
  return(p)
}
