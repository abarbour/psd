###
###  Default method for qualcon, the quality control (whitening) system
###
qualcon.default <-function(x, critvar=10) {
  ###
  # PORT of RLP's qualcon.m
  # with mods by abarbour
  # Dec 2011
  #
  # porting:  4 Jan 2012
  # testing:  4 Jan 2012
  #           noise: rnorm(n,sd=0.2)
  #           sig:   .001*(1:n) + triang(n)+3*sin(pi*1:n/180+pi/4)+3*sin(2*pi*1:n/180+pi/4)
  #           temporary outliers are easily found
  #              x[100:120] <- 5.4
  #              x[1000:1020] <- 5.5
  #              
  ###
  #
  #  Searches for data errors indicated by the whitened series
  #  falling outside fixed limits.
  #  Plots +-100 points of x and whitened x around error section.
  #  badpoints <- number error regions detected
  #
  ##
  ## Args:	
  ##
  ## Returns:	
  ##
  ## TODO(abarbour):	
  ##
  lx <- length(x)
  x.w <- whiten(x)
  wcrit <- critvar * sqrt(var(x.w))
  par(ask=TRUE)
  toret <- show.white(wcrit, 1:lx, x, x.w, pltlabs=FALSE)
  Icrit <- which(abs(x.w) > wcrit)
  lIc <- length(Icrit)
  if (lIc > 0){
    badpoints <- 0
    i2 <- 0
    for ( j  in  1:length(Icrit) ) {
      if ( Icrit[j] >= i2 ){
        i1 <- Icrit[j]-100 
        i2 <- Icrit[j]+100
        k1 <- max(1, i1)  
        k2 <- min(i2, lx)
        k1k2 <- k1:k2
        show.whiten(wcrit, k1k2, x, x.w, limsc=2.8)
        badpoints <- 1 + badpoints
      }
    }
  }
  return(invisible(toret))
} 
# end qualcon.default
###