##
##  Default method for riedsid, the Riedel & Sidorenko taper optimization
##
.riedsid.default <- function(spec, ntaper) {
  ###
  # PORT of RLP's riedsid.m
  # abarbour
  # Dec 2011
  #
  # porting:  Jan 3, 2012
  # testing:
  ###
  #
  #  Estimates optimal number of tapers at each frequency of
  #  given spec, based on Riedel-Sidorenko MSE recipe and other
  #  tweaks due to RLP.  
  #  ntaper is the vector of taper lengths in the previous iteration.
  #
  ##
  ## Args:	
  ##
  ## Returns:	
  ##
  ## TODO(abarbour):
  ##

  message("\t\ttaper optimization")
  eps <- 1e-78  #  A small number to protect against zeros
  
  spec <- as.vector(spec)
  
  nf <- envAssignGet("num_freqs", length(spec))
  
  Ones <- ones(nf) # was rowvec
  Zeros <- zeros(nf) # was rowvec
  
  if (length(ntaper)==1) { 
    ntap <- ntaper*Ones
  } else {
    ntap <- ntaper
  }
  # find the minimum by column for nf/2, 
  nspan <<- as.matrix(apply(cbind(nf*Ones/2, 7*ntap/5),1,min))
  
  # rowvec:
  #nspan <- t( t( round( apply(rbind(0.5*nf*Ones, 1.4*ntap), 2, min) ) ) )
  
  #  Create log spec, and pad to handle begnning and end values
  nadd <- 1 + max(nspan)
  # [spec(nadd:-1:2); spec; spec(nf-1:-1:nf-nadd)]
  Y <- log(eps + c(spec[nadd:2], spec, spec[(nf-1):(nf-nadd)]))
  #  R <- spec"/spec <- Y" + (Y')^2  2nd form preferred for consistent smoothing
  #  
  d2Y <- t(Zeros)
  dY <- d2Y
  for (  j  in  1:nf ) {
    j1 <- j-nspan[j]+nadd-1 
    j2 <- j+nspan[j]+nadd-1 
    #  Over an interval proportional to taper length, fit a least
    #  squares quadratic to Y to estimate smoothed 1st, 2nd derivs
    jr <- j1:j2
    u <- jr - (j1+j2)/2
    u2 <- u^2
    L <- j2-j1+1 
    L2 <- L^2
    uzero <- (L2 - 1)/12
    # first deriv
    dY[j,1] <- u %*% Y[jr] * (12/(L*(L2 - 1)))
    # second deriv
    d2Y[j,1] <- (u2 - uzero) %*% Y[jr] * (360/(L*(L2 - 1)*(L2 - 4)))
  }
  #
  #  Riedel-Sidorenko recipe (eq 13): kopt <- (12*abs(spec ./ d2spec)).^0.4 but
  #  parabolic weighting in psdcore requires: (480)^0.2*abs(spec./d2spec)^0.4
  #
  #  Original form:  kopt <- 3.428*abs(spec ./ d2spec).^0.4
  kopt <- round( 3.428 / abs(eps + d2Y + dY^2)^0.4 )
  #  Curb run-away growth of kopt due to zeros of spec'' limits
  #  slopes to be < 1 in magnitude, preventing greedy averaging:
  #  Scan forward and create new series where slopes <= 1
  state<-0
  for ( j  in  2:nf ) {
    if (state == 0) {
       slope <- kopt[j]-kopt[j-1]
       if (slope >= 1 ) {
         state <- 1
         kopt[j] <- kopt[j-1]+1
       }
    } else {
       if (kopt[j] >= kopt[j-1]+1) {
         kopt[j] <- kopt[j-1]+1
       } else {
         state <- 0
       }
    }
  }
  #  Scan backward to bound slopes >= -1
  state <- 0
  for ( j  in  nf:2 ) {
    if (state == 0) {
       slope <- kopt[j-1]-kopt[j]
       if (slope >= 1) {
         state <- 1
         kopt[j-1] <- kopt[j]+1
       }
    } else {
       if (kopt[j-1] >= kopt[j]+1) {
         kopt[j-1] <- kopt[j]+1
       } else {
         state <- 0
       }
    }
  }
  ##  Never average over more than the spec length!
  ## colMeans or rowMeans [ ]
  kopt <- t( apply( rbind(t(kopt), Ones*round(nf/2)), 2, min) )
  ##
  return(invisible(kopt))
} 
# end riedsid.default
###