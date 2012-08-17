##
##  Default method for riedsid, the Riedel & Sidorenko taper optimization
##
.riedsid.default <- function(spec, ntaper, 
                             restrict.deriv=c("slope",
                                              "loess",
                                              "friedman.super",
                                              "none")) {
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

  ##message("\t\ttaper optimization")
  ##
  eps <- 1e-78  #  A small number to protect against zeros
  
  spec <- as.vector(spec)
  
  nf <- envAssignGet("num_freqs", length(spec))
  
  Ones <- ones(nf) # was rowvec, now col
  Zeros <- zeros(nf) # was rowvec, now col
  
  if (length(ntaper)==1) { 
    ntap <- ntaper*Ones
  } else {
    ntap <- ntaper
  }
  # find the minimum by column for 1/2 nf, 7/5 ntap
  # ones is colvec, ntap is colvec
  minmat <- cbind(nf*Ones/2, 7*ntap/5)
  # margin==1 produces a colvec of rowwise minimums
  nspan <- as.matrix(apply(minmat, MARGIN=1, FUN=base::min))
  # when it was rowvec:
  #nspan <- t( t( round( apply(rbind(0.5*nf*Ones, 1.4*ntap), 2, min) ) ) )
  
  #  Create log spec, and pad to handle begnning and end values
  nadd <- 1 + max(nspan)
  # [spec(nadd:-1:2); spec; spec(nf-1:-1:nf-nadd)]
  Y <- log(eps + c(spec[nadd:2], spec, spec[(nf-1):(nf-nadd)]))
  #  R <- spec"/spec <- Y" + (Y')^2  2nd form preferred for consistent smoothing
  dY <- d2Y <- Zeros # was row: t(Zeros)
  nadd1 <- nadd-1
  for (  j  in  1:nf ) {
    j1 <- j-nspan[j]+nadd1 
    j2 <- j+nspan[j]+nadd1 
    #  Over an interval proportional to taper length, fit a least
    #  squares quadratic to Y to estimate smoothed 1st, 2nd derivs
    jr <- j1:j2
    u <- jr - (j1 + j2)/2 # rowvec 
    u2 <- u*u             # rowvec
    L <- j2-j1+1          # constant
    L2 <- L*L             # constant
    LL2 <- L*L2           # constant
    LL2L <- LL2 - L       # constant
    uzero <- (L2 - 1)/12  # constant
    # first deriv
    dY[j] <- u %*% Y[jr] * 12 / LL2L
    # second deriv
    d2Y[j] <- (u2 - uzero) %*% Y[jr] * 360/LL2L/(L2-4)
  }
  #
  #  Riedel-Sidorenko recipe (eq 13): kopt <- (12*abs(spec ./ d2spec)).^0.4 but
  #  parabolic weighting in psdcore requires: (480)^0.2*abs(spec./d2spec)^0.4
  #
  #  Original form:  kopt <- 3.428*abs(spec ./ d2spec).^0.4
  kopt <- round( 3.428 / abs(eps + d2Y + dY*dY)^0.4 )
  #
  restrict.deriv <- match.arg(restrict.deriv)
  do.restrict <- ifelse(restrict.deriv=="none",FALSE,TRUE)
  if (do.restrict){
    ##message(restrict.deriv)
    kseq <- 1:nf
    #  Curb run-away growth of kopt due to zeros of spec'' limits
    if (restrict.deriv=="slope"){
      #  slopes to be < 1 in magnitude, preventing greedy averaging:
      #  Scan forward and create new series where slopes <= 1
      derivs <- splineGrad(kopt, kseq, plot.derivs=F)
      slopes <- derivs$dydx
      state<-0
      for ( j  in  2:nf ) {
        slope <- slopes[j]
        if (state == 0) {
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
    } else if (restrict.deriv=="loess"){
      # Loess regression smoothing
      kopt <- as.matrix(predict(
        stats::loess(y ~ x, data.frame(x=kseq,y=kopt), span=.08, degree=0)),
                        data.frame(x=kseq))
    } else if (restrict.deriv=="friedman.super"){
      # Friedman's super smoother.
      # Meh.
      kopt <- as.matrix(stats::supsmu(kseq,kopt,, span=.3, bass=2)$y)
    } else{
      # nothing
      do.restrict <- FALSE
    }
  }
  if (!do.restrict){
    warning("taper optimization is unbound")
  }
  ##  Never average over more than the spec length!
  ## kopt is colvec, Ones is colvec
  minmat2 <- cbind(kopt, Ones*round(nf/2))
  # margin==1 priduce colvec of rowwise minimums
  kopt.bound <- as.matrix( apply(minmat2, MARGIN=1, FUN=base::min) )
  ##
  return(invisible(kopt.bound))
} 
# end riedsid.default
###