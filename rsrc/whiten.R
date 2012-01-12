whiten <- function(x) {
  ###
  # PORT of RLP's riedsid.m
  # abarbour
  # Dec 2011
  #
  # porting: Jan 3, 2012
  # testing: 
  #     Jan 3, 2012: 
  #        appears to work fine with 'filter' but not as well with 'open'
  #        noise<-rnorm(n,sd=0.2)
  #        signal<-.001*(1:n) + triang(n)+3*sin(pi*1:n/180+pi/4)+3*sin(2*pi*1:n/180+pi/4)
  ###
  #
  #  Applies a short convolution filter to time series x that
  #  flattens its spectrum (making it more nearly white).
  #  The resultant approximately white-noise time series makes
  #  nonstationarity and error spikes much more prominent.  
  #  The new series can be used for preliminary data screening.
  # 
  #  Method:  solves Yule-Walker equations using autocorrelation
  #  functions based on a sine-multitaper spectrum.
  #
  ##
  ## Args:	
  ##
  ## Returns:	
  ##
  ## TODO(abarbour):
  ##
  ## Prelims
  nord <- 4
  ntap <- 10
  nx <- length(x)
  if (nx < 10){ 
    print('Series too short to whiten') 
    return
  }
  ##
  #  Calculate psd with multitapers (fixed num tapers)
  spec <- psdcore(x, ntap, plotpsd=FALSE)
  nf <- length(spec)        
  #  Fourier transform the psd to get the autocovariance function
  autoc <- Re(fft( c(spec[1]/2, spec[2:nf-1], matrix(0,nf-1,1)) ))
  autoc.n <- autoc / autoc[1]
  #  Form the Yule-Walker equations and solve them
  Toe <- toeplitz(autoc.n[1:nord])
  coeffs <- solve(Toe, autoc.n[2:(nord+1)])
  #  Convolve filter over original time series
  # 'filter' == 'same' (matlab), confirm [ ]
  whitex <- convolve(x-mean(x), c(-1,coeffs), type="filter")
  ##
  return(invisible(whitex))
}
# end whiten
