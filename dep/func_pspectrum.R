..todep.pspectrum.default <- function(x, 
                                      fsamp=1, 
                                      tapcap=1e3, 
                                      ntapinit=10, 
                                      niter=3, 
                                      ndec=1,
                                      units=c("time","signal"),
                                      prewhiten=TRUE,
                                      plotpsd=TRUE, 
                                      ylims=c(.07,3e4), xlims=c(0,0.5),
                                      devmode=FALSE) {
  ###
  # PORT of RLP's pspectrum.m
  # abarbour
  # Dec 2011
  #
  # porting: 10 Jan 2012
  #   [fixed 3-Feb-12] adaptive estimation limited to decimation==1 at the moment
  # testing:
  ###
  #
  #  Adaptive multitaper estimator of power spectral density (psd) of
  #  the stationary time series x .
  #
  #  fsamp is the sampling frequency = 1/(sampling interval).  If fsamp
  #  is absent, use fsamp=1.
  #
  #  psd of length nf gives spectrum at nf evenly spaced frequencies: 
  #     f<-[ 0, df, 2*df, ... (nf-1)*df]' ('), where nf = 1+ n/2, and n=length(x),
  #     and df=1/T.  If n is odd, x is truncated by 1.
  #
  # -------  Tuning parameters -------------
  #   tapcap=maximum number of tapers allowed per freq then uncertainty
  #   of estimates >= psd/sqrt(Cap).
  ##
  ## Args:  
  ##
  ## Returns:	
  ##
  ## TODO(abarbour):	
  ##    [ ] fix frequency vector when decimation is > 1 (only a plotting issue really)
  ##    [ ] structured return
  ##    [ ] type <- match.arg(type)
  ##
  # --- setup the environment ---
  rlp_initEnv(refresh=TRUE)
  
  # Cap the number of tapers to prevent runaway
  Cap <- abs(tapcap)
  if (Cap == 0 | Cap > 1e5){ Cap <- 1e3 }
  
  #  Number of refinement iterations usually <= 5
  Niter <- rlp_envAssignGet("num_iter", abs(niter))
  
  ntapinit <- rlp_envAssignGet("init_tap", abs(ntapinit))
  
  #   ndec: number of actual psd calculations is n/ndec
  #   the rest are filled in with interpolation.  Sampling in
  #   frequency is variable to accommodate spectral shape
  lx <- length(x)
  if (lx < 10000) ndec <- 1
  #
  PSDFUN <- psdcore
  if (devmode) {
    warning("operating in development mode")
    PSDFUN <- .devpsdcore
  }
  #            -----------------
  #  Get pilot estimate of psd with fixed number of tapers and no decimation
  message(sprintf("\t>>>> Pilot estimation with\t%i\ttapers",ntapinit))
  psd <- PSDFUN(x, ntaper=ntapinit, 
                ndecimate=1, 
                plotpsd=FALSE, xlims=xlims, as.spec=TRUE)
  plot(psd)
  ## add this to psdcore return?
  nf <- nrow(psd$freq)
  ##nf <- length(psd)
  rlp_envAssign("num_freqs", nf)
  Ones <<- ones(nf)  # row vec of ones
  ntaper <<- ntapinit * Ones
  
  if ( plotpsd && Niter > 0 ){
    require(RColorBrewer, quietly=TRUE, warn.conflicts=FALSE)
    pal <- brewer.pal(8, "Paired")
  }
  
  if ( Niter > 0){
    message("\t\t>>>> Adaptive spectrum refinement:")
    for ( iterate in  1:Niter ) {
      message(sprintf("\t\t\t>>>> taper optimization round\t%02i",iterate))
      ## calculate optimal tapers
      kopt <- riedsid(psd$spec, ntaper)
      stopifnot(exists('kopt'))
      ## contrain the total number of tapers
      ntaper <- kopt
      ntaper[ntaper>Cap] <- Cap
      psd <- PSDFUN(x, ntaper=ntaper, ndecimate=ndec, 
                    plotpsd=FALSE, plotcolor=pal[iterate], as.spec=TRUE)
      plot(psd, col=pal[iterate], add=T)
    }
  }
  ##  Scale to physical units and provide frequency vector
  ## [ ] ?
  psd$spec <- psd$spec/fsamp
  psd$freq <- psd$freq*fsamp
  psd$sample.rate <- fsamp
  ##
  return(invisible(psd))
} 
# end pspectrum.default
###
