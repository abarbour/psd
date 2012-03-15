###
###  Default method for pspectrum, the main function used for
###  adaptive estimation
###
pspectrum.default <- function(x, 
                      fsamp=1, 
                      tapcap=1e3, 
                      ntapinit=10, 
                      niter=2, 
                      ndec=1,
                      units=c("time","signal"),
                      plotpsd=TRUE, 
                      ylims=c(.07,3e4), xlims=c(0,0.5)) {
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
  initEnv(refresh=TRUE)
  
  # Cap the number of tapers to prevent runaway
  Cap <- abs(tapcap)
  if (Cap == 0 || Cap > 1e5){ Cap <- 1e3 }

  #  Number of refinement iterations usually <= 5
  envAssign("num_iter", abs(niter))
  Niter <- envGet("num_iter")
  
  #   ndec: number of actual psd calculations is n/ndec
  #   the rest are filled in with interpolation.  Sampling in
  #   frequency is variable to accommodate spectral shape
  lx <- length(x)
  if (lx < 10000) ndec <- 1
  #
  #            -----------------
  #  Get pilot estimate of psd with fixed number of tapers and no decimation
  msg <- sprintf("\t>>>> Pilot estimation with\t%i\ttapers\n",ntapinit)
  cat(msg)
  psd <- psdcore(x, ntaper=ntapinit, ndecimate=1, plotpsd=plotpsd, xlims=xlims)
  envAssign("num_freqs", length(psd))
  nf <- envGet("num_freqs")
  Ones <- ones(1, nf)  # row vec
  ntaper <- ntapinit * Ones
  
  if ( plotpsd && Niter > 0 ){
    require(RColorBrewer, quietly=TRUE, warn.conflicts=FALSE)
    pal <- brewer.pal(8, "Paired")
  }
  
  if ( Niter > 0){
    cat("\t\t>>>> Adaptive spectrum refinement:\n")
    for ( iterate in  1:Niter ) {
      cat(sprintf("\t\t\t>>>> taper optimization round\t%02i\n",iterate))
      kopt <- riedsid(psd, ntaper) # riedsid resets nf
      # choose the minimum between Cap and koopt
      ntaper <- t(as.matrix(apply(rbind(ones(1,nf)*Cap, kopt),2,min)))
      psd <- psdcore(x, ntaper=ntaper, ndecimate=ndec, 
                     plotpsd=plotpsd, plotcolor=pal[iterate])
    }
  }
  #  Scale to physical units and provide frequency vector
  psd <- psd/fsamp
  nf <- envGet("num_freqs")
  f <- t(t( seq.int(0, fsamp/2, length.out=nf) ))

  # psd class? [ ]
  psd.df <- data.frame(f=f, psd=psd, ntaper=t(ntaper))
  #psd.df <- list(freq=f, psd=psd, ntap=t(ntaper), call=match.call())
  # for method print to show call
  # move to psd method [ ]
#   cat("\t>>>> Results summary:\n")
#   print(summary(psd.df))
  return(invisible(psd.df))
  # return a structure.  From acf():
  #acf.out <- structure(.Data = list(acf = acf, type = type, 
  #                     n.used = sampleT, lag = lag, series = series, 
  #                     snames = colnames(x)), class = "acf")
  # from spec.pgram:
  #  spg.out <- list(freq = freq, spec = spec, coh = coh, phase = phase,
  #    kernel = kernel, df = df, bandwidth = bandwidth, n.used = N, 
  #    orig.n = N0, series = series, snames = colnames(x), 
  #    method = ifelse(!is.null(kernel), "Smoothed Periodogram", 
  #    "Raw Periodogram"), taper = taper, pad = pad, detrend = detrend, 
  #    demean = demean)
  #  class(spg.out) <- "spec"
  #
  #  any like this?
  #  .NotYetImplemented()
} 
# end pspectrum.default
###