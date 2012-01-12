pspectrum <- function(x, 
                      fsamp=1, 
                      tapcap=1e3, 
                      ntapinit=10, 
                      niter=5, 
                      ndec=1,
                      units=c("time","signal"),
                      plot=TRUE, ylims=c(.07,3e4)) {
  ###
  # PORT of RLP's pspectrum.m
  # abarbour
  # Dec 2011
  #
  # porting: 10 Jan 2012
  #   [ ] adaptive estimation limited to decimation==1 at the moment
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
  Cap <- abs(tapcap)
  if (Cap == 0 || Cap > 1000){ Cap <- 1000 }
  #  Niter<-number of refinement iterations usually <= 5
  Niter <<- abs(niter)
  if (Niter == 0){ Niter <<- 1 }
  #   ndec: number of actual psd calculations is n/ndec
  #   the rest are filled in with interpolation.  Sampling in
  #   frequency is variable to accommodate spectral shape
  lx <- length(x)
  if (lx < 10000){
    ndec <- 1
  }
  # --- env
  #   psdenv <- new.env(parent=baseenv())
  #
  #            -----------------
  #  Get pilot estimate of psd with fixed number of tapers and no decimation
  psd <- psdcore(x, ntaper=ntapinit, ndecimate=1, plot=plot)
  nf <<- length(psd)
  ones <- matrix(1,1,nf)  # row vec
  #  Iterative refinement of spectrum 
  ntaper <- ntapinit * ones
  cat("\t>>>> Adaptive estimation:\n")
  for ( iterate in  1:Niter ) {
    cat(sprintf("\t\t>>>> taper optimization round\t%02i\n",iterate))
    kopt <- riedsid(psd, ntaper)
    # riedsid resets nf
    ntaper <- t(as.matrix(apply(rbind(matrix(1,1,nf)*Cap, kopt),2,min)))
    psd <- psdcore(x, ntaper=ntaper, ndecimate=ndec, plot=FALSE)
  }
  #  Scale to physical units and provide frequency vector
  psd <- psd/fsamp
  nf <<- length(psd)
  f <- t(t(seq(0, fsamp/2, length.out=nf)))
  if (plot==TRUE) {
    lims <- ylims
    par(las=1)
    plot(f[2:nf], psd[2:nf], 
         main="Adaptive Sine-multitaper PSD Estimation",
         sub=sprintf("%i iterations",Niter),
         log="y",
         ylab=sprintf("PSD rel. 1 %s**2 * N * dT",units[2]),
         ylim=lims, 
         xlab=sprintf("Freq. in 1/N/dT (dT in %s)",units[1]), 
         xlim=c(0,fsamp/2), xaxs="i",
         type="s")
  }
  # psd class? [ ]
  psd.df <- data.frame(f=f, psd=psd, ntaper=t(ntaper))
  cat("\t>>>> Results summary:\n")
  print(summary(psd.df))
  return(invisible(psd.df))
} # end pspectrum
