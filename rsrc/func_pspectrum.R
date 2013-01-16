###
###  Default method for pspectrum, the main function used for
###  adaptive estimation
###
pilot_spec <- function(...) UseMethod("pilot_spec")
pilot_spec.default <- function(x, x.frequency=1, ntap=10, ...){
  stopifnot(length(ntap)==1)
  # initial spectrum:
  Pspec <- psdcore(X.d=x, X.frq=x.frequency, ntaper=ntap, 
                   ndecimate=1L, demean=TRUE, detrend=TRUE, as.spec=TRUE) #, ...)
  num_frq <- length(Pspec$freq)
  num_tap <- length(Pspec$taper)
  stopifnot(num_tap <= num_frq)
  Ptap <- Pspec$taper
  # generate a series, if necessary
  if (num_tap < num_frq) Ptap <- rep.int(Ptap[1], num_frq)
  # return taper object
  Pspec$taper <- as.taper(Ptap, setspan=TRUE)
  #
  return(Pspec)
}

new_histlist <- function(adapt_stages){
  histlist <- vector("list", 3) # freq, list-tap, list-psd
  names(histlist) <- c("freq", "stg_kopt", "stg_psd")
  num_pos <- 1 + adapt_stages # pilot + adapts
  histlist[[2]] <- histlist[[3]] <- vector("list", adapt_stages+1)
  rlp_envAssignGet("histlist",histlist)
}

update_adapt_history <- function(stage, ntap, psd, freq=NULL){
  histlist <- rlp_envGet("histlist")
  # stage==0 <--> index==1
  stg_ind <- stage+1
  nulfrq <- is.null(freq)
  if (!nulfrq) histlist$freq <- freq
  histlist$stg_kopt[[stg_ind]] <- ntap
  histlist$stg_psd[[stg_ind]] <- psd
  if (is.null(histlist$freq) & stage>0) warning("freqs absent despite non-pilot stage update")
  rlp_envAssignGet("histlist",histlist)
}

adapt_message <- function(stage){
  stopifnot(stage>=0)
  if (stage==0){stage <- paste(stage,"(pilot)")}
  message(sprintf("Stage  %s  estimation", stage))
}

pspectrum <- function(x, x.frqsamp, ntap_pilot=10, niter=2, verbose=TRUE, ...) UseMethod("pspectrum")
pspectrum.default <- function(x, x.frqsamp, ntap_pilot=10, niter=2, verbose=TRUE, ...){
  for (stage in 0:niter){
    if (verbose) adapt_message(stage)
    if (stage==0){
      # --- setup the environment ---
      rlp_initEnv(refresh=TRUE,verbose=verbose)
      # --- pilot spec ---
      Pspec <- pilot_spec(x=x, x.frequency=x.frqsamp, ntap=ntap_pilot)
      # --- history ---
      save_hist <- ifelse(niter < 10, TRUE, FALSE)
      if (save_hist){
        new_histlist(niter)
        update_adapt_history(0, Pspec$taper, Pspec$spec, Pspec$freq)
      }
      x <- 0 # to prevent passing orig data back/forth
    }
    ## calculate optimal tapers
    kopt <- riedsid(Pspec$spec, Pspec$taper)
    stopifnot(exists('kopt'))
    ## reapply to spectrum
    Pspec <- psdcore(X.d=x, ntaper=kopt)
    ## update history
    if (save_hist) update_adapt_history(stage, Pspec$taper, Pspec$spec)
  }
  return(Pspec)
}

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
