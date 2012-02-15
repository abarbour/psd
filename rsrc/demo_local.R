##setwd("/Users/abarbour/kook.processing/R/dev/packages/rlpSpec/rsrc")
run.demo <- function(){
  ## Some basic tests
  ##  to demonstrate the adaptive estimation process
  ##
  ## Args:	
  ##
  ## Returns:	
  ##
  ## TODO(abarbour):
  ##
  #setwd("~/kook.processing/R/dev/packages/rlpSpec/rsrc")
  # orig port loc: matlab.port/psd/test.rfiles")
  ##
  ## methods
  ##
  source("rsrc/usemethods.R")
  source("rsrc/S3_summary.R")
  source("rsrc/S3_print.R")
  source("rsrc/S3_plot.R")
  ##
  ## ported source
  ##
  source("rsrc/func_psdcore.R")
  source("rsrc/func_pspectrum.R")
  source("rsrc/func_qualcon.R")
  source("rsrc/func_riedsid.R")
  source("rsrc/func_supps.R")
  source("rsrc/func_whiten.R")
  ##
  ## packages
  ##
  library(signal)
  # for triang
  ##
  ##  the meat and potatoes
  ##
  set.seed(282)
  n<-1256*10
  ##
  noise<-rnorm(n,sd=0.2)
  signal<-.001*(1:n) + triang(n)+3*sin(pi*1:n/180+pi/4)+3*sin(2*pi*1:n/180+pi/4)
  x<-signal+noise
  x<-signal
  ##
  ## Begin Adaptive estimation and plotting
  ##
  psd.x <<- pspectrum(x)
  plot.psd(psd.x)
  ##
  ## now some real data...
  ##  mag.dir <- "/Users/abarbour/kook.processing/R/dev/packages/rlpSpec/data"
  ##  mag <- cbind(
  ##    read.table(paste(mag.dir,"mag.raw",sep="/"),
  ##		  colClasses="numeric", col.names="raw"),
  ##	  read.table(paste(mag.dir,"mag.clean",sep="/"), 
  ##		  colClasses="numeric", col.names="clean"))
  load("data/mag.rda")
  plot(mag$raw, type="s", main="Magnetometer data")
  text(200, 40, "raw series")
  lines(mag$clean-50, type="s", col="red")
  text(350, -100, "cleaned series (-50)", col="red")
  ##  raw: magnetometer data with outliers
  ##  clean: the same magnetometer data, but cleaned of outliers
  ##
  ##  This demonstrates the importance of having 'clean' data for the purposes of spectrum
  ##  estimation, since they can badly bias low frequencies (and in most cases not as many 
  ##  tapers need to be applied locally to smooth the spectrum).
  ##
  psd.mag.raw <<- pspectrum(mag$raw)
  plot.psd(psd.mag.raw)
  ##
  ## since it's clear we have outliers, let's run qualcon
  # 3 Feb-12 (r.2.14.1 throwing error - use colMean)
  #   qualcon(psd.mag.raw)
  ##
  ## remove outlilers, and redo spectrum estimation
  psd.mag.cln <<- pspectrum(mag$clean)
  plot.psd(psd.mag.cln)
  ##
  ## the difference is not negligible
  ##
} 
# end run.demo
