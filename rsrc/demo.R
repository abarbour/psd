setwd("/Users/abarbour/kook.processing/R/dev/packages/rlpSpec/rsrc")
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
  setwd("~/kook.processing/R/dev/packages/rlpSpec/rsrc")
  # orig port loc: matlab.port/psd/test.rfiles")
  ##
  ## ported source
  ##
  source("funcs2.R")
  source("psdcore.R")
  source("pspectrum.R")
  source("qualcon.R")
  source("riedsid.R")
  source("whiten.R")
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
  #
  ##  mag.dir <- "/Users/abarbour/kook.processing/R/dev/packages/rlpSpec/data"
  ##  mag <- cbind(
  ##	  read.table(paste(mag.dir,"mag.raw",sep="/"),
  ##		  colClasses="numeric", col.names="raw"),
  ##	  read.table(paste(mag.dir,"mag.clean",sep="/"), 
  ##		  colClasses="numeric", col.names="clean"))
  load("../data/mag.rda")
  ##
  ## adaptive estimation and plotting
  ##
  ##  raw: magnetometer data with outliers
  ##  clean: the same magnetometer data but cleaned of outliers
  ##
  ##  This demonstrates the importance of having 'clean' data for the purposes of spectrum
  ##  estimation, since they can badly bias low frequencies.
  ##
  psd.mag.raw <<- pspectrum(mag$raw)
  plot.psd(psd.mag.raw)
  ##
  ## since it's clear we have outliers, let's run qualcon
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
