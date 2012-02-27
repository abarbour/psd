##
##  Prepare the raw text files here for use in rlpSpec
##  Feb 2012
##
## for the NA section interp/impute/etc.
##
source('/Users/abarbour/kook.processing/R/dev/timetasks/merge/funcs.R')
##
dtStrp <- function(x){
  #"POSIXct" no good for colClasses, so manually strptime
  strptime(x, format="%Y-%m-%dT%H:%M:%OS", tz="UTC")
}
##
readBSM <- function(fi, nalevel=1e4){
  toret <- read.table(fi, h=TRUE, 
                      colClasses=c("character", rep("numeric",4)), 
                      na.strings="99999", 
                      col.names=c("dt","CH0","CH1","CH2","CH3"),
                      stringsAsFactors=FALSE)
  toret$CH0[toret$CH0>nalevel] <- NA
  toret$CH1[toret$CH1>nalevel] <- NA
  toret$CH2[toret$CH2>nalevel] <- NA
  toret$CH3[toret$CH3>nalevel] <- NA
  print(summary(toret))
  return(invisible(toret))
}
##
procBSM <- function(dat){
  dat$dt <- dtStrp(dat$dt)
  dat$CH0 <- naOP(dat$CH0, quiet=FALSE)
  dat$CH1 <- naOP(dat$CH1, quiet=FALSE)
  dat$CH2 <- naOP(dat$CH2, quiet=FALSE)
  dat$CH3 <- naOP(dat$CH3, quiet=FALSE)
  return(invisible(dat))
}
##
library(RColorBrewer)
pal <- brewer.pal(n=4,name="Paired")
##
########
########
##
## a few hours of 20Hz BSM data, linearized
## Varian, and Pinon Flat Observatory
##
bsm <- list(varian="data/bsm/B073.ALL_20.l.txt.gz",
            pinyon="data/bsm/B084.ALL_20.l.txt.gz")
bsm$varian_dat <- procBSM(readBSM(bsm$varian))
bsm$pinyon_dat <- procBSM(readBSM(bsm$pinyon))
save(bsm, file="bsm.rda")

# plot(dat$dt, dat$CH0,type="s", ylim=15*c(-1,0), col=pal[1])
# lines(dat$dt, dat$CH1,type="s", col=pal[2])
# lines(dat$dt, dat$CH2,type="s", col=pal[3])
# lines(dat$dt, dat$CH3,type="s", col=pal[4])
# 
# psbsm<-pspectrum(dat$CH0, niter=2)
# rpsbsm<-spectrum(dat$CH0, pad=1, taper=0.2, spans=c(50,50), log="dB")
# 
# plot(log10(20*rpsbsm$freq),10*log10(rpsbsm$spec),type="s",ylim=40*c(-1,1),xlim=c(-3,1))
# lines(log10(20*psbsm$f),10*log10(psbsm$psd),type="s",col="red")                                                          


