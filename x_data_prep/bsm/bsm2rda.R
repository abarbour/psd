##
##  Prepare the raw text files here for use in rlpSpec
##  Feb 2012
##
## for the NA section interp/impute/etc.
##
source('/Users/abarbour/kook.processing/R/dev/timetasks/merge/funcs.R')
##
dtStrp <- function(x){
  # "POSIXct" no good for colClasses, so manually strptime
  # (ensures fractional seconds are accounted for
  strptime(x, format="%Y-%m-%dT%H:%M:%OS", tz="UTC")
}
##
readBSM <- function(fi, nalevel=1e4){
  # read in a ascii table (gethf_20sps.csh, bottle_merge.py, bottle.py)
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
  # first convert the timestamp, making sure fractional seconds
  # are accounted for
  dat$dt <- dtStrp(dat$dt)
  # (from funcs.R) deal with regions filled with NAs (interp, impute, etc.)
  dat$CH0 <- naOP(dat$CH0, quiet=FALSE)
  dat$CH1 <- naOP(dat$CH1, quiet=FALSE)
  dat$CH2 <- naOP(dat$CH2, quiet=FALSE)
  dat$CH3 <- naOP(dat$CH3, quiet=FALSE)
  return(invisible(dat))
}
##
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

