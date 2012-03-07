##
##  Convert Scarborough data to R format
##
setwd("~/kook.processing/R/dev/packages/rlpSpec/data/scarb")
scarb <- read.table("chan5.dat", h=F, nrow=1e5, col.names="Ez")
delt = 0.016  # sec
sps = 1/delt	# [62.5] Hz
#The record was begun on May 27, 2009 at 6:40 (GMT or local?)
tsSt<-as.POSIXct(origin=0, x=strptime('2009-05-27 06:40:00', '%Y-%m-%d %H:%M:%S'), tz="GMT")
scarbts <- ts(data=scarb$Ez, start=0, frequency=sps)
scarbl <- list(ts=scarbts, tsst=tsSt, sps=sps, dt=delt, component="Ez")
# write file
save(scarbl, file="./scarb.Rda")
##