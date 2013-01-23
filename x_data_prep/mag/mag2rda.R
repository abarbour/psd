(owd<-getwd())
mag.dir <- "."

## Raw intensity
mr <- read.table(paste(mag.dir,"mag.raw",sep="/"), colClasses="numeric", col.names="raw")
## Edited intensity
mc <- read.table(paste(mag.dir,"mag.clean",sep="/"), colClasses="numeric", col.names="clean")
## Difference
mdiff <- zapsmall(mc$clean - mr$raw)
## Distance
km <- 0:(length(mr)-1)

## Data.frame
mag <- data.frame(km=km, raw=mr, clean=mc, mdiff=mdiff)

## Edits (include?)
#mag$edit <- FALSE
#mag$edit[abs(mag$diff) > 0] <- TRUE

## Rdata
magsat <- mag
save(magsat, file = "magsat.rda")
##
setwd(owd)
getwd()
###
