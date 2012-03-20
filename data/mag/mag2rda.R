(owd<-getwd())
#setwd("data/mag")
mag.dir <- "."
mag <- cbind(
	read.table(paste(mag.dir,"mag.raw",sep="/"),
		colClasses="numeric", col.names="raw"),
	read.table(paste(mag.dir,"mag.clean",sep="/"),
		colClasses="numeric", col.names="clean")
	)
mag$diff <- mag$clean - mag$raw
mag$diff[abs(mag$diff) < 1e-6 ] <- 0
###
save(mag, file = "mag.rda")
###
setwd(owd)
getwd()
###
