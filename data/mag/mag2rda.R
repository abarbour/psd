owd<-getwd()
setwd("data")
mag.dir <- "."
mag <- cbind(
	read.table(paste(mag.dir,"mag.raw",sep="/"),
		colClasses="numeric", col.names="raw"),
	read.table(paste(mag.dir,"mag.clean",sep="/"),
		colClasses="numeric", col.names="clean")
	)
plot(mag$raw, type="s", main="Magnetometer data")
text(200, 40, "raw series")
lines(mag$clean-50, type="s", col="red")
text(350, -100, "cleaned series (-50)", col="red")
save(mag, file = "mag.rda")
###
setwd(owd)
getwd()
###
