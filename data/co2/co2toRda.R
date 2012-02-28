##
##  Builtin only has until 1997, this is through 2003
##
library(zoo)
# read in url of (latest?) co2 data
co2raw<-read.table(
	url("ftp://cdiac.esd.ornl.gov/pub/maunaloa-co2/maunaloa.co2"),
	h=T,nrows=46,comment.char="*", na.string="-99.99")
names(co2raw) <- c('yr',1:12,"an","anfit")
# reshape and exclude annuals
co2m <- matrix(t(co2raw[2:13]), ncol=1)
# interpolate, and create timeseries
co2ts <- ts(na.approx(zoo(co2m)),start=c(1958,3), freq=12)
co2tsd <- decompose(co2ts)
save(co2tsd, file="co2tsd.rda")

psc<-pspectrum(na.approx(co2tsd$x),fsamp=12)
psct<-pspectrum(na.approx(co2tsd$trend),fsamp=12)
pscs<-pspectrum(na.approx(co2tsd$seasonal),fsamp=12)
pscr<-pspectrum(na.approx(co2tsd$random),fsamp=12)

plot(pscs$f,log10(pscs$psd), type="l", ylim=7*c(-1,1))
lines(pscr$f,log10(pscr$psd),col="red")
lines(psct$f,log10(psct$psd),col="blue")
lines(psc$f,log10(psc$psd),col="green")
