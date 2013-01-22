##
##  Builtin only has until 1997, but this is the full monthly record
##
## read in url of co2 data (regularly updated, full series)
co2raw<-read.table(
  url("ftp://ftp.cmdl.noaa.gov/ccg/co2/trends/co2_mm_mlo.txt"),
	h=F, na.string="-99.99")
names(co2raw) <- c('yr',"m","yf","co2","co2interp","trnd","days")
frq <- 12
## convert to timeseries
co2ts <- ts(co2raw$co2interp, start=c(1958,3), freq=frq)
## save raw (unnecessary really)
save(co2ts, file="./co2ts_full.rda")
## timeseries decomposition
co2tsd <- decompose(co2ts)
save(co2tsd, file="./co2tsd_full.rda")