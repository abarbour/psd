rm(list=ls())
#library(ellipse)
library(zoo)
library(reshape2)
library(ggplot2)
TODT <- function(dts){strptime(dts, tz="UTC", format="%Y-%m-%dT%H:%M:%OS")}
#
#
# Data appear to be nearly 10 days of 1Hz data, which gives nearly 3/4 million
# points.  This can take a while to load.
# gzcat Tohoku/B084_B/B084.txt.gz | awk '{print NF}' | uniq -c 
# 777602 13
#
#
#http://stackoverflow.com/questions/4756989/how-to-load-data-quickly-into-r
# subsequently:
# http://cran.r-project.org/web/packages/saves/index.html
#
#        Column 1        DateTime                UTC
#1        Column 2        Eee+Enn                 microstrain
# 2       Column 3        Eee+Enn_tide            microstrain
#  3      Column 4        Eee+Enn_baro            microstrain
#   4     Column 5        Eee-Enn                 microstrain
#    5    Column 6        Eee-Enn_tide            microstrain
#     6   Column 7        Eee-Enn_baro            microstrain
#      7  Column 8        2Ene                    microstrain
#       8 Column 9        2Ene_tide               microstrain
#9        Column 10       2Ene_baro               microstrain
# 10       Column 11       BaroPressure            millibar
#   11     Column 12       PorePressure            millibar
#        Column 13       Version                 generation date
#        Column 14       TiltX                   microradian
#        Column 15       TiltY                   microradian
#
# DateTime(UTC)
#1  Eee+Enn(ms)
#2  Eee+Enn_tide(ms)
#3  Eee+Enn_baro(ms)
#4  Eee-Enn(ms)
#5  Eee-Enn_tide(ms)
#6  Eee-Enn_baro(ms)
#7  2Ene(ms)
#8  2Ene_tide(ms)
#9  2Ene_baro(ms)
#10 BaroPressure(mbar)
#11 PorePressure(mbar)
#Version
#
#2011-03-11T00:00:00 <- start
#2011-03-11T05:46:24 <- origin of Mw9 Tohoku-oki Eq, 2011
#http://earthquake.usgs.gov/earthquakes/eqinthenews/2011/usc0001xgp/
#
fi <- "Tohoku/B084_B/B084.txt.gz"
fi <- "Tohoku/B084_B/B084.tensor.LAB.20110311_Tohoku.txt.gz"
eqorigin <- TODT("2011-03-11T05:46:24 UTC")
# P and S travel times from the iasp 91 tau-p model
taup <- c(712.520, 1301.851)
# backaz, dist, depth, in km
# Honshu back azim is -52.4
eqd <- c(backaz=-52.4, geod=8613, depth=30) 
#
#
# read, and exclude version field
system.time(dat <- read.table(fi, 
		nrows=16e3,
		skip=14e3,
		header=TRUE,
		colClasses=c("character",rep("numeric",11),"NULL") ,
		col.names=c("Dts",
			"areal",
			"areal_interp",
			"areal.tide",
			"areal.baro",
			"gamma1",
			"gamma1_interp",
			"gamma1.tide",
			"gamma1.baro",
			"gamma2",
			"gamma2_interp",
			"gamma2.tide",
			"gamma2.baro",
			"pressure.atm",
			"pressure.pore",
			"version")
		))
# NA values
dat[dat==999999] <- NA
# Set some attributes:
#	units:
attr(dat,"units") <- list(tzone="UTC", 
			  strain="microstrain",
			  pressure="mBar",
			  distance="km",
			  traveltime="seconds")
#	earthquake info:
attr(dat,"iasp") <- list(geodetic=eqd,
			 model="91 tau-p model",
			 tau_P_S=taup,
			 origin=eqorigin)
summary(dat)

# Add some helpful information:
dat$Dt <- TODT(dat$Dts)
dat$Origin.secs <- as.numeric(dat$Dt - eqorigin)
dat$epoch <- "preseismic"
dat$epoch[dat$Origin.secs>=(taup[1]-1)] <- "seismic"
summary(dat)

# Fill NA forward/backwards to create consistent file
dat[,2:12] <- na.locf(na.locf(zoo(dat[,2:12]), na.rm=FALSE), fromLast=TRUE, na.rm=FALSE)
summary(dat)

# Rename, and save locally
Tohoku <- dat
save(Tohoku, file="Tohoku.rda", compress = 'xz')
rm(Tohoku)

# For plotting, melt and facet
datm <- melt(dat[,c("areal","gamma1","gamma2","pressure.pore","Origin.secs","epoch")],
	     id.vars=c("Origin.secs","epoch"))
p <- ggplot(datm, aes(x=Origin.secs,y=value,group=variable)) +
     geom_line() +
     facet_grid(variable~epoch,scales="free") +
     geom_vline(xintercept=taup, linetype="dotted") +
     theme_bw()
print(p)
ggsave("Tohoku.pdf",p)
