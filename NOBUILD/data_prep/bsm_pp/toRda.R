rm(list=ls())

library(ellipse)
library(psd)

opar <- par(no.readonly = TRUE)
PSD <- function(x) psdcore(x, 1, 10)

stanames <- paste0("B0",c(81,87,84,88))
Stas <- data.frame(stas=as.character(stanames), ppsd.scale=(1e9*c(.05,.25,1,10))**2)

Honshu <- "2011-03-11 05:46:24" #UTC
Honshu.ct <- as.POSIXct(strptime(Honshu, tz="UTC", "%Y-%m-%d %H:%M:%S"))
#
POSDF <- function(x,
	name=deparse(substitute(x)),
	nterms=length(x),
	sec.samp=1, sec.forward=0, start.time=Honshu.ct){
	#Time sequence
	tseq <- seq(from=Honshu.ct-sec.forward, by=sec.samp, length.out=nterms)
	posdf <- tseq
	return(posdf)
}
for (sta in Stas$stas){
  fis <- paste(sta,"ppar","long",sep=".")
  dat <- read.table(fis, header=FALSE)
  dat <- transform(dat, V1=1*V1, V2=100*V2)
  dE <- dat$V1
  dP <- dat$V2
  Honshu <- "2011-03-11 05:46:24" #UTC
  Honshu.ct <- as.numeric(as.POSIXct(strptime(Honshu, tz="UTC", "%Y-%m-%d %H:%M:%S")))
  Ets <- ts(dE, start = Honshu.ct)
  Pts <- ts(dP, start = Honshu.ct) #+44)
  #
  names(dat) <- c("areal","pressure")
  attr(dat, "station") <- sta
  dat$trnd <- seq_len(nrow(dat))
  EPlm <- lm(pressure ~ areal + trnd, dat)
  plot(dat[,1:2], main=sta)
  plot(ellipse(EPlm, which=c('areal','trnd'), level=0.95), type='l', col="red")
  points(EPlm$coefficients['areal'], EPlm$coefficients['trnd'])
  lines(ellipse(EPlm, which=c('areal','trnd'), level=0.90), col="red",lty=3)
  EPpsd <- Epsd <- PSD(dE)
  Ppsd <- PSD(dP)
  tmpdf <- Stas[Stas==sta,]
  Ppsd$spec <- Ppsd$spec/tmpdf$ppsd.scale
  EPpsd$spec <- EPpsd$spec/Ppsd$spec
  plot(EPpsd, log="dB", main=sprintf("%s PP_sc/BSM strain/%g Pa)", sta, 1/sqrt(tmpdf$ppsd.scale)))
  plot(Epsd, log="dB", ylim=c(-220,-110),
	  col="red", main=sprintf("%s PP_sc (%g, blue) and BSM (red)", sta, 1/sqrt(tmpdf$ppsd.scale)))
  plot(Ppsd, log="dB", add=TRUE, col="blue")
  #
  plot(ts.union(Ets,Pts), main=sta)
}

par(opar)
