rm(list=ls())

library(ellipse)
library(psd)
library(reshape2)
library(ggplot2)

bsmnm <- read.table("bsmnm.txt",header=TRUE)

opar <- par(no.readonly = TRUE)
PSD <- function(x) psdcore(x, 1, 20)

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
#
#
STAFUN <- function(sta){
  fis <- paste(sta,"ppar","long",sep=".")
  dat <- read.table(fis, header=FALSE)
  dat <- transform(dat, V1=1*V1, V2=100*V2)
  dE <- dat$V1
  dP <- dat$V2
  Honshu <- "2011-03-11 05:46:24" #UTC
  Honshu.ct <- as.numeric(as.POSIXct(strptime(Honshu, tz="UTC", "%Y-%m-%d %H:%M:%S")))
  Ets <- ts(dE, start = Honshu.ct)
  Pts <- ts(dP, start = Honshu.ct) #+44)
  plot(ts.union(Ets,Pts), main=sta)
  plot(dP ~ dE, main=sta)
  #
  names(dat) <- c("areal","pressure")
  attr(dat, "station") <- sta
  dat$trnd <- seq_len(nrow(dat))
  EPlm <- lm(pressure + 1 ~ areal + trnd + 1, dat)
  #plot(EPlm)
  EPpsd <- Epsd <- PSD(dE)
  Ppsd <- PSD(dP)
  tmpdf <- Stas[Stas==sta,]
  Ppsd$spec <- Ppsd$spec/tmpdf$ppsd.scale
  EPpsd$spec <- EPpsd$spec/Ppsd$spec
  plot(Epsd, log="dB", ylim=c(-210,-110),
	  col="red", main=sprintf("%s PP_sc (%g, blue) and BSM (red)", sta, 1/sqrt(tmpdf$ppsd.scale)))
  plot(Ppsd, log="dB", add=TRUE, col="blue")
  #
  lines(P50 ~ freq, bsmnm, type="s")
  lines(P10 ~ freq, bsmnm, lty=3, type="s")
  #
  #
  plot(EPpsd, log="dB", 
	  main=sprintf("%s PP_sc/BSM strain/%g Pa)", sta, 1/sqrt(tmpdf$ppsd.scale)),
	  ylim=60*c(-1,1))
  Rpsd <- PSD(residuals(EPlm))
  plot(Rpsd, log="dB", add=TRUE, col="grey")
  #
  return(EPlm)
}

ELLPLT <- function(sta,mods){
  lmod <- get(sta, mods)
  ell <- ellipse(lmod, which=c('areal','trnd'), level=0.95)
  plot(ell, type='l',col="red", main=sta)
  points(lmod$coefficients['areal'], lmod$coefficients['trnd'])
  lines(ellipse(lmod, which=c('areal','trnd'), level=0.90), col="red",lty=3)
  return(ell[,"areal"])
}

pdf("honshu.pdf")
print(res1 <- sapply(Stas$stas, FUN=STAFUN, simplify=FALSE))
par(mfrow=c(2,2))
print(res2 <- sapply(Stas$stas, FUN=ELLPLT, mods=res1, simplify=TRUE))

res2m <- melt(res2)
res2m$value <- log10(abs(res2m$value))
res2m <- transform(res2m, Var2=reorder(Var2, value))
g <- ggplot(res2m, aes(y=value, x=Var2, group=Var2))
g+stat_summary(fun.y = mean, fun.ymin = min, fun.ymax = max,
	geom="pointrange",aes(colour=Var2))+coord_flip()

dev.off()
#par(opar)
