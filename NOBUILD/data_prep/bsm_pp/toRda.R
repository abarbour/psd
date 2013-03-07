rm(list=ls())

library(ellipse)
library(psd)
library(reshape2)
library(ggplot2)

bsmnm <- read.table("bsmnm.txt",header=TRUE)

opar <- par(no.readonly = TRUE)

## PSD calc
PSD <- function(x) psdcore(x, 1, 30)
PSDP <- function(x) pilot_spec(x, 1, 15)
PSDA <- function(x) pspectrum(x, niter=1)

stanames <- paste0("B0",c(81,87,84,88))
flt <- read.csv("/Users/abarbour/kook.processing/research/BSMCoseismic/OffsetEstimation/studies/bsm/fltgeods.min.txt")
Stas <- data.frame(sta=as.character(stanames),
		   ppsd.scale=(1e9*c(.05,.25,1,10))**2,
		   pspecadd=c(F,T,T,T))
Stas <- merge(Stas, flt, by="sta")
Stas$cols <- as.numeric(factor(Stas$sta))


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
PPPSD <- function(sta){
  fis <- paste(sta,"ppar","long",sep=".")
  dat <- read.table(fis, header=FALSE)
  tmpdf <- Stas[Stas==sta,]
  padd <- tmpdf$pspecadd
  pcol <- tmpdf$cols
  # Pressure in hPa, strain in strain
  dat <- transform(dat, V1=1*V1, V2=100*V2)
  #dE <- dat$V1
  dP <- dat$V2
  Ppsd <- PSDP(dP)
  Ppsd$freq <- log10(Ppsd$freq)
  par(las=1)
  plot(Ppsd, log="dB", add=padd, col=pcol, ylim=c(-10,80), yaxs="i")
  text(.35,18+5*pcol,sta, col=pcol)
  return(Ppsd)
}
#
STAFUN <- function(sta){
  fis <- paste(sta,"ppar","long",sep=".")
  dat <- read.table(fis, header=FALSE)
  # Pressure in hPa, strain in strain
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

pdf("ppspec.pdf")
resp <- sapply(Stas$sta, FUN=PPPSD, simplify=FALSE)
abline(h=0,lty=3)
dev.off()
pdf("honshu.pdf")
print(res1 <- sapply(Stas$sta, FUN=STAFUN, simplify=FALSE))
par(mfrow=c(2,2))
print(res2 <- sapply(Stas$sta, FUN=ELLPLT, mods=res1, simplify=TRUE))
dev.off()

res2m <- melt(res2)
names(res2m) <- c("Pos","sta","value")
res2m <- merge(res2m, flt, by="sta")
res2m$value <- log10(abs(res2m$value))
res2m <- transform(res2m, sta=reorder(sta, geodkm))

#g <- ggplot(res2m, aes(y=value, x=geodkm, group=sta))
#g <- ggplot(res2m, aes(y=value, x=1/sqrt(geodkm), group=sta))
#g+stat_summary(fun.y = mean, fun.ymin = min, fun.ymax = max, geom="pointrange",aes(colour=sta))
VP <- function(bulkmod, nu=0.3, dens=2650){
	rhoVp2 <- 3*bulkmod*(1-nu)/(1+nu)
	Vp2 <- rhoVp2/dens
	return(sqrt(Vp2))
}
res2m$Vp <- VP(10**res2m$value)

g <- ggplot(res2m, aes(y=value, x=sqrt(geodkm)))+geom_jitter(alpha=0.1)
p <- g + xlim(1,4) + ylim(6.5,10) +
  xlab("sqrt( SJF-normal distance, km )")+ ylab("Effective bulk modulus, log10 Pa")+
  geom_text(data=Stas,aes(label=sta, y=7),size=3,angle=90)+
  scale_colour_discrete(guide="none")
ggsave("bulkmod.pdf", p)

#g <- ggplot(res2m, aes(y=2*Vp/1e3, x=(geodkm)))+geom_jitter(alpha=0.1)
#p <- g + 
#  # xlim(1,4) + 
#  # ylim(6.5,10) +
#  xlab("sqrt( SJF-normal distance, km )")+ ylab("Effective Vp at nu=0.25, rho=2.6; km/s")+
#  #geom_text(data=Stas,aes(label=sta, y=7),size=3,angle=90)+
#  scale_colour_discrete(guide="none")
#ggsave("vp.pdf", p)

