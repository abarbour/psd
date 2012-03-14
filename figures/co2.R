##
##  Creates figure(s) in text for CO2 data
##
setwd("~/kook.processing/R/dev/packages/rlpSpec/figures")
source('funcload.R')
load("../data/co2/co2tsd.rda")
library(zoo)

# Adaptive multitaper estimation: rlpSpec::pspectrum
frq<-12
deseas <- na.approx(co2tsd$x) - na.approx(co2tsd$seasonal)
psc <-pspectrum(na.approx(co2tsd$x), fsamp=frq, plotpsd=F, niter=4)
pscs<-pspectrum(deseas, fsamp=frq, plotpsd=F, niter=4)

# R built-in spectrum estimation: spec.pgram
pad<-1
tap<-0.2
rpsc <-spectrum(na.approx(co2tsd$x), pad=pad, taper=tap, plot=F)
rpscs<-spectrum(deseas, pad=pad, taper=tap, plot=F)

# Previously publish spectra (Thoning et al 1989)
thonig <- read.table("../data/co2/thoning_spec_1989.txt",h=T)

# create a composite data frame with:
# spec.pgram spectra
df <- data.frame(f=rpsc$freq, psd=rpsc$spec, src="spec.pgram", decomp="Full")
df <- rbind(df, data.frame(f=rpscs$freq, psd=rpscs$spec, src="spec.pgram", decomp="Seasonally Adjusted"))
# multitaper spectra
df <- rbind(df, data.frame(f=psc$f, psd=psc$psd, src="rlpSpec", decomp="Full"))
df <- rbind(df, data.frame(f=pscs$f, psd=pscs$psd, src="rlpSpec", decomp="Seasonally Adjusted"))
# and Thonig spectra for reference 
df <- rbind(df, data.frame(f=thonig$freq_cpy, psd=thonig$spec_ppmco2, src="", decomp="Full"))
df <- rbind(df, data.frame(f=thonig$freq_cpy, psd=thonig$spec_ppmco2, src="", decomp="Seasonally Adjusted"))

## Now start plotting

library(ggplot2)
library(RColorBrewer)
pal <- brewer.pal(n=8,"Paired")

g <- ggplot(df, aes(x=f, y=10*log10(psd), colour=src, size=src))
p <- g+geom_line()+facet_grid(.~decomp)
(p + 
  scale_x_continuous("Frequency, Cycles per Year", limits=c(0,5.5), expand=c(0,0))+
  scale_y_continuous("Power, dB (ppm^2 * year)",
                     breaks=seq(-50,22,by=10), limits=c(-50,22), expand=c(0,0))+
   scale_colour_manual("Spectrum\nestimator", values=c(pal[1:2],"black"), 
                       breaks=c("spec.pgram","rlpSpec",""),
                       labels=c("spec.pgram\n","rlpSpec\n","Thoning et al\n(1976-1985)"))+
   scale_size_manual(values=c(0.9,1.0,0.5), guide="none")+
   #+ scale_linetype_manual(values=factor(c(0,0,1)), legend=F)
   opts(title="Mauna Loa CO2 concentration: Power spectra, 1959 - 2012") + theme_bw()
 )
##
ggsave("./co2.pdf", height=3.2, width=7)
##
g <- ggplot(df, aes(x=log2(f), y=10*log10(psd), colour=src, size=src))
p <- g+geom_line()+facet_grid(.~decomp)
(p + 
  scale_x_continuous("Frequency, Cycles-per-Year", 
                     limits=log2(c(.3,5.5)), 
                     breaks=-1:2, labels=c("1/2",1,2,4),
                     expand=c(0,0))+
  scale_y_continuous("Power, dB (ppm^2 * year)",
                     breaks=seq(-50,22,by=10), limits=c(-50,22), expand=c(0,0))+
  scale_colour_manual("Spectrum\nestimator", values=c(pal[1:2],"black"), 
                      breaks=c("spec.pgram","rlpSpec",""),
                      labels=c("spec.pgram\n","rlpSpec\n","Thoning et al\n(1976-1985)"))+
  scale_size_manual(values=c(0.9,1.0,0.5), guide="none")+
 #+ scale_linetype_manual(values=factor(c(0,0,1)), legend=F)
 opts(title="Mauna Loa monthly CO2: Power spectra, 1959 - 2012") + theme_bw()
 )
##
ggsave("./co2log.pdf", height=3.2, width=7)
##