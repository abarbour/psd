##
##  Replicates figure in text
##

load("../data/co2/co2tsd.rda")

library(zoo)

# Adaptive multitaper estimation: rlpSpec::pspectrum
frq<-12
psc <-pspectrum(na.approx(co2tsd$x),        fsamp=frq, plotpsd=F)
psct<-pspectrum(na.approx(co2tsd$trend),    fsamp=frq, plotpsd=F)
pscs<-pspectrum(na.approx(co2tsd$seasonal), fsamp=frq, plotpsd=F)
pscr<-pspectrum(na.approx(co2tsd$random),   fsamp=frq, plotpsd=F)

# R built-in spectrum estimation: spec.pgram
pad<-1
tap<-0.2
rpsc <-spectrum(na.approx(co2tsd$x),        pad=pad, taper=tap, plot=F)
rpsct<-spectrum(na.approx(co2tsd$trend),    pad=pad, taper=tap, plot=F)
rpscs<-spectrum(na.approx(co2tsd$seasonal), pad=pad, taper=tap, plot=F)
rpscr<-spectrum(na.approx(co2tsd$random),   pad=pad, taper=tap, plot=F)

# Previously publish spectra (Thoning et al 1989)
thonig <- read.table("../data/co2/thoning_spec_1989.txt",h=T)

# create a composite data frame with:
# spec.pgram spectra
df <- data.frame(f=rpsc$freq, psd=rpsc$spec, src="spec.pgram", decomp="Full")
df <- rbind(df, data.frame(f=rpsct$freq, psd=rpsct$spec, src="spec.pgram", decomp="Trend"))
df <- rbind(df, data.frame(f=rpscs$freq, psd=rpscs$spec, src="spec.pgram", decomp="Seasonal"))
df <- rbind(df, data.frame(f=rpscr$freq, psd=rpscr$spec, src="spec.pgram", decomp="Residual"))
# multitaper spectra
df <- rbind(df, data.frame(f=psc$f, psd=psc$psd, src="rlpSpec", decomp="Full"))
df <- rbind(df, data.frame(f=psct$f, psd=psct$psd, src="rlpSpec", decomp="Trend"))
df <- rbind(df, data.frame(f=pscs$f, psd=pscs$psd, src="rlpSpec", decomp="Seasonal"))
df <- rbind(df, data.frame(f=pscr$f, psd=pscr$psd, src="rlpSpec", decomp="Residual"))
# and Thonig spectra for reference 
df <- rbind(df, data.frame(f=thonig$freq_cpy, psd=thonig$spec_ppmco2, src="", decomp="Full"))
df <- rbind(df, data.frame(f=thonig$freq_cpy, psd=thonig$spec_ppmco2, src="", decomp="Trend"))
df <- rbind(df, data.frame(f=thonig$freq_cpy, psd=thonig$spec_ppmco2, src="", decomp="Seasonal"))
df <- rbind(df, data.frame(f=thonig$freq_cpy, psd=thonig$spec_ppmco2, src="", decomp="Residual"))

## Now start plotting

library(ggplot2)
library(RColorBrewer)
pal <- brewer.pal(n=8,"Paired")

g <- ggplot(df, aes(x=f, y=20*log10(psd), colour=src, size=src))
p <- g+geom_line()+facet_grid(.~decomp)
(p + 
  scale_x_continuous("Frequency, Cycles per Year", limits=c(0,5.5), expand=c(0,0))+
  scale_y_continuous("CO2 Power, dB rel ppm^2 * year", 
                     breaks=seq(-180,30,by=30),limits=c(-200,44), expand=c(0,0))+
   scale_colour_manual("Spectrum\nestimator", values=c(pal[1:2],"black"), 
                       breaks=c("spec.pgram","rlpSpec",""),
                       labels=c("spec.pgram","rlpSpec","Thoning et al\n(1976-1985)"))+
   scale_size_manual(values=c(0.9,1.0,0.5), legend=F)+
   #+ scale_linetype_manual(values=factor(c(0,0,1)), legend=F)
   opts(title="Mauna Loa CO2 Concentration Power Spectra: 1959 - 2012")
  #+ theme_bw()
 )
ggsave("./co2.pdf", height=3.2, width=10)