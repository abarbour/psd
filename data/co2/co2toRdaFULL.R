##
##  Builtin only has until 1997, this is the full monthly record
##
# read in url of co2 data (regularly updated, full series)
co2raw<-read.table(
  url("ftp://ftp.cmdl.noaa.gov/ccg/co2/trends/co2_mm_mlo.txt"),
	h=F, na.string="-99.99")
names(co2raw) <- c('yr',"m","yf","co2","co2interp","trnd","days")
co2ts <- ts(co2raw$co2interp, start=c(1958,3), freq=12)
save(co2ts, file="data/co2/co2ts_full.rda")
co2tsd <- decompose(co2ts)
save(co2tsd, file="data/co2/co2tsd_full.rda")

library(zoo)

psc<-pspectrum(na.approx(co2tsd$x),fsamp=12)
psct<-pspectrum(na.approx(co2tsd$trend),fsamp=12)
pscs<-pspectrum(na.approx(co2tsd$seasonal),fsamp=12)
pscr<-pspectrum(na.approx(co2tsd$random),fsamp=12)

# plot(co2tsd)
# plot(pscs$f,log10(pscs$psd), type="l", ylim=7*c(-1.1,0.5), xlim=c(0,5.5), xaxs="i")
# lines(pscr$f,log10(pscr$psd),col="red")
# lines(psct$f,log10(psct$psd),col="blue")
# lines(psc$f,log10(psc$psd),col="green")

rpsc<-spectrum(na.approx(co2tsd$x),pad=1,taper=0.2,plot=F)
rpsct<-spectrum(na.approx(co2tsd$trend),pad=1,taper=0.2,plot=F)
rpscs<-spectrum(na.approx(co2tsd$seasonal),pad=1,taper=0.2,plot=F)
rpscr<-spectrum(na.approx(co2tsd$random),pad=1,taper=0.2,plot=F)

thonig <- read.table("data/co2/thonig_spec_1989.txt",h=T)

library(ggplot2)
library(RColorBrewer)
pal <- brewer.pal(n=8,"Paired")

# spec.pgram spectra
df <- data.frame(f=rpsc$freq, psd=rpsc$spec, src="spec.pgram", decomp="Full")
df <- rbind(df, data.frame(f=rpsct$freq, psd=rpsct$spec, src="spec.pgram", decomp="Trend"))
df <- rbind(df, data.frame(f=rpscs$freq, psd=rpscs$spec, src="spec.pgram", decomp="Seasonal"))
df <- rbind(df, data.frame(f=rpscr$freq, psd=rpscr$spec, src="spec.pgram", decomp="Residual"))
# now add the MT spectra
df <- rbind(df, data.frame(f=psc$f, psd=psc$psd, src="rlpSpec", decomp="Full"))
df <- rbind(df, data.frame(f=psct$f, psd=psct$psd, src="rlpSpec", decomp="Trend"))
df <- rbind(df, data.frame(f=pscs$f, psd=pscs$psd, src="rlpSpec", decomp="Seasonal"))
df <- rbind(df, data.frame(f=pscr$f, psd=pscr$psd, src="rlpSpec", decomp="Residual"))
# and Thonig spectra for reference
df <- rbind(df, data.frame(f=thonig$freq_cpy, psd=thonig$spec_ppmco2, src="", decomp="Full"))
df <- rbind(df, data.frame(f=thonig$freq_cpy, psd=thonig$spec_ppmco2, src="", decomp="Trend"))
df <- rbind(df, data.frame(f=thonig$freq_cpy, psd=thonig$spec_ppmco2, src="", decomp="Seasonal"))
df <- rbind(df, data.frame(f=thonig$freq_cpy, psd=thonig$spec_ppmco2, src="", decomp="Residual"))

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
ggsave("figures/co2.pdf", height=3.2, width=10)




