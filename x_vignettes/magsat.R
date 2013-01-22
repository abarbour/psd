rm(list=ls())
##
##load("mag.rda")
data(magsat)
stopifnot(exists("magsat"))
##
subset(magsat, abs(mdiff)>0)
#    km    raw   clean     mdiff
#    403  0  209.1 -3.6355 -212.7355
#    717  0 -248.7 -9.7775  238.9225
##
plt_Mag01 <- function(){
  plot(raw ~ km, mag, type="s",
       main="MAGSAT Airborne-magnetometer data",
       xaxs="i", 
       xlab="Along-path distance, km", 
       yaxs="i", ylim=c(-300,225), 
       ylab="Field strength, nT",
       lwd=1)
  mtext("Horizontal-component field")
  text(200, 100, "raw series")
  lines(clean-50 ~ km, mag, type="s", col="red", lwd=1)
  text(350, -150, "cleaned series (-50)", col="red")
  lines(40*edit-320 ~ km, mag, type="s", col="blue", lwd=2)
}
plt_Mag01()
##
pdf("./magsat.pdf",height=3.3)
plt_Mag01()
dev.off()
##
