#vertical axis rotation of [areal gamma1 gamm2] referenced to X degrees (compass
#dir):
#
#       [areal (invariant), c*gamma1 + s*gamma2, -s*gamma1 + c*gamma2]
#where
#       c = cos(2*rotation)
#       s = sin(2*rotation)
library(psd)
data(Tohoku) #load("Tohoku.rda")
library(signal)
#
# Back azimuth, in radians
bazd <- attributes(Tohoku)$iasp$geodetic["backaz"][[1]]
baz <- bazd * pi/180
#
# shears components
#areal <- as.matrix(Tohoku[,c("areal")])
#shears <- as.matrix(Tohoku[,c("gamma1","gamma2")])
# remove mean values
#areal <- apply(areal, FUN=function(x)x-mean(x), 2)
#shears <- apply(shears, FUN=function(x)x-mean(x), 2)
#
bf <- butter(5,0.001,"high")  # 5th order butterworth filter 0.001*Nyquist
PROCDAT <- function(dat, fw=bf){
  require(signal)
  return(filtfilt(fw, dat))
}
DEMEAN <- function(mat) apply(mat, FUN=function(x)x-mean(x), 2)
Porep <- DEMEAN(as.matrix(Tohoku[,c("pressure.pore")]))
Areal <- DEMEAN(as.matrix(Tohoku[,c("areal")]))
Shears <- DEMEAN(as.matrix(Tohoku[,c("gamma1","gamma2")]))
#Shears[,2] <- Shears[,2]/2
#
SHEARROT <- function(backaz, shears){
  # 1) put shears into a complex vector
  shears.c <- complex(real=shears[,"gamma1"], imaginary=shears[,"gamma2"])
  # Positive increments in argument correspont to counter-clockwise
  # rotations; hence, rotation into backaz (cw in North-0) should
  # be -1*backaz
  # 2) modify argument by adding rotation
  backaz <- backaz*2
  shears.r <- complex(modulus = Mod(shears.c), argument = Arg(shears.c) - backaz)
  return(shears.r)
}
AREALS <- function(backaz.deg, areal, shears.){
  message(backaz.deg-bazd)
  backaz. <- backaz.deg * pi / 180
  shears.r <- SHEARROT(backaz., shears.)
  # 3) rotate shears in radial/transverse extension
  Eh <- Re(shears.r) + areal
  Er <- Im(shears.r)# + areal
  # 4) regress equivalent areal strains
  #areal_r <- Eh + Er
  #Efit <- lm(areal_r ~ areal)
  # intercept is implicit, and the first coeff
  # second coeff is areal_r to areal fit
  #arc <- coefficients(Efit)[[2]]
  #resnorm <- norm(res <- as.matrix(residuals(Efit)))
  #
  Efit <- lm(Eh ~ Er)
  fc <- coefficients(Efit)
  afc <- fc[[2]]
  message(afc)
  #  res <- (Eh+Er)
  resnorm <- norm(res <- as.matrix(residuals(Efit)))
  #  resnorm
  #  afc - 1
  #1-log10(afc)
  #1-log10(resnorm)
  sqrt(var(res))
}
FUN <- function(az) AREALS(az, areal=Areal, shears.=Shears)
#Opt <- optim(baz*180/pi, FUN, gr=NULL, method="L-BFGS-B", lower=-360, upper=360)
#Opt <- optim(bazd, FUN, gr=NULL, method="L-BFGS-B", lower=-90, upper=90)
bazd <- 360+bazd
#bazd <- 307
OPT <- function(){
  Opt <- optim(bazd, FUN, gr=NULL, method="L-BFGS-B", lower=0, upper=360)
  ##Opt <- optim(0, FUN, gr=NULL, method="L-BFGS-B", lower=-90, upper=90)
  bazopt <- Opt$par[1]
}
if (!exists("bazopt")) bazopt <- OPT()

shears.r <- SHEARROT(bazd*pi/180, Shears) # does not minimize transverse extension
shears.r <- SHEARROT(bazopt*pi/180, Shears)
Eh <- Re(shears.r) + Areal/2
Er <- Im(shears.r) + Areal/2
shears.sjf <- SHEARROT((225)*pi/180, Shears)
g1.fp <- Re(shears.sjf)
g2.fp <- Im(shears.sjf)

pdf("Tohoku_rotated.pdf",h=11,w=8.5)
#par(mfcol=c(2,1), pty="m")
par(pty="m", oma=rep(1,4))

main <- sprintf("B084-Tohoku strains")
delts <- 0.55*rev(seq_len(11)-4)
plot(g1.fp+delts[1], type="l", lty=1, ylim=range(delts)*1.1,
	ylab="microstrain", xlab="arbitrary time, seconds",
	main=main, col="dark red")
mtext(sprintf("g1, g2 rotated by %.01f to min areal-misfit (%.01f rel eq backaz)",bazopt,bazopt-bazd))
lines(g2.fp+delts[2], col="dark red")
#lines(log10(abs(Er/Eh))+delts[3], col="green")
lines(Eh+delts[3])
lines(Er+delts[4])
Areal_r <- Eh + Er
Efit <- lm(Areal_r ~ Areal)
lines(Areal_r+delts[5])
lines(Areal+delts[6], col="dark red")
lines(residuals(Efit)+delts[7], col="dark grey")
Pfit <- lm(Porep ~ Areal)
psc <- coefficients(Pfit)[2]
#lines(Porep/psc+delts[8], col="light blue")
Pfit_r <- lm(Porep ~ Areal_r)
psc_r <- coefficients(Pfit_r)[2]
#lines(Porep/psc_r+delts[8], col="blue")

# Volume strain
#vertical strain for Poissons ratio nu
#nu <- 0.25
#Ez <- -1*nu*as.matrix(uEz)/(1-nu)
#Ekk <- Areal+Ez
#uEz <- residuals(Efit)
Pfit <- lm(Porep ~ Areal)
sc <- coefficients(Pfit)
psc <- sc[2]
Pfit_k <- lm(Porep ~ Areal_r)
sc_k <- coefficients(Pfit_k)
psc_k <- sc_k[2]
lines((Pkk <- Porep/psc_k)+delts[8], col="blue")
lines(residuals(Pfit_k)/psc_k+delts[9], col="dark grey")
lines((Pkk <- Porep/psc)+delts[10], col="blue")
text(5000, delts[10]*1.10, sprintf("e/e ratio of  nu=%.02f  for  nu/(1-2*nu)",(R<-psc_k/psc)/(1+2*R)),cex=0.8)
lines(residuals(Pfit)/psc+delts[11], col="dark grey")
#
text(5000, delts+0.1,
     c("g1","g2",
	     "ett (hoop)",
	     "err (radial)", 
	     "ett+err",
	     "areal", 
	     "res ett+err ~ areal",
	     sprintf("pp ~ ett+err * %.02g 1/Pa",1e-6/psc_k),
	     "res pp ~ ett + err",
	     sprintf("pp ~ areal * %.02g 1/Pa",1e-6/psc),
	     "res pp ~ areal"),
     cex=0.9)
dev.off()
