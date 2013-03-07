#vertical axis rotation of [areal gamma1 gamm2] referenced to X degrees (compass
#dir):
#
#       [areal (invariant), c*gamma1 + s*gamma2, -s*gamma1 + c*gamma2]
#where
#       c = cos(2*rotation)
#       s = sin(2*rotation)
load("Tohoku.rda")
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
  shears.r <- complex(modulus = Mod(shears.c), argument = Arg(shears.c) - backaz)
  return(shears.r)
}
AREALS <- function(backaz.deg, areal, shears.){
  backaz. <- backaz.deg * pi / 180
  shears.r <- SHEARROT(backaz., shears.)
  # 3) rotated shears
  Eh <- Re(shears.r)
  Er <- Im(shears.r)
  # 4) regress equivalent areal strains
  areal_r <- Eh + Er
  Efit <- lm(areal_r ~ areal)
  # intercept is implicit, and the first coeff
  # second coeff is areal_r to areal fit
  arc <- coefficients(Efit)[[2]]
  resnorm <- norm(res <- as.matrix(residuals(Efit)))
  message(backaz.deg-bazd)
  sqrt(var(res))
}
FUN <- function(az) AREALS(az, areal=Areal, shears.=Shears)
#Opt <- optim(baz*180/pi, FUN, gr=NULL, method="L-BFGS-B", lower=-360, upper=360)
Opt <- optim(bazd, FUN, gr=NULL, method="L-BFGS-B", lower=-90, upper=90)
bazopt <- Opt$par[1]
bazopt

shears.r <- SHEARROT(bazopt*pi/180, Shears)
Eh <- Re(shears.r)
Er <- Im(shears.r)

#par(mfcol=c(2,1), pty="m")
par(pty="m")

delts <- 0.4*rev(seq_len(8)-4)
plot(Shears[,1]+delts[1], type="l", lty=1, ylim=range(delts)*1.2)
lines(Shears[,2]+delts[2], col="red")
lines(Eh+delts[3])
lines(Er+delts[4], col="red")
lines(areal+delts[5], col="black")
areal_r <- Eh + Er
Efit <- lm(areal_r ~ areal)
summary(Efit)

lines(areal_r+delts[6], col="red")
lines(residuals(Efit)+delts[7], col="grey")
text(2000, delts+0.1,
     c("g1","g2","ett hoop", "err radial", 
	     "areal", "ett+err", "res ett+err ~ areal","pore pressure (sc)"))

Pfit <- lm(Porep ~ areal)
psc <- coefficients(Pfit)[2]
Pfit_r <- lm(Porep ~ areal_r)
psc_r <- coefficients(Pfit_r)[2]
lines(Porep/psc_r+delts[8], col="blue")
#lines(Porep/psc+delts[8], col="light blue")
#
