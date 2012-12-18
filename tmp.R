# x <- 1:10
# 
# #http://stackoverflow.com/questions/11445702/ema-computation-using-filter-function-in-r
# 
# # Since filter computes
# # y[n] = x[n] + alpha * y[n-1]
# # you need to rescale the result.
# 
# f <- function(x,win) {
#   alpha <- 2/(win+1)
#   filter(x, 1-alpha, method="recursive", side=1, init=x[1]/alpha)*alpha
# }
# x <- 1:10
# k <- 3
# # getEMA2(x,k)
# f(x,k) # identical
# 
# #http://stackoverflow.com/questions/9265553/r-filter-a-vector-using-a-function
# #
# isGoodNumber <- function(X) 
# {
#  if (X==5) return(TRUE) else return(FALSE)
# }
# isGoodNumber = Vectorize(isGoodNumber)
# v<-c(1,2,3,4,5,5,5,5)
# isGoodNumber(v)
# # equiv to unlist(lapply(X=v,FUN=isGoodNumber))
# v[isGoodNumber(v)]
# #*or*
# Filter( isGoodNumber, v)
#
#
# derivs <- splineGrad(kopt, kseq, plot.derivs=F)
# slopes <- derivs$dydx # size nf

kopt <- c(1:10,12,14,16,20,30,20,16,14,12,10:1)
slopes <- c(rep(1,10),rep(2,2),0,rep(-2,2),rep(-1,10))
nf <- length(kopt)
x <- 1:nf
# derivs <- splineGrad(x, kopt, TRUE)
# slopes <- (derivs$dydx)
plot(x, slopes, type="h", lwd=4)

dkm <- 2
dk.max <- dkm*sign(slopes)
(dat <- data.frame(slope=(slopes[2:nf]), # dk
                   ntaps=kopt[2:nf], # k
                   tosub=FALSE,
                   ntaps.max=dk.max[2:(nf)],
                   ntaps.subs=kopt[1:(nf-1)] # val to sub
)) 
dat$ntaps2 <- dat$ntaps
dat$tosub[(abs(dat$slope)) >= dkm] <- TRUE
dat
dat$ntaps2[dat$tosub] <- dat$ntaps.subs[dat$tosub]
dat
plot(x, kopt, type="h")
lines(x[2:nf], dat$ntaps, col="blue", lwd=2)
lines(x[2:nf], dat$ntaps2, col="red", lwd=2)

# 'slope' is essentially the state of slope, e.g. is it g.t. 1?
# [ ] make the max slope settable
#
# set boolean vector of TRUE/FALSE whether slope is greater than max.slope

# state<-0
# for ( j  in  2:nf ) {
#   slope <- slopes[j]
#   if (state == 0) {
#     if (slope >= 1 ) {
#       state <- 1
#       kopt[j] <- kopt[j-1]+1
#     }
#   } else {
#     if (kopt[j] >= kopt[j-1]+1) {
#       kopt[j] <- kopt[j-1]+1
#     } else {
#       state <- 0
#     }
#   }
# }
library(Peaks)
kopt <- c(rep(1,20),1:10,12,14,16,20,150,20,16,14,16,18,20,30,20,18,16,14,12,10:1,rep(1,10))
plot(kopt, type="l", col="red", ylim=1.2*c(0.1,max(kopt)), log="y")
MC.probs <- seq(1,5,by=0.5)
smarks <- sapply(X=MC.probs, FUN=function(x){
  MC.win <- round(length(kopt)*x)
  smark <- SpectrumSmoothMarkov(kopt, MC.win)
  # poor form!
  toret <- sapply(X=smark, FUN=function(x) max(1,x))
  lines(toret, lwd=x)
  toret
  }, simplify=TRUE)
# points(kopt.adj <- sapply(X=kopt + floor(round(kopt-smark)), FUN=function(x) max(1,x)),pch=3,col="blue")
lines(round(rowMeans(cbind(kopt,smarks))),col="green",lwd=4)
lines(kopt, col="red",lwd=4)
