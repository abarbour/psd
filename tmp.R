x <- 1:10

#http://stackoverflow.com/questions/11445702/ema-computation-using-filter-function-in-r

# Since filter computes
# y[n] = x[n] + alpha * y[n-1]
# you need to rescale the result.

f <- function(x,win) {
  alpha <- 2/(win+1)
  filter(x, 1-alpha, method="recursive", side=1, init=x[1]/alpha)*alpha
}
x <- 1:10
k <- 3
# getEMA2(x,k)
f(x,k) # identical

#http://stackoverflow.com/questions/9265553/r-filter-a-vector-using-a-function

#
# derivs <- splineGrad(kopt, kseq, plot.derivs=F)
# slopes <- derivs$dydx
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
