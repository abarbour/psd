
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
