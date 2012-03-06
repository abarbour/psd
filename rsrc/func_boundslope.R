##
##  Default method to forward and backwards constrain the slope of
##  the second derivative of the spectrum
##
boundslope.default <- function(kopt){
  ##
  ## Args:  
  ##
  ## Returns:	
  ##
  ## TODO(abarbour):
  ##
  #
  #  Curb run-away growth of kopt due to zeros of d2psd (S'') 
  #  limits slopes to be < 1 in magnitude, preventing greedy averaging:
  #
  #
  # ?lapply: lapply returns a list of the same length as X, each element of 
  # which is the result of applying FUN to the corresponding element of X.
  #
  # lapply(X, FUN, ...)
  #X    a vector (atomic or list) or an expression object. Other objects 
  #     (including classed objects) will be coerced by base::as.list.
  #FUN	the function to be applied to each element of X: see ‘Details’. In the 
  #     case of functions like +, %*%, the function name must be backquoted or quoted.
  #...	optional arguments to FUN.
  #
  nK <- length(kopt)
  Aspan <- 2:nK
  kopt.A <- kopt[Aspan]
  Bspan <- 1:(nK-1)
  kopt.B <- kopt[Bspan]
  # first difference
  kopt.dK <- diff(kopt, lag=1)
  df <- data.frame(x=Aspan, slopes=kopt.dK, dir="init")
  # where slopes are <= 1
  states <- rep(0,(nK-1))
  dKinds <- kopt.dK <= 1
  kopt.A[dKinds] <- kopt.B[dKinds]+1
  kopt[Aspan] <- kopt.A
  # first difference
  kopt.dK <- diff(kopt, lag=1)
  df <- rbind(df, data.frame(x=Aspan, slopes=kopt.dK, dir="fw"))
  ##
  ## now reverse direction
  ##
  rAspan <- nK:2
  kopt.A <- kopt[rAspan]
  Bspan <- (nK-1):1
  kopt.B <- kopt[Bspan]
  # first difference
  kopt.dK <- diff(rev(kopt), lag=1)
  # where slopes are >= -1
  states <- rep(0,(nK-1))
  dKinds <- kopt.dK >= -1
  kopt.A[dKinds] <- kopt.B[dKinds]+1
  kopt[rAspan] <- kopt.A
  ##
  # first difference
  kopt.dK <- diff(kopt, lag=1)
  df <- rbind(df, data.frame(x=Aspan, slopes=kopt.dK, dir="bw"))
  library(ggplot2)
  print(ggplot(df, aes(x=x, y=slopes, colour=dir))+ geom_line())
  ##
  return(invisible(kopt))
}


dostate <- function(Kopt,tspan,constraint){
  # rev [] and kopt []
  Aspan <- tspan
  kopt.A <- Kopt[Aspan]
  Bspan <- Aspan - 1
  kopt.B <- Kopt[Bspan]
  # first difference
  kopt.dK <- diff(Kopt, lag=1)
  states <- rep(0, length(Aspan))
  # find indices for the constraint
  dKinds <- kopt.dK >= constraint
  # replace
  kopt.A[dKinds] <- kopt.B[dKinds] + 1
  Kopt[Aspan] <- kopt.A
  # return
  return(Kopt)
}