qualcon <-function(x, critvar=10) {
  ###
  # PORT of RLP's qualcon.m
  # with mods by abarbour
  # Dec 2011
  #
  # porting:  4 Jan 2012
  # testing:  4 Jan 2012
  #           noise: rnorm(n,sd=0.2)
  #           sig:   .001*(1:n) + triang(n)+3*sin(pi*1:n/180+pi/4)+3*sin(2*pi*1:n/180+pi/4)
  #           temporary outliers are easily found
  #              x[100:120] <- 5.4
  #              x[1000:1020] <- 5.5
  #              
  ###
  #
  #  Searches for data errors indicated by the whitened series
  #  falling outside fixed limits.
  #  Plots +-100 points of x and whitened x around error section.
  #  badpoints <- number error regions detected
  #
  lx <- length(x)
  x.w <- whiten(x)
  wcrit <- critvar * sqrt(var(x.w))
  par(ask=TRUE)
  toret <- show.white(wcrit, 1:lx, x, x.w, pltlabs=FALSE)
  Icrit <- which(abs(x.w) > wcrit)
  lIc <- length(Icrit)
  if (lIc > 0){
    badpoints <- 0
    i2 <- 0
    for ( j  in  1:length(Icrit) ) {
      if ( Icrit[j] >= i2 ){
        i1 <- Icrit[j]-100 
        i2 <- Icrit[j]+100
        k1 <- max(1, i1)  
        k2 <- min(i2, lx)
        k1k2 <- k1:k2
        show.white(wcrit, k1k2, x, x.w, limsc=2.8)
        badpoints <- 1 + badpoints
      }
    }
  }
  return(invisible(toret))
} # end qualcon

show.white <-function(wcrit, k1k2, x, x.w, limsc=2.5, pltlabs=TRUE) {
  #  Display the series and its filtered version
  t1 <- k1k2[1]
  te <- k1k2[length(k1k2)]
  xplt <- x[k1k2]
  xplt <- xplt - median(xplt)
  xplt.w <- x.w[k1k2]
  x.wcriti <- which(abs(xplt.w) > wcrit) + t1 - 1
  x.wcrit <- x.w[x.wcriti]
  #
  plot(k1k2, xplt, type="s", col="gray", 
       ylim=c(-1*limsc*wcrit, limsc*wcrit),
       ylab="value",
       xlab="term",
       main="Outlier Inspector")
  
  abline(h=wcrit*c(1,-1), lty=3, col="black")
  abline(h=wcrit*c(2,-2), lty=3, col="red")
  if (pltlabs==TRUE){
    text(t1+25, -0.35*wcrit, "demeaned F(Xa:Xb)", pos=1, cex=0.8, col="dark gray")
    text(t1+17, -0.60*wcrit, "AR prew. F(X)", pos=1, cex=0.8, col="blue")
    text(t1+20, -1*wcrit, sprintf("threshold (%.3f)\ncrossings (+)",wcrit), pos=1, cex=0.8, col="black")
    text(t1+10, -2*wcrit, "2x thresh", pos=3, cex=0.8, col="red")
  }
  if (length(x.wcriti)>0){
    points(x.wcriti, x.wcrit, col="black", pch=3, cex=0.8)
  }
  lines(k1k2, xplt.w, type="s", col="blue")
  # display some info
  #   terms
  cat('Display of terms: ', t1, te, "\n")
  #   outliers
  if (length(x.wcriti>0)){
    cat("Threshold crossing indices:\n")
    toret <- t(rbind(x.wcriti,x.wcrit))
    print(toret)
    return(invisible(toret))
  } else {
    cat(sprintf("No critical outliers found for threshold of  %.03f\n",wcrit))
  }
} # end show.white
