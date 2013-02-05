###
###  S3 (or S4) methods for PLOTTING
###
plot.psd <- function(psd.df, logx=TRUE, xlabel=NULL, ylabel=NULL, niter=NULL, showplot=TRUE,...){
  ##
  ## Plot the results of the PSD estimation
  ##
  ## Args:  
  ##
  ## Returns:  
  ##
  ## TODO(abarbour): 
  ## [ ] convert this to a method for class 'psd'
  ## [ ] MOD: http://learnr.wordpress.com/2009/05/26/ggplot2-two-or-more-plots-sharing-the-same-legend/
  ##
  require(sfsmisc, quietly=TRUE, warn.conflicts=FALSE)  # for nice labels
  require(ggplot2, quietly=TRUE, warn.conflicts=FALSE)  # for plot engine
  
  dims <- dim.data.frame(psd.df)
  nrow <- dims[1]
  nvar <- dims[2]
  pltdf <- psd.df[2:(nrow-1),]
#   if (is.null(niter) && exists("Niter")){niter <- Niter}
  if (exists("num_iter",envir=psdenv) & is.null(niter)){niter <- envGet("num_iter")}
  # percent spectral uncertainty from  number of tapers
  # sigma^2 ~ 6 S^2 / 5 K
#   A*(1-(seA-1)
  pltdf$sigma2 <- (pltdf$psd**2)*1.2/pltdf$ntaper
  
  ## plot setup and label breaks
  pltdf$x <- pltdf$f
  #
  if (is.null(ylabel)){ylabel <- "PSD, log10 units**2 * N * dT"}
  if (is.null(xlabel)){xlabel <- "Frequency, 1/N/dT"}
  if (logx){
    xlabel <- sprintf("%s log10",xlabel)
    pltdf$x <- log10(pltdf$x)
  }
  g <- ggplot(pltdf, aes(x=x, y=psd))#, group=src))
  atY <- 10^seq(-4, 10, by=1)
  atYL <- sfsmisc::axTexpr(2, at=atY, drop.1=TRUE)
  # plot grobs
  p <- g +
    # std err
    geom_ribbon(size=0.25, colour="black", aes(ymax=(1+sqrt(sigma2)/psd), ymin=(1-sqrt(sigma2)/psd), fill="sig")) +
    # tapers
    geom_ribbon(size=0.25, colour="black", aes(ymax=ntaper, ymin=envGet("init_tap"), fill="tap")) +
    # psd
    geom_path(size=0.5)+
    scale_x_continuous(xlabel)+ #, expand=c(0,0))+
    scale_y_log10(ylabel, breaks=atY, labels=atYL, expand=c(0,0))+
    scale_fill_discrete("",
                        breaks=c("sig","tap"), 
                        labels=c("Uncertainty: +- Sigma", "Optim. tapers: Rel. pilot spec."))+
    theme_bw()+
    opts(title=sprintf("Adaptive Sine-multitaper PSD Estimation\n%i iterations",niter),
         legend.position=c(0.83, 0.90))
  if (showplot){print(p)}
  return(invisible(p))
}
# end plot.psd
##
##
##
# plot.qualcon or plot.whiten from show.whiten in suppfuncs [ ]
show.whiten <- function(wcrit, k1k2, x, x.w, limsc=2.5, pltlabs=TRUE) {
  ##  Display the series and its filtered version, interactively
  ##
  ## Args:  
  ##
  ## Returns:	
  ##
  ## TODO(abarbour):	
  ##     a) interative outlier identification or incorporate other package?
  ##     b) make this a plot method
  ##
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
} 
# end show.whiten
#
#