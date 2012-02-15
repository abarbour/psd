###
###  S3 (or S4) methods for printing, plotting, etc
###
##
## summaries
##
summary.qual <- function(object, ...){
  ###
  res <- list(call=object$call)
  class(res) <- "summary.qual"
  res
}
summary.psd <- function(object, ...) {
  # from doc example
  #   se <- sqrt(diag(object$vcov)) 
  #   tval <- coef(object) / se
  #   TAB <- cbind(Estimate = coef(object), StdErr = se,
  #                t.value = tval, p.value = 2*pt(-abs(tval), df=object$df))
  #   res <- list(call=object$call, coefficients=TAB)
  res <- list(call=object$call)
  class(res) <- "summary.psd" 
  res
}
##
## printing
##
print.summary.qual <- function(x, ...){
  ## develop
  cat("Call:\n") 
  print(x$call) 
  cat("\n")
}
print.summary.psd <- function(x, ...) {
  ## develop
  cat("Call:\n") 
  print(x$call) 
  cat("\n")
  #   printCoefmat(x$coefficients, P.value=TRUE, has.Pvalue=TRUE)
}
print.psd <- function(psd, ...){
  ##
  ## Args:  
  ##
  ## Returns:	
  ##
  ## TODO(abarbour):
  ##
  cat("Call:\n") 
  print(x$call)
  cat("\n%%%%% PSD SUMMARY %%%%%\n")
  summary(psd)
}
##
## plotting
##
plot.psd <- function(psd.df, niter=NULL, ...){
  ## Plot the results of the PSD estimation
  ##
  ## Args:  
  ##
  ## Returns:	
  ##
  ## TODO(abarbour): convert this to a method [ ]
  ##
  require(ggplot2, quietly=TRUE, warn.conflicts=FALSE)
  require(gridExtra, quietly=TRUE, warn.conflicts=FALSE)
  
  dims <- dim.data.frame(psd.df)
  nrow <- dims[1]
  nvar <- dims[2]
  pltdf <- psd.df[2:(nrow-1),]
  
  optset <- function(title="", legendpos="none"){
    #https://kohske.wordpress.com/2010/12/25/drawing-on-full-region-in-ggplot2/
    opts.full <- opts( 
      title=title,
      plot.title = theme_text(size=14, lineheight=.1, face="italic", hjust=0),
      legend.background = theme_rect(colour="gray"),
      legend.position = legendpos,
      panel.background = theme_blank(),
      #     panel.grid.major = theme_blank(),
      panel.grid.minor = theme_blank(),
      panel.margin = unit(0,"null"),
      plot.margin = rep(unit(0,"null"),4),
      axis.ticks = theme_blank(),
      axis.text.x = theme_blank(),
      axis.text.y = theme_blank(),
      #     axis.title.x = theme_blank(),
      #     axis.title.y = theme_blank(),
      axis.ticks.length = unit(0,"null"),
      axis.ticks.margin = unit(0,"null")
      )
  }
  if (is.null(niter) && exists("Niter")){niter <- Niter}
  # plot the spectrum in dB, colored by ntapers
  p1 <- ggplot(pltdf, aes(x=f, y=20*log10(psd)))+
    geom_step(size=1, aes(colour=log10(ntaper)))+
    scale_colour_gradient("Tapers applied\n(log10)",low = "black", high = "red")+
    #,breaks=(1:3), labels=10**(1:3))+
    scale_x_log10("Frequency, log10 1/N/dT", minor_breaks=NA)+
    scale_y_continuous("PSD, dB rel. units**2 * N * dT", minor_breaks=NA)+
    optset("A: Tapered spectrum", c(0.18,0.50))
  # the number of tapers, colored by spectral levels
  p2 <- ggplot(pltdf, aes(x=f, y=ntaper))+
    geom_step(size=1, aes(colour=20*log10(psd)))+
    scale_colour_gradient("PSD, dB",low = "black", high = "red")+
    scale_y_log10("Tapers, log10", breaks=10**(0:3), labels=10**(0:3))+
    scale_x_log10("Frequency, log10 1/N/dT", minor_breaks=NA)+
    optset("B: Number of tapers applied", c(0.15,0.60))
  
  #   print(p1)
  plts <- grid.arrange(p1, p2, 
                       main=sprintf("Adaptive Sine-multitaper PSD Estimation\n%i iterations",niter))
}
# end plot.psd
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
###