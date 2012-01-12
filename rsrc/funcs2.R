#
mod <- function(x,y){
  x1<-trunc(trunc(x/y)*y)
  z<-trunc(x)-x1
  z
}
#
plot.psd <- function(psd.df, niter=NULL){
  require(ggplot2, quietly=T)
  require(gridExtra, quietly=T)
#   require(RColorBrewer)
  
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
    scale_colour_gradient("Tapers applied",low = "black", high = "red")+
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
