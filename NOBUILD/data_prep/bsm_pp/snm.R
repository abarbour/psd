SNM <- function(molten=TRUE, verbose=TRUE){
  bsm <- read.table("http://www.seismosoc.org/publications/BSSA_html/bssa_101-5/2011062-esupp/BSSA-D-11-00062_S1_noisemod-bsm.txt", header=TRUE)
  bsm$meter.type <- "PBO borehole strainmeters (Anza, CA)"
  lsm <- read.table( "http://www.seismosoc.org/publications/BSSA_html/bssa_101-5/2011062-esupp/BSSA-D-11-00062_S2_noisemod-lsm.txt", header=TRUE)
  lsm$meter.type <- "UCSD laser strainmeters"
  nm <- rbind(bsm, lsm)
  attr(nm, "source.doi") <- "10.1785/0120110062"
  if (molten){
    require(reshape2)
    idv <- c("meter.type","freq")
    nm <- melt(nm, id.vars=idv)
    attr(nm,"molten.vars") <- idv
  }
  if (verbose){
    on.exit(
    message(cat("A ggplot2-style plot would be something like this:\n
    nm <- SNM(molten=TRUE)
    ggplot(nm)+
    geom_path(aes(x=log10(freq), y=value, linetype=meter.type, colour=variable))\n")))
    return(invisible(nm))
  } else {
    return(nm)
  }
}

bsm <- read.table("bsmnm.txt",header=TRUE)
bsm$meter.type <- "PBO borehole strainmeters (Anza, CA)"
lsm <- read.table("lsmnm.txt",header=TRUE)
lsm$meter.type <- "UCSD laser strainmeters"
snm <- rbind(bsm, lsm)
attr(snm, "source.doi") <- "10.1785/0120110062"
attr(snm, "generator") <- structure(SNM)
save(snm, file="strainnoise.rda")
