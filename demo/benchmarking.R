library(multitaper)
library(rbenchmark)
library(plyr)
library(ggplot2)
library(reshape2)
# refresh latest core functionality
source('~/nute.processing/development/rlpSpec/rsrc/func_psdcore.R')
PSD <- .psdcore.default
#
set.seed(1234)
# differences in processing time reduce with orders of magnitude increases
run.bench <- function(nd, nt=8, reps=10){
  print(nd)
  initEnv(refresh=TRUE)
  X.d <- arima.sim(list(order = c(1,1,0), ar = 0.9),n=nd)
  X.p <- prewhiten(X.d, 10, plot=FALSE, verbose=FALSE)
  ntaps <- rep.int(nt,nd)
  # not really any improvement in speed for byte-compiled version:
  #library(compiler)
  #PSDc <- cmpfun(.psdcore.default)
  psdbench <- benchmark(PSD(X.d, ntaper=nt, plot=FALSE, force_calc=TRUE),
                        PSD(X.d, ntaper=ntaps, plot=FALSE, force_calc=TRUE),
                        multitaper::spec.mtm(X.d,k=nt,plot=FALSE),
                        replications=reps)
  psdbench$num_terms <- nd
  psdbench$num_taps <- nt
  psdbench$test <- c("rlpSpec::psdcore","rlpSpec::psdcore (vec.)","multitaper::spec.mtm")
  return(psdbench)
}
#
nds <- round(10**seq.int(from=1,to=4.0,by=0.1))
allbench <- lapply(X=nds, FUN=function(x) run.bench(nd=x))
save(allbench, file="~/Google Drive/PUB/GEOKOOK/rlpspec_allbench.Rda")

allbench.df <- plyr::ldply(allbench)
allbench.df.drp <- subset(allbench.df, 
                          select = c(test, num_terms, user.self, sys.self, elapsed, relative) 
                          )
allbench.df.mlt <- melt(allbench.df.drp, id.vars=c("test","num_terms"))
head(tmpd<-plyr::ddply(allbench.df.mlt,
                       .(variable,num_terms),
                  summarise, # Note: can get confused with Hmisc::summarize
                  summary="medians",
                  value=mean_cl_normal(value)[1,1]))
tests<-unique(allbench.df$test)
allmeds <- ldply(lapply(X=tests, FUN=function(x,df=tmpd){df$test <- x; return(df)}))


g <- ggplot(data=allbench.df.mlt, aes(x=log10(num_terms), y=log2(value), colour=test, group=test)) + 
  scale_colour_discrete(guide="none") + geom_point() 
g2 <- g + facet_grid(variable~test, scales="free_y")
g2 + geom_path(colour="black", data=allmeds, aes(group=test))

# Profiling
do.prof <- function(){
  #initEnv(refresh=TRUE)
  nd <- 1e2 # differences in processing time reduce with orders of magnitude increases
  X.d <- arima.sim(list(order = c(1,1,0), ar = 0.9),n=nd)
  nt <- 8
  Rprof()
  for (i in 1:1000)
    mypsd <- PSD(X.d, ntaper=nt, plot=FALSE, force_calc=TRUE)
  Rprof(NULL)
  prof <- summaryRprof()
  head(prof$by.self)
  head(prof$by.total)
}

