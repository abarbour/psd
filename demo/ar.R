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
  return(psdbench)
}
#
nds <- round(10**seq.int(from=1,to=4.0,by=0.1))
allbench <- lapply(X=nds, FUN=function(x) run.bench(nd=x))
allbench.df <- plyr::ldply(allbench)
allbench.df.drp <- subset(allbench.df, 
                          select = c(test, num_terms, user.self, sys.self, elapsed, relative) 
                          )
allbench.df.mlt <- melt(allbench.df.drp, id.vars=c("test","num_terms"))
#Using test as id variables

# plot the benchmark data
stopifnot(exists("allbench.df"))
g <- ggplot(allbench.df, aes(x=log10(num_terms), colour=test)) + 
  scale_colour_discrete(guide="none")
g + geom_line(aes(y=log2(user.self)))
g + geom_line(aes(y=log2(sys.self)))
g + geom_line(aes(y=log2(elapsed)))
g + geom_point(size=4,aes(y=log2(relative)))

g <- ggplot(allbench.df.mlt, aes(x=log10(num_terms), 
                                 y=log2(value), 
                                 colour=test,
                                 group=num_terms
                                 )
            ) + 
              scale_colour_discrete(guide="none") +
              stat_summary(fun.data="median_hilow", 
                           #colour="dark grey", 
                           geom="crossbar", 
                           width=0.2)
(p1 <- g + geom_line(aes(group=test)))

p1 + facet_grid(variable~test,scales="free")

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

