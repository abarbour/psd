#http://www.bioconductor.org/packages/release/bioc/vignettes/Rgraphviz/inst/doc/Rgraphviz.R
#source("http://bioconductor.org/biocLite.R")
library(Rgraphviz)
Z <- c(
  "pspectrum",  #ps
  "pilot_spec", #p_s
  "riedsid", 	#r
  	"parabolic_weights", 	#p_w
  "psdcore",	#p
  	"prewhiten",		#pw
  	"constrain_tapers", 	#ct
	"minspan",		#m
  "spectral_properties"		#sp
)
##
plot(g1<-randomGraph(Z, 1:4, 0.3))
##
# niter -- num adapts
# ** -- check dependency [ ]
#
#	ps	p_s	r	p_w	p	pw	ct	m	sp
#ps	0	1	niter	0	niter	0	0	0	0
#p_s	1	0	0	0	1	0	0	0	0
#r	niter	0	0	1	0	0	**	**	0
#p_w	0	0	1	0	0	0	0	0	0
#p	niter	1	0	0	0	1	1	**	0
#pw	0	0	0	0	1	0	0	0	0
#ct	0	0	**	0	1	0	0	**	0
#m	0	0	**	0	**	0	**	0	0
#sp	0	0	0	0	0	0	0	0	0

#?graphNEL
set.seed(123)
V <- LETTERS[1:4]
nl <- length(Z) #4
edL <- vector("list", length=nl) #=4
names(edL) <- Z #V
for(i in 1:nl) #:4
  edL[[i]] <- list(edges=(nl+1)-i, weights=1) #runif(1)) #=5-i
gR <- graphNEL(nodes=Z, edgeL=edL, edgemode="dir") #=V
#edges(gR)
#edgeWeights(gR)
#plot(gR)
##
