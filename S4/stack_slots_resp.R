## setClass
## A simple class with two slots
track <- setClass("track", representation(x="numeric", y="numeric"))

## setAs
setAs("track", "numeric", function(from) from@y)

t1 <- new("track", x=1:20, y=(1:20)^2)

as(t1, "numeric")

psd4as <- asS4(psd3)
slotNames(psd4as)
TOSPECS4 <- function(from){
  to <- specS4(from$freq, 
               freq = from$freq,
               spec = from$spec,
               coh = from$coh,
               phase = from$phase,
               kernel = from$kernel,
               snames = from$snames)
  return(to)
}

psd3 <- spectrum(rnorm(1e3), plot=FALSE)
tmpp <- psd3
class(tmpp) <- "spec2"
setClass('spec',prototype=specS4())
setOldClass("spec2", S4Class="specS4")

setAs("spec", "specS4", function(from, to){
  from <- nonull.spec(from)
  new(to, 
      freq = from$freq,
      spec = from$spec,
      coh = from$coh,
      phase = from$phase,
      kernel = from$kernel,
      snames = from$snames)
})

as(psd3, 'specS4')
as(tmpp, "specS4")

## The next example shows:
##  1. A virtual class to define setAs for several classes at once.
##  2. as() using inherited information

setClass("ca", representation(a = "character", id = "numeric"))

setClass("cb", representation(b = "character", id = "numeric"))

setClass("id")
setIs("ca", "id")
setIs("cb", "id")


setAs("id", "numeric", function(from) from@id)

CA <- new("ca", a = "A", id = 1)
CB <- new("cb", b = "B", id = 2)

setAs("cb", "ca", function(from, to )new(to, a=from@b, id = from@id))

as(CB, "numeric")