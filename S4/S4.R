
psd3 <- spectrum(rnorm(1e3), plot=FALSE)

#tmpp <- psd3
#class(tmpp) <- "spec2"

specS4 <- setClass("specS4",
  representation=representation(freq="numeric",
                                spec="numeric",
                                coh="numeric",
                                phase="numeric",
                                kernel="numeric",
                                df="numeric",
                                bandwidth="numeric",
                                n.used="numeric",
                                orig.n="numeric",
                                series="character",
                                snames="character",
                                method="character",
                                taper="numeric",
                                pad="numeric",
                                detrend="logical",
                                demean="logical"),
                   prototype = prototype(coh=numeric(0),
                                         phase=numeric(0),
                                         kernel=numeric(0),
                                         df=Inf,
                                         snames="",
                                         detrend=FALSE,
                                         demean=FALSE),
#                    validity = function(object) {
#                      #if (is.null(object@coh))
#                      # force validity
#                      TRUE
#                    }
  )


setClass('spec',prototype=specS4())
setAs("spec", "specS4", function(from, to){
  new(to, 
      freq = from$freq,
      spec = from$spec,
      coh = as.numeric(from$coh),
      phase = as.numeric(from$phase),
      kernel = as.numeric(from$kernel),
      snames = as.character(from$snames))
})
as(psd3, 'specS4')
