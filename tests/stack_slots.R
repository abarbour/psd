psd3 <- spectrum(rnorm(1e3), plot=FALSE)
summary(psd3)
# Length Class  Mode     
# freq      500    -none- numeric  
# spec      500    -none- numeric  
# coh         0    -none- NULL     
# phase       0    -none- NULL     
# kernel      0    -none- NULL     
# df          1    -none- numeric  
# bandwidth   1    -none- numeric  
# n.used      1    -none- numeric  
# orig.n      1    -none- numeric  
# series      1    -none- character
# snames      0    -none- NULL     
# method      1    -none- character
# taper       1    -none- numeric  
# pad         1    -none- numeric  
# detrend     1    -none- logical  
# demean      1    -none- logical 
class(unclass(psd3))
# [1] "list"
is.object(psd3) & !isS4(psd3)
# [1] TRUE
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

psd4 <- specS4()
class(unclass(psd4))
# [1] "S4"
validObject(psd4)
# [1] TRUE

spec_check <- function(X) UseMethod("spec_check")
spec_check.NULL <- function(X){
  message("spec check NULL")
  #X[X==NULL] <- numeric(0)
  return(numeric(0))
}
spec_check.default <- function(X){
  return(X)
}

##
##
nonull <- function(psd) UseMethod("nonull")
nonull.spec <- function(psd){
  stopifnot(inherits(psd, 'spec', FALSE))
  # spec.pgram may return NULL for these:
  psd$coh <- as.numeric(psd$coh)
  psd$phase <- as.numeric(psd$phase)
  psd$kernel <- as.numeric(psd$kernel)
  psd$snames <- as.character(psd$snames)
  return(psd)
}

as.specS4 <- function(psd) UseMethod("as.specS4")
as.specS4.spec <- function(psd){
  stopifnot(inherits(psd, 'spec', FALSE))
  psd <- nonull.spec(psd)
  S4spec <- specS4()
  S4spec@freq <- psd$freq
  S4spec@spec <- psd$spec
  S4spec@coh <- psd$coh
  S4spec@phase <- psd$phase
  S4spec@kernel <- psd$kernel
  S4spec@snames <- psd$snames
  S4spec@df <- psd$df
  S4spec@bandwidth <- psd$bandwidth
  S4spec@n.used <- psd$n.used
  S4spec@orig.n <- psd$orig.n
  S4spec@series <- psd$series
  S4spec@snames <- psd$snames
  S4spec@method <- psd$method
  S4spec@taper <- psd$taper
  S4spec@pad <- psd$pad
  S4spec@detrend <- psd$detrend
  S4spec@demean <- psd$demean
  #... is there a better way to do this?
  #
  return(S4spec)
}
as.specS4.specS4 <- function(psd){
  stopifnot(inherits(psd, 'specS4', FALSE))
  message("already specS4")
  return(invisible(psd))
}

psd4c1 <- as.specS4(psd4)
psd4c2 <- as.specS4(psd3)

# setClass("Matrix",
#          representation(Dim = "integer", Dimnames = "list", "VIRTUAL"),
#          prototype = prototype(Dim = integer(2), Dimnames = list(NULL,NULL)),
#          validity = function(object) {
#            Dim <- object@Dim
#            if (length(Dim) != 2)
#              return("Dim slot must be of length 2")
#            if (any(Dim < 0))
#              return("Dim slot must contain non-negative values")
#            Dn <- object@Dimnames
#            if (!is.list(Dn) || length(Dn) != 2)
#              return("'Dimnames' slot must be list of length 2")
#            lDn <- sapply(Dn, length)
#            if (lDn[1] > 0 && lDn[1] != Dim[1])
#              return("length(Dimnames[[1]])' must match Dim[1]")
#            if (lDn[2] > 0 && lDn[2] != Dim[2])
#              return("length(Dimnames[[2]])' must match Dim[2]")
#            ## 'else'	ok :
#            TRUE
#          })

