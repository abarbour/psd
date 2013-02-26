#RDEX#\dontrun{
require(psd)
##
## Objects with class 'tapers'
##
is.tapers(as.tapers(1))
is.tapers(as.tapers(1:10))
is.tapers(as.tapers(matrix(1:10,ncol=1)))
as.tapers(list(x=1:10,y=1:30)) # note dimensions
as.tapers(x<-data.frame(x=1:10,y=10:19))
as.tapers(x, min_taper=3, max_taper=10)
# class 'character' is in-coercible; raise error
try(as.tapers(c("a","b")), silent=TRUE)
#RDEX#}
