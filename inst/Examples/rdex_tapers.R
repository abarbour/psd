#' dontrun{
##
## Objects with class 'taper'
##
is.taper(as.taper(1))
is.taper(as.taper(1:10))
is.taper(as.taper(matrix(1:10,ncol=1)))
as.taper(list(x=1:10,y=1:30)) # note dimensions
as.taper(x<-data.frame(x=1:10,y=10:19))
as.taper(x, min_taper=3, max_taper=10)
# class 'character' is in-coercible; raise error
try(as.taper(c("a","b")))
#' }
