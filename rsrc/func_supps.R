###
###  Default methods for supplementary functions 
###    (to make porting from Matlab to R tractable and supplement)
###
# TODO(abarbour):
#	ONES
#	ZEROS
##
mod.default <- function(x,y){
  ## modulo division
  ## R %/% requires strict consideration of order of operations whereas this is
  ## function internal and thus less prone to error, but perhaps less efficient?
  ##
  ## Args:	x	val 1
  ##		y	val 2
  ##
  ## Returns:	modulo division of x vs y
  ##
  x1 <- trunc(trunc(x/y)*y)
  z <- trunc(x)-x1
  z
}
# end mod
###