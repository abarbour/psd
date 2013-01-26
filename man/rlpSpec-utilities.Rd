\docType{package}
\name{rlpSpec-utilities}
\alias{as.rowvec}
\alias{char2envir}
\alias{colvec}
\alias{dB}
\alias{decibels}
\alias{envir2char}
\alias{is.spec}
\alias{is.taper}
\alias{mod}
\alias{mod.default}
\alias{modulo}
\alias{na_mat}
\alias{na_mat.default}
\alias{ones}
\alias{ones.default}
\alias{rlpSpec-utilities}
\alias{rlpSpec-utilities-package}
\alias{rowvec}
\alias{vardiff}
\alias{vector_reshape}
\alias{vector_reshape.default}
\alias{zeros}
\alias{zeros.default}
\title{Various utility functions.}
\usage{
  vardiff(Xd)

  dB(Rat)

  char2envir(envchar)

  envir2char(envir)

  vector_reshape(x,
    vec.shape = c("horizontal", "vertical"))

  \method{vector_reshape}{default} (x,
    vec.shape = c("horizontal", "vertical"))

  colvec(x)

  rowvec(x)

  is.spec(Obj)

  is.taper(Obj)

  na_mat(nrow, ncol = 1)

  \method{na_mat}{default} (nrow, ncol = 1)

  zeros(nrow)

  \method{zeros}{default} (nrow)

  ones(nrow)

  \method{ones}{default} (nrow)

  mod(X, Y)

  \method{mod}{default} (X, Y)
}
\arguments{
  \item{Xd}{object to difference}

  \item{Rat}{numeric; A ratio to convert to decibels
  (\code{dB}).}

  \item{envchar}{An object with class 'character'.}

  \item{envir}{An object of class 'environment'.}

  \item{x}{An object to reshape (\code{vector_reshape}).}

  \item{vec.shape}{choice between horizontally-long or
  vertically-long vector.}

  \item{Obj}{An object to test for class inheritance.}

  \item{nrow}{integer; the number of rows to create.}

  \item{ncol}{integer; the number of columns to create
  (default 1).}

  \item{X}{numeric; the "numerator" of the modulo division}

  \item{Y}{numeric; the "denominator" of the modulo
  division}
}
\value{
  numeric

  numeric

  \code{char2envir} returns the result of evaluating the
  object: an environment object; \code{envir2char} returns
  the result of deparsing the environment name: a character
  string.

  \code{vector_reshape} returns a "reshaped" vector,
  meaning it has had it's dimensions changes so that it has
  either one row (if \code{vec.shape=="horizontal"}), or
  one column (\code{"vertical"}).

  \code{is.spec} and \code{is.taper} both return logicals
  about whether or not the object does have class 'spec' or
  'taper', respectively

  \code{na_mat} returns a matrix of dimensions
  \code{(nrow,ncol)} with \code{NA} values, the
  representation of which is set by \code{NA_real_}

  For \code{zeros} or \code{ones} respectively, a matrix
  vector with \code{nrow} zeros or ones.

  \code{mod} returns the result of a modulo division, which
  is equivalent to \code{(X) \%\% (Y)}.
}
\description{
  Various utility functions are used:

  \code{vardiff} reeturns the variance of the first
  difference of the series.

  \code{dB} returns an object converted to decibels.

  \code{char2envir} converts a character string of an
  environment name to an evaluated name; whereas,
  \code{envir2char} converts an environment name to a
  character string.

  \code{vector_reshape} reshapes a vector into another
  vector.

  \code{colvec} returns the object as a vertically long
  vector; whereas \code{rowvec} returns the object as a
  horizontally long vector.

  \code{is.spec} reports whether an object has class S3
  class 'spec', as would one returned by, for example,
  \code{spectrum}.

  \code{is.taper} reports whether an object has S3 class
  'taper', as would one returned by, for example,
  \code{\link{as.taper}}.

  \code{na_mat} populates a matrix of specified dimensions
  with \code{NA} values.

  \code{zeros} populate a column-wise matrix with zeros;
  whereas, \code{ones} populates a column-wise matrix with
  ones

  \code{mod} finds the modulo division of X and Y.
}
\details{
  Decibels are defined as \deqn{10 \log{}_{10} x}.

  \code{colvec, rowvec} are simple wrapper functions to
  \code{vector_reshape}.

  Modulo division has higher order-of-operations ranking
  than other arithmetic operations; hence, \code{x + 1 \%\%
  y} is equivalent to \code{x + (1 \%\% y)} which can
  produce confusing results. \code{mod} is simply a series
  of \code{trunc} commands which reduces the chance for
  unintentionally erroneous results.
}
\note{
  \code{char2envir} ensures the \code{envchar} object is a
  character, so that something is not unintentionally
  evaluated; \code{envir2char} simply deparses the object
  name.

  The performance of \code{mod} has not been tested against
  the \code{\%\%} arithmetic method -- it may or may not be
  slower for large numeric vectors.
}
\examples{
db(1) # signal is equal <--> zero dB
##
## Evaluate character strings
print(.rlpenv)
char2envir(.rlpenv)
char2envir("some nonexistent environment") # error
# and environment objects:
print(.GlobalEnv)
envir2char(.GlobalEnv)
envir2char(.rlpSpecEnv)
char2envir(some_nonexistent_environment) # error
##
##
## Check for spec object:
# quick power spectral density
x <- rnorm(1e2, sd=10)
psd <- spectrum(x, plot=FALSE)
# return is class 'spec'
is.spec(psd) # TRUE
#
# but the underlying structure is just a list
psd <- unclass(psd)
is.spec(psd) # FALSE
##
## check for taper
is.taper(x)
is.taper(as.taper(x))
##
## matrix and vector creation:
# NA matrix
nd <- 5
na_mat(nd)
na_mat(nd,nd-1)
# zeros
zeros(nd)
zeroes(nd)
# and ones
ones(nd)
##
##
# \\dontrun{
##
## modulo division
##
x <- 1:10
mc1a <- mod(1,2)
mc2a <- mod(1+x,2)
mc1b <- 1 \%\% 2
mc2b <- 1 + x \%\% 2
mc2c <- (1 + x) \%\% 2
all.equal(mc1a, mc1b) # TRUE
all.equal(mc2a, mc2b) # "Mean absolute difference: 2"
all.equal(mc2a, mc2c) # TRUE
##
# }
}
\author{
  A.J. Barbour <andy.barbour@gmail.com>
}
\references{
  For \code{\link{mod}}: see Peter Dalgaard's explanation
  of the non-bug (#14771) I raised (instead I should've
  asked it on R-help):
  \url{https://bugs.r-project.org/bugzilla3/show_bug.cgi?id=14771\#c2}
}
\seealso{
  \code{\link{rlpSpec-package}}

  \code{\link{rlpSpec-package}}

  \code{\link{rlpSpec-package}}

  \code{\link{as.taper}}
}
\keyword{S3methods}
\keyword{arithmetic-operations}
\keyword{environment}
\keyword{inherits}
\keyword{is}
\keyword{matrix-creation}
\keyword{matrix-manipulation}
\keyword{methods}
\keyword{modulo-division}
\keyword{normalization}
\keyword{utilities}
\keyword{vector-creation}
\keyword{vector-manipulation}

