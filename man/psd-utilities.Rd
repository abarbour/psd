\name{psd-utilities}
\alias{as.rowvec}
\alias{char2envir}
\alias{colvec}
\alias{dB}
\alias{db}
\alias{decibels}
\alias{envir2char}
\alias{is.spec}
\alias{is.tapers}
\alias{mod}
\alias{mod.default}
\alias{modulo}
\alias{na_mat}
\alias{na_mat.default}
\alias{ones}
\alias{ones.default}
\alias{psd-utilities}
\alias{rowvec}
\alias{varddiff}
\alias{vardiff}
\alias{vector_reshape}
\alias{vector_reshape.default}
\alias{zeros}
\alias{zeros.default}
\title{Various utility functions.}
\usage{
vardiff(Xd, double.diff = FALSE)

varddiff(Xd)

dB(Rat, invert = FALSE, pos.only = TRUE, is.power = FALSE)

char2envir(envchar)

envir2char(envir)

vector_reshape(x, vec.shape = c("horizontal", "vertical"))

\method{vector_reshape}{default}(x, vec.shape = c("horizontal", "vertical"))

colvec(x)

rowvec(x)

is.spec(Obj)

is.tapers(Obj)

na_mat(nrow, ncol = 1)

\method{na_mat}{default}(nrow, ncol = 1)

zeros(nrow)

\method{zeros}{default}(nrow)

ones(nrow)

\method{ones}{default}(nrow)

mod(X, Y)

\method{mod}{default}(X, Y)
}
\arguments{
  \item{Xd}{object to difference}

  \item{double.diff}{logical; should the double difference
  be used instead?}

  \item{Rat}{numeric; A ratio to convert to decibels
  (\code{dB}).}

  \item{invert}{logical; assumes \code{Rat} is already in
  decibels, so return ratio}

  \item{pos.only}{logical; if \code{invert=FALSE}, sets
  negative or zero values to NA}

  \item{is.power}{logical; should the factor of 2 be
  included in the decibel calculation?}

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
\code{char2envir} returns the result of evaluating the
object: an environment object; \code{envir2char} returns
the result of deparsing the environment name: a character
string.

\code{vector_reshape} returns a "reshaped" vector, meaning
it has had it's dimensions changes so that it has either
one row (if \code{vec.shape=="horizontal"}), or one column
(\code{"vertical"}).

\code{is.spec} and \code{is.tapers} both return logicals
about whether or not the object does have class 'spec' or
'tapers', respectively

\code{na_mat} returns a matrix of dimensions
\code{(nrow,ncol)} with \code{NA} values, the
representation of which is set by \code{NA_real_}

For \code{zeros} or \code{ones} respectively, a matrix
vector with \code{nrow} zeros or ones.

\code{mod} returns the result of a modulo division, which
is equivalent to \code{(X) \%\% (Y)}.
}
\description{
\emph{The various utility functions are:}

\code{vardiff} returns the variance of the first (or
second) difference of the series. \code{varddiff} is a
convenience wrapper to return variance for the second
difference.

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

\code{is.spec} reports whether an object has class S3 class
'spec', as would one returned by, for example,
\code{spectrum}.

\code{is.tapers} reports whether an object has S3 class
'tapers', as would one returned by, for example,
\code{\link{as.tapers}}.

\code{na_mat} populates a matrix of specified dimensions
with \code{NA} values.

\code{zeros} populate a column-wise matrix with zeros;
whereas, \code{ones} populates a column-wise matrix with
ones.  \emph{Note that \code{nrow} is enforced to be at
least 1 for both functions.}

\code{mod} finds the modulo division of X and Y.
}
\details{
Decibels are defined as \eqn{10 \log{}_{10}
\frac{X_1}{X_2}}, unless \code{is.power=TRUE} in which
\eqn{\mathrm{db} X^2 \equiv 20 \log{}_{10} X^2}

\code{colvec, rowvec} are simple wrapper functions to
\code{vector_reshape}.

Modulo division has higher order-of-operations ranking than
other arithmetic operations; hence, \code{x + 1 \%\% y} is
equivalent to \code{x + (1 \%\% y)} which can produce
confusing results. \code{mod} is simply a series of
\code{trunc} commands which reduces the chance for
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
#RDEX#\dontrun{
require(psd)
##
## Various utilities
##
set.seed(1234)
X <- rnorm(1e2)
##
## Matrix and vector creation:
##
# NA matrix
nd <- 5
na_mat(nd)
na_mat(nd,nd-1)
# zeros
zeros(nd)
# and ones
ones(nd)
##
## Check for tapers object:
##
is.tapers(X)
is.tapers(as.tapers(X))
##
## Check for spec object:
##
PSD <- spectrum(X, plot=FALSE)
# return is class 'spec'
is.spec(PSD) # TRUE
# but the underlying structure is just a list
PSD <- unclass(PSD)
is.spec(PSD) # FALSE
##
## Environment character strings
##
print(envname <- get_psd_env_name())
print(envir <- char2envir(envname))
try(char2envir("some nonexistent environment"), silent=TRUE) # error
# and environment objects:
print(.GlobalEnv)
envir2char(.GlobalEnv)
envir2char(envir)
try(envir2char(some_nonexistent_environment), silent=TRUE) # error
##
## decibels
##
dB(1) # signal is equal <--> zero dB
sig <- 1e-10
all.equal(sig, dB(dB(sig), invert=TRUE))
pow <- sig**2
all.equal(pow, dB(dB(sig, is.power=TRUE), invert=TRUE, is.power=TRUE))
## 
## Variance of difference series
##
vardiff(X)
varddiff(X)
all.equal(vardiff(X, TRUE), varddiff(X))
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
#RDEX#}
}
\author{
A.J. Barbour <andy.barbour@gmail.com>
}
\references{
For \code{\link{mod}}: see Peter Dalgaard's explanation of
the non-bug (#14771) I raised (instead I should've asked it
on R-help):
\url{https://bugs.r-project.org/bugzilla3/show_bug.cgi?id=14771\#c2}
}
\seealso{
\code{\link{psd-package}}, \code{\link{as.tapers}}
}
\keyword{S3methods}
\keyword{arithmetic-operations}
\keyword{decibel}
\keyword{environment}
\keyword{first-difference}
\keyword{inherits}
\keyword{is}
\keyword{matrix-creation}
\keyword{matrix-manipulation}
\keyword{methods}
\keyword{modulo-division}
\keyword{normalization}
\keyword{utilities}
\keyword{variance}
\keyword{vector-creation}
\keyword{vector-manipulation}

