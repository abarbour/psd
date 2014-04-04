\name{riedsid}
\alias{riedsid}
\alias{riedsid.default}
\alias{riedsid.spec}
\title{Constrained, optimal tapers using the Riedel & Sidorenko--Parker method.}
\usage{
riedsid(PSD, ntaper, tapseq = NULL, Deriv.method = c("local_qls", "spg"),
  Local.loss = c("Optim", "Less", "More"), constrained = TRUE,
  c.method = NULL, verbose = TRUE, ...)

\method{riedsid}{spec}(PSD, ...)

\method{riedsid}{default}(PSD, ntaper, tapseq = NULL,
  Deriv.method = c("local_qls", "spg"), Local.loss = c("Optim", "Less",
  "More"), constrained = TRUE, c.method = NULL, verbose = TRUE, ...)
}
\arguments{
  \item{PSD}{vector or class 'spec'; the spectral values
  used to optimize taper numbers}

  \item{ntaper}{scalar or vector; number of tapers to apply
  optimization}

  \item{tapseq}{vector; representing positions or
  frequencies (same length as \code{PSD})}

  \item{Deriv.method}{character string; choice of gradient
  estimation method}

  \item{Local.loss}{string; sets how sensitive the spectral
  derivatives are}

  \item{constrained}{logical; should the taper constraints
  be applied to the optimum tapers?}

  \item{c.method}{string; constraint method to use if
  \code{constrained=TRUE}}

  \item{verbose}{logical; should messages be printed?}

  \item{...}{optional argments passed to
  \code{\link{constrain_tapers}}}
}
\value{
Object with class 'tapers'.
}
\description{
Estimates the optimal number of tapers at each frequency of
given PSD, using a modified Riedel-Sidorenko MSE recipe
(RS-RLP).
}
\details{
The optimization is as follows. First, weighted derivatives
of the input PSD are computed. Using those derivates the
optimal number of tapers is found through the RS-RLP
formulation. Constraints are then placed on the practicable
number of tapers.

\subsection{Taper constraints}{ The parameter
\code{c.method} provides an option to change the method of
taper constraints.  A description of each may be found in
the documentation for \code{\link{constrain_tapers}}.

Once can use \code{constrained=FALSE} to turn off all taper
constraints; this could lead to strange behavior though. }

\subsection{Spectral derivatives}{ The parameter
\code{Deriv.method} determines which method is used to
estimate derivatives. \itemize{ \item{\code{"local_qls"}}{
(\strong{default}) uses quadratic weighting and local
least-squares estimation; then, \code{Local.loss} can alter
slightly the weighting to make the derivatives more or less
succeptible to changes in spectral values. Can be slower
than \code{"spg"}.} \item{\code{"spg"}}{ uses
\code{\link{splineGrad}}; then, additional arguments may be
passed to control the smoothness of the derivatives (e.g
\code{spar} in \code{smooth.spline}).} } }
}
\section{Warning}{
  The \code{"spg"} can become numerically unstable, and
  it's not clean when it will be preferred to the
  \code{"local_qls"} method other than for efficiency's
  sake.
}
\examples{
#RDEX#\dontrun{
require(psd)
##
## Riedel-Sidorenko--Parker taper optimization
##
set.seed(1234)
# some params
nd <- 512 # num data
ntap <- 10 # num tapers
nrm <- 40 # sharpness of the peaks rel 2*variance
#
# create a pseudo spectrum
# with broad peaks
riex <- rnorm(nd) + nrm*abs(cos(pi*(x<-0:(nd-1))/180) + 1.2)
riex <- riex + 8*nrm*dcauchy(x, nd/3)
riex <- riex + 5*nrm*dnorm(x, nd/2)
# flat regions
riex[riex<25] <- 25
ried <- dB(riex, invert=TRUE)
#
# optimize tapers
rtap <- riedsid(riex, ntaper=ntap)
#
# plot
op <- par(no.readonly = TRUE)
par(mfrow=c(2,1), mar=rep(1.3,4), mai=rep(0.6,4))
# ... the mock spectrum
plot(riex, type="h", xaxs="i", ylim=c(0,200)) 
# ... the optimal tapers
plot(rtap, log="y") 
# original tapers:
lines(as.tapers(rep.int(ntap,nd)), col="blue")
par(op)
#RDEX#}
}
\author{
A.J. Barbour <andy.barbour@gmail.com> adapted original by
R.L. Parker.
}
\seealso{
\code{\link{constrain_tapers}}, \code{\link{psdcore}},
\code{smooth.spline}
}
\keyword{riedel-sidorenko}
\keyword{tapers}
\keyword{tapers-constraints}

