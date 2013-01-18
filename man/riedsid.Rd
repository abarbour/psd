\name{riedsid}
\alias{riedsid}
\alias{riedsid.default}
\alias{riedsid.spec}
\title{Constrained, optimal tapers using the Riedel & Sidorenko--Parker method.}
\usage{
  riedsid(psd, ntaper, tapseq = NULL, constrained = TRUE,
    c.method = NULL, ...)

  riedsid.spec(pspec, ...)

  riedsid.default(psd, ntaper, tapseq = NULL,
    constrained = TRUE, c.method = NULL, ...)
}
\arguments{
  \item{psd}{vector; the spectral values used to optimize
  taper numbers}

  \item{pspec}{object with class 'spec'}

  \item{ntaper}{scalar or vector; number of tapers to apply
  optimization}

  \item{tapseq}{vector; representing positions or
  frequencies (same length as psd)}

  \item{constrained}{logical; should the taper constraints
  be applied to the optimum tapers?}

  \item{c.method}{string; constraint method to use if
  \code{constrained=TRUE}}

  \item{...}{optional argments passed to
  \code{constrain_tapers}}
}
\description{
  Estimates the optimal number of tapers at each frequency
  of given PSD; the optimal value chosen is based on the
  Riedel-Sidorenko MSE recipe, and modified by R.L. Parker.
}
\details{
  First, weighted derivatives of the input PSD are
  computed. Using those derivates the optimal number of
  tapers is found through the RS-RLP formulation. The
  default state returns constrained optimal-tapers using
  \code{\link{constrain_tapers}}, but this can be turned
  off with \code{constrained=FALSE}.
}
\examples{
#'
#' \\dontrun{
#'
##
## Riedel-Sidorenko--Parker taper optimization
##
set.seed(1234)
## some params
nd <- 1e3 # num data
ntap <- 10 # num tapers
nrm <- 40 # sharpness of the peaks rel 2*variance
#
## create a pseudo spectrum
riex <- rnorm(nd) + nrm*abs(cos(pi*(0:(nd-1))/180) + 1.2)
#
## optimize tapers
rtap <- riedsid(riex, ntaper=ntap)
# ... again, but do not constrain
rtap2 <- riedsid(riex, ntaper=ntap, constrained=FALSE)
#
## plot
op <- par(no.readonly = TRUE)
par(mfrow=c(2,1), mar=rep(1.3,4), mai=rep(0.6,4))
# ... the mock spectrum
plot(riex, type="h", xaxs="i") 
# ... the optimal tapers
# unconstrained: note logarithmic peaks from noise in spectrum
# even when variance is 20 times smaller than peak amplitudes
plot(rtap2, log="y") 
# constrained (red curve): much better behaved
lines(rtap) 
# original tapers:
lines(as.taper(rep.int(ntap,nd)), col="blue")
par(op)
##
#' }
}
\author{
  A.J. Barbour <andy.barbour@gmail.com> ported original by
  R.L. Parker.
}
\seealso{
  \code{\link{constrain_tapers}}, \code{\link{psdcore}}
}
\keyword{riedel-sidorenko}
\keyword{taper}
\keyword{taper-constraints}

