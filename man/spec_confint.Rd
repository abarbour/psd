\name{spec_confint}
\alias{spec_confint}
\alias{spec_confint.default}
\alias{spec_confint.spec}
\alias{spec_confint.tapers}
\title{Multitaper PSD confidence intervals.}
\usage{
  spec_confint(dof, p = 0.95, as.db = FALSE)

  \method{spec_confint}{spec} (dof, p = 0.95,
    as.db = FALSE)

  \method{spec_confint}{tapers} (dof, p = 0.95,
    as.db = FALSE)

  \method{spec_confint}{default} (dof, p = 0.95,
    as.db = FALSE)
}
\arguments{
  \item{dof}{numeric; the degrees of freedom \eqn{\nu}}

  \item{p}{numeric; the coverage probability \eqn{p}, bound
  within \eqn{[0,1)}}

  \item{as.db}{logical; should the values be returned as
  decibels?}
}
\value{
  A \code{data.frame} with the following properties (and
  names): \itemize{ \item{\code{lower}: Based on upper tail
  probabilities (\eqn{p})} \item{\code{upper}: Based on
  lower tail probabilities (\eqn{1-p})}
  \item{\code{median}: Based on lower tail probabilities
  (\eqn{p=0.5})} \item{\code{approx}: Approximation based
  on \eqn{1/\sqrt(\nu - 1)}.} }
}
\description{
  Multitaper PSD confidence intervals.
}
\details{
  The errors are estimated from the number of degrees of
  freedom \eqn{\nu} by evaluating the
  \eqn{\chi_{p,\nu}^{2}(\nu,\nu)} distribution for an
  optional coverage probability \eqn{p} (defaulting to
  \eqn{p=0.95}). Additionally, the \eqn{p=0.5} values and
  an approximation from \eqn{1/\sqrt{\nu - 1}} are
  returned.

  A more sophisticated (and complicated) approach would be
  to estimate via jack-knifing (Prieto et al 2007), but
  this is not yet made available.

  Additive uncertainties \eqn{\delta S} are returned, such
  that the spectrum with confidence interval is \eqn{S \pm
  \delta S}.
}
\examples{
\dontrun{
require(psd)
##
## Confidence intervals from taper numbers
##
sp <- spectral_properties(as.tapers(1:50), p=0.95, db.ci=TRUE)
par(las=1)
plot(stderr.chi.upper ~ taper, sp, type="s",
       ylim=c(-10,20), yaxs="i", xaxs="i",
       xlab=expression("number of tapers ("* nu/2 *")"), ylab="dB",
       main="Spectral uncertainties")
mtext("(additive factor)", line=.3)
lines(stderr.chi.lower ~ taper, sp, type="s")
lines(stderr.chi.median ~ taper, sp, type="s", lwd=2)
lines(stderr.chi.approx ~ taper, sp, type="s", col="red",lwd=2)
# to reach 3 db width confidence interval at p=.95
abline(v=33, lty=3)
legend("topright",
        c(expression("Based on "* chi^2 *"(p,"*nu*") and (1-p,"*nu*")"),
          expression(""* chi^2 *"(p=0.5,"*nu*")"),
          "approximation"),
lwd=c(1,3,3), col=c("black","black","red"), bg="white")
##
}
}
\author{
  A.J. Barbour <andy.barbour@gmail.com>, modified from the
  \code{spec.ci} function inside \code{stats::plot.spec}.
}
\seealso{
  \code{\link{spectral_properties}},
  \code{\link{psd-package}}, \code{plot.spec},
  \code{\link{dB}}
}
\keyword{degrees-of-freedom}
\keyword{properties}
\keyword{tapers}
\keyword{uncertainty}

