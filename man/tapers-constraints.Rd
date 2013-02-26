\name{tapers-constraints}
\alias{constrain_taper_friedman_smooth}
\alias{constrain_taper_loess_smooth}
\alias{constrain_taper_markov_chain}
\alias{constrain_taper_simple_slope}
\alias{constrain_tapers}
\alias{constrain_tapers.tapers}
\alias{ctap_friedman}
\alias{ctap_friedman.tapers}
\alias{ctap_loess}
\alias{ctap_loess.tapers}
\alias{ctap_markov}
\alias{ctap_markov.tapers}
\alias{ctap_simple}
\alias{ctap_simple.tapers}
\alias{minspan}
\alias{minspan.tapers}
\alias{tapers-constraints}
\title{Taper constraint methods.}
\usage{
  minspan(tapvec, ...)

  \method{minspan}{tapers} (tapvec, ...)

  constrain_tapers(tapvec, tapseq = NULL,
    constraint.method = c("simple.slope", "markov.chain", "loess.smooth", "friedman.smooth", "none"),
    verbose = TRUE, ...)

  \method{constrain_tapers}{tapers} (tapvec, tapseq = NULL,
    constraint.method = c("simple.slope", "markov.chain", "loess.smooth", "friedman.smooth", "none"),
    verbose = TRUE, ...)

  ctap_simple(tapvec, tapseq = NA, maxslope = 1, ...)

  \method{ctap_simple}{tapers} (tapvec, tapseq = NA,
    maxslope = 1, ...)

  ctap_markov(tapvec, tapseq = NA,
    chain.width = round(5 * length(tapvec)), ...)

  \method{ctap_markov}{tapers} (tapvec, tapseq = NA,
    chain.width = round(5 * length(tapvec)),
    normalize = TRUE, ...)

  ctap_loess(tapvec, tapseq = NULL, loess.span = 0.3,
    loess.degree = 1, verbose = TRUE, ...)

  \method{ctap_loess}{tapers} (tapvec, tapseq = NULL,
    loess.span = 0.3, loess.degree = 1, verbose = TRUE,
    ...)

  ctap_friedman(tapvec, tapseq = NULL, smoo.span = 0.3,
    smoo.bass = 2, verbose = TRUE, ...)

  \method{ctap_friedman}{tapers} (tapvec, tapseq = NULL,
    smoo.span = 0.3, smoo.bass = 2, verbose = TRUE, ...)
}
\arguments{
  \item{tapvec}{\code{'tapers'} object; the number of
  tapers at each frequency}

  \item{tapseq}{vector; positions or frequencies --
  necessary for smoother methods}

  \item{constraint.method}{character; method to use for
  constraints on tapers numbers}

  \item{verbose}{logical; should warnings and messages be
  given?}

  \item{maxslope}{integer; constrain based on this maximum
  first difference}

  \item{chain.width}{scalar; the width the MC should
  consider for the change probability}

  \item{normalize}{logical; should the refined tapers be
  normalized?}

  \item{loess.span}{scalar; the span used in \code{loess}}

  \item{loess.degree}{scalar; the polynomial degree}

  \item{smoo.span}{scalar; fraction of the observations in
  the span of the running lines smoother}

  \item{smoo.bass}{scalar; controls the smoothness of the
  fitted curve}

  \item{...}{optional arguments (unused)}
}
\value{
  An object with class \code{'tapers'}.
}
\description{
  In the Riedel-Sidorenko recipe, the number of optimal
  tapers at each frequency is strongly dependent on the
  first and second derivatives of the spectrum. It is
  crucial to enforce constraints on the number of actual
  tapers applied; this is because the derivatives of
  "noisy" series can be bogus.

  \code{\link{minspan}} sets the maximum span a tapers
  object may have, which is necessary because it would be
  nonsense to have more tapers than the length of the
  series.

  \code{\link{constrain_tapers}} refines the number of
  tapers; the method by which it does this is chosen with
  the \code{constraint.method} parameter. See
  \strong{Constraint methods} section for descriptions of
  each method. Below is a summary of the function
  associated with each \code{constraint.method}: \itemize{
  \item \code{'simple.slope'} uses
  \code{\link{ctap_simple}} \item \code{'markov.chain'}
  uses \code{\link{ctap_markov}} \item
  \code{'loess.smooth'} uses \code{\link{ctap_loess}} \item
  \code{'friedman.smooth'} uses \code{\link{ctap_friedman}}
  \item \code{'none'} returns unbounded tapers. }
}
\details{
  \code{\link{minspan}} bounds the number of tapers to
  within the minimum of either the maximum number of tapers
  found in the object, or the half-length of the series.
}
\section{Details of Constraint Methods}{
  \subsection{via first differencing (default)}{
  \code{\link{ctap_simple}} is the default, and preferred
  constraint method. The algortihm uses first-differencing
  to modify the number of tapers in the previous position.
  Effectively, the constraint is based on a causal,
  1st-order Finite Impulse-response Filter (FIR) which
  makes the method sensitive to rapid changes in the number
  of tapers; naturally, smoother spectra tend to produce
  less fluctuation in taper numbers, which makes this well
  suited for adaptive processing.

  This produces, generally, the most stable results,
  meaning repeatedly running the constraint will not change
  values other than on the first execution; the same cannot
  be said for the other methods.

  In pure-R this algorithm can be very slow; however, here
  we have included it as dynamically loaded c-code so it it
  reasonably fast. }

  \subsection{via Markov Chain}{ \code{\link{ctap_markov}}
  uses a Markov Chain, based on the theory of quantum-well
  probability chains, which are commonly used in gamma-ray
  spectroscopy.

  The main function behind this method is from Morhac
  (2008): \code{SpectrumSmoothMarkov}. This calculates the
  probability that the number of tapers would have changed
  (from it's previous value); it is very fast.  Details of
  the theory behind this algorithm may be found in Morhac
  (2008) and Silagadze (1996). }

  \subsection{via LOESS smoothing}{
  \code{\link{ctap_loess}} uses \code{loess} to smooth the
  taper vector; is can be very slow thanks to quadratic
  scaling. }

  \subsection{via Friedman super-smoothing}{
  \code{\link{ctap_friedman}} uses \code{supsmu}, the
  Friedman super-smoother. }
}

\section{Warning}{
  \code{\link{ctap_markov}} can produce "unstable" results
  in the sense that for successive application on taper
  vectors, even modest sized serially-correlated peaks
  tends to sharpen; hence, this method should be used with
  care, unless the intention is to specifically enhance
  peaks.  The \code{'chain.width'} parameter controls the
  broadness of the a priori distribution.  As a rule of
  thumb: the smaller the parameter is, the shorter the
  tails become.

  \code{\link{ctap_loess}} results tend to be strongly
  dependent on the tuning parameters given to \code{loess}
  (for obvious reasons); hence, some effort should be given
  to understand their effect, and/or re-tuning them if
  needed.

  \code{\link{ctap_friedman}} results are generally poor in
  my opinion; hence, the method may be removed in future
  releases.
}
\examples{
\dontrun{
require(psd)
##
## Taper constraint procedures
##
data(magnet)
X <- magnet$clean
##
## spectrum, then riedsid
kopt <- riedsid(PSD <- psdcore(X, ntaper=10, refresh=TRUE))
kopt.loess  <- riedsid(PSD, c.method="loess.smooth")
kopt.super  <- riedsid(PSD, c.method="friedman.smooth")
kopt.markov <- riedsid(PSD, c.method="markov.chain")
#
plot(kopt, log="y", ylim =c(.1, 3e2))
lines(kopt.super, log="y", col="red")
lines(kopt.loess, log="y", col="green")
lines(kopt.markov, log="y", col="orange")
##
##
## To compare all the methods at once:
demo("ctap")
##
}
}
\author{
  A.J. Barbour <andy.barbour@gmail.com> and R.L.Parker. AJB
  adapted some of RLP's original code, and wrote the main
  function in \code{\link{ctap_simple}} for dynamic loading
  C-code. The main function used by
  \code{\link{ctap_markov}} is from Morhac (2008).
}
\references{
  Morhac, M. (2008), Peaks: Peaks, \emph{R package},
  \strong{version 0.2}

  Silagadze, Z.K. (1996), A new algorithm for automatic
  photopeak searches, \emph{Nucl. Instrum. Meth. A},
  \strong{376} 451,
  \url{http://arxiv.org/abs/hep-ex/9506013}
}
\seealso{
  \code{\link{splineGrad}}, \code{\link{riedsid}}
}
\keyword{riedel-sidorenko}
\keyword{tapers}
\keyword{tapers-constraints}

