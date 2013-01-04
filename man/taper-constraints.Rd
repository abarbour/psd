\docType{package}
\name{taper-constraints}
\alias{constrain_taper_friedman_smooth}
\alias{constrain_taper_loess_smooth}
\alias{constrain_taper_markov_chain}
\alias{constrain_taper_simple_slope}
\alias{constrain_tapers}
\alias{constrain_tapers.taper}
\alias{ctap_friedman}
\alias{ctap_friedman.taper}
\alias{ctap_loess}
\alias{ctap_loess.taper}
\alias{ctap_markov}
\alias{ctap_markov.taper}
\alias{ctap_simple}
\alias{ctap_simple.taper}
\alias{minspan}
\alias{minspan.taper}
\alias{taper-constraints}
\alias{taper-constraints-package}
\title{Taper constraint methods.}
\usage{
  minspan(tapvec, nf = length(tapvec), ...)

  minspan.taper(tapvec, nf = length(tapvec), ...)

  constrain_tapers(tapvec, tapseq = NULL,
    constraint.method = c("simple.slope", "markov.chain", "loess.smooth", "friedman.smooth", "none"),
    min.tapers = 1, verbose = TRUE, ...)

  constrain_tapers.taper(tapvec, tapseq = NULL,
    constraint.method = c("simple.slope", "markov.chain", "loess.smooth", "friedman.smooth", "none"),
    min.tapers = 1, verbose = TRUE, ...)

  ctap_simple(tapvec, tapseq = NA, maxslope = 1, ...)

  ctap_simple.taper(tapvec, tapseq = NA, maxslope = 1, ...)

  ctap_markov(tapvec, tapseq = NA,
    chain.width = round(5 * length(tapvec)), ...)

  ctap_markov.taper(tapvec, tapseq = NA,
    chain.width = round(5 * length(tapvec)), ...)

  ctap_loess(tapvec, tapseq = NULL, loess.span = 0.3,
    loess.degree = 1, verbose = TRUE, ...)

  ctap_loess.taper(tapvec, tapseq = NULL, loess.span = 0.3,
    loess.degree = 1, verbose = TRUE, ...)

  ctap_friedman(tapvec, tapseq = NULL, smoo.span = 0.3,
    smoo.bass = 2, verbose = TRUE, ...)

  ctap_friedman.taper(tapvec, tapseq = NULL,
    smoo.span = 0.3, smoo.bass = 2, verbose = TRUE, ...)
}
\arguments{
  \item{tapvec}{'taper' object; the number of tapers at
  each frequency}

  \item{nf}{scalar; number of positions or frequencies}

  \item{...}{(unused) optional argments}

  \item{tapvec}{vector; the number of tapers at each
  frequency}

  \item{tapseq}{vector; positions or frequencies --
  necessary for smoother methods}

  \item{constraint.method}{method to use for constraints on
  tapers numbers}

  \item{min.tapers}{integer; the minimum number of tapers}

  \item{verbose}{logical; should warnings and messages be
  given?}

  \item{...}{optional argments passed to the
  constraint.method}

  \item{tapvec}{'taper' object; the number of tapers at
  each frequency}

  \item{tapseq}{(unused) vector; positions or frequencies
  -- necessary for smoother methods}

  \item{maxslope}{(\code{'simple.slope'}) integer;
  constrain based on this maximum first difference}

  \item{...}{(unused) optional argments}

  \item{tapvec}{'taper' object; the number of tapers at
  each frequency}

  \item{tapseq}{(unused) vector; positions or frequencies
  -- necessary for smoother methods}

  \item{chain.width}{scalar; the width the MC should
  consider for the change probability}

  \item{...}{(unused) optional argments}

  \item{tapvec}{'taper' object; the number of tapers at
  each frequency}

  \item{tapseq}{vector; positions or frequencies --
  necessary for smoother methods}

  \item{loess.span}{scalar; the span used in \code{loess}}

  \item{loess.degree}{scalar; the polynomial degree}

  \item{verbose}{logical; should warnings and messages be
  given?}

  \item{...}{optional argments; unused}

  \item{tapvec}{'taper' object; the number of tapers at
  each frequency}

  \item{tapseq}{vector; positions or frequencies --
  necessary for smoother methods}

  \item{smoo.span}{scalar; fraction of the observations in
  the span of the running lines smoother}

  \item{smoo.bass}{scalar; controls the smoothness of the
  fitted curve}

  \item{verbose}{logical; should warnings and messages be
  given?}

  \item{...}{optional argments; unused}
}
\value{
  an object with class 'taper', with a constrained taper
  numbers

  an object with class 'taper', with a constrained taper
  numbers

  an object with class 'taper', with a constrained taper
  numbers

  an object with class 'taper', with a constrained taper
  numbers

  an object with class 'taper', with a constrained taper
  numbers
}
\description{
  Taper constraint methods.

  Find the max span of taper object

  Apply constraints on the number of tapers

  Constrain tapers with first differencing.

  Constrain tapers with a Markov Chain, based on
  quantum-well probability.

  Constrain tapers with the Loess smoother.

  Constrain tapers with the Friedman 'super-smoother'.
}
\details{
  minimize the number of tapers: min(nf/2, 7*tapvec/5)

  Refines the number of tapers; the method by which it does
  this may be chosen by the user.

  This is the default, and preferred constraint method. The
  algortihm uses first-differencing to modify the number of
  tapers in the previous position.  Effectively, the
  constraint is based on a causal, 1st-order Finite
  Impulse-response Filter (FIR) which makes the method
  sensitive to rapid changes in the number of tapers;
  naturally, smoother spectra tend to produce less
  fluctuation in taper numbers, which makes this well
  suited for adaptive processing.

  In pure-R this algorithm can be very slow; however, here
  we have included it as dynamically loaded c-code so it it
  reasonably fast.

  This method uses a Markov Chain, based on the theory of
  quantum-well probability, commonly used in gamma-ray
  spectroscopy.

  The main function behind this method is
  \link[Peaks]{SpectrumSmoothMarkov} which calculates the
  probability that the number of tapers would have changed
  (from it's previous value); it is very fast.  Details of
  the theory behind this algorithm may be found in Morhac
  (2008) and Silagadze (1996).

  This method uses \code{stats::loess} to smooth the taper
  vector; is can be very slow thanks to quadratic scaling.

  This method uses the base Fiedman super-smoother.
}
\note{
  The value of \code{nf} should represent the total number
  of tapers from the original spectrum; the default assumes
  the taper object represents that size.

  The results obtained by \strong{\code{'simple.slope'}}
  are generally the most stable, meaning repeatedly running
  the constraint will not change values other than on the
  first execution; the same cannot be said for the other
  methods.

  This algorithm can produce "unstable" results in the
  sense that for successive application on taper vectors,
  even modest sized serially-correlated peaks tends to
  sharpen; hence, this method should be used with care,
  unless the intention is to specifically enhance peaks.
  The \code{'chain.width'} parameter controls the broadness
  of the a priori distribution.  As a rule of thumb: the
  smaller the parameter is, the shorter the tails become.

  The results obtained can be strongly dependent on the
  \code{loess} tuning parameters; hence, some effort should
  be given to understand their effect, and/or re-tuning the
  default parameters.

  The results obtained by \strong{\code{'friedman.smooth'}}
  are generally poor; hence, the method may be removed in
  future releases.
}
\examples{
\dontrun{
ntap <- as.taper(1:10)
par(mfrow=c(2,1))
plot(ntap)
plot(minspan(ntap))
par(mfrow=c(1,1))
}
\dontrun{
ntap <- as.taper(1:10)
ctap_simple(ntap)
}
\dontrun{
ntap <- as.taper(1:10)
ctap_markov(ntap)
}
\dontrun{
ntap <- as.taper(1:10)
ctap_loess(ntap)
}
\dontrun{
ntap <- as.taper(1:10)
ctap_friedman(ntap)
}
}
\author{
  Andrew Barbour <andy.barbour@gmail.com>

  Andrew Barbour <andy.barbour@gmail.com> ported original
  by R.L.Parker.

  Andrew Barbour <andy.barbour@gmail.com> ported original
  by R.L.Parker to R and C.

  Andrew Barbour <andy.barbour@gmail.com> ported original
  by R.L.Parker to R and C.

  Andrew Barbour <andy.barbour@gmail.com> adapted the
  algorithm from Morhac (2008) for use here.

  Andrew Barbour <andy.barbour@gmail.com>

  Andrew Barbour <andy.barbour@gmail.com>
}
\references{
  Morhac, M. (2008), Peaks: Peaks, \emph{R package},
  \strong{version 0.2}

  Silagadze, Z.K. (1996), A new algorithm for automatic
  photopeak searches, \emph{Nucl. Instrum. Meth. A},
  \strong{376} 451.

  \url{http://arxiv.org/abs/hep-ex/9506013}
}
\seealso{
  \code{\link{constrain_tapers}}

  \code{\link{ctap_simple}}, \code{\link{ctap_markov}},
  \code{\link{riedsid}}

  \code{\link{constrain_tapers}}

  \code{\link{constrain_tapers}},
  \code{\link[Peaks]{SpectrumSmoothMarkov}}

  \code{\link{constrain_tapers}}, \code{\link{loess}}

  \code{\link{constrain_tapers}}, \code{\link{supsmu}}
}
\keyword{taper}
\keyword{taper-constraints}

