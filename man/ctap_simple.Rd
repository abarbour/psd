\name{ctap_simple}
\alias{constrain_taper_simple_slope}
\alias{ctap_simple}
\alias{ctap_simple.taper}
\title{ctap_simple}
\usage{
  ctap_simple(tapvec, tapseq = NA, maxslope = 1, ...)

  ctap_simple.taper(tapvec, tapseq = NA, maxslope = 1, ...)
}
\arguments{
  \item{tapvec}{'taper' object; the number of tapers at
  each frequency}

  \item{tapseq}{(unused) vector; positions or frequencies
  -- necessary for smoother methods}

  \item{maxslope}{(\code{'simple.slope'}) integer;
  constrain based on this maximum first difference}

  \item{...}{(unused) optional argments}
}
\value{
  an object with class 'taper', with a constrained taper
  numbers
}
\description{
  Constrain tapers with first differencing.
}
\details{
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
}
\note{
  The results obtained by \strong{\code{'simple.slope'}}
  are generally the most stable, meaning repeatedly running
  the constraint will not change values other than on the
  first execution; the same cannot be said for the other
  methods.
}
\examples{
\dontrun{
ntap <- as.taper(1:10)
ctap_simple(ntap)
}
}
\seealso{
  \code{\link{constrain_tapers}}
}

