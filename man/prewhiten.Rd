\name{prewhiten}
\alias{prewhiten}
\alias{prewhiten.default}
\alias{prewhiten.ts}
\title{Prewhiten a series.}
\usage{
  prewhiten(tser, AR.max = 0L, detrend = TRUE,
    demean = TRUE, plot = TRUE, verbose = TRUE)

  prewhiten.default(tser, x.fsamp = 1, x.start = c(1, 1),
    ...)

  prewhiten.ts(tser, AR.max = 0L, detrend = TRUE,
    demean = TRUE, plot = TRUE, verbose = TRUE)
}
\arguments{
  \item{tser}{vector; An object to prewhiten.}

  \item{AR.max}{numeric; the maximum AR order to fit.}

  \item{detrend}{logical; Should a trend (and mean) be
  removed?}

  \item{demean}{logical; Should a mean value be removed?}

  \item{plot}{logical; Should the results be plotted?}

  \item{verbose}{logical; Should messages be printed?}

  \item{x.fsamp}{sampling frequency (for non \code{ts}
  objects)}

  \item{x.start}{start time of observations (for non
  \code{ts} objects)}

  \item{...}{variables passed to \code{\link{prewhiten.ts}}
  (for non \code{ts} objects)}
}
\description{
  Remove (optionally) mean, trend, and Auto Regressive (AR)
  model from the original series.
}
\section{Mean and trend (\code{AR.max < 1})}{
  Power spectral density estimates can become badly biased
  (especially at lower frequencies) if a signal of the form
  \eqn{f(x) = A x + B} is not removed from the series. If
  \code{detrend=TRUE} a model of this form is removed over
  the entire series using a linear least-squares estimator;
  in this case a mean value is removed regardless of the
  logical state of \code{demean}. To remove \emph{only} a
  mean value, set \code{detrend=FALSE} and (obviously)
  \code{demean=TRUE}.
}

\section{Auto Regressive (AR) innovations (\code{AR.max > 0})}{
  When an autoregressive model is removed from a
  non-stationary series, the residuals are known as
  'innovations', and may be stationary (or very-nearly
  stationary). This functions fits an AR model [order at
  least 1, but up to and including AR(\code{AR.max})] to
  the series by solving the Yule-Walker equations; however,
  AIC is used to estimate the highest significant order,
  which means that higher-order components may not
  necessarily be fit.

  The resulting innovations can be used to better estimate
  the stationary component of the original signal, or even
  in an interactive editing process to manually identify
  outliers.

  A quick way to determine whether this may be needed for
  the series is to run \code{acf} on it, and see if
  significant auto-correlations are found.
}

\section{NA values}{
  \code{NA} values are allowed.  If present, the function
  \code{\link[zoo]{na.locf}}, which stands for "Last
  Observation Carried Forward", is used to impute them with
  real numbers.
}
\author{
  Andrew Barbour <andy.barbour@gmail.com>
}
\seealso{
  \code{\link{psdcore}}, \code{\link{pspectrum}}
}
\keyword{S3methods}
\keyword{autoregressive-innovations}
\keyword{prewhiten}
\keyword{timeseries}

