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
\section{NA values}{
  \code{NA} values are allowed.  If present, the function
  \code{\link[zoo]{na.locf}}, which stands for "Last
  Observation Carried Forward", is used to impute them with
  real numbers.
}

