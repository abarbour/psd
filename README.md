# psd 

Adaptive, sine multitaper power spectral density estimation for R

by Andrew J Barbour, Jonathan Kennel, and Robert L Parker

[![](https://www.r-pkg.org/badges/version-last-release/psd?color=green)](https://cran.r-project.org/package=psd) [![Travis Build Status](https://travis-ci.org/abarbour/psd.svg?branch=master)](https://travis-ci.org/abarbour/psd) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/abarbour/psd?branch=master&svg=true)](https://ci.appveyor.com/project/abarbour/psd) [![Code Coverage](https://codecov.io/github/abarbour/psd/coverage.svg?branch=master)](https://codecov.io/github/abarbour/psd?branch=master)  [![Downloads](https://cranlogs.r-pkg.org/badges/psd)](https://www.r-pkg.org/pkg/psd) [![License](https://img.shields.io/badge/license-GPL-lightgrey.svg)](https://www.gnu.org/licenses/gpl-2.0.html) [![Citation](https://img.shields.io/badge/published-CAGEO-red.svg)](https://doi.org/10.1016/j.cageo.2013.09.015)

### **_Latest News_**

_As of version 2.0, one can calculate the multivariate PSD ("cross spectrum") between two signals._

## Description

This is an [R](https://www.r-project.org)
package for computing univariate power spectral density
estimates with little or no tuning effort.
We employ sine multitapers, allowing the number to vary with frequency
in order to reduce mean square error, the sum of squared bias and
variance, at each point.  The approximate criterion of
[Riedel and Sidorenko (1995)](https://doi.org/10.1109/78.365298)
is modified to prevent runaway averaging that otherwise occurs when
the curvature of the spectrum goes to zero.  An iterative procedure
refines the number of tapers employed at each frequency.  The resultant
power spectra possess significantly lower variances 
than those of traditional, non-adaptive estimators.  The sine tapers also provide
useful spectral leakage suppression.  Resolution and uncertainty can
be estimated from the number of degrees of freedom (twice the number
of tapers).

This technique is particularly suited to long time series, because
it demands only one numerical Fourier transform, and requires no
costly additional computation of taper functions, like the Slepian
functions.  It also avoids the degradation of the low-frequency
performance associated with record segmentation 
in Welch's method.
Above all, the adaptive process relieves the user of the need to set
a tuning parameter, such as time-bandwidth product or segment length,
that fixes frequency resolution for the entire frequency interval; instead
it provides frequency-dependent spectral resolution tailored to the
shape of the spectrum itself.

`psd` elegantly handles
spectra with large dynamic range and mixed-bandwidth features|features
typically found in geophysical datasets.  

## How to Cite

Bob and Andy have a [paper in Computers & Geosciences][1]
to accompany this software [(download a pdf, 1MB)][pdf]; it describes the theory behind
the estimation process, and how we apply it in practice.
If you find `psd` useful in your research, we kindly request
you cite our paper. See also:

    citation("psd")

## Getting Started

You can to install the package and it's dependencies
with [CRAN](https://cran.r-project.org/package=psd)
(from within the `R` environment):

    install.packages("psd")

then load the package library

    library(psd)

We have included a dataset to play with, namely `Tohoku`, which represents
recordings of high-frequency borehole strainmeter data during
teleseismic waves from the 2011 Mw 9.0 Tohoku 
earthquake ([original data source](https://goo.gl/Gx7Ww)).
Access and inspect these data with:

    data(Tohoku)
    print(str(Tohoku))

The 'preseismic' data has interesting spectral features, so we
subset it, and analyze the areal strain (the change in borehole
diameter):

    Dat <- subset(Tohoku, epoch=="preseismic")
    Areal <- ts(Dat$areal)

For the purposes of improving the accuracy of the spectrum, we remove a linear trend:

    Dat <- prewhiten(Areal, plot=FALSE)

Now we can calculate the adaptive PSD:

    mtpsd <- pspectrum(Dat[['prew_lm']], plot=TRUE)
    print(class(mtpsd))

In the previous example the `plot=TRUE` flag produces a comparison with a basic periodogram, but
we can also visualize the spectrum with builtin plotting methods:

    plot(mtpsd, log="dB")

The spectral uncertainty can be easily calculated:

    sprop <- spectral_properties(mtpsd)
    with(sprop, {
        plot(taper/max(taper), type="h", ylim=c(0,2), col="dark grey") 
        lines(stderr.chi.lower)
        lines(stderr.chi.upper)
    })

### Installing the Development Version

Should you wish to install the development version
of this software, the [remotes][2] library
will be useful:

    library(remotes)
    install_github("abarbour/psd")

[1]: https://doi.org/10.1016/j.cageo.2013.09.015
[2]: https://cran.r-project.org/package=remotes
[pdf]: https://github.com/abarbour/psd/raw/master/paper/2014.barbour_parker.official.CAGEO3272.pdf
