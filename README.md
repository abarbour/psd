#psd: Adaptive, sine multitaper power spectral density estimation in R#

##by Andrew J Barbour and Robert L Parker##

###Description###
This is an `R`
package for computing univariate power spectral density
estimates with little or no tuning effort.
We employ sine multitapers, allowing the number to vary with frequency
in order to reduce mean square error, the sum of squared bias and
variance, at each point.  The approximate criterion of
[Riedel and Sidorenko (1995)](http://dx.doi.org/10.1109/78.365298)
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

##Getting Started##

Firstly you'll need to install the package and it's dependencies
from [CRAN](http://cran.r-project.org/web/packages/psd/)
(from within the `R` environment):

    install.packages("psd", dependencies=TRUE)

then load the package library

    library(psd)

We have included a dataset to play with, named `Tohoku`, which represents
recordings of
high-frequency borehole strainmeter data during
teleseismic waves from the 2011 Mw 9.0 Tohoku 
earthquake ([source](http://goo.gl/Gx7Ww)):

    data(Tohoku)
    print(str(Tohoku))

The 'preseismic' data has interesting spectral features, so
subset it, and use the areal strain:

    Dat <- subset(Tohoku, epoch=="preseismic")
    Areal <- ts(Dat$areal)

Remove a linear trend:

    Dat <- prewhiten(Areal, plot=FALSE)

And calculate the adaptive PSD:

    mtpsd <- pspectrum(Dat$prew_lm, plot=TRUE)
    print(class(mtpsd))

Visualize the spectrum with builtin methods:

    plot(mtpsd, log="dB")

and the spectral uncertainty:

    sprop <- spectral_properties(mtpsd)
    Ntap <- sprop$taper/max(sprop$taper)
    plot(Ntap, type="h", ylim=c(0,2), col="dark grey") 
    lines(sprop$stderr.chi.lower)
    lines(sprop$stderr.chi.upper)

