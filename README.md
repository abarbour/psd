README for rlpSpec

Creation date:  Wed Jan 11 14:27:31 PST 2012

data/           datasets 
rsrc/           rfunctions

This work is a port of Robert Parker's (Matlab version of) PSD*

From the site:

"""
	P.S.D. stands for the Power Spectral Density, often simply
    referred to as the power spectrum, a property of stationary
    stochastic time series that usually provides remarkable insights
    into the physical processes in question. Loosely speaking, the
    PSD gives the energy density in each frequency band; it is a
    statistic of the underlying random process, and so it can only
    be estimated from a given realization, that is, a given finite
    collection of data. There are probably dozens of different
    estimation algorithms for the job. The program psd provides a
    small number of easily used methods, with the intent of providing
    a good spectral estimate without requiring the user to set too
    many tuning parameters.

	Specifically, psd is based on the sine
    multitaper approach of Riedel and Sidorenko to spectral estimation.
    All spectral estimators smooth the spectrum in order to reduce
    variance; the S-R method applies smoothing that varies with
    frequency to minimize the sum of variance and bias at each
    point.
"""

* http://igppweb.ucsd.edu/~parker/Software/index.html
