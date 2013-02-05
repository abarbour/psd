#!/bin/sh
set -e
#
# generate flow diagrams:
# http://yuml.me/diagram/scruffy/class/samples
#
# %3E ">"
# %3C "<"
# %20 " "
##
YUML="http://yuml.me/diagram/scruffy/class/"
YUML="http://yuml.me/diagram/class/"

## call notes
cat << XXX | yuml.fawk | sh
[note: PSPECTRUM controls the adaptive process. It repeats NITER times.{bg:darkorange}]
[note: With each iteration RIEDSID optimizes the number of tapers depending on the shape of the current spectrum.{bg:cornsilk}]
XXX
mv yuml.png yuml_n.png

## call graph
cat << XXX | yuml.fawk | sh
[PSPECTRUM]-.-%3E[RIEDSID]
[RIEDSID]-.-%3E[PSDCORE]
[PSDCORE]NITER times%20-.-%3E[PSPECTRUM]
[PILOT_SPEC]%3C-%3E[PSDCORE]
[PILOT_SPEC]%3C-%3E[PSPECTRUM]
[CONSTRAIN_TAPERS]%3C-%3E[RIEDSID]
XXX
mv yuml.png yuml_d.png

open yuml_?.png
