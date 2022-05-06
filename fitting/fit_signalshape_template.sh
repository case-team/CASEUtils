#!/bin/bash

# Example case: we have MC signal samples at masses of 1500, 2500 and 3500 GeV
# H5s must are expected to have been created from CASE h5_maker.py script
# --dcbModel option causes the fit model to be a double crystal ball shape instead of the gaussian core + single crystal ball previously used
python fit_signalshapes.py -M 1500 2500 3500 -i /path/to/signal_m_1500.h5 /path/to/signal_m_2500.h5 /path/to/signal_m_3500.h5 --dcbModel -o /output/dir/

# using the fits, interpolate into different mass points of your choice (given in --masses option)
# "signal_model" currently only supported when set to "graviton"
# The input (-i option) shoud point to output "full_fit.root" file in -o directory of fit_signalshapes.py script
# Example here: interpolate into masses for a scan from 1500 GeV to 3500 GeV using 100 GeV steps
# In the output directory, ROOT files containing the interpolated parameters can be found, which can be used with the normal "dijetfit.py" script. Make sure to also choose the
# --dcb-model option there if you did so when fitting the signalshapes.
python interpolation.py -s signal_model -i /output/dir/full_fit.root --masses 1500 1600 1700 1800 1900 2000 2100 2200 2300 2400 2500 2600 2700 2800 2900 3000 3100 3200 3300 3400 3500 -o /some/output/directory


