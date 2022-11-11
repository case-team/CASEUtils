# Dijet fit for CASE analysis


Based on previous VAE [version](https://github.com/jngadiub/VAEDijetFit)

### Prerequisites

Higgs combine tools is needed, either standalone or from cmssw. See [here](http://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/) for instructions to get latest version.

### How To

The main fitting code is `dijetfit.py`.
It takes as input an h5 file with a branch of 'mjj', the dijet mass of events
to be fit. 
One also must supply a signal shape, a ROOT file that was outputted from
`fit_signalshapes.py`  or `interpolation.py`.
There are several other options for the fit as well, which are relatively
self-explanatory.
The fit automatically performs an F-test to determine if 2,3 or 4 parameter QCD
fit function should be used before doing the signal + background fit.
The outputs of the fit are saved in a json.

`fit_signalshapes.py` takes as inputs lists of h5's (all which have an `mjj` branch) and 
a corresponding list of signal massses. By default the signal shape is only fit
in the range of +/- 20\% of the signal mass (this can be changed with an option).
By default it uses a single crystal ball with a Gaussian core but we have been
generally using the double crystall ball fit (`--dcb-model` option). 

`interpolation.py` is used to interpolate signal shapes to different dijet
masses.

Example commands to produce the interpolated signal shapes:

```
python fit_signalshapes.py -M 2000 3000 5000 --inputFiles XToYYprime_M2000_MY80_MYprime170_sig_shape.h5 XToYYprime_M3000_MY80_MYprime170_sig_shape.h5 XToYYprime_M5000_MY80_MYprime170_sig_shape.h5 --dcb-model --fitRange 0.4 -o XToYY_shapes/

 python interpolation.py -i XToYY_shapes/full_fit.root -s diboson -o XToYY_shapes --masses 1600 1700 1800 1900 2100 2200 2300 2400 2500 2600 2700 2800 2900 3100 3200 3300 3400 3500 3600 3700 3800 3900 4000 4100 4200 4300 4400 4500 4600 4700 4800 4900 5100 5200 5300 5400 5500 5600 5700 5800 5900
```

