import matplotlib.pyplot as plt
import numpy as np
import h5py
import json
import os
import optparse

from dijetfit_counting import dijetfit

parser = optparse.OptionParser()
parser.add_option("-i","--inputFile",dest="inputFile",default='fit_inputs/no_selection_03p.h5',help="input h5 file")
parser.add_option("-o","--outputFile",dest="outputFile",default='test_mjj_scan.json',help="Where to store the output (json file)")
parser.add_option("-p","--plotDir",dest="plotDir",default='scan_plots/',help="Where to store the fitting plots")

(options,args) = parser.parse_args()

## TODO find elegant way to make this command line input
mjj_fit_points = [1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500, 3600, 3700, 3800, 3900, 4000]

if not os.path.exists(options.plotDir):
    os.mkdir(options.plotDir)

class FitOptions:
    def __init__(self, mjj):
        self.mass = mjj
        self.inputFile = options.inputFile
        self.plotDir = os.path.join(options.plotDir, "mass_"+str(mjj))

results = {}

for mjj_point in mjj_fit_points:
    fit_options = FitOptions(mjj_point)
    results[mjj_point] = dijetfit(fit_options)

with open(options.outputFile, 'w') as f:
    json.dump(results, f)
