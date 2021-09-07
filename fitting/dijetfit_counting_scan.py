import matplotlib.pyplot as plt
import numpy as np
import h5py
import json
import os
import optparse
import yaml
from dijetfit_counting import dijetfit

parser = optparse.OptionParser()
parser.add_option("-i","--inputFile",dest="inputFile",default='example_scan_input.yml',help="input YAML file mapping mass values (in GeV) to input hdf5 files.")
parser.add_option("-o","--outputFile",dest="outputFile",default='test_mjj_scan.json',help="Where to store the output (json file)")
parser.add_option("-p","--plotDir",dest="plotDir",default='scan_plots/',help="Where to store the fitting plots")
parser.add_option("-b", "--blinded", dest="blinded", action="store_true", default=False, help="Blinding the signal region for the fit.")
(options,args) = parser.parse_args()


with open(options.inputFile, 'r') as stream:
    input_files = yaml.safe_load(stream)

print "mass values =", input_files.keys()

if not os.path.exists(options.plotDir):
    os.mkdir(options.plotDir)

class FitOptions:
    def __init__(self, mjj, input_file, blinded):
        self.mass = mjj
        self.inputFile = input_file
        self.plotDir = os.path.join(options.plotDir, "mass_"+str(mjj))
        self.blinded = blinded

results = {}

for mjj_point in input_files.keys():
    print "m =", mjj_point
    fit_options = FitOptions(mjj_point, input_files[mjj_point], options.blinded)
    results[mjj_point] = dijetfit(fit_options)

with open(options.outputFile, 'w') as f:
    json.dump(results, f)
