from Sampler import BlackBox, MaxMultiSampler
from optparse import OptionParser
import os
from os.path import join, dirname, exists
import numpy as np
import json

parser = OptionParser()
parser.add_option("-s", "--seed", type=int, default=123, help="Random seed")
parser.add_option("-i", "--input", type=str,
                  default="/eos/cms/store/group/phys_b2g/CASE/h5_files/UL/",
                  help="Input base directory.")
parser.add_option("-o", "--output", type=str, default="./mixed_years/",
                  help="Output directory.")
parser.add_option("-d", "--subdir", type=str, default="new/",
                  help=("Subdirectory after the year-specific folder "
                        "(<input base dir>/<year dir>/<subdir>/<actual h5 input samples>)."
                        "Provide empty string if subdir should be omitted."))
parser.add_option("-c", "--config", type=str,
                  default=join(dirname(__file__), "sample_signal_config.json"),
                  help="JSON config file containing zears and relative lumies.")
options, args = parser.parse_args()


base_dir = options.input
out_dir = options.output
subdir = options.subdir
with open(options.config) as json_file:
    config = json.load(json_file)
    years = config["years"]
    lumies = config["lumies"]
np.random.seed(options.seed)

assert len(years) > 1
for year in years:
    assert year in lumies.keys()

if not exists(out_dir):
    os.makedirs(out_dir)

h5_files = {}
for year in years:
    h5_files[year] = {x for x in os.listdir(join(base_dir, year, subdir))
                      if x.endswith(".h5")}

# Find signals present in all years
signal_file_list = list(set.intersection(*h5_files.values()))
signal_file_list.sort()

for signal_name in signal_file_list:
    print("\nworking on: "+signal_name)
    sig_list = [join(base_dir, year, subdir, signal_name)
                for year in years]
    lumi_list = [lumies[year] for year in years]
    multi_sampler = MaxMultiSampler(sig_list, lumi_list)

    BB = BlackBox(multi_sampler.get_samplers(), keys=[], nBatches=1)
    f_out_name = join(out_dir, signal_name.replace(".h5", ""))
    BB.writeOut(f_out_name)
