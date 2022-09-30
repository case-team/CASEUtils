from H5_maker import * 
from ttbar_h5_maker import *
import glob


parser = OptionParser()
parser.add_option("-f", "--flag", dest = "flag", default = -1234, type=int, help="Flag to label what type of process this is (QCD, ttbar, signal, etc)")
parser.add_option("--sys", default = False, action = 'store_true', help="Add additional info the h5's for systematics")
parser.add_option("--top_ptrw", default = False, action = 'store_true', help="Include ttbar top pt reweighting factors")
parser.add_option("--ttbar", default = False, action = 'store_true', help="Semi leptonic ttbar version of h5 maker (different preselection)")
parser.add_option("--sample_type", default = "MC", help="MC or data")
parser.add_option("-i", "--input", dest = "fin", default = '', help="Input file name")
parser.add_option("-o", "--output", dest = "fout", default = 'test.h5', help="Output file name")
parser.add_option("-j", "--json", default = '', help="Json file name")
parser.add_option("-y", "--year", type=int, default = 2016, help="Year the sample corresponds to")
parser.add_option("-n", "--nEvents",  type=int, default = -1, help="Maximum number of events to output (-1 to run over whole file)")

options, args = parser.parse_args()

if(options.flag == -1234):
    print("No --flag option set. You must specify what type of process this is! \n" )
    exit(1)

if(".root" in options.fin): input_files = [options.fin]
else: input_files = glob.glob(options.fin + "*.root")

if(options.ttbar):
    NanoReader_TTbar(options.flag, inputFileNames = input_files, outputFileName = options.fout, json = options.json, year = str(options.year), 
        nEventsMax = options.nEvents, include_systematics = options.sys, do_top_ptrw = options.top_ptrw, sampleType = options.sample_type)
else:

    NanoReader(options.flag, inputFileNames = input_files, outputFileName = options.fout, json = options.json, year = options.year, 
        nEventsMax = options.nEvents, include_systematics = options.sys, do_top_ptrw = options.top_ptrw)

