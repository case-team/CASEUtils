from H5_maker import * 
from ttbar_h5_maker import *
from gen_ttbar_h5_maker import *
from gen_h5_maker import *
from utils import *
import glob


parser = OptionParser()
parser.add_option("-f", "--flag", dest = "flag", default = -1234, type=int, help="Flag to label what type of process this is (QCD, ttbar, signal, etc)")
parser.add_option("--sys", default = False, action = 'store_true', help="Add additional info the h5's for systematics")
parser.add_option("--top_ptrw", default = False, action = 'store_true', help="Include ttbar top pt reweighting factors")
parser.add_option("--ttbar", default = False, action = 'store_true', help="Semi leptonic ttbar version of h5 maker (different preselection)")
parser.add_option("--tW", default = False, action = 'store_true', help="tW sample")
parser.add_option("--herwig", default = False, action = 'store_true', help="Semi leptonic ttbar version of h5 maker (different preselection)")
parser.add_option("--gen", default = False, action = 'store_true', help="Gen level ttbar version of h5 maker (different preselection)")
parser.add_option("--sample_type", default = "MC", help="MC or data")
parser.add_option("-i", "--input", dest = "fin", default = '', help="Input file name")
parser.add_option("-o", "--output", dest = "fout", default = 'test.h5', help="Output file name")
parser.add_option("-j", "--json", default = '', help="Json file name")
parser.add_option("-y", "--year", type=str, default = "2016", help="Year the sample corresponds to")
parser.add_option("--nJobs", type =int, default =0,  help="Year the sample corresponds to")
parser.add_option("--iJob", type =int, default =0,  help="Year the sample corresponds to")
parser.add_option("-n", "--nEvents",  type=int, default = -1, help="Maximum number of events to output (-1 to run over whole file)")

options, args = parser.parse_args()

if(options.flag == -1234):
    print("No --flag option set. You must specify what type of process this is! \n" )
    exit(1)

if(".root" in options.fin): input_files = [options.fin]
else: input_files = get_file_list(options.fin)

if(options.nJobs > 0):
    flist_new = [input_files[i] for i in range(len(input_files)) if (i % options.nJobs) == options.iJob]
    input_files = flist_new

print(input_files)


if(options.gen and options.ttbar):
    NanoReader_Gen_TTbar(options.flag, inputFileNames = input_files, outputFileName = options.fout, json = options.json, year = options.year, 
        nEventsMax = options.nEvents, include_systematics = options.sys, do_top_ptrw = options.top_ptrw, sampleType = options.sample_type, herwig = options.herwig)

elif(options.gen):
    NanoReader_Gen(options.flag, inputFileNames = input_files, outputFileName = options.fout, json = options.json, year = options.year, 
        nEventsMax = options.nEvents, include_systematics = options.sys, do_top_ptrw = options.top_ptrw, sampleType = options.sample_type, herwig = options.herwig)
elif(options.ttbar):
    NanoReader_TTbar(options.flag, inputFileNames = input_files, outputFileName = options.fout, json = options.json, year = options.year, 
        nEventsMax = options.nEvents, include_systematics = options.sys, do_top_ptrw = options.top_ptrw, sampleType = options.sample_type, tW = options.tW, herwig = options.herwig)
else:

    NanoReader(options.flag, inputFileNames = input_files, outputFileName = options.fout, json = options.json, year = options.year, 
        nEventsMax = options.nEvents, include_systematics = options.sys, do_top_ptrw = options.top_ptrw)

