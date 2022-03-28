from H5_maker import * 


parser = OptionParser()
parser.add_option("-f", "--flag", dest = "flag", default = -1234, type=int, help="Flag to label what type of process this is (QCD, ttbar, signal, etc)")
parser.add_option("--sys", default = False, action = 'store_true', help="Add additional info the h5's for systematics")
parser.add_option("--top_ptrw", default = False, action = 'store_true', help="Include ttbar top pt reweighting factors")
parser.add_option("-i", "--input", dest = "fin", default = '', help="Input file name")
parser.add_option("-o", "--output", dest = "fout", default = 'test.h5', help="Output file name")
parser.add_option("-j", "--json", default = '', help="Json file name")
parser.add_option("-y", "--year", type=int, default = 2016, help="Year the sample corresponds to")
parser.add_option("-n", "--nEvents",  type=int, default = -1, help="Maximum number of events to output (-1 to run over whole file)")

options, args = parser.parse_args()

if(options.flag == -1234):
    print("No --flag option set. You must specify what type of process this is! \n" )
    exit(1)

NanoReader(options.flag, inputFileNames = [options.fin], outputFileName = options.fout, json = options.json, year = options.year, 
        nEventsMax = options.nEvents, include_systematics = options.sys, do_top_ptrw = options.top_ptrw)

