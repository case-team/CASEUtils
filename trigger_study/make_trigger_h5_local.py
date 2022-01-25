import os
from trigger_H5_maker import NanoReader
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-i", "--input", dest="fin", default='', help="Input file name (text file containing list of root files or single root file)")
parser.add_option("-o", "--output", dest="fout", default='test.h5', help="Output file name")
parser.add_option("-y", "--year", dest="year", type=int, default=2016, help="Year the sample corresponds to")
parser.add_option("-n", "--nEvents", dest="nEvents",  type=int, default=-1, help="Maximum number of events to output (-1 to run over whole file)")
parser.add_option("-f", "--filenumber", dest="file_number", type=int, default=-1, help="The number of the file in the given infile list to run on.")

options, args = parser.parse_args()

if options.year==2016:
    json = "jsons/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt"
elif options.year==2017:
    json = "jsons/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt"
elif options.year==2018:
    json = "jsons/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt"
else:
    raise RuntimeError("Year unknown!")

if options.fin.endswith(".txt"):
    print "input is a text file of input files:", options.fin
    clean_file_list = []
    raw_file_list = open(options.fin, 'r').readlines()
    for entry in raw_file_list:
        if entry.startswith("#"):
            continue
        filepath = entry.replace('\n','')
        clean_file_list.append(filepath)
    if options.file_number==-1:
        print "running full input file list..."
        input_files = clean_file_list
    else:
        print "running full input file nr {} from list...".format(options.file_number)
        input_files = [clean_file_list[options.file_number]]

elif options.fin.endswith(".root"):
    print "running on single provided input file..."
    input_files = [options.fin]

else:
    raise RuntimeError("Input file format unknown!!")

if os.path.isdir(options.fout):
    outfile = os.path.join(options.fout, "ntuple{}.h5".format("_"+str(options.file_number) if options.file_number!=-1 else "_full"))
else:
    outfile = options.fout

NanoReader(inputFileNames=input_files, outputFileName=outfile, json=json, year=options.year, nEventsMax=options.nEvents)

