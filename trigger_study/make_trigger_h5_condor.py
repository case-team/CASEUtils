#! /usr/bin/env python

###
### Macro for submitting trigger ntuple production to HTCondor.
###

import sys
import os 
from commands import getoutput
from argparse import ArgumentParser
import json
from math import ceil

EXE_PATH = os.path.dirname(os.path.abspath(__file__))
CONDOR_LOG_PATH = os.path.join(EXE_PATH, "condor_logs")
NANO_TOOLS_PATH = "." # overwritten by cmd args
CERTIFICATE_PATH = "" # overwritten by cmd args

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('-q', '--queue',   dest='queue',
        choices=['espresso', 'microcentury', 'longlunch', 'workday', 'tomorrow', 'testmatch', 'nextweek'],
        type=str, default='tomorrow', action='store', help="select queue for submission" )
    parser.add_argument('-o', '--outdir',  dest='outdir', type=str,
        default='/eos/user/m/msommerh/CASE/trigger_samples', action='store', help="set out file location." )
    parser.add_argument('-y', '--year', dest='year', type=str, default='2016',
        action='store', help="set year." )
    parser.add_argument('-rs', '--resubmit',   dest='resubmit', type=str, action='store',
        default="", help="Indicate file containing titles of the sample for resubmission." )
    parser.add_argument('-rf', '--resubmit_file',   dest='resubmit_file', type=int, action='store',
        default=-1, help="Indicate output file number for single resubmission." )
    parser.add_argument('-n', '--nano_tools_path',  dest='nano_tools_path', type=str,
        default='.', action='store',
        help="path to the src directory of a CMSSW environemnt where the NanoAOD tools are installed" )
    parser.add_argument('-c', '--certificate',  dest='certificate', type=str,
        default='/tmp/x509up_u110866', action='store',
        help="path to GRID certificate" )
    args = parser.parse_args()

    #jobflavour = 'espresso' #max 30min
    #jobflavour = 'microcentury' #max 1h
    #jobflavour = 'longlunch' #max 2h
    #jobflavour = 'workday' #max 8h
    #jobflavour = 'tomorrow' #max 1d
    #jobflavour = 'testmatch' #max 3d

    NANO_TOOLS_PATH = args.nano_tools_path
    CERTIFICATE_PATH = args.certificate  

else:
  args = None

def getFileListDAS(dataset):
        """Get list of files from DAS."""
        dataset  = dataset.replace('__','/')
        instance = 'prod/global'
        if 'USER' in dataset:
            instance = 'prod/phys03'
        cmd = 'dasgoclient -query="file dataset={} instance={}"'.format(dataset, instance)
        cmd_out  = getoutput( cmd )
        tmpList  = cmd_out.split(os.linesep)
        filelist = [ ]
        for line in tmpList:
          if '.root' in line:
            #filelist.append("root://xrootd-cms.infn.it/"+line)
            filelist.append("root://cms-xrd-global.cern.ch/"+line)
        return filelist


def submitJobs(title, infiles, outdir, jobflavour):
    path = EXE_PATH
    outpath = os.path.join(outdir, title)
    if not os.path.exists(outpath):
        os.makedirs(outpath)
        print "Directory "+outpath+" created."
    else:
        print "Directory "+outpath+" already exists."

    if not os.path.exists(CONDOR_LOG_PATH):
        os.makedirs(CONDOR_LOG_PATH)
        print "Directory "+CONDOR_LOG_PATH+" created."
    else:
        print "Directory "+CONDOR_LOG_PATH+" already exists."

    workdir = os.path.join(CONDOR_LOG_PATH, "tmp_"+title)
    if not os.path.exists(workdir):
        os.makedirs(workdir)
        print "Directory "+workdir+" created."
    else:
        print "Directory "+workdir+" already exists."
    os.chdir(workdir)

    #write executable file for submission
    with open('job.sh', 'w') as fout:
        fout.write("#!/bin/sh\n")
        fout.write("echo\n")
        fout.write("echo\n")
        fout.write("echo 'START---------------'\n")
        fout.write("echo 'WORKDIR ' ${PWD}\n")

        fout.write("cd "+str(NANO_TOOLS_PATH)+"\n")
        fout.write("export SCRAM_ARCH=slc6_amd64_gcc700\n" )
        fout.write("eval `scram runtime -sh`\n")
        fout.write("echo 'cmssw release = ' $CMSSW_BASE\n")
        #fout.write("source ./setupEnv.sh\n")
        fout.write("NANOPATH=\"$( dirname \"$(readlink -f -- \"$0\")\" )\"\n")
        fout.write("export PYTHONPATH=$NANOPATH:$PYTHONPATH\n")
 
        fout.write("cd "+str(EXE_PATH)+"\n")

        fout.write("export X509_USER_PROXY="+CERTIFICATE_PATH+"\n")
        fout.write("use_x509userproxy=true\n")

        fout.write("########## input arguments ##########\n")
        fout.write("file_nr=$1\n")
        fout.write("#####################################\n")
        fout.write("echo\n")
        fout.write("python make_trigger_h5_local.py -i {} -o {} -y {} -f {}\n".format(infiles, outpath, args.year, args.resubmit_file if args.resubmit_file!=-1 else r"$file_nr"))
        fout.write("echo 'STOP---------------'\n")
        fout.write("echo\n")
        fout.write("echo\n")

    #submit job
    os.system("chmod 755 job.sh")
    if args.resubmit_file==-1:
        os.system("mv job.sh "+title+".sh")
        makeSubmitFileCondor(title+".sh", title, jobflavour, path+"/"+infiles)
    else:
        os.system("mv job.sh "+title+"_"+str(args.resubmit_file)+".sh")
        makeSubmitFileCondor(title+"_"+str(args.resubmit_file)+".sh", title, jobflavour, path+"/"+infiles)
    os.system("condor_submit submit.sub")
    print "job submitted"
    os.chdir(path)

def makeSubmitFileCondor(exe, jobname, jobflavour, infiles):
    print "make options file for condor job submission"
    submitfile = open("submit.sub", "w")
    submitfile.write("executable            = "+exe+"\n")
    submitfile.write("arguments             = $(ProcId)\n")
    submitfile.write("output                = "+jobname+".$(ClusterId).$(ProcId).out\n")
    submitfile.write("error                 = "+jobname+".$(ClusterId).$(ProcId).err\n")
    submitfile.write("log                   = "+jobname+".$(ClusterId).log\n")
    submitfile.write('+JobFlavour           = "'+jobflavour+'"\n')
    if args.resubmit_file==-1:
        file_list = []
        file_content = open(infiles, 'r').readlines()
        for entry in file_content:
            if not entry.startswith("#"):
                filepath = entry.replace('\n','')
                file_list.append(filepath)
        nJobs = int(len(file_list))
        submitfile.write("queue {}".format(nJobs))
    else:
        submitfile.write("queue") 
    submitfile.close()

def main():

        outdir = args.outdir
        
        if args.resubmit_file != -1: print "only submitting output file nr", args.resubmit_file

        ## load data sets from file
        if args.year in ['2016','2017','2018']:
                data_set_file = 'samples/samples_{}_{}.json'.format("SingleMuon", args.year)
        else:
                print "Unknown year. Abort submission!!"
                sys.exit()
    
        with open(data_set_file, 'r') as json_file:
                data_sets = json.load(json_file)
    
        if args.resubmit != "":
                rs_titles = open(args.resubmit, 'r').readlines()
                for n, entry in enumerate(rs_titles):
                        rs_titles[n] = entry.replace("\n", "")
                print "resubmitting the following samples:", rs_titles
    
        for title in data_sets.keys():
    
                if args.resubmit != "":
                        if title not in rs_titles: continue 
    
                data_set = data_sets[title]
                infiles = "filelists/"+title+".txt"
    
                ## create filelist from DAS
                if not os.path.exists("filelists"):
                    os.makedirs("filelists")                    
                txtfile = open(infiles, "w")
                txtfile.write("# created from {}\n".format(data_set))
                filelist = getFileListDAS(data_set)
                for entry in filelist:
                      txtfile.write(entry+"\n")
                txtfile.close()
    
                ## submit job
                submitJobs(title, infiles, outdir, args.queue) 
        print
        print
        print "your jobs:"
        os.system("condor_q")
        userName=os.environ['USER']
        print
        print 'Done submitting jobs!'
        print

if __name__ == "__main__":
    print
    main()
    print "done"
