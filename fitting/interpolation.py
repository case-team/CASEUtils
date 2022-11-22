import ROOT
from array import array
import os, sys, re, optparse, pickle, shutil, json, copy
from Utils import returnString
import argparse
from os.path import join
ROOT.gROOT.SetBatch(True)


parser = argparse.ArgumentParser(description='Interpolate double crystal ball template from signal template fits of masses where we have MC samples available.')
parser.add_argument("-s", "--signal-model", dest="signal_model", required=True, help="Choose which signal model is used for interpolation.")
parser.add_argument("-i", "--input", dest="in_file", required=True, help="Input ROOT file containing the fitted values of the signal template fits.")
parser.add_argument("--masses", dest="masses", nargs="+", type=float, help="List of target masses (in GeV) for which parameter template JSON files should be created", default=None)
parser.add_argument("-o","--output",dest="output",help="Output directory", default='.')
parser.add_argument("-m","--min",dest="min",type=float, help="minimum x",default=0)
parser.add_argument("-M","--max",dest="max",type=float, help="maximum x",default=0)


args = parser.parse_args()
### use same settings as in diboson search
all_graphs = {"diboson" : "mean:spline,sigma:spline,alpha:pol1,sign:pol1,alpha2:pol2,sign2:pol2"}


#define output dictionary
rootFile = ROOT.TFile(args.in_file)
graphStr = all_graphs[args.signal_model].split(',')
parameterization = {}
mass_templates = {}
if args.masses is not None:
    mass_templates = dict([(mass, {}) for mass in args.masses])

ff=ROOT.TFile(join(args.output, args.signal_model + "_interpolation.root"),"RECREATE")
ff.cd()
print("graphStr: ", graphStr)
for string in graphStr:
    comps =string.split(':')
    graph=rootFile.Get(comps[0])
    #import pdb; pdb.set_trace()

    if "pol" in comps[1]:
        func=ROOT.TF1(comps[0]+"_func",comps[1],0,13000)
        #if (comps[0]=="sign") and (args.signal_model=="graviton"):
        #    func.SetParameter(0,0)
        #func=ROOT.TF1(comps[0]+"_func","[0]-[1]*x" ,0,13000)
    elif  comps[1]=="llog":
        func=ROOT.TF1(comps[0]+"_func","[0]+[1]*log(x)",1,13000)
        func.SetParameters(1,1)
    elif  "laur" in comps[1]:
        order=int(comps[1].split("laur")[1])
        st='0'
        for i in range(0,order):
            st=st+"+["+str(i)+"]"+"/x^"+str(i)
        print('Laurent String: ', st)
        func=ROOT.TF1(comps[0]+"_func",st,1,13000)
        for i in range(0,order):
            func.SetParameter(i,0)
    elif comps[1]=="sqrt":
        func = ROOT.TF1(comps[0]+"_func","[0]+[1]*sqrt(x)",1,13000)
        st= "work in progress"
    elif comps[1]=="1/sqrt":
        func = ROOT.TF1(comps[0]+"_func","[0]+[1]/sqrt(x)",1,13000)
        st = "work in progress"
    elif comps[1]=="spline":
        print("fit spline")
        func = ROOT.TSpline3(comps[0],graph)

    print("function: ", func)
    if comps[1]=="spline":
        print("Spline used, no need to fit it")
    else:
        print(comps[1])
        print('fit function ' + func.GetName())
        graph.Fit(func, "", "", args.min, args.max)
        graph.Fit(func, "", "", args.min, args.max)
        graph.Fit(func, "", "", args.min, args.max)
    
    if args.masses is not None:
        for mass in args.masses:
            mass_templates[mass][comps[0]] = func.Eval(mass)
            
        
    parameterization[comps[0]]=returnString(func,comps[1])

    c = ROOT.TCanvas()
    c.SetRightMargin(0.2)
    graph.Draw()
    func.SetLineColor(ROOT.kRed)
    func.Draw("lcsame")
    name = args.signal_model + "_interpolation_" + comps[0]
    c.SaveAs(join(args.output, name + ".png"))
    c.SaveAs(join(args.output, name + ".pdf"))
    c.SaveAs(join(args.output, name+".C"))
    graph.Write(comps[0])
    func.Write(comps[0]+"_func")
ff.Close()
f = open(join(args.output, args.signal_model + "_interpolation.json"), "w")
json.dump(parameterization, f)
f.close()
if args.masses is not None:
    for mass in args.masses:
        
        graphs = {}
        outfile_name = join(
            args.output,
            "{}_interpolation_M{}.root".format(args.signal_model, mass))

        outfile = ROOT.TFile(outfile_name, "RECREATE")
        
        for key, value in mass_templates[mass].iteritems():
            graphs[key] = ROOT.TGraphErrors()
            graphs[key].SetPoint(0, mass, value)

        outfile.cd()
        for key, value in mass_templates[mass].iteritems():
            graphs[key].Write(key)
        outfile.Close()

        with open(join(args.output, "{}_interpolation_M{}.json".format(args.signal_model, mass)),"w") as f:
            json.dump(mass_templates[mass], f)

print("Done!")
