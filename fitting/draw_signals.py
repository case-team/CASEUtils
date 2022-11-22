import ROOT
import json
from numpy import random
from array import array
import os,sys,commands
from Fitter import *
from DataCardMaker import *
ROOT.gROOT.SetBatch(True)

f_shape = ""


f_sig1 = "XToYY_shapes/sig_fit_3000.root"
f_sig2 = "XToYY_TNT_shapes/sig_fit_3000.root"
f_sig3 = "XToYY_TNT_lowsig_shapes/sig_fit_3000.root"

var = ROOT.RooRealVar("mjj", "mjj", 2500., 3500.)
card = DataCardMaker('sig_test')
card.addDCBSignalShape('sig1', 'mjj', f_sig1, {'CMS_scale_j': 1.0}, {'CMS_res_j': 1.0})
card.addDCBSignalShape('sig2', 'mjj', f_sig2, {'CMS_scale_j': 1.0}, {'CMS_res_j': 1.0})
card.addDCBSignalShape('sig3', 'mjj', f_sig3, {'CMS_scale_j': 1.0}, {'CMS_res_j': 1.0})

frame = var.frame()
#dcb1 = get_signal_shape(f_sig)
c1 =ROOT.TCanvas("c1","",800,800)
card.w.Print()
card.w.ls()
sig1 = card.w.pdf('sig1_JJ_sig_test')
sig2 = card.w.pdf('sig2_JJ_sig_test')
sig3 = card.w.pdf('sig3_JJ_sig_test')

sig1.plotOn(frame, ROOT.RooFit.LineColor(ROOT.kRed + 1), ROOT.RooFit.Name("No Cut") )
sig2.plotOn(frame, ROOT.RooFit.LineColor(ROOT.kBlue ), ROOT.RooFit.Name("TNT Cut (5 #sigma inj.)") )
sig3.plotOn(frame, ROOT.RooFit.LineColor(ROOT.kGreen ), ROOT.RooFit.Name("TNT Cut (0 #sigma inj.)") )
frame.SetTitle("")

frame.Draw()
legend = ROOT.TLegend(0.3, 0.3)
l1 = ROOT.TLine(0,0,0,0)
l2 = ROOT.TLine(0,0,0,0)
l3 = ROOT.TLine(0,0,0,0)
l1.SetLineColor(ROOT.kRed + 1)
l2.SetLineColor(ROOT.kBlue)
l3.SetLineColor(ROOT.kGreen)
legend.AddEntry(l1,"No Cut",  "l") 
legend.AddEntry(l2, "TNT Cut (5 #sigma)", "l") 
legend.AddEntry(l3, "TNT Cut (0 #sigma)", "l") 
legend.Draw()
c1.Print("sig_fit_test.png")

