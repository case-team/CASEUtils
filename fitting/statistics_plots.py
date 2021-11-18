#! /usr/bin/env python

import os, sys, getopt
import copy
import time
import math
import json
import optparse
from array import array
from ROOT import gROOT, gRandom, gStyle
from ROOT import TFile, TTree, TCut, TH1F, TH2F, TGraph, TGraph2D, TGraphErrors, TGraphAsymmErrors, TLine
from ROOT import TStyle, TCanvas, TPad
from ROOT import TLegend, TLatex, TText, TColor

## TODO remove hardcoded lower limit in exclusion limits plots

## currently saved in the counting experiment fit results
#results['chi2']
#results['ndof']
#results['fit_prob']
#results['nPars_QCD']
#results['signif']
#results['p_value']
#results['obs_lim_events']
#results['exp_lim_events']
#results['exp_lim_plus1sigma_events']
#results['exp_lim_minus1sigma_events']
#results['exp_lim_plus2sigma_events']
#results['exp_lim_minus2sigma_events']

usage = "usage: %prog [options]"
parser = optparse.OptionParser(usage)
parser.add_option("-i", "--inputFile", action="store", type="string", dest="inputFile", default="test_mjj_scan.json")
parser.add_option("-p", "--plotDir", action="store", type="string", dest="plotDir", default="test_scan_plots/")

(options, args) = parser.parse_args()
 
gStyle.SetOptStat(0)

gROOT.SetBatch(True)

input_file = options.inputFile
output_dir = options.plotDir

if not os.path.exists(output_dir):
    os.mkdir(output_dir)

def statistics_plots():

    with open(input_file) as f:
        fit_results = json.load(f)

    mass = list(fit_results.keys())
    mass.sort()
    print "mass =", mass
    
    limits_obs = [fit_results[m_label]['obs_lim_events'] for m_label in mass]
    limits_exp_plus2sigma = [fit_results[m_label]['exp_lim_plus2sigma_events'] for m_label in mass]
    limits_exp_minus2sigma = [fit_results[m_label]['exp_lim_minus2sigma_events'] for m_label in mass]

    min_lim = min(min(limits_obs), min(limits_exp_minus2sigma))
    max_lim = max(max(limits_obs), max(limits_exp_plus2sigma))

    Obs0s = TGraph()
    Exp0s = TGraph()
    Exp1s = TGraphAsymmErrors()
    Exp2s = TGraphAsymmErrors()
    Sign = TGraph()
    pVal = TGraph()

    for i, m_label in enumerate(mass):
        m = float(m_label)
        n = Exp0s.GetN()
        Obs0s.SetPoint(n, m, fit_results[m_label]['obs_lim_events'])
        Exp0s.SetPoint(n, m, fit_results[m_label]['exp_lim_events'])
        Exp1s.SetPoint(n, m, fit_results[m_label]['exp_lim_events'])
        Exp1s.SetPointError(n, 0., 0. ,(fit_results[m_label]['exp_lim_events']-fit_results[m_label]['exp_lim_minus1sigma_events']), (fit_results[m_label]['exp_lim_plus1sigma_events']-fit_results[m_label]['exp_lim_events']))
        Exp2s.SetPoint(n, m, fit_results[m_label]['exp_lim_events'])
        Exp2s.SetPointError(n, 0., 0., (fit_results[m_label]['exp_lim_events']-fit_results[m_label]['exp_lim_minus2sigma_events']), (fit_results[m_label]['exp_lim_plus2sigma_events']-fit_results[m_label]['exp_lim_events']))
        Sign.SetPoint(n, m, fit_results[m_label]['signif'])
        pVal.SetPoint(n, m, fit_results[m_label]['p_value'])

    Exp2s.SetLineWidth(2)
    Exp2s.SetLineStyle(1)
    Obs0s.SetLineWidth(3)
    Obs0s.SetMarkerStyle(0)
    Obs0s.SetLineColor(1)
    Exp0s.SetLineStyle(2)
    Exp0s.SetLineWidth(3)
    Exp1s.SetFillColor(417) #kGreen+1
    Exp1s.SetLineColor(417) #kGreen+1
    Exp2s.SetFillColor(800) #kOrange
    Exp2s.SetLineColor(800) #kOrange
    Exp2s.GetXaxis().SetTitle("m_{jj} (GeV)")
    Exp2s.GetXaxis().SetTitleSize(Exp2s.GetXaxis().GetTitleSize()*1.25)
    Exp2s.GetXaxis().SetNoExponent(True)
    Exp2s.GetXaxis().SetMoreLogLabels(True)
    Exp2s.GetYaxis().SetTitle("Signal Events")
    Exp2s.GetYaxis().SetTitleOffset(1.5)
    Exp2s.GetYaxis().SetNoExponent(True)
    Exp2s.GetYaxis().SetMoreLogLabels()

    Sign.SetLineWidth(2)
    Sign.SetLineColor(629)
    Sign.GetXaxis().SetTitle("m_{jj} (GeV)")
    Sign.GetXaxis().SetTitleSize(Sign.GetXaxis().GetTitleSize()*1.25)
    Sign.GetYaxis().SetTitle("Significance")

    pVal.SetLineWidth(2)
    pVal.SetLineColor(629)
    pVal.GetXaxis().SetTitle("m_{jj} (GeV)")
    pVal.GetXaxis().SetTitleSize(pVal.GetXaxis().GetTitleSize()*1.25)
    pVal.GetYaxis().SetTitle("local p-Value")

    # ---------- Exclusion Limits ----------
    c1 = TCanvas("c1", "Exclusion Limits", 800, 600)
    c1.cd()
    #SetPad(c1.GetPad(0))
    c1.GetPad(0).SetTopMargin(0.06)
    c1.GetPad(0).SetRightMargin(0.05)
    c1.GetPad(0).SetLeftMargin(0.12)
    c1.GetPad(0).SetTicks(1, 1)
    #c1.GetPad(0).SetGridx()
    #c1.GetPad(0).SetGridy()
    c1.GetPad(0).SetLogy()
    Exp2s.Draw("A3")
    Exp1s.Draw("SAME, 3")
    Exp0s.Draw("SAME, L")
    Obs0s.Draw("SAME, L")
    #setHistStyle(Exp2s)
    Exp2s.GetXaxis().SetTitleSize(0.050)
    Exp2s.GetYaxis().SetTitleSize(0.050)
    Exp2s.GetXaxis().SetLabelSize(0.045)
    Exp2s.GetYaxis().SetLabelSize(0.045)
    Exp2s.GetXaxis().SetTitleOffset(0.90)
    Exp2s.GetYaxis().SetTitleOffset(1.25)
    Exp2s.GetYaxis().SetMoreLogLabels(True)
    Exp2s.GetYaxis().SetNoExponent(True)
    Exp2s.GetYaxis().SetRangeUser(20, max_lim*1.1) ## FIXME hardcoded lower limit FIXME
    Exp2s.GetXaxis().SetRangeUser(float(mass[0])-50, float(mass[-1])+50)

    # legend
    top = 0.9
    nitems = 4

    leg = TLegend(0.63, top-nitems*0.3/5., 0.96, top)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0) #1001
    leg.SetFillColor(0)
    leg.SetHeader("95% CL upper limits")
    leg.AddEntry(Obs0s, "Observed", "l")
    leg.AddEntry(Exp0s, "Expected", "l")
    leg.AddEntry(Exp1s, "#pm 1 std. deviation", "f")
    leg.AddEntry(Exp2s, "#pm 2 std. deviation", "f")
    leg.Draw()
    latex = TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.045)
    latex.SetTextFont(42)

    leg2 = TLegend(0.12, 0.225-2*0.25/5., 0.65, 0.225)
    leg2.SetBorderSize(0)
    leg2.SetFillStyle(0) #1001
    leg2.SetFillColor(0)
    c1.GetPad(0).RedrawAxis()

    leg2.Draw()
    Obs0s.Draw("SAME, L")
    c1.GetPad(0).Update()

    if not gROOT.IsBatch(): raw_input("Press Enter to continue...")

    c1.Print(os.path.join(output_dir, "exclusion_limits.pdf"))
    c1.Print(os.path.join(output_dir, "exclusion_limits.png"))

    # ---------- Significance ----------
    c2 = TCanvas("c2", "Significance", 800, 600)
    c2.cd()
    c2.GetPad(0).SetTopMargin(0.06)
    c2.GetPad(0).SetRightMargin(0.05)
    c2.GetPad(0).SetTicks(1, 1)
    c2.GetPad(0).SetGridx()
    c2.GetPad(0).SetGridy()
    Sign.GetYaxis().SetRangeUser(0., 5.)
    Sign.GetXaxis().SetRangeUser(float(mass[0])-50, float(mass[-1])+50)
    Sign.Draw("AL3")
    c2.Print(os.path.join(output_dir, "significance.pdf"))
    c2.Print(os.path.join(output_dir, "significance.png"))

    # ---------- p-Value ----------
    c3 = TCanvas("c3", "p-Value", 800, 600)
    c3.cd()
    c3.GetPad(0).SetTopMargin(0.06)
    c3.GetPad(0).SetRightMargin(0.05)
    c3.GetPad(0).SetTicks(1, 1)
    c3.GetPad(0).SetGridx()
    c3.GetPad(0).SetGridy()
    c3.GetPad(0).SetLogy()
    pVal.Draw("AL3")
    pVal.GetYaxis().SetRangeUser(2.e-7, 0.5)

    ci = [1., 0.317310508, 0.045500264, 0.002699796, 0.00006334, 0.000000573303, 0.000000001973]
    line = TLine()
    line.SetLineColor(922)
    line.SetLineStyle(7)
    text = TLatex()
    text.SetTextColor(922)
    text.SetTextSize(0.025)
    text.SetTextAlign(12)
    for i in range(1, len(ci)-1):
        line.DrawLine(pVal.GetXaxis().GetXmin(), ci[i]/2, pVal.GetXaxis().GetXmax(), ci[i]/2);
        text.DrawLatex(pVal.GetXaxis().GetXmax()*1.01, ci[i]/2, "%d #sigma" % i);

    c3.Print(os.path.join(output_dir, "p_values.pdf"))
    c3.Print(os.path.join(output_dir, "p_values.png"))


if __name__ == "__main__":
    statistics_plots()
