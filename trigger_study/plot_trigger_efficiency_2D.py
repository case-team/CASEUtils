#! /usr/bin/env python

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
#from __future__ import unicode_literals

import numpy as np
import os
import h5py
from ROOT import gROOT, gStyle, TH1F, TH2F, TGraph, TGraphAsymmErrors, TCanvas, TLegend, TLatex, Double, TFile
from root_numpy import fill_hist
from argparse import ArgumentParser
from array import array
from collections import OrderedDict

from plotting_utils import drawCMS

gROOT.SetBatch(True)
gStyle.SetOptStat(0)

PRESELECTIONS = {
    "pt_min": 300.,
    "pt_max": -1,
    "eta_max": 2.5,
    "m_min": 0.,
    "ID": 2, # put 1 for loose, 2 for tight, 4 for tight leptonveto
    "dEta_min": 0.,
    "dEta_max": 1.3,
    "leptonVeto": 1. # 1: no selection, 0: filter out vetoed events
    }

VARIABLES = { # index of variable for jet1 and jet2 respectively
    "pt": (0,4),
    "eta": (1,5),
    "m": (2,6),
    "ID": (3,7),
    "dEta": 8,
    "leptonVeto": 9
    }

ID_TITLES = {
    0: "No", # currently not passing ntuple preselection
    1: "Loose", # currently not passing ntuple preselection
    2: "Tight",
    4: "Tight LeptonVeto"
    }

LUMIS = {
    "2016": 35920.,
    "2017": 41530.,
    "2018": 59740.,
    "run2": 137190.
    }

TRIGGERS = { # index by which they are stored in the ntuples
    "HLT_Mu50": 1,
    "HLT_IsoMu27": 2,
    "HLT_PFHT800": 3,
    "HLT_PFHT900": 4,
    "HLT_PFHT1050": 5,
    "HLT_PFJet450": 6,
    "HLT_AK8PFJet450": 7,
    "HLT_PFJet500": 8,
    "HLT_AK8PFHT700_TrimR0p1PT0p03Mass50": 9,
    "HLT_AK8PFHT750_TrimMass50": 10,
    "HLT_AK8PFHT800_TrimMass50": 11,
    "HLT_AK8PFJet360_TrimMass30": 12,
    "HLT_AK8PFJet400_TrimMass30": 13,
    "HLT_PFHT650_WideJetMJJ900DEtaJJ1p5": 14
    }

REFERENCE_TRIGGERS = [
    "HLT_Mu50",
    "HLT_IsoMu27"
    ]

# CASE AN choice
TRIGGER_YEARS = { # can define era-specific triggers here. X stands for all eras except the separately specified ones
    "2016": { "X": ["HLT_PFHT800", "HLT_PFJet450"],
              "H": ["HLT_PFHT900", "HLT_PFJet450", "HLT_AK8PFJet450"]},
    "2017": { "X": ["HLT_PFHT1050", "HLT_PFJet500"]},
    "2018": { "X": ["HLT_PFHT1050", "HLT_PFJet500"]}
    }

# Diboson analysis choice
DIBOSON_TRIGGERS = { # can define era-specific triggers here. X stands for all eras except the separately specified ones
    "2016": { "X": ["HLT_PFHT800", "HLT_PFJet450", "HLT_PFHT650_WideJetMJJ900DEtaJJ1p5",
                    "HLT_AK8PFHT700_TrimR0p1PT0p03Mass50", "HLT_AK8PFJet360_TrimMass30"],
              "H": ["HLT_PFHT900", "HLT_PFJet450", "HLT_AK8PFJet450", "HLT_PFHT650_WideJetMJJ900DEtaJJ1p5",
                    "HLT_AK8PFHT700_TrimR0p1PT0p03Mass50", "HLT_AK8PFJet360_TrimMass30"]},
    "2017": { "B": ["HLT_PFHT1050", "HLT_PFJet500"],
              "C": ["HLT_PFHT1050", "HLT_PFJet500", "HLT_AK8PFHT750_TrimMass50", "HLT_AK8PFJet360_TrimMass30"],
              "D": ["HLT_PFHT1050", "HLT_PFJet500", "HLT_AK8PFHT750_TrimMass50", "HLT_AK8PFJet360_TrimMass30"],
              "X": ["HLT_PFHT1050", "HLT_PFJet500", "HLT_AK8PFHT800_TrimMass50", "HLT_AK8PFJet400_TrimMass30"]},
    "2018": { "X": ["HLT_PFHT1050", "HLT_PFJet500", "HLT_AK8PFHT800_TrimMass50", "HLT_AK8PFJet400_TrimMass30" ]}
    }


BINNING = [526, 565, 606, 649, 693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546]#, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010]#, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6099]

def trigger_hist_fill_2D(hist_all, hist_pass, ntuple_file_info, x_var="mjj", y_var="m1"):

    ntuple_file, (_, year, era) = ntuple_file_info
    try:
        current_file = h5py.File(ntuple_file, "r")
    except IOError:
        print("file corrupted:", ntuple_file, "!!!")
        return
    trigger_vars = current_file["trigger_variables"][()]
    preselection_vars = current_file["preselection_variables"][()]

    if trigger_vars.shape[0] !=preselection_vars.shape[0]:
        print("shape mismatch in file", ntuple_file_info, "!!!")
        return

    # Apply preselections
    pt_mask = (preselection_vars[:,VARIABLES["pt"][0]]>PRESELECTIONS["pt_min"]) \
        & (preselection_vars[:,VARIABLES["pt"][1]]>PRESELECTIONS["pt_min"])
    eta_mask = (np.abs(preselection_vars[:,VARIABLES["eta"][0]])<PRESELECTIONS["eta_max"])\
        & (np.abs(preselection_vars[:,VARIABLES["eta"][1]])<PRESELECTIONS["eta_max"])
    m_mask = (preselection_vars[:,VARIABLES["m"][0]]>=PRESELECTIONS["m_min"])\
        & (preselection_vars[:,VARIABLES["m"][1]]>=PRESELECTIONS["m_min"])
    ID_mask = (preselection_vars[:,VARIABLES["ID"][0]] >= PRESELECTIONS["ID"]) \
        & (preselection_vars[:,VARIABLES["ID"][1]]>=PRESELECTIONS["ID"])
    dEta_mask = (preselection_vars[:,VARIABLES["dEta"]]<PRESELECTIONS["dEta_max"])
    lepton_mask = (preselection_vars[:,VARIABLES["leptonVeto"]]<=PRESELECTIONS["leptonVeto"])
    if PRESELECTIONS["pt_max"]!=-1:
        pt_mask = (pt_mask) & (preselection_vars[:,VARIABLES["pt"][0]]<=PRESELECTIONS["pt_max"]) \
            & (preselection_vars[:,VARIABLES["pt"][1]]<=PRESELECTIONS["pt_max"])
    if PRESELECTIONS["dEta_min"]!=0:
        dEta_mask = (dEta_mask) & (preselection_vars[:,VARIABLES["dEta"]]>=PRESELECTIONS["dEta_min"])
    reference_trigger_mask = pt_mask!=pt_mask
    for trigger in REFERENCE_TRIGGERS:
        reference_trigger_mask = (reference_trigger_mask)|(trigger_vars[:,TRIGGERS[trigger]]==1)
    preselection_mask = pt_mask & eta_mask & m_mask & ID_mask & dEta_mask & lepton_mask & reference_trigger_mask

    # Put together year/era-wise triggers
    trigger_mask = pt_mask!=pt_mask
    if era in TRIGGER_YEARS[year].keys():
        current_triggers = TRIGGER_YEARS[year][era]
    else:
        current_triggers = TRIGGER_YEARS[year]["X"]
    for trigger in current_triggers:
        trigger_mask = (trigger_mask)|(trigger_vars[:,TRIGGERS[trigger]]==1)

    if x_var != "mjj":
        if x_var.endswith("1"):
            var_index_x = VARIABLES[x_var.replace("1","")][0]
        elif x_var.endswith("2"):
            var_index_x = VARIABLES[x_var.replace("2","")][1]
        else:
            var_index_x = VARIABLES[x_var][0]
        x_vals = preselection_vars[:,var_index_x]
    else:
        x_vals = trigger_vars[:,0]
   
    if y_var != "mjj":
        if y_var.endswith("1"):
            var_index_y = VARIABLES[y_var.replace("1","")][0]
        elif y_var.endswith("2"):
            var_index_y = VARIABLES[y_var.replace("2","")][1]
        else:
            var_index_y = VARIABLES[y_var][0]
        y_vals = preselection_vars[:,var_index_y]
    else:
        y_vals = trigger_vars[:,0]
 
    fill_hist(hist_all, np.vstack((x_vals, y_vals)).T[preselection_mask])
    fill_hist(hist_pass, np.vstack((x_vals, y_vals)).T[(preselection_mask)&(trigger_mask)])


def print_hist_content(hist):
    n_bins = hist.GetNbinsX()
    bin_list = [hist.GetBinLowEdge(0)]
    content_list = []
    for i in range(n_bins+1):
        bin_list.append(hist.GetBinLowEdge(i+1))
        content_list.append(hist.GetBinContent(i))
    print("bins =", bin_list)
    print("content =", content_list, "\n")
    return bin_list, content_list


def print_threshold(graph, efficiency_threshold):
    n_points = graph.GetN()
    for i in range(n_points):
        x_val = Double()
        y_val = Double()
        graph.GetPoint(i, x_val, y_val)
        if y_val > efficiency_threshold:
            print("An efficiency of {}% > {}% is reached at {}".format(
                y_val*100, efficiency_threshold*100., x_val))
            return x_val, y_val


def trigger_efficiency_line_2D(data_path, base_dataset, year, output_file, binning_x, binning_y, x_var="mjj", y_var="m1"):

    print("preselections:", PRESELECTIONS)
    print("\ntriggers:", TRIGGER_YEARS)

    bin_array_x = array("d", binning_x)
    bin_array_y = array("d", binning_y)
    hist_pass = TH2F("pass_{}_{}".format(base_dataset, year),
        "pass_{}_{}".format(base_dataset, year), len(bin_array_x)-1,
        bin_array_x, len(bin_array_y)-1, bin_array_y)
    hist_all = TH2F("all_{}_{}".format(base_dataset, year),
        "all_{}_{}".format(base_dataset, year), len(bin_array_x)-1,
        bin_array_x, len(bin_array_y)-1, bin_array_y)
   
    # Get list of ntuple file paths
    file_list = []
    for data_path_content in os.listdir(data_path):
        candidate_dir = os.path.join(data_path, data_path_content)
        decoded_content = data_path_content.split("_")
        if os.path.isdir(candidate_dir) and decoded_content[0]==base_dataset and decoded_content[1]==str(year):
            for candidate_dir_content in os.listdir(candidate_dir):
                if candidate_dir_content.endswith(".h5"):
                    file_list.append((os.path.join(candidate_dir, candidate_dir_content), decoded_content))
    
    # Fill each file content into the histograms
    for i, ntuple_file in enumerate(file_list):
        #if i>50:
        #    print("\n------------------------------- breaking here -------------------------------\n")
        #    break
        if i%100.==0:
            print("filling in file number", i)
        trigger_hist_fill_2D(hist_all, hist_pass, ntuple_file, x_var=x_var, y_var=y_var)

    # Make efficiency graph
    #hist_pass.Sumw2()
    #hist_all.Sumw2()
    eff = hist_pass.Clone("efficiency")
    eff.Divide(hist_all)
    #print_hist_content(hist_all)
    return eff, (hist_pass, hist_all)


def combined_efficiency_2D(hist_pass_list, hist_all_list):
    assert len(hist_pass_list)==len(hist_all_list)
    
    for i in range(len(hist_pass_list)):
        if i==0:
            hist_pass = hist_pass_list[i].Clone("pass_combined")
            hist_all = hist_all_list[i].Clone("all_combined")
        else:
            hist_pass.Add(hist_pass_list[i])
            hist_all.Add(hist_all_list[i])

    eff = hist_pass.Clone("efficiency")
    eff.Divide(hist_all)
    return eff, (hist_pass, hist_all)


def trigger_efficiency_plot_2D(line_dict, year, output_file, binning_x, binning_y, x_var="mjj", y_var="m1"):
 
    c1 = TCanvas("c1", "Trigger Efficiency", 800, 700)
    c1.cd(1)

    if x_var != "mjj":
        if x_var =="pt1":
            x_label = "leading jet p_{T} (GeV)"
        elif x_var =="pt2":
            x_label = "subleading jet p_{T} (GeV)"
        elif x_var =="m1":
            x_label = "subleading jet mass (GeV)"
        elif x_var =="m2":
            x_label = "subleading jet mass (GeV)"
    else:
        x_label = "m_{jj} (GeV)"

    if y_var != "mjj":
        if y_var =="pt1":
            y_label = "leading jet p_{T} (GeV)"
        elif y_var =="pt2":
            y_label = "subleading jet p_{T} (GeV)"
        elif y_var =="m1":
            y_label = "leading jet mass (GeV)"
        elif y_var =="m2":
            y_label = "subleading jet mass (GeV)"
    else:
        y_label = "m_{jj} (GeV)"

    for i, line_title in enumerate(line_dict.keys()):
        line_dict["graph"].SetContour(50)
        line_dict["graph"].Draw("COLZ")
        line_dict["graph"].SetTitle(";{};{};{}".format(x_label, y_label, "trigger efficiency"))
        line_dict["graph"].SetMinimum(0.)
        line_dict["graph"].SetMaximum(1.)

        line_dict["graph"].GetXaxis().SetTitleSize(0.045)
        line_dict["graph"].GetYaxis().SetTitleSize(0.045)
        line_dict["graph"].GetYaxis().SetTitleOffset(1.1)
        line_dict["graph"].GetXaxis().SetTitleOffset(1.05)
        line_dict["graph"].GetZaxis().SetTitleOffset(1.1)
        line_dict["graph"].GetXaxis().SetLimits(binning_x[0], binning_x[-1])
        line_dict["graph"].GetYaxis().SetLimits(binning_y[0], binning_y[-1])

    c1.SetTopMargin(0.05)
    c1.SetRightMargin(0.15)

    text = TLatex()
    text.SetTextColor(1)
    text.SetTextFont(42)
    text.SetTextAlign(11)
    text.SetTextSize(0.045)
    ID_string = "PF jet %s ID" % (ID_TITLES[PRESELECTIONS["ID"]])
    eta_string = "jet #eta < %.1f" % (PRESELECTIONS["eta_max"])
    pt_string = "jet p_{T} > %.0f GeV" % (PRESELECTIONS["pt_min"])
    deta_string = "|#Delta#eta| < %.1f" % (PRESELECTIONS["dEta_max"])
    if PRESELECTIONS["dEta_min"]!=0:
        deta_string =  "%.1f < "%(PRESELECTIONS["dEta_min"]) + deta_string
    if PRESELECTIONS["pt_max"]!=-1:
        pt_string = "%.0f < jet p_{T} < %.0f GeV" % (PRESELECTIONS["pt_min"], PRESELECTIONS["pt_max"])
    if PRESELECTIONS["m_min"] > 0:
        m_string = "jet m > %.0f GeV" % (PRESELECTIONS["m_min"])
        text.DrawLatexNDC(0.565, 0.3, "#splitline{#splitline{#splitline{%s}{%s}}{#splitline{%s}{%s}}}{%s}" % (ID_string, eta_string, pt_string, deta_string, m_string))
    else:
        text.DrawLatexNDC(0.565, 0.3, "#splitline{#splitline{%s}{%s}}{#splitline{%s}{%s}}" % (ID_string, eta_string, pt_string, deta_string))

    text.Draw("SAME")

    drawCMS(LUMIS[year], "Preliminary", year=year)

    c1.Print(output_file)


def trigger_efficiency_2D(data_path, base_dataset, year, output_file,
    binning_x, binning_y, x_var="mjj", y_var="m1"):

    if year=="run2":
        print("producing all three years and overlaying them...")
        graph_2016, (hist_pass_2016, hist_all_2016) = trigger_efficiency_line_2D(
            data_path, base_dataset, "2016", output_file, 
            binning_x, binning_y, x_var=x_var, y_var=y_var)
        graph_2017, (hist_pass_2017, hist_all_2017) = trigger_efficiency_line_2D(
            data_path, base_dataset, "2017", output_file,
            binning_x, binning_y, x_var=x_var, y_var=y_var)
        graph_2018, (hist_pass_2018, hist_all_2018) = trigger_efficiency_line_2D(
            data_path, base_dataset, "2018", output_file,
            binning_x, binning_y, x_var=x_var, y_var=y_var)

        hist_pass_list = [hist_pass_2016, hist_pass_2017, hist_pass_2018]
        hist_all_list = [hist_all_2016, hist_all_2017, hist_all_2018]
        graph_run2, (hist_pass_combined, hist_all_combined) = combined_efficiency_2D(hist_pass_list, hist_all_list)
        line_dict = {"graph": graph_run2, "color": 1,"label": "combined (137 fb^{-1})"}

    else:
        graph, (hist_pass_combined, hist_all_combined) = trigger_efficiency_line_2D(
            data_path, base_dataset, year, output_file, binning_x, binning_y,
            x_var=x_var, y_var=y_var)
        line_dict = {"graph": graph, "color": 2, "label": None}

    trigger_efficiency_plot_2D(line_dict, year, output_file, binning_x, binning_y, x_var=x_var, y_var=y_var)



if __name__ == "__main__": 
    parser = ArgumentParser()
    parser.add_argument("-i", "--input_path",  dest="data_path", type=str,
        default="/eos/user/m/msommerh/CASE/trigger_samples", action='store',
        help="Path to the directory containing all ntuples (in subdirectories with primary dataset name and year).")
    parser.add_argument("-o", "--output_file",  dest="output_file", type=str,
        default="trigger_efficiency.pdf", action='store',
        help="Name of output image file including extension.")
    parser.add_argument("-y", "--year",  dest="year", type=str,
        default="2016", action='store', choices=["2016", "2017", "2018", "run2"],
        help="Data taking year (2016, 2017, 2018, run2).")
    parser.add_argument("-d", "--dataset",  dest="base_dataset", type=str,
        default="SingleMuon", action='store',
        help="Name primary dataset (e.g. SingleMuon, JetHT).")
    parser.add_argument("--diboson",  dest="diboson", default=False, action='store_true',
        help="Use diboson triggers instead of CASE ones.")
    parser.add_argument("--eta",  dest="eta", type=float,
        default=-1, action='store',
        help="Jet eta max for preselection.")
    parser.add_argument("--pt",  dest="pt", type=float,
        default=-1, action='store',
        help="Jet pt min for preselection.")
    parser.add_argument("--pt-max",  dest="pt_max", type=float,
        default=-1, action='store',
        help="Jet pt max for preselection.")
    parser.add_argument("--m",  dest="m", type=float,
        default=-1, action='store',
        help="Jet mass min for preselection.")
    parser.add_argument("--deta",  dest="deta", type=float,
        default=-1, action='store',
        help="Dijet deta max for preselection.")
    parser.add_argument("--deta-min",  dest="deta_min", type=float,
        default=-1, action='store',
        help="Dijet deta min for preselection.")
    parser.add_argument("--jet-id",  dest="jet_id", type=float,
        default=-1, action='store',
        help="Jet ID for preselection.")
    parser.add_argument("--lepton-veto",  dest="lepton_veto", default=False, action='store_true',
        help="Apply a lepton veto.")
    parser.add_argument("--x_variable",  dest="x_variable", type=str,
        default="mjj", action='store', choices=["mjj", "pt1", "pt2", "m1", "m2"],
        help="X-axis variable with respect to which the efficiency is drawn.")
    parser.add_argument("--y_variable",  dest="y_variable", type=str,
        default="m1", action='store', choices=["mjj", "pt1", "pt2", "m1", "m2"],
        help="Y-axis variable with respect to which the efficiency is drawn.")
    args = parser.parse_args()

    # Adjust preselections from args
    if args.eta!=-1:
        PRESELECTIONS["eta_max"] = args.eta
    if args.pt!=-1:
        PRESELECTIONS["pt_min"] = args.pt
    if args.pt_max!=-1:
        PRESELECTIONS["pt_max"] = args.pt_max
    if args.deta!=-1:
        PRESELECTIONS["dEta_max"] = args.deta
    if args.deta_min!=-1:
        PRESELECTIONS["dEta_min"] = args.deta_min
    if args.m!=-1:
        PRESELECTIONS["m_min"] = args.m
    if args.jet_id!=-1:
        PRESELECTIONS["ID"] = args.jet_id
    if args.lepton_veto:
        PRESELECTIONS["leptonVeto"] = 0.

    if args.diboson:
        TRIGGER_YEARS = DIBOSON_TRIGGERS

    if args.x_variable != "mjj":
        if "pt" in args.x_variable:
            binning_x = [200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300,
                         310, 320, 330, 340, 350, 360, 370, 380, 390, 400, 410,
                         420, 430, 440, 450, 460, 470, 480, 490, 500, 520, 540,
                         560, 580, 600, 650, 700, 750, 800, 900, 1000]
        elif "m" in args.x_variable:
            binning_x = [0,10, 20, 30, 40, 50, 60, 70, 80, 90, 100,
                         110, 120, 130, 140, 150, 160, 170, 180, 190, 200,
                         210, 220, 230, 240, 250, 260, 270, 280, 290, 300,
                         310, 320, 330, 340, 350, 360, 370, 380, 390, 400,
                         410, 420, 430, 440, 450, 460, 470, 480, 490, 500]
    else:
        binning_x = BINNING                       

    if args.y_variable != "mjj":
        if "pt" in args.y_variable:
            binning_y = [200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300,
                         310, 320, 330, 340, 350, 360, 370, 380, 390, 400, 410,
                         420, 430, 440, 450, 460, 470, 480, 490, 500, 520, 540,
                         560, 580, 600, 650, 700, 750, 800, 900, 1000]
        elif "m" in args.y_variable:
            binning_y = [0,10, 20, 30, 40, 50, 60, 70, 80, 90, 100,
                         110, 120, 130, 140, 150, 160, 170, 180, 190, 200,
                         210, 220, 230, 240, 250, 260, 270, 280, 290, 300,
                         310, 320, 330, 340, 350, 360, 370, 380, 390, 400,
                         410, 420, 430, 440, 450, 460, 470, 480, 490, 500]
    else:
        binning_y = BINNING                       


    trigger_efficiency_2D(args.data_path , args.base_dataset, args.year,
        args.output_file, binning_x, binning_y,
        x_var=args.x_variable, y_var=args.y_variable)
