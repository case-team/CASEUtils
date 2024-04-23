import h5py
import os
import ROOT
import numpy as np
import sys
from parse import parse
from utils import *




def get_num_events(dname, max_files = -1, dir_pre = ""):
    f_list = get_file_list(dname)

    if(max_files > 0):
        f_list = f_list[:max_files]
    elif(max_files == -10):
        #match files exactly
        f_pre_list = get_file_list(dir_pre)
        f_list_matched = []
        for f in f_pre_list:
            _,tag,_ = parse("{}/nano_mc2018_{}_{}.root", f)
            for f1 in f_list:
                if (tag in f1): 
                    f_list_matched.append(f1)
                    break



    num_evts = 0
    num_total = 0

    print(len(f_list))
    for fname in f_list:
        f = ROOT.TFile.Open(fname, "r")
        t = f.Get("Events")
        num_evts += t.GetEntries()
        

        #before skimming / preselection
        t_runs = f.Get("Runs")
        for entry in t_runs:
            num_total += entry.genEventCount

        f.Close()
    return num_evts, num_total


def redo_preselection_eff(f_input,  fname_preselected = "", xsec = 1.0, max_files = -1):

    print(f_input)
    n_total = 0

    if(".root" in fname_preselected):
        f_pre = ROOT.TFile.Open(fname_preselected, "r")
        t_pre = f_pre.Get("Events")
        n_preselected  = t_pre.GetEntries()

        t_runs = f_pre.Get("Runs")
        t_runs.GetEntry(0)
        n_total = t_runs.ReadBranch("genEventCount")

    else:

        n_preselected, n_total = get_num_events(fname_preselected)



    with h5py.File(f_input, "a") as f:
        print(list(f.keys()))
        n_saved = f['event_info'].shape[0]
        print("Saved %i" % n_saved)
        if(n_saved <1): return
        raw_preselection_eff = n_saved / n_preselected
        weight_corr_factor =   f['preselection_eff'][0]/ raw_preselection_eff

        corr_preselection_eff = (n_saved / n_total) * weight_corr_factor
        print("Skim eff: %.4f" % (n_preselected / n_total))
        print('preselection eff: raw : %4f, orig : %.4f, corrected : %.4f'  % (raw_preselection_eff, f['preselection_eff'][0], corr_preselection_eff))

        if('preselection_eff_corr' in f.keys()):
            f['preselection_eff_corr'][0] = corr_preselection_eff
        else:
            f.create_dataset("preselection_eff_corr", data=np.array([corr_preselection_eff]))

        if (xsec > 0.):
            if('sys_weights' in list(f.keys())):
                gen_weights = f['sys_weights'][:,0]
            else:
                gen_weights = f['event_info'][:,3]

            rw_factor = xsec * 1000. * corr_preselection_eff / np.sum(gen_weights)

            norm_weights  = (gen_weights * rw_factor).reshape(-1)
            print("Total weight is %s" % np.sum(norm_weights))

            if('norm_weights' in f.keys()):
                del f['norm_weights']
            f.create_dataset("norm_weights", chunks = True, data= norm_weights, maxshape = None)

if(__name__ == "__main__"):
    print(sys.argv)
    redo_preselection_eff(sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]), int(sys.argv[5]))

