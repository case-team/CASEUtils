import h5py
import numpy as np
import ROOT
import sys


#Using CMSSW Wrapped fj contrib modules
#https://github.com/cms-jet/NanoAODJMARTools/tree/master
ROOT.gSystem.Load("libPhysicsToolsNanoAODJMARTools.so")


def compute_taus(pf_cands):

    pfCandsVec = ROOT.vector("TLorentzVector")()
    ncands = 0
    for i,c in enumerate(pf_cands):
        if(c[3] > 0.0001):
            v = ROOT.TLorentzVector(c[0], c[1], c[2], c[3])
            pfCandsVec.push_back(v)
            ncands+=1

    nSubJ = ROOT.NsubjettinessWrapper( 1, 0.8, 0, 6 ) ## beta, cone size, measureDef 0=Normalize, axesDef 6=onepass_kt_axes

    maxTau = 6

    nsub1 = nSubJ.getTau(maxTau, pfCandsVec)
    return nsub1, ncands

#fname = "Lund_output_files_gen_mar13/TT_powheg.h5"
fnames = sys.argv[1:]
print(fnames)


for fname in fnames:
    print("Opening " + fname)
    f = h5py.File(fname, "a")
    if('jet1_extraInfo' in list(f.keys())):
        print("Already made, skipping")
        continue

    cands = f['jet1_PFCands'][()]

    outshape = (cands.shape[0], 7)
    info = np.zeros(outshape)

    for i,cand, in enumerate(cands):
        taus, ncands = compute_taus(cand)

        #match format of reco files
        info[i] = [taus[0], taus[1], taus[2], taus[3], taus[4], taus[5], ncands]

    f.create_dataset("jet1_extraInfo", data=info)

    if('jet2_PFCands' in f.keys()):
        cands = f['jet2_PFCands'][()]

        outshape = (cands.shape[0], 7)
        info = np.zeros(outshape)

        for i,cand, in enumerate(cands):
            taus, ncands = compute_taus(cand)

            #match format of reco files
            info[i] = [taus[0], taus[1], taus[2], taus[3], taus[4], taus[5], ncands]

        f.create_dataset("jet2_extraInfo", data=info)


    f.close()

