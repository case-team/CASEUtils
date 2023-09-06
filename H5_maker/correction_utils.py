import awkward as ak
import h5py
import numpy as np
import correctionlib

import ROOT
from ROOT import TLorentzVector, TFile
import numpy as np
import h5py
from optparse import OptionParser
import sys
import utils


from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import *
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.tools import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.JetSysColl import JetSysColl, JetSysObj
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import eventLoop
from PhysicsTools.NanoAODTools.postprocessing.framework.preskimming import preSkim

#pdg ID
top_ID = 6
W_ID = 24
B_ID = 5
MAXLEP_ID = 16

#from cris https://github.com/farakiko/boostedhiggs/blob/main/boostedhiggs/corrections.py#L231
lepton_corrections = {
    "trigger": {
        "muon": {  # For Mu50 (| TkMu50 )
            "2016APV": "NUM_Mu50_or_TkMu50_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose",
            "2016": "NUM_Mu50_or_TkMu50_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose",
            "2017": "NUM_Mu50_or_OldMu100_or_TkMu100_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose",
            "2018": "NUM_Mu50_or_OldMu100_or_TkMu100_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose",
        },
    },
    "id": {
        "muon": {
            "2016APV": "NUM_TightID_DEN_TrackerMuons",
            "2016": "NUM_TightID_DEN_TrackerMuons",
            "2017": "NUM_TightID_DEN_TrackerMuons",
            "2018": "NUM_TightID_DEN_TrackerMuons",
        },
    },
    "iso": {
        "muon": {
            "2016APV": "NUM_TightRelIso_DEN_TightIDandIPCut",
            "2016": "NUM_TightRelIso_DEN_TightIDandIPCut",
            "2017": "NUM_TightRelIso_DEN_TightIDandIPCut",
            "2018": "NUM_TightRelIso_DEN_TightIDandIPCut",
        },
    },
}

#def build_lumimask(filename):
#    from coffea.lumi_tools import LumiMask
#    with importlib.resources.path("boostedhiggs.data", filename) as path:
#        return LumiMask(path)
#lumi_masks = {
#    "2016": build_lumimask("Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt"),
#    "2017": build_lumimask("Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt"),
#    "2018": build_lumimask("Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt"),
#}


"""
CorrectionLib files are available from: /cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration - synced daily
"""
pog_correction_path = "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/"
pog_jsons = {
    "muon": ["MUO", "muon_Z.json.gz"],
    "electron": ["EGM", "electron.json.gz"],
    "pileup": ["LUM", "puWeights.json.gz"],
    "jet": ["JME", "jmar.json.gz"],
    "btag": ["BTV", "btagging.json.gz"],
}

def ang_dist(phi1, phi2):
    dphi = phi1 - phi2
    if(dphi < -math.pi):
        dphi += 2.* math.pi
    if(dphi > math.pi):
        dphi -= 2.*math.pi
    return dphi

def deltaR(o1, o2):
    #o1_vec = ROOT.TLorentzVector()
    #o2_vec = ROOT.TLorentzVector()
    #o1_vec.SetPtEtaPhiM(o1.pt, o1.eta, o1.phi, 0.)
    #o2_vec.SetPtEtaPhiM(o2.pt, o2.eta, o2.phi, 0.)
    #return o1.DeltaR(o2)
    return ((o1.eta - o2.eta)**2 + ang_dist(o1.phi, o2.phi)**2)**(0.5)

def get_UL_year(year):
    if year == "2016":
        year = "2016postVFP"
    elif year == "2016APV":
        year = "2016preVFP"
    return f"{year}_UL"


def get_pog_json(obj, year):
    try:
        pog_json = pog_jsons[obj]
    except:
        print(f'No json for {obj}')
    year = get_UL_year(year)
    return f"{pog_correction_path}POG/{pog_json[0]}/{year}/{pog_json[1]}"
    # os.system(f"cp {pog_correction_path}POG/{pog_json[0]}/{year}/{pog_json[1]} boostedhiggs/data/POG_{pog_json[0]}_{year}_{pog_json[1]}")
    # fname = ""
    # with importlib.resources.path("boostedhiggs.data", f"POG_{pog_json[0]}_{year}_{pog_json[1]}") as filename:
    #     fname = str(filename)
    # print(fname)
    # return fname


def get_puID_SF(jet, year):
    if(jet.pt > 50): return 1.0,1.0,1.0

    cset = correctionlib.CorrectionSet.from_file(get_pog_json("jet", year))
    map_name = "PUJetID_eff"
    wp = "L"
    nominal = cset[map_name].evaluate(abs(jet.eta), jet.pt, "nom", wp)
    up = cset[map_name].evaluate(abs(jet.eta), jet.pt, "up", wp)
    down = cset[map_name].evaluate(abs(jet.eta), jet.pt, "down", wp)

    return nominal,up,down



def get_bjet_SF(jet, year, cset = None, sample = "deepJet_comb", wp = "M"):

    ul_year = get_UL_year(year)
    if(cset is None): cset = correctionlib.CorrectionSet.from_file(get_pog_json("btag", year))
    if(jet.hadronFlavour >= 4): #charm and b
        flavor = int(jet.hadronFlavour)
        key = sample
    else: 
        flavor = 0
        key = sample.replace("comb", "incl")



    nominal = cset[key].evaluate("central", "M", jet.hadronFlavour, abs(jet.eta), jet.pt)
    up = cset[key].evaluate("up", "M", jet.hadronFlavour, abs(jet.eta), jet.pt)
    down = cset[key].evaluate("down", "M", jet.hadronFlavour, abs(jet.eta), jet.pt)
    return nominal, up, down



def get_lepton_weights(lepton, year, lepton_type="muon"):
    ul_year = get_UL_year(year)
    if lepton_type == "electron":
        ul_year = ul_year.replace('_UL', '')
    cset = correctionlib.CorrectionSet.from_file(get_pog_json(lepton_type, year))

    def set_isothreshold(corr, value, lepton_pt, lepton_type):
        iso_threshold = {
            "muon": 55.,
            "electron": 120.,
        }[lepton_type]
        if corr == "trigger_iso":
            value[lepton_pt > iso_threshold] = 1.
        elif corr == "trigger_noniso":
            value[lepton_pt < iso_threshold] = 1.
        elif corr == "isolation":
            value[lepton_pt > iso_threshold] = 1.
        return value

    def get_clip(lep_pt, lep_eta, lepton_type,corr=None):
        clip_pt = [0., 2000]
        clip_eta = [-2.4999, 2.4999]
        if lepton_type == "electron":
            clip_pt = [10.0, 499.999]
            if corr == "reco":
                clip_pt = [20.1, 499.999]
        elif lepton_type == "muon":
            clip_pt = [30., 1000.]
            clip_eta = [0., 2.3999]
            if corr == "trigger_noniso":
                clip_pt = [52., 1000.]
        lepton_pt = np.clip(lep_pt, clip_pt[0], clip_pt[1])
        lepton_eta = np.clip(lep_eta, clip_eta[0], clip_eta[1])
        return lepton_pt,lepton_eta

    lep_pt = np.array(lepton.pt)
    lep_eta = np.array(lepton.eta)
    if lepton_type=="muon": lep_eta = np.abs(lep_eta)

    values = {}
    values["nominal"] = 1.
    for corr,corrDict in lepton_corrections.items():
        if lepton_type not in corrDict.keys():
            continue
        if year not in corrDict[lepton_type].keys():
            continue
        json_map_name = corrDict[lepton_type][year]

        lepton_pt,lepton_eta = get_clip(lep_pt, lep_eta, lepton_type, corr)

        nom = cset[json_map_name].evaluate(ul_year, lepton_eta, lepton_pt, "sf")
        values["nominal"] *= nom


        values[corr + "_up"] = cset[json_map_name].evaluate(ul_year, lepton_eta, lepton_pt, "systup") / nom
        values[corr + "_down"] = cset[json_map_name].evaluate(ul_year, lepton_eta, lepton_pt, "systdown") / nom

        #for key, val in values.items():
            ## restrict values to 1 for some SFs if we are above/below the ISO threshold 
            # Don't understand, comment out
            #values[key] = set_isothreshold(corr, val, np.array(ak.fill_none(lepton.pt, 0.)), lepton_type)

    return values

def get_pdf_weight(inTree):
    weights = np.array(list(inTree.readBranch("LHEPdfWeight")))
    weights = weights[1:-2]
    unc = np.std(weights)
    return 1. + unc, 1. - unc

def get_pileup_weight(year, nPU):
    """
    Should be able to do something similar to lepton weight but w pileup
    e.g. see here: https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/LUMI_puWeights_Run2_UL/
    """
    cset = correctionlib.CorrectionSet.from_file(get_pog_json("pileup", year))
    if('APV' in year ): year = '2016'

    year_to_corr = {'2016': 'Collisions16_UltraLegacy_goldenJSON',
                    '2017': 'Collisions17_UltraLegacy_goldenJSON',
                    '2018': 'Collisions18_UltraLegacy_goldenJSON',
                    }

    values = {}

    nom = cset[year_to_corr[year]].evaluate(nPU, "nominal")
    up = cset[year_to_corr[year]].evaluate(nPU, "up")/ nom
    down  = cset[year_to_corr[year]].evaluate(nPU, "down") / nom

    # add weights (for now only the nominal weight)
    return nom, up, down



def isFinal(statusFlag):
    #check if isLastCopy flag is set
    mask = 1 << 13 #13th bit of status flag
    return (statusFlag & mask) != 0


def get_top_ptrw(event, top = None, anti_top = None):


    if(top is None or anti_top is None): #refind gen particles

        GenPartsColl = Collection(event, "GenPart")

        for genPart in GenPartsColl:
            if(abs(genPart.pdgId) == top_ID and isFinal(genPart.statusFlags)):
                if(genPart.pdgId > 0): 
                    if(top is None): top = genPart
                    else: print("WARNING : Extra top ? ")
                else: 
                    if(anti_top is None): anti_top = genPart
                    else: print("WARNING : Extra antitop ? ")

        if(top is None or anti_top is None):
            print("Couldnt find ttbar pair !")
            return 1.0, 1.0, 1.0
    
    alpha = 0.0615
    beta = 0.0005

    top_pt = min(top.pt, 500.)
    anti_top_pt = min(anti_top.pt, 500.)

    #factors to correct normalization of up and down variations
    up_norm_factor = 0.941
    down_norm_factor = 1.064

    nom = np.exp(alpha - beta * top_pt) * np.exp(alpha - beta * anti_top_pt)
    up = np.exp(alpha - 1.5 * beta * top_pt) * np.exp(alpha - 1.5 * beta * anti_top_pt) / nom / up_norm_factor
    down = np.exp(alpha - 0.5 * beta * top_pt) * np.exp(alpha - 0.5 * beta * anti_top_pt) / nom /down_norm_factor

    #print(top_pt, anti_top_pt, nom)

    return nom, up, down

def get_parent_top(coll, p):
    #find top quark at start of decay chain
    if(p.genPartIdxMother < 0): return None
    mother = coll[p.genPartIdxMother]
    if(abs(mother.pdgId) == top_ID): return p.genPartIdxMother
    return get_parent_top(coll, mother)



def get_ttbar_gen_parts(event, ak8_jet):

    GenPartsColl = Collection(event, "GenPart")

    top = anti_top = W = anti_W = fermion = anti_fermion = b_quark = None


    for genPart in GenPartsColl:
        #tops
        if(abs(genPart.pdgId) == top_ID and isFinal(genPart.statusFlags)):
            if(genPart.pdgId > 0): 
                if(top is None): top = genPart
                else: print("WARNING : Extra top ? ")
            else: 
                if(anti_top is None): anti_top = genPart
                else: print("WARNING : Extra antitop ? ")
        m = genPart.genPartIdxMother
        #W's
        if(abs(genPart.pdgId) == W_ID and isFinal(genPart.statusFlags)):
            if(genPart.pdgId > 0): 
                if(W is None): W = genPart
                else: print("WARNING : Extra W ? ")
            else: 
                if(anti_W is None): anti_W = genPart
                else: print("WARNING : Extra anti W ? ")

    if(top is None or anti_top is None or W is None or anti_W is None):
        print("Couldnt find top or W: ")
        print(top, anti_top, W, anti_W)
    else:
        close_W, close_top = (W,top) if (deltaR(W, ak8_jet) < deltaR(anti_W, ak8_jet)) else (anti_W,anti_top)



    for genPart in GenPartsColl:
        #quarks or leptons from W decay
        m = genPart.genPartIdxMother
        if(abs(genPart.pdgId) <= MAXLEP_ID and m > 0 and GenPartsColl[m] is close_W):
            if(genPart.pdgId > 0): 
                if(fermion is None): fermion = genPart
                else: print("WARNING : Extra quark ? ")
            else: 
                if(anti_fermion is None): anti_fermion = genPart
                else: print("WARNING : Extra anti quark ? ")

        #find b quark from top
        #if(abs(genPart.pdgId) == B_ID): print("b mother:", GenPartsColl[m].pdgId)
        if(abs(genPart.pdgId) == B_ID and GenPartsColl[m] is close_top):
            if(b_quark is None): b_quark = genPart
            else: print("WARNING : Extra quark ? ")



    return top, anti_top, W, anti_W, fermion, anti_fermion, b_quark

def check_matching(jet, f1, f2, b_quark):
    #check if quarks are inside ak8 jet
    #0 = no matching, 1 = W_matched, 2 = top_matched

    f1_in = f1 is not None and abs(f1.pdgId) <= B_ID and deltaR(jet,f1) < 0.8
    f2_in = f2 is not None and abs(f2.pdgId) <= B_ID and deltaR(jet,f2) < 0.8
    b_in = b_quark is not None and deltaR(jet,b_quark) < 0.8

    W_match = f1_in and f2_in
    top_match = W_match and b_in

    if(top_match): return 2
    elif(W_match): return 1
    else: return 0



