# Read nanoAOD with PF constituents (aka pancakes), apply a pre-selection and output to an H5 file format
# adjusted from the CASE H5_maker.py to store trigger-specific variables

import ROOT
from ROOT import TLorentzVector, TFile
import numpy as np
import h5py
from optparse import OptionParser

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import *
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.tools import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.JetSysColl import JetSysColl, JetSysObj
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import eventLoop
from PhysicsTools.NanoAODTools.postprocessing.framework.preskimming import preSkim

triggers = [
    "HLT_Mu50", "HLT_IsoMu27", # reference triggers
    "HLT_PFHT800", "HLT_PFHT900", "HLT_PFHT1050",
    "HLT_PFJet450", "HLT_AK8PFJet450", "HLT_PFJet500",
    "HLT_AK8PFHT700_TrimR0p1PT0p03Mass50", "HLT_AK8PFHT750_TrimMass50", "HLT_AK8PFHT800_TrimMass50",
    "HLT_AK8PFJet360_TrimMass30", "HLT_AK8PFJet400_TrimMass30",
    "HLT_PFHT650_WideJetMJJ900DEtaJJ1p5" 
    ]

## pre-preselections:
MIN_JET_PT = 200.
MIN_MJJ = 500.
MAX_ETA = 2.5


def append_h5(f, name, data):
    prev_size = f[name].shape[0]
    f[name].resize(( prev_size + data.shape[0]), axis=0)
    f[name][prev_size:] = data


def deltaPhi(phi1, phi2):
    dphi = phi1-phi2
    while dphi>=np.pi:
        dphi -= 2*np.pi
    while dphi<-np.pi:
        dphi += 2*np.pi
    return dphi


def deltaR(object1, object2):
    return (((object1.eta-object2.eta)**2) + ((deltaPhi(object1.phi, object2.phi))**2))**0.5


class Outputer:
    def __init__(self, outputFileName="out.root", batch_size = 5000):
        self.batch_size = batch_size
        self.output_name = outputFileName
        self.first_write = False
        self.idx = 0
        self.nBatch = 0
        self.reset()

    def reset(self):
        self.idx = 0
        self.trigger_variables = np.zeros((self.batch_size, 1+len(triggers)), dtype=np.float16)
        self.preselection_variables = np.zeros((self.batch_size, 10), dtype=np.float16)

    def fill_event(self, inTree, jet1, jet2, mjj, pass_trigger, lepton_veto):

        d_eta = abs(jet1.eta - jet2.eta)

        trigger_variables = [mjj]+pass_trigger
        preselection_variables = [jet1.pt, jet1.eta, jet1.msoftdrop, jet1.jetId, jet2.pt, jet2.eta, jet2.msoftdrop, jet2.jetId, d_eta, lepton_veto]

        self.trigger_variables[self.idx] = np.array(trigger_variables, dtype=np.float16)
        self.preselection_variables[self.idx] = np.array(preselection_variables, dtype=np.float16)

        self.idx +=1
        if(self.idx % self.batch_size == 0): self.write_out()

    def write_out(self):
        self.idx = 0
        print("Writing out batch %i \n" % self.nBatch)
        self.nBatch += 1
        write_size = self.trigger_variables.shape[0]

        if(not self.first_write):
            self.first_write = True
            print("First write, creating dataset with name %s \n" % self.output_name)
            with h5py.File(self.output_name, "w") as f:
                f.create_dataset("trigger_variables", data=self.trigger_variables, chunks=True, maxshape=(None, self.trigger_variables.shape[1]))
                f.create_dataset("preselection_variables", data=self.preselection_variables, chunks=True, maxshape=(None, self.preselection_variables.shape[1]))

        else:
            with h5py.File(self.output_name, "a") as f:
                append_h5(f,'trigger_variables',self.trigger_variables)
                append_h5(f,'preselection_variables',self.preselection_variables)
        self.reset()

    def final_write_out(self):
        if(self.idx < self.batch_size):
            print("Last batch only filled %i events, shortening arrays \n" % self.idx)
            self.trigger_variables = self.trigger_variables[:self.idx] 
            self.preselection_variables = self.preselection_variables[:self.idx]

        self.write_out()


def NanoReader(inputFileNames=["in.root"], outputFileName="out.root", json='', year=2016, nEventsMax=-1):

    filters = ["Flag_goodVertices",
            "Flag_globalTightHalo2016Filter",
            "Flag_eeBadScFilter", 
            "Flag_HBHENoiseFilter",
            "Flag_HBHENoiseIsoFilter",
            "Flag_ecalBadCalibFilter",
            "Flag_EcalDeadCellTriggerPrimitiveFilter",
            "Flag_BadChargedCandidateFilter",
    ]

    if(year == 2016): filters.append("Flag_CSCTightHaloFilter")

    if(year == 2016):
        ht_triggers = ["HLT_PFHT800", "HLT_PFHT900"]
        pt_triggers = ["HLT_PFJet450", "HLT_AK8PFJet450"]
        ht_sub_triggers = ["HLT_AK8PFHT700_TrimR0p1PT0p03Mass50"]
        pt_sub_triggers = ["HLT_AK8PFJet360_TrimMass30"]
    elif(year == 2017):
        ht_triggers = ["HLT_PFHT1050"]
        pt_triggers = ["HLT_PFJet500"]
        ht_sub_triggers = ["HLT_AK8PFHT750_TrimMass50", "HLT_AK8PFHT800_TrimMass50"]
        pt_sub_triggers = ["HLT_AK8PFJet360_TrimMass30", "HLT_AK8PFJet400_TrimMass30"]
    elif(year == 2018):
        ht_triggers = ["HLT_PFHT1050"]
        pt_triggers = ["HLT_PFJet500"]
        ht_sub_triggers = ["HLT_AK8PFHT800_TrimMass50"]
        pt_sub_triggers = ["HLT_AK8PFJet400_TrimMass30"]
    else:
        print("Invalid year option of %i. Year must be 2016, 2017, or 2018! \n" % year)
        exit(1)

    mjj_cut = MIN_MJJ
    
    nFiles = len(inputFileNames)
    print("Will run over %i files and output to %s" % (nFiles, outputFileName))
    count = 0
    saved = 0

#----------------- Begin loop over files ---------------------------------

    for fileName in inputFileNames:

        print("Opening file %s" % fileName)

        inputFile = TFile.Open(fileName)
        if(not inputFile): #check for null pointer
            print("Unable to open file %s, exting \n" % fileName)
            return 1

        #get input tree
        inTree = inputFile.Get("Events")

        # pre-skimming
        if(json != ''):
            elist,jsonFilter = preSkim(inTree, json)

            #number of events to be processed 
            nTotal = elist.GetN() if elist else inTree.GetEntries()
            
            print('Pre-select %d entries out of %s '%(nTotal,inTree.GetEntries()))
            inTree= InputTree(inTree, elist) 
        else:
            nTotal = inTree.GetEntries()
            inTree= InputTree(inTree) 
            print('Running over %i entries \n' % nTotal)

        out = Outputer(outputFileName)

        # Grab event tree from nanoAOD
        eventBranch = inTree.GetBranch('event')
        treeEntries = eventBranch.GetEntries()

# -------- Begin Loop over tree-------------------------------------

        entries = inTree.entries
        for entry in xrange(entries):

            if count % 10000 == 0 :
                print('--------- Processing Event ' + str(count) +'   -- percent complete ' + str(100*count/nTotal/nFiles) + '% -- ')

            count +=1
            # Grab the event
            event = Event(inTree, entry)

            passTrigger = [False for trig in triggers]
            passFilter = True
            for fil in filters:
                passFilter = passFilter and inTree.readBranch(fil)
            if(not passFilter): continue
            for j, trig in enumerate(triggers):
                try:
                    passTrigger_current = inTree.readBranch(trig)
                except RuntimeError:
                    passTrigger_current = False
                passTrigger[j] = passTrigger_current

            AK8Jets = Collection(event, "FatJet")

            jet_min_pt = MIN_JET_PT
            jet1 = jet2 = jet3 =  None
        
            for jet in AK8Jets:
                if((jet.jetId & 2 == 2) and jet.pt > jet_min_pt and abs(jet.eta) < MAX_ETA):
                    if(jet1 == None or jet.pt > jet1.pt):
                        jet3 = jet2
                        jet2 = jet1
                        jet1 = jet
                    elif(jet2 == None or jet.pt > jet2.pt):
                        jet3 = jet2
                        jet2 = jet
                    elif(jet3 == None or jet.pt > jet3.pt):
                        jet3 = jet

            if(jet1 == None or jet2 == None): continue

            #Order jets so jet1 is always the higher mass one
            if(jet1.msoftdrop < jet2.msoftdrop):
                temp = jet1
                jet1 = jet2
                jet2 = temp

            j1_4vec = ROOT.Math.PtEtaPhiMVector(jet1.pt, jet1.eta, jet1.phi, jet1.msoftdrop)
            j2_4vec = ROOT.Math.PtEtaPhiMVector(jet2.pt, jet2.eta, jet2.phi, jet2.msoftdrop)

            # lepton veto
            lepton_veto = False
            electrons = Collection(event, "Electron")
            for electron in electrons:
                if electron.pt < 35 or np.abs(electron.eta) > 2.5 or not electron.cutBased_HEEP:
                    continue
                if deltaR(jet1, electron) < 0.8:
                    lepton_veto = True
                elif deltaR(jet2, electron) < 0.8:
                    lepton_veto = True
                if lepton_veto:
                    break

            muons = Collection(event, "Muon")
            for muon in muons:
                if muon.pt < 30 or np.abs(muon.eta) > 2.4 or not (muon.highPtId&2)==2 or muon.pfRelIso04_all >= 0.05:
                    continue
                if deltaR(jet1, muon) < 0.8:
                    lepton_veto = True
                elif deltaR(jet2, muon) < 0.8:
                    lepton_veto = True
                if lepton_veto:
                    break

            dijet = j1_4vec + j2_4vec
            mjj = dijet.M()

            if(mjj< mjj_cut): continue

            saved+=1
            out.fill_event(inTree, jet1, jet2, mjj, passTrigger, lepton_veto)
            if(nEventsMax > 0 and saved >= nEventsMax): break
# -------- End Loop over tree-------------------------------------
# -------- End Loop over files-------------------------------------

    efficiency = float(saved)/count
    out.final_write_out()
    print("Done. Selected %i events. Selection efficiency is %.3f \n" % (saved, efficiency))
    print("Outputed to %s" % outputFileName)
    return saved


