# Read nanoAOD with PF constituents (aka pancakes), apply a pre-selection and output to an H5 file format

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






def nPFCounter(index, event):
    count = 0
    jet_indices = event.FatJetPFCands_jetIdx
    length = len(event.FatJetPFCands_jetIdx)
    for i in range(length):
        if jet_indices[i] == index:
            count += 1
    
    return count

def inHEMRegion(jet, year):
    if(year == 2018):
        return jet.eta > -3.2 and jet.eta < -1.3 and jet.phi  > -1.57 and jet.phi < -0.87
    else: return False

class Outputer:
    def __init__(self, outputFileName="out.root", batch_size = 5000, truth_label = 0, sample_type="MC", sort_pfcands=True, year = -1):
        self.batch_size = batch_size
        self.output_name = outputFileName
        self.sample_type = sample_type
        self.first_write = False
        self.truth_label = np.array([[truth_label]]*batch_size, dtype=np.int8)
        self.idx = 0
        self.nBatch = 0
        self.n_pf_cands = 100 #how many PF candidates to save (max)
        self.sort_pfcands = sort_pfcands
        self.year = year
        self.reset()

    def reset(self):
        self.idx = 0
        self.jet1_PFCands = np.zeros((self.batch_size, self.n_pf_cands,4), dtype=np.float16)
        self.jet2_PFCands = np.zeros((self.batch_size, self.n_pf_cands, 4), dtype=np.float16)
        self.jet1_extraInfo = np.zeros((self.batch_size, 7), dtype=np.float32)
        self.jet2_extraInfo = np.zeros((self.batch_size, 7), dtype=np.float32)
        self.jet_kinematics = np.zeros((self.batch_size, 14), dtype=np.float32)
        self.event_info = np.zeros((self.batch_size, 8), dtype=np.float32)


    def is_leptonic_decay(self, event):
    #leptonic decays = generator level lepton coming from a particle with mass > 20 GeV
    #not perfect but works pretty well testing on Wkk and W' grav samples

        nGenParts =event.nGenPart

        GenPartsColl = Collection(event, "GenPart")

        for genPart in GenPartsColl:
            m = genPart.genPartIdxMother
            if(m < 0 or abs(m) > nGenParts):
                continue
            elif abs(genPart.pdgId) >= 11 and abs(genPart.pdgId) < 18 and GenPartsColl[m].mass > 20: 
                return True

        return False

    def get_pfcands_sorted(self, pfcands):
        
        
        pfcands_pt = np.sqrt(pfcands[:, 0]**2 + pfcands[:, 1]**2)
        sorted_idx = np.flip(np.argsort(pfcands_pt))
        pfcands = pfcands[sorted_idx]
        
        return pfcands.astype(np.float16)
            
    
    def fill_event(self, inTree, event, jet1, jet2, jet3, PFCands, subjets, mjj, num_jets):

        if self.sample_type == "data":
            genWeight = 1
            leptonic_decay = False
        else:
            genWeight = inTree.readBranch('genWeight')
            leptonic_decay = self.is_leptonic_decay(event)
        
        MET = inTree.readBranch('MET_pt')
        MET_phi = inTree.readBranch('MET_phi')
        eventNum = inTree.readBranch('event')
        run = inTree.readBranch('run')

        event_info = [eventNum, MET, MET_phi, genWeight, leptonic_decay, run, self.year, num_jets]

        d_eta = abs(jet1.eta - jet2.eta)


        jet_kinematics = [mjj, d_eta, jet1.pt, jet1.eta, jet1.phi, jet1.msoftdrop, jet2.pt, jet2.eta, jet2.phi, jet2.msoftdrop]

        if(jet3 != None):
            jet_kinematics.extend([jet3.pt, jet3.eta, jet3.phi, jet3.msoftdrop])
        else:
            jet_kinematics.extend([0., 0., 0., 0.])

        
        #maximum deepcsv from top 2 subjets of the fatjet
        jet1_btag = jet2_btag = -1.
        
        if(jet1.subJetIdx1 >= 0):
            jet1_btag = subjets[jet1.subJetIdx1].btagDeepB
        if(jet1.subJetIdx2 >= 0):
            jet1_btag = max(jet1_btag, subjets[jet1.subJetIdx2].btagDeepB)

        if(jet2.subJetIdx1 >= 0):
            jet2_btag = subjets[jet2.subJetIdx1].btagDeepB
        if(jet2.subJetIdx2 >= 0):
            jet2_btag = max(jet2_btag, subjets[jet2.subJetIdx2].btagDeepB)

        jet1_extraInfo = [jet1.tau1, jet1.tau2, jet1.tau3, jet1.tau4, jet1.lsf3, jet1_btag, jet1.nPFConstituents]
        jet2_extraInfo = [jet2.tau1, jet2.tau2, jet2.tau3, jet2.tau4, jet2.lsf3, jet2_btag, jet2.nPFConstituents]
        #print(jet1.PFConstituents_Start, jet1.PFConstituents_Start + jet1.nPFConstituents, jet2.PFConstituents_Start, jet2.PFConstituents_Start + jet2.nPFConstituents)

        j1_nPF = min(self.n_pf_cands, jet1.nPFConstituents)
        j2_nPF = min(self.n_pf_cands, jet2.nPFConstituents)
        range1 = range(jet1.PFConstituents_Start, jet1.PFConstituents_Start + j1_nPF, 1)
        range2 = range(jet2.PFConstituents_Start, jet2.PFConstituents_Start + j2_nPF, 1)
        jet1_PFCands = []
        jet2_PFCands = []
        for idx in range1:
            cand = ROOT.Math.PtEtaPhiMVector(PFCands[idx].pt, PFCands[idx].eta, PFCands[idx].phi, PFCands[idx].mass)
            jet1_PFCands.append([cand.Px(), cand.Py(), cand.Pz(), cand.E()])

        for idx in range2:
            cand = ROOT.Math.PtEtaPhiMVector(PFCands[idx].pt, PFCands[idx].eta, PFCands[idx].phi, PFCands[idx].mass)
            jet2_PFCands.append([cand.Px(), cand.Py(), cand.Pz(), cand.E()])


        self.event_info[self.idx] = np.array(event_info, dtype=np.float32)
        self.jet_kinematics[self.idx] = np.array(jet_kinematics, dtype = np.float32)
        self.jet1_extraInfo[self.idx] = np.array(jet1_extraInfo, dtype = np.float32)
        self.jet2_extraInfo[self.idx] = np.array(jet2_extraInfo, dtype = np.float32)
        
        # sort PFCands by pt
        if self.sort_pfcands:
            self.jet1_PFCands[self.idx,:jet1.nPFConstituents] = self.get_pfcands_sorted(np.array(jet1_PFCands, dtype = np.float32))
            self.jet2_PFCands[self.idx,:jet2.nPFConstituents] = self.get_pfcands_sorted(np.array(jet2_PFCands, dtype = np.float32))
        else:
            self.jet1_PFCands[self.idx,:jet1.nPFConstituents] = np.array(jet1_PFCands, dtype = np.float16)
            self.jet2_PFCands[self.idx,:jet2.nPFConstituents] = np.array(jet2_PFCands, dtype = np.float16)

        self.idx +=1
        if(self.idx % self.batch_size == 0): self.write_out()


    def write_out(self):
        self.idx = 0
        print("Writing out batch %i \n" % self.nBatch)
        self.nBatch += 1
        write_size = self.event_info.shape[0]
        truth_label_write = self.truth_label[:write_size]

        if(not self.first_write):
            self.first_write = True
            print("First write, creating dataset with name %s \n" % self.output_name)
            with h5py.File(self.output_name, "w") as f:
                f.create_dataset("truth_label", data=truth_label_write, chunks = True, maxshape=(None,1))
                f.create_dataset("event_info", data=self.event_info, chunks = True, maxshape=(None, self.event_info.shape[1]))
                f.create_dataset("jet_kinematics", data=self.jet_kinematics, chunks = True, maxshape=(None, self.jet_kinematics.shape[1]))
                f.create_dataset("jet1_extraInfo", data=self.jet1_extraInfo, chunks = True, maxshape=(None, self.jet1_extraInfo.shape[1]))
                f.create_dataset("jet2_extraInfo", data=self.jet2_extraInfo, chunks = True, maxshape=(None, self.jet2_extraInfo.shape[1]))
                f.create_dataset("jet1_PFCands", data=self.jet1_PFCands, chunks = True, maxshape=(None, self.jet1_PFCands.shape[1], 4), compression='gzip')
                f.create_dataset("jet2_PFCands", data=self.jet2_PFCands, chunks = True, maxshape=(None, self.jet2_PFCands.shape[1], 4), compression='gzip')

        else:
            with h5py.File(self.output_name, "a") as f:
                utils.append_h5(f,'truth_label',truth_label_write)
                utils.append_h5(f,'event_info',self.event_info)
                utils.append_h5(f,'jet_kinematics',self.jet_kinematics)
                utils.append_h5(f,'jet1_extraInfo',self.jet1_extraInfo)
                utils.append_h5(f,'jet2_extraInfo',self.jet2_extraInfo)
                utils.append_h5(f,'jet1_PFCands',self.jet1_PFCands)
                utils.append_h5(f,'jet2_PFCands',self.jet2_PFCands)

        self.reset()

    def final_write_out(self, eff):
        if(self.idx < self.batch_size):
            print("Last batch only filled %i events, shortening arrays \n" % self.idx)
            self.jet1_PFCands = self.jet1_PFCands[:self.idx]
            self.jet2_PFCands = self.jet2_PFCands[:self.idx]
            self.jet1_extraInfo = self.jet1_extraInfo[:self.idx]
            self.jet2_extraInfo = self.jet2_extraInfo[:self.idx]
            self.jet_kinematics = self.jet_kinematics[:self.idx] 
            self.event_info = self.event_info[:self.idx]

        self.write_out()
        with h5py.File(self.output_name, "a") as f:
            f.create_dataset("preselection_eff", data=np.array([eff]))


def NanoReader(process_flag, inputFileNames=["in.root"], outputFileName="out.root", json = '', year = 2016, nEventsMax = -1, sampleType = "MC", sort_pfcands=True):
    
    if not ((sampleType == "MC") or (sampleType=="data")):
        print("Error! sampleType needs to be set to either data or MC! Please set correct option and retry.")
        sys.exit()
    
    #Applying standard MET filters: https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2#Analysis_Recommendations_for_ana
    filters = ["Flag_goodVertices",
    "Flag_globalSuperTightHalo2016Filter",
    "Flag_HBHENoiseFilter",
    "Flag_HBHENoiseIsoFilter",
    "Flag_EcalDeadCellTriggerPrimitiveFilter",
    "Flag_BadPFMuonFilter",
    "Flag_BadPFMuonDzFilter",
    "Flag_eeBadScFilter", 
    ]
    if((year == 2016) or (year == 2016.5)): filters.append("Flag_CSCTightHaloFilter")
    if(year == 2017 or year == 2018): filters.append("Flag_ecalBadCalibFilter")

    if(year == 2016):
        triggers = ["HLT_PFHT800", "HLT_PFJet450"]
    #2016.5 corresponds to run 2016H, which needs specific triggers
    elif (year == 2016.5):
        triggers = ["HLT_PFHT900", "HLT_PFJet450", "HLT_AK8PFJet450"]
    elif(year == 2017):
        triggers = ["HLT_PFHT1050", "HLT_AK8PFJet500"]
    elif(year == 2018):
        triggers = ["HLT_PFHT1050", "HLT_AK8PFJet500"]
    else:
        print("Invalid year option of %i. Year must be 2016, 2017, or 2018! \n" % year)
        exit(1)

        triggers = [
                'HLT_PFHT780',
                'HLT_PFHT890',
                'HLT_PFHT1050',
                'HLT_PFJet500',
                'HLT_AK8PFJet450',
                'HLT_AK8PFJet500'
                ]

    mjj_cut = 1200.
    
    nFiles = len(inputFileNames)
    print("Will run over %i files and output to %s with truth label %i" % (nFiles, outputFileName, process_flag))
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

        out = Outputer(outputFileName, truth_label =  process_flag, sample_type=sampleType, sort_pfcands=sort_pfcands, year = year)


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



            
            passTrigger = False
            passFilter = True
            for fil in filters:
                passFilter = passFilter and inTree.readBranch(fil)
            if(not passFilter): continue
            
            # Apply triggers only to data and MC
            for trig in triggers:
                passTrigger = passTrigger or inTree.readBranch(trig)

            if(not passTrigger): continue



            PFCands = Collection(event, "FatJetPFCands")
                        
            try:
                mytest = event.FatJetPFCands_eta
            except:
                FullPFCands = Collection(event, "PFCands")
                for cand in PFCands:
                    cand.eta = FullPFCands[cand.pFCandsIdx].eta
                    cand.phi = FullPFCands[cand.pFCandsIdx].phi
                    cand.mass = FullPFCands[cand.pFCandsIdx].mass
            
            AK8Jets = Collection(event, "FatJet")
            #MuonsCol = Collection(event, "Muon")
            #ElectronsCol = Collection(event, "Electron")
            #PhotonsCol = Collection(event, "Photon")
            subjets = Collection(event, "SubJet")

            jet_min_pt = 300
            #keep 2 jets with pt > 200, tight id
            jet1 = jet2 = jet3 =  None
        
            pf_conts_start = 0 #keep track of indices for PF candidates
            jet_index = 0
            num_jets = 0
            for jet in AK8Jets:
                #jetId : bit1 = loose, bit2 = tight, bit3 = tightLepVeto
                #want tight id
                if((jet.jetId & 2 == 2) and abs(jet.eta) < 2.5):
                    jet.PFConstituents_Start = pf_conts_start
                    if(jet.pt > 50): num_jets+=1
                    if((jet1 == None or jet.pt > jet1.pt) and jet.pt > jet_min_pt):
                        jet3 = jet2
                        jet2 = jet1
                        jet1 = jet
                    elif((jet2 == None or jet.pt > jet2.pt) and jet.pt > jet_min_pt):
                        jet3 = jet2
                        jet2 = jet
                    elif(jet3 == None or jet.pt > jet3.pt):
                        jet3 = jet
                
                
                try:
                    pf_conts_start += jet.nPFConstituents # try/except needed because FatJet_nPFConstituents isn't stored when using PFNano
                except:
                    jet.nPFConstituents = nPFCounter(jet_index, event)
                    pf_conts_start += jet.nPFConstituents
                
                jet_index += 1
            
            
            

            if(jet1 == None or jet2 == None): continue


            if(inHEMRegion(jet1, year) or inHEMRegion(jet2, year)): continue

            #Order jets so jet1 is always the higher mass one
            if(jet1.msoftdrop < jet2.msoftdrop):
                temp = jet1
                jet1 = jet2
                jet2 = temp

            j1_4vec = ROOT.Math.PtEtaPhiMVector(jet1.pt, jet1.eta, jet1.phi, jet1.msoftdrop)
            j2_4vec = ROOT.Math.PtEtaPhiMVector(jet2.pt, jet2.eta, jet2.phi, jet2.msoftdrop)

            dijet = j1_4vec + j2_4vec
            mjj = dijet.M()

            if(mjj< mjj_cut): continue
        

            saved+=1
            out.fill_event(inTree, event, jet1, jet2, jet3, PFCands, subjets, mjj, num_jets)
            if(nEventsMax > 0 and saved >= nEventsMax): break
# -------- End Loop over tree-------------------------------------
# -------- End Loop over files-------------------------------------

    efficiency = float(saved)/count
    out.final_write_out(efficiency)
    print("Done. Selected %i events. Selection efficiency is %.3f \n" % (saved, efficiency))
    print("Outputed to %s" % outputFileName)
    return saved


