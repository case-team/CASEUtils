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


sys_weights_map = {
        'nom_weight' : 0,
        'pdf_up' : 1,
        'pdf_down': 2,
        'prefire_up': 3,
        'prefire_down' : 4,
        'pileup_up' : 5 ,
        'pileup_down' : 6,
        'btag_up' : 7,
        'btag_down' : 8,
        'PS_ISR_up' : 9,
        'PS_ISR_down' : 10,
        'PS_FSR_up' : 11,
        'PS_FSR_down' : 12,
        'F_up' : 13,
        'F_down' : 14,
        'R_up' : 15,
        'R_down' : 16,
        'RF_up' : 17,
        'RF_down' : 18
        }

JME_vars_map = {
        'pt_JES_up' : 0,
        'm_JES_up' : 1,
        'pt_JES_down' : 2,
        'm_JES_down' : 3,
        'pt_JER_up' : 4,
        'm_JER_up' : 5,
        'pt_JER_down' : 6,
        'm_JER_down' : 7,
        'm_JMS_up' : 8,
        'm_JMS_down' : 9,
        'm_JMR_up' : 10,
        'm_JMR_down' : 11
        }






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

def get_branch_mean(inTree, branch_name):

    inTree.Draw(branch_name)
    temp = ROOT.gPad.GetPrimitive("htemp")
    return temp.GetMean()



class Outputer:
    def __init__(self, outputFileName="out.root", batch_size = 5000, truth_label = 0, sample_type="MC", sort_pfcands = False, include_systematics = True, year = 2017):

        self.batch_size = batch_size
        self.output_name = outputFileName
        self.sample_type = sample_type
        self.first_write = False
        self.truth_label = np.array([[truth_label]]*batch_size, dtype=np.int8)
        self.idx = 0
        self.nBatch = 0
        self.n_pf_cands = 100 #how many PF candidates to save (max)
        self.include_systematics = include_systematics
        self.n_SVs = 10 #how many SVs candidates to save (max)
        self.sort_pfcands = sort_pfcands
        self.year = year
        self.reset()

    def reset(self):
        self.idx = 0
        self.jet1_PFCands = np.zeros((self.batch_size, self.n_pf_cands,4), dtype=np.float16)
        self.jet2_PFCands = np.zeros((self.batch_size, self.n_pf_cands, 4), dtype=np.float16)
        self.jet1_SVs = np.zeros((self.batch_size, self.n_SVs,6), dtype=np.float32)
        self.jet2_SVs = np.zeros((self.batch_size, self.n_SVs, 6), dtype=np.float32)
        self.jet1_extraInfo = np.zeros((self.batch_size, 7), dtype=np.float32)
        self.jet2_extraInfo = np.zeros((self.batch_size, 7), dtype=np.float32)
        self.jet_kinematics = np.zeros((self.batch_size, 14), dtype=np.float32)
        self.event_info = np.zeros((self.batch_size, 8), dtype=np.float32)
        self.sys_weights = np.zeros((self.batch_size, 19), dtype=np.float32)
        self.jet1_JME_vars = np.zeros((self.batch_size, 12), dtype=np.float32)
        self.jet2_JME_vars = np.zeros((self.batch_size, 12), dtype=np.float32)


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

    def get_weight_avgs(self, inTree):
        #avg across whole sample (before preselection), allows effect on preselection efficiency to be accounted for
        ROOT.gROOT.SetBatch(True)

        self.avg_weights = dict()

        sys_branch_names = ["Pileup__nom", "Pileup__up", "Pileup__down", 
                "lead_sjbtag_corr__nom",  "lead_sjbtag_corr__up",  "lead_sjbtag_corr__down",  "sublead_sjbtag_corr__nom",  "sublead_sjbtag_corr__up",  "sublead_sjbtag_corr__down",  
                "Pdfweight__up", "Pdfweight__down", 
                "PSWeight[0]",  "PSWeight[1]", "PSWeight[2]", "PSWeight[3]",
                ]

        #read number of scale variations
        event = Event(inTree, 0)
        nLHEScale = inTree.readBranch("nLHEScaleWeight")
        for i in range(nLHEScale): sys_branch_names.append("LHEScaleWeight[%i]" % i)
        if(self.year == 2016 or self.year == 2017):
            sys_branch_names += ['Prefire__nom', 'Prefire__up', 'Prefire__down']

        print("Avg. weights: ")
        for sys_branch in sys_branch_names:
            sys_mean = get_branch_mean(inTree, sys_branch)
            print(sys_branch, sys_mean)
            self.avg_weights[sys_branch] = sys_mean

    def normed_weight(self, tree, branch):
        return tree.readBranch(branch) / self.avg_weights[branch]
            
    
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
        SVs = Collection(event, 'FatJetSVs')

        event_info = [eventNum, MET, MET_phi, genWeight, leptonic_decay, run, self.year, num_jets]

        jet1.pt_corr = jet1.pt
        jet2.pt_corr = jet2.pt

        jet1.msoftdrop_corr = jet1.msoftdrop
        jet2.msoftdrop_corr = jet2.msoftdrop


        sys_weights = []
        jet1_JME_vars = []
        jet2_JME_vars = []
        if(self.include_systematics):
            #These branches are not in the regualr PFNano/Pancakes should have been added by TIMBER

            #PDF's
            pdf_up = self.normed_weight(inTree, 'Pdfweight__up')
            pdf_down = self.normed_weight(inTree, 'Pdfweight__down')
            
            #Prefire
            if(self.year == 2016 or self.year == 2017):
                prefire_nom = self.normed_weight(inTree, 'Prefire__nom')
                prefire_up = self.normed_weight(inTree, 'Prefire__up') / prefire_nom
                prefire_down = self.normed_weight(inTree, 'Prefire__down') / prefire_nom
            else:
                prefire_nom = prefire_up = prefire_down = 1.0

            #Pileup
            pileup_nom = self.normed_weight(inTree, 'Pileup__nom')
            pileup_up = self.normed_weight(inTree, 'Pileup__up') / pileup_nom
            pileup_down = self.normed_weight(inTree, 'Pileup__down') / pileup_nom

            #btag TODO divide out average weight so normalization unchanged?
            btag_nom =  self.normed_weight(inTree, 'lead_sjbtag_corr__nom') * self.normed_weight(inTree, 'sublead_sjbtag_corr__nom')
            btag_up =  self.normed_weight(inTree, 'lead_sjbtag_corr__up') * self.normed_weight(inTree, 'sublead_sjbtag_corr__up') / btag_nom
            btag_down =  self.normed_weight(inTree, 'lead_sjbtag_corr__down') * self.normed_weight(inTree, 'sublead_sjbtag_corr__down') / btag_nom

            #PS weights
            #Older samples don't have
            nPS = inTree.readBranch("nPSWeight")
            PS_ISR_up = PS_ISR_down = PS_FSR_up = PS_FSR_down  = 1.0
            if(nPS > 1):
                #order https://cms-nanoaod-integration.web.cern.ch/integration/cms-swCMSSW_10_6_X/mc106Xul17_doc.html
                PS_weights = inTree.readBranch("PSWeight")
                PS_ISR_up = PS_weights[0] / self.avg_weights['PSWeight[0]']
                PS_FSR_up = PS_weights[1] / self.avg_weights['PSWeight[1]']
                PS_ISR_down = PS_weights[2] / self.avg_weights['PSWeight[2]']
                PS_FSR_down = PS_weights[3] / self.avg_weights['PSWeight[3]']


            #Renorm / Fac weights

            nScale = inTree.readBranch("nLHEScaleWeight")
            RF_down = R_down = F_down = F_up = R_up = RF_up = 1.
            if(nScale == 9):
                #order https://cms-nanoaod-integration.web.cern.ch/integration/cms-swCMSSW_10_6_X/mc106Xul17_doc.html
                scale_weights = inTree.readBranch("LHEScaleWeight")
                
                RF_down = scale_weights[0] / self.avg_weights['LHEScaleWeight[0]']
                R_down = scale_weights[1] / self.avg_weights['LHEScaleWeight[1]']
                F_down = scale_weights[3] / self.avg_weights['LHEScaleWeight[3]']
                F_up = scale_weights[5] / self.avg_weights['LHEScaleWeight[5]']
                R_up = scale_weights[7] / self.avg_weights['LHEScaleWeight[7]']
                RF_up = scale_weights[8] / self.avg_weights['LHEScaleWeight[8]']
            elif(nScale == 8):
                #some files only have 8 weights
                scale_weights = inTree.readBranch("LHEScaleWeight")
                
                RF_down = scale_weights[0] / self.avg_weights['LHEScaleWeight[0]']
                R_down = scale_weights[1] / self.avg_weights['LHEScaleWeight[1]']
                F_down = scale_weights[3] / self.avg_weights['LHEScaleWeight[3]']
                F_up = scale_weights[4] / self.avg_weights['LHEScaleWeight[4]']
                R_up = scale_weights[6] / self.avg_weights['LHEScaleWeight[6]']
                RF_up = scale_weights[7] / self.avg_weights['LHEScaleWeight[7]']




            gen_weight = prefire_nom * pileup_nom * btag_nom
            sys_weights = [gen_weight, pdf_up, pdf_down, prefire_up, prefire_down, pileup_up, pileup_down, btag_up, btag_down, 
            PS_ISR_up, PS_ISR_down, PS_FSR_up, PS_FSR_down, F_up, F_down, R_up, R_down, RF_up, RF_down]

            test = event.lead_sjbtag_corr__vec

            #clip extreme variations
            self.sys_weights[self.idx] = np.clip(np.array(sys_weights, dtype=np.float32), 1e-3, 1e3)

            self.jet1_JME_vars[self.idx] = jet1.JME_vars
            self.jet2_JME_vars[self.idx] = jet2.JME_vars


            






        d_eta = abs(jet1.eta - jet2.eta)




        jet_kinematics = [mjj, d_eta, jet1.pt_corr, jet1.eta, jet1.phi, jet1.msoftdrop_corr, jet2.pt_corr, jet2.eta, jet2.phi, jet2.msoftdrop_corr]

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

        #SV's
        jet1_SVs = []
        jet2_SVs = []

        for SV in SVs:
            SV_vec = [SV.mass, SV.pt, SV.ntracks, SV.normchi2, SV.dxysig, SV.d3dsig]
            if(SV.jetIdx == jet1.idx):
                jet1_SVs.append(SV_vec)
            elif(SV.jetIdx == jet2.idx):
                jet2_SVs.append(SV_vec)

        j1_nSVs = min(len(jet1_SVs), self.n_SVs)
        j2_nSVs = min(len(jet2_SVs), self.n_SVs)

        jet1_SVs = jet1_SVs[:j1_nSVs]
        jet2_SVs = jet2_SVs[:j2_nSVs]

        if(j1_nSVs > 0): self.jet1_SVs[self.idx, :j1_nSVs] = np.array(jet1_SVs, dtype = np.float32)
        if(j2_nSVs > 0): self.jet2_SVs[self.idx, :j2_nSVs] = np.array(jet2_SVs, dtype = np.float32)

        #print(self.jet2_SVs[self.idx])



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
                f.create_dataset("jet1_SVs", data=self.jet1_SVs, chunks = True, maxshape=(None, self.n_SVs, 6), compression='gzip')
                f.create_dataset("jet2_SVs", data=self.jet2_SVs, chunks = True, maxshape=(None, self.n_SVs, 6), compression='gzip')
                if(self.include_systematics):
                    f.create_dataset("sys_weights", data=self.sys_weights, chunks = True, maxshape=(None, self.sys_weights.shape[1]))
                    f.create_dataset("jet1_JME_vars", data=self.jet1_JME_vars, chunks = True, maxshape=(None, self.jet1_JME_vars.shape[1]))
                    f.create_dataset("jet2_JME_vars", data=self.jet2_JME_vars, chunks = True, maxshape=(None, self.jet2_JME_vars.shape[1]))

        else:
            with h5py.File(self.output_name, "a") as f:
                utils.append_h5(f,'truth_label',truth_label_write)
                utils.append_h5(f,'event_info',self.event_info)
                utils.append_h5(f,'jet_kinematics',self.jet_kinematics)
                utils.append_h5(f,'jet1_extraInfo',self.jet1_extraInfo)
                utils.append_h5(f,'jet2_extraInfo',self.jet2_extraInfo)
                utils.append_h5(f,'jet1_PFCands',self.jet1_PFCands)
                utils.append_h5(f,'jet2_PFCands',self.jet2_PFCands)
                utils.append_h5(f,'jet1_SVs',self.jet1_SVs)
                utils.append_h5(f,'jet2_SVs',self.jet2_SVs)
                if(self.include_systematics):
                    utils.append_h5(f,'sys_weights',self.sys_weights)
                    utils.append_h5(f,'jet1_JME_vars',self.jet1_JME_vars)
                    utils.append_h5(f,'jet2_JME_vars',self.jet2_JME_vars)

        self.reset()

    def final_write_out(self, eff, eff_JES_up = -1., eff_JES_down = -1., eff_JER_up = -1., eff_JER_down = -1.):
        if(self.idx < self.batch_size):
            print("Last batch only filled %i events, shortening arrays \n" % self.idx)
            self.jet1_PFCands = self.jet1_PFCands[:self.idx]
            self.jet2_PFCands = self.jet2_PFCands[:self.idx]
            self.jet1_extraInfo = self.jet1_extraInfo[:self.idx]
            self.jet2_extraInfo = self.jet2_extraInfo[:self.idx]
            self.jet_kinematics = self.jet_kinematics[:self.idx] 
            self.event_info = self.event_info[:self.idx]
            self.jet1_SVs = self.jet1_SVs[:self.idx]
            self.jet2_SVs = self.jet2_SVs[:self.idx]
            if(self.include_systematics):
                self.sys_weights = self.sys_weights[:self.idx]
                self.jet1_JME_vars = self.jet1_JME_vars[:self.idx]
                self.jet2_JME_vars = self.jet2_JME_vars[:self.idx]

        self.write_out()
        self.preselection_eff = eff
        with h5py.File(self.output_name, "a") as f:
            f.create_dataset("preselection_eff", data=np.array([eff]))
            if(self.include_systematics):
                f.create_dataset("preselection_eff_JES_up", data=np.array([eff_JES_up]))
                f.create_dataset("preselection_eff_JES_down", data=np.array([eff_JES_down]))
                f.create_dataset("preselection_eff_JER_up", data=np.array([eff_JER_up]))
                f.create_dataset("preselection_eff_JER_down", data=np.array([eff_JER_down]))

    def add_d_eta_eff(self, d_eta_cut = 1.3):
        with h5py.File(self.output_name, "a") as f:
            d_eta = f['jet_kinematics'][:, 1]
            d_eta_mask = d_eta < d_eta_cut
            if(not self.include_systematics):
                d_eta_eff = np.mean(d_eta_mask)
            else:
                weight_pass = np.sum(f['sys_weights'][:,0][d_eta_mask])
                weight_tot = np.sum(f['sys_weights'][:,0])
                d_eta_eff = weight_pass / weight_tot

            print("Delta eta cut (< %.2f) eff is %.3f " % (d_eta_cut, d_eta_eff))
            f.create_dataset("d_eta_eff", data=np.array([d_eta_eff]))

    def normalize_sys_weights(self):
        if(self.include_systematics):
            with h5py.File(self.output_name, "a") as f:
                cur_weights = f['sys_weights'][:]
                normed_weights = np.copy(cur_weights[:])
                weight_avg = np.mean(cur_weights, axis = 0)
                print(weight_avg)

                #renormalize so nominal weight avgs to 1, change preselection eff
                nom_weight_avg = weight_avg[0]
                f['preselection_eff'][0] *= nom_weight_avg
                f['sys_weights'][:,0] /= nom_weight_avg

                self.preselection_eff = f['preselection_eff'][0] 
                

                f['preselection_eff_JES_up'][0] *= nom_weight_avg
                f['preselection_eff_JES_down'][0] *= nom_weight_avg
                f['preselection_eff_JER_up'][0] *= nom_weight_avg
                f['preselection_eff_JER_down'][0] *= nom_weight_avg


                #for key,idx in sys_weights_dict.items():
                
        







def NanoReader(process_flag, inputFileNames=["in.root"], outputFileName="out.root", json = '', year = 2016, nEventsMax = -1, sampleType = "MC", 
        sort_pfcands=True,  include_systematics = False):
    
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

    n_JES_up = n_JES_down = n_JER_up = n_JER_down = 0


#----------------- Begin loop over files ---------------------------------

    for fileName in inputFileNames:

        print("Opening file %s" % fileName)

        inputFile = TFile.Open(fileName)
        if(not inputFile): #check for null pointer
            print("Unable to open file %s, exting \n" % fileName)
            return 1

        #get input tree
        TTree = inputFile.Get("Events")

        # pre-skimming
        if(json != ''):
            elist,jsonFilter = preSkim(TTree, json)

            #number of events to be processed 
            nTotal = elist.GetN() if elist else TTree.GetEntries()
            
            print('Pre-select %d entries out of %s '%(nTotal,TTree.GetEntries()))


            inTree= InputTree(TTree, elist) 
        else:
            nTotal = TTree.GetEntries()
            inTree= InputTree(TTree) 
            print('Running over %i entries \n' % nTotal)

        out = Outputer(outputFileName, truth_label =  process_flag, sample_type=sampleType, sort_pfcands=sort_pfcands, include_systematics = include_systematics, year = year)

        if(include_systematics):
            out.get_weight_avgs(inTree)


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
            if(not passFilter): 
                continue
            
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
            AK4Jets = Collection(event, "Jet")
            subjets = Collection(event, "SubJet")



                


            jet_min_pt = 300
            #keep 2 jets with pt > 200, tight id
            jet1 = jet2 = jet3 =  None
        
            pf_conts_start = 0 #keep track of indices for PF candidates
            jet_index = 0
            num_jets = 0
            for i,jet in enumerate(AK8Jets):
                jet.idx = i
                #jetId : bit1 = loose, bit2 = tight, bit3 = tightLepVeto
                #want tight id
                #select highest two pt jets, keep track of 3rd
                if((jet.jetId & 2 == 2) and abs(jet.eta) < 2.5):
                    jet.PFConstituents_Start = pf_conts_start
                    if(jet.pt > 50): num_jets+=1
                    if((jet1 == None or jet.pt > jet1.pt and jet.pt > 50.)):
                        jet3 = jet2
                        jet2 = jet1
                        jet1 = jet
                    elif((jet2 == None or jet.pt > jet2.pt and jet.pt > 50.)):
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


            jet1_4vec = ROOT.Math.PtEtaPhiMVector(jet1.pt, jet1.eta, jet1.phi, jet1.msoftdrop)
            jet2_4vec = ROOT.Math.PtEtaPhiMVector(jet2.pt, jet2.eta, jet2.phi, jet2.msoftdrop)


            dijet = jet1_4vec + jet2_4vec
            mjj = dijet.M()



            if(include_systematics):
                #JME Corrections 
                dijet_idx1  = inTree.readBranch("DijetIdx1")
                dijet_idx2  = inTree.readBranch("DijetIdx2")

                if(dijet_idx1 != jet1.idx or dijet_idx2 != jet2.idx):
                    print("Dijet indices from TIMBER and this selector don't match!")
                    print("TIMBER : %i %i. This : %i %i " %(dijet_idx1, dijet_idx2, jet1.idx, jet2.idx))
                    print(jet1.msoftdrop, jet2.msoftdrop)
                    sys.exit(1)


                #nominal
                jet1.pt_corr = inTree.readBranch("FatJet1_pt_corr")
                jet2.pt_corr = inTree.readBranch("FatJet2_pt_corr")

                jet1.msoftdrop_corr = inTree.readBranch("FatJet1_msoftdrop_corr")
                jet2.msoftdrop_corr = inTree.readBranch("FatJet2_msoftdrop_corr")


                #systematics
                jet1_pt_JES_up = inTree.readBranch("FatJet1_pt_JES_up")
                jet2_pt_JES_up = inTree.readBranch("FatJet2_pt_JES_up")
                jet1_msoftdrop_JES_up = inTree.readBranch("FatJet1_msoftdrop_JES_up")
                jet2_msoftdrop_JES_up = inTree.readBranch("FatJet2_msoftdrop_JES_up")

                jet1_pt_JES_down = inTree.readBranch("FatJet1_pt_JES_down")
                jet2_pt_JES_down = inTree.readBranch("FatJet2_pt_JES_down")
                jet1_msoftdrop_JES_down = inTree.readBranch("FatJet1_msoftdrop_JES_down")
                jet2_msoftdrop_JES_down = inTree.readBranch("FatJet2_msoftdrop_JES_down")

                jet1_pt_JER_up = inTree.readBranch("FatJet1_pt_JER_up")
                jet2_pt_JER_up = inTree.readBranch("FatJet2_pt_JER_up")
                jet1_msoftdrop_JER_up = inTree.readBranch("FatJet1_msoftdrop_JER_up")
                jet2_msoftdrop_JER_up = inTree.readBranch("FatJet2_msoftdrop_JER_up")

                jet1_pt_JER_down = inTree.readBranch("FatJet1_pt_JER_down")
                jet2_pt_JER_down = inTree.readBranch("FatJet2_pt_JER_down")
                jet1_msoftdrop_JER_down = inTree.readBranch("FatJet1_msoftdrop_JER_down")
                jet2_msoftdrop_JER_down = inTree.readBranch("FatJet2_msoftdrop_JER_down")

                jet1_msoftdrop_JMS_up = inTree.readBranch("FatJet1_msoftdrop_JMS_up")
                jet2_msoftdrop_JMS_up = inTree.readBranch("FatJet2_msoftdrop_JMS_up")

                jet1_msoftdrop_JMS_down = inTree.readBranch("FatJet1_msoftdrop_JMS_down")
                jet2_msoftdrop_JMS_down = inTree.readBranch("FatJet2_msoftdrop_JMS_down")

                jet1_msoftdrop_JMR_up = inTree.readBranch("FatJet1_msoftdrop_JMR_up")
                jet2_msoftdrop_JMR_up = inTree.readBranch("FatJet2_msoftdrop_JMR_up")

                jet1_msoftdrop_JMR_down = inTree.readBranch("FatJet1_msoftdrop_JMR_down")
                jet2_msoftdrop_JMR_down = inTree.readBranch("FatJet2_msoftdrop_JMR_down")



                jet1.JME_vars = [jet1_pt_JES_up, jet1_msoftdrop_JES_up, jet1_pt_JES_down, jet1_msoftdrop_JES_down, 
                               jet1_pt_JER_up, jet1_msoftdrop_JER_up, jet1_pt_JER_down, jet1_msoftdrop_JER_down,
                               jet1_msoftdrop_JMS_up, jet1_msoftdrop_JMS_down, jet1_msoftdrop_JMR_up, jet1_msoftdrop_JMR_down]

                jet2.JME_vars = [jet2_pt_JES_up, jet2_msoftdrop_JES_up, jet2_pt_JES_down, jet2_msoftdrop_JES_down, 
                               jet2_pt_JER_up, jet2_msoftdrop_JER_up, jet2_pt_JER_down, jet2_msoftdrop_JER_down,
                               jet2_msoftdrop_JMS_up, jet2_msoftdrop_JMS_down, jet2_msoftdrop_JMR_up, jet2_msoftdrop_JMR_down]



                #compute mjj for JES, JER variations
                jet1_4vec_JES_up = ROOT.Math.PtEtaPhiMVector(jet1_pt_JES_up, jet1.eta, jet1.phi, jet1.msoftdrop_corr)
                jet2_4vec_JES_up = ROOT.Math.PtEtaPhiMVector(jet2_pt_JES_up, jet2.eta, jet2.phi, jet2.msoftdrop_corr)
                mjj_JES_up = (jet1_4vec_JES_up + jet2_4vec_JES_up).M()

                jet1_4vec_JES_down = ROOT.Math.PtEtaPhiMVector(jet1_pt_JES_down, jet1.eta, jet1.phi, jet1.msoftdrop_corr)
                jet2_4vec_JES_down = ROOT.Math.PtEtaPhiMVector(jet2_pt_JES_down, jet2.eta, jet2.phi, jet2.msoftdrop_corr)
                mjj_JES_down = (jet1_4vec_JES_down + jet2_4vec_JES_down).M()

                jet1_4vec_JER_up = ROOT.Math.PtEtaPhiMVector(jet1_pt_JER_up, jet1.eta, jet1.phi, jet1.msoftdrop_corr)
                jet2_4vec_JER_up = ROOT.Math.PtEtaPhiMVector(jet2_pt_JER_up, jet2.eta, jet2.phi, jet2.msoftdrop_corr)
                mjj_JER_up = (jet1_4vec_JER_up + jet2_4vec_JER_up).M()

                jet1_4vec_JER_down = ROOT.Math.PtEtaPhiMVector(jet1_pt_JER_down, jet1.eta, jet1.phi, jet1.msoftdrop_corr)
                jet2_4vec_JER_down = ROOT.Math.PtEtaPhiMVector(jet2_pt_JER_down, jet2.eta, jet2.phi, jet2.msoftdrop_corr)
                mjj_JER_down = (jet1_4vec_JER_down + jet2_4vec_JER_down).M()


                #compute if modified 4-vecs pass pre-selection
                if(mjj_JES_up >= mjj_cut and jet1_pt_JES_up >= jet_min_pt and jet2_pt_JES_up >= jet_min_pt): n_JES_up +=1
                if(mjj_JES_down >= mjj_cut and jet1_pt_JES_down >= jet_min_pt and jet2_pt_JES_down >= jet_min_pt): n_JES_down +=1
                if(mjj_JER_up >= mjj_cut and jet1_pt_JER_up >= jet_min_pt and jet2_pt_JER_up >= jet_min_pt): n_JER_up +=1
                if(mjj_JER_down >= mjj_cut and jet1_pt_JER_down >= jet_min_pt and jet2_pt_JER_down >= jet_min_pt): n_JER_down +=1
                

                if(jet2.msoftdrop_corr > jet1.msoftdrop_corr):
                    #if corrected mass of jet2 is now larger than jet1, swap
                    temp = jet1
                    jet1 = jet2
                    jet2 = temp


                #recompute mjj for nominal case
                jet1_4vec = ROOT.Math.PtEtaPhiMVector(jet1.pt_corr, jet1.eta, jet1.phi, jet1.msoftdrop_corr)
                jet2_4vec = ROOT.Math.PtEtaPhiMVector(jet2.pt_corr, jet2.eta, jet2.phi, jet2.msoftdrop_corr)

                #print(jet1.msoftdrop, jet1.msoftdrop_corr, jet2.msoftdrop, jet2.msoftdrop_corr)



                dijet = jet1_4vec + jet2_4vec
                mjj = dijet.M()



            if(mjj< mjj_cut or jet1.pt < jet_min_pt or jet2.pt < jet_min_pt): continue
        

            saved+=1
            out.fill_event(inTree, event, jet1, jet2, jet3, PFCands, subjets, mjj, num_jets)
            if(nEventsMax > 0 and saved >= nEventsMax): break
# -------- End Loop over tree-------------------------------------
# -------- End Loop over files-------------------------------------

    efficiency = float(saved)/count
    efficiency_JES_up = efficiency_JES_down = efficiency_JER_up = efficiency_JER_down = -1.

    if(include_systematics):
        efficiency_JES_up = float(n_JES_up)/count
        efficiency_JES_down = float(n_JES_down)/count
        efficiency_JER_up = float(n_JER_up)/count
        efficiency_JER_down = float(n_JER_down)/count


    out.final_write_out(efficiency, efficiency_JES_up, efficiency_JES_down, efficiency_JER_up, efficiency_JER_down)
    out.normalize_sys_weights()
    out.add_d_eta_eff()
    print("Done. Selected %i events. Selection efficiency is %.3f \n" % (saved, out.preselection_eff))
    if(include_systematics):
        with h5py.File(out.output_name, "r") as f:
                print("Eff JES_up %.3f " % f['preselection_eff_JES_up'][0] )
                print("Eff JES_down %.3f " % f['preselection_eff_JES_down'][0] )
                print("Eff JER_up %.3f " % f['preselection_eff_JER_up'][0])
                print("Eff JER_down %.3f " % f['preselection_eff_JER_down'][0] )

    print("Outputed to %s" % outputFileName)
    return saved


