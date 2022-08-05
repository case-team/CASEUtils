# Read nanoAOD with PF constituents (aka pancakes), apply a pre-selection and output to an H5 file format
from H5_maker import *
from array import array
import correctionlib



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

def rel_pt(mu, jet):
    mu_vec = ROOT.TLorentzVector()
    jet_vec = ROOT.TLorentzVector()

    mu_vec.SetPtEtaPhiM(mu.pt, mu.eta, mu.phi, 0.)
    mu_vec.SetPz(0.) # pt vec only (?)
    jet_vec.SetPtEtaPhiM(jet.pt, jet.eta, jet.phi, jet.mass)
    return mu_vec.Perp(jet_vec.Vect())





class Outputer_TTbar(Outputer):
    def __init__(self, outputFileName="out.root", batch_size = 5000, truth_label = 0, sample_type="MC", 
            sort_pfcands = False, include_systematics = True, do_top_ptrw = False, year = 2017):

        self.batch_size = batch_size
        self.output_name = outputFileName
        self.sample_type = sample_type
        self.first_write = False
        self.truth_label = np.array([[truth_label]]*batch_size, dtype=np.int8)
        self.idx = 0
        self.nBatch = 0
        self.n_pf_cands = 100 #how many PF candidates to save (max)
        self.include_systematics = include_systematics
        self.do_top_ptrw = do_top_ptrw
        self.top_weights = []
        self.sort_pfcands = sort_pfcands
        self.year = year
        self.reset()

    def reset(self):
        self.idx = 0
        self.jet1_PFCands = np.zeros((self.batch_size, self.n_pf_cands,4), dtype=np.float16)
        self.jet1_extraInfo = np.zeros((self.batch_size, 9), dtype=np.float32)
        self.jet_kinematics = np.zeros((self.batch_size, 4), dtype=np.float32)
        self.btag_jet_info = np.zeros((self.batch_size, 5), dtype=np.float32)
        self.mu_info = np.zeros((self.batch_size, 4), dtype=np.float32)
        self.event_info = np.zeros((self.batch_size, 6), dtype=np.float32)
        self.sys_weights = np.zeros((self.batch_size, 21), dtype=np.float32)
        self.jet1_JME_vars = np.zeros((self.batch_size, 12), dtype=np.float32)


    
    def fill_event(self, inTree, event, jet1, sel_mu, btag_jet):
        #jet1 is ak8 jet

        if self.sample_type == "data":
            genWeight = 1
        else:
            genWeight = inTree.readBranch('genWeight')
        
        MET = inTree.readBranch('MET_pt')
        MET_phi = inTree.readBranch('MET_phi')
        eventNum = inTree.readBranch('event')
        run = inTree.readBranch('run')
        subjets = Collection(event, "SubJet")
        PFCands = list(Collection(event, "PFCands"))
        PFCandsIdxs = list(Collection(event, "FatJetToPFCands"))

        #try:
        #    mytest = event.FatJetPFCands_eta
        #except:
        #    FullPFCands = Collection(event, "PFCands")
        #    for cand in PFCands:
        #        cand.eta = FullPFCands[cand.pFCandsIdx].eta
        #        cand.phi = FullPFCands[cand.pFCandsIdx].phi
        #        cand.mass = FullPFCands[cand.pFCandsIdx].mass

        event_info = [eventNum, MET, MET_phi, genWeight, run, self.year]

        jet1.pt_corr = jet1.pt

        jet1.msoftdrop_corr = jet1.msoftdrop


        sys_weights = []
        jet1_JME_vars = []
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

            top_ptrw_nom = top_ptrw_up = top_ptrw_down = 1.0
            if(self.do_top_ptrw):
                top_ptrw_nom = inTree.readBranch("TptReweight__nom")
                top_ptrw_up = inTree.readBranch("TptReweight__up") / top_ptrw_nom
                top_ptrw_down = inTree.readBranch("TptReweight__down") / top_ptrw_nom
                self.top_weights.append(top_ptrw_nom)






            gen_weight = prefire_nom * pileup_nom * btag_nom * top_ptrw_nom
            sys_weights = [gen_weight, pdf_up, pdf_down, prefire_up, prefire_down, pileup_up, pileup_down, btag_up, btag_down, 
            PS_ISR_up, PS_ISR_down, PS_FSR_up, PS_FSR_down, F_up, F_down, R_up, R_down, RF_up, RF_down, top_ptrw_up, top_ptrw_down]

            test = event.lead_sjbtag_corr__vec

            #clip extreme variations
            self.sys_weights[self.idx] = np.clip(np.array(sys_weights, dtype=np.float32), 1e-3, 1e3)

            self.jet1_JME_vars[self.idx] = jet1.JME_vars


            

        jet_kinematics = [jet1.pt_corr, jet1.eta, jet1.phi, jet1.msoftdrop_corr]
        btag_jet_info = [btag_jet.pt, btag_jet.eta, btag_jet.phi, btag_jet.mass, btag_jet.btagDeepB]
        mu_info = [sel_mu.pt, sel_mu.eta, sel_mu.phi, sel_mu.charge]

        
        #maximum deepcsv from top 2 subjets of the fatjet
        jet1_btag = -1.
        
        if(jet1.subJetIdx1 >= 0):
            jet1_btag = subjets[jet1.subJetIdx1].btagDeepB
        if(jet1.subJetIdx2 >= 0):
            jet1_btag = max(jet1_btag, subjets[jet1.subJetIdx2].btagDeepB)

        jet1_extraInfo = [jet1.tau1, jet1.tau2, jet1.tau3, jet1.tau4, jet1.lsf3, jet1_btag, jet1.nPFConstituents, jet1.deepTagMD_H4qvsQCD, jet1.deepTagMD_WvsQCD]

        j1_nPF = min(self.n_pf_cands, jet1.nPFConstituents)
        range1 = PFCandsIdxs[jet1.pf_cands_start : jet1.pf_cands_start + j1_nPF] # indices of pf cands

        jet1_PFCands = []
        for i,conv in enumerate(range1):
            idx = conv.candIdx
            if(i > j1_nPF): break
            cand = ROOT.Math.PtEtaPhiMVector(PFCands[idx].pt, PFCands[idx].eta, PFCands[idx].phi, PFCands[idx].mass)
            jet1_PFCands.append([cand.Px(), cand.Py(), cand.Pz(), cand.E()])


        #SV's

        self.event_info[self.idx] = np.array(event_info, dtype=np.float32)
        self.jet_kinematics[self.idx] = np.array(jet_kinematics, dtype = np.float32)
        self.jet1_extraInfo[self.idx] = np.array(jet1_extraInfo, dtype = np.float32)
        self.mu_info[self.idx] = np.array(mu_info, dtype = np.float32)
        self.btag_jet_info[self.idx] = np.array(btag_jet_info, dtype = np.float32)
        
        # sort PFCands by pt
        if self.sort_pfcands:
            self.jet1_PFCands[self.idx,:jet1.nPFConstituents] = self.get_pfcands_sorted(np.array(jet1_PFCands, dtype = np.float32))
        else:
            self.jet1_PFCands[self.idx,:jet1.nPFConstituents] = np.array(jet1_PFCands, dtype = np.float16)

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
                f.create_dataset("btag_jet_info", data=self.btag_jet_info, chunks = True, maxshape=(None, self.btag_jet_info.shape[1]))
                f.create_dataset("mu_info", data=self.mu_info, chunks = True, maxshape=(None, self.mu_info.shape[1]))
                f.create_dataset("jet1_extraInfo", data=self.jet1_extraInfo, chunks = True, maxshape=(None, self.jet1_extraInfo.shape[1]))
                f.create_dataset("jet1_PFCands", data=self.jet1_PFCands, chunks = True, maxshape=(None, self.jet1_PFCands.shape[1], 4), compression='gzip')
                if(self.include_systematics):
                    f.create_dataset("sys_weights", data=self.sys_weights, chunks = True, maxshape=(None, self.sys_weights.shape[1]))
                    f.create_dataset("jet1_JME_vars", data=self.jet1_JME_vars, chunks = True, maxshape=(None, self.jet1_JME_vars.shape[1]))

        else:
            with h5py.File(self.output_name, "a") as f:
                utils.append_h5(f,'truth_label',truth_label_write)
                utils.append_h5(f,'event_info',self.event_info)
                utils.append_h5(f,'jet_kinematics',self.jet_kinematics)
                utils.append_h5(f,'jet1_extraInfo',self.jet1_extraInfo)
                utils.append_h5(f,'jet1_PFCands',self.jet1_PFCands)
                utils.append_h5(f, 'btag_jet_info', self.btag_jet_info)
                utils.append_h5(f, 'mu_info', self.mu_info)
                if(self.include_systematics):
                    utils.append_h5(f,'sys_weights',self.sys_weights)
                    utils.append_h5(f,'jet1_JME_vars',self.jet1_JME_vars)

        self.reset()

    def final_write_out(self, eff, eff_JES_up = -1., eff_JES_down = -1., eff_JER_up = -1., eff_JER_down = -1.):
        if(self.idx < self.batch_size):
            print("Last batch only filled %i events, shortening arrays \n" % self.idx)
            self.jet1_PFCands = self.jet1_PFCands[:self.idx]
            self.jet1_extraInfo = self.jet1_extraInfo[:self.idx]
            self.jet_kinematics = self.jet_kinematics[:self.idx] 
            self.btag_jet_info = self.btag_jet_info[:self.idx]
            self.mu_info = self.mu_info[:self.idx]
            self.event_info = self.event_info[:self.idx]
            if(self.include_systematics):
                self.sys_weights = self.sys_weights[:self.idx]
                self.jet1_JME_vars = self.jet1_JME_vars[:self.idx]

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

                if(self.do_top_ptrw):
                    avg_top_ptrw = np.mean(self.top_weights)
                    f.create_dataset("top_ptrw_avg", data=np.array([avg_top_ptrw]))




def NanoReader_TTbar(process_flag, inputFileNames=["in.root"], outputFileName="out.root", json = '', year = 2016, nEventsMax = -1, sampleType = "MC", 
        sort_pfcands=True,  include_systematics = False, do_top_ptrw = False):
    
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
    #"Flag_BadPFMuonDzFilter",
    "Flag_eeBadScFilter", 
    ]
    if((year == 2016) or (year == 2016.5)): filters.append("Flag_CSCTightHaloFilter")
    if(year == 2017 or year == 2018): filters.append("Flag_ecalBadCalibFilter")

    triggers = ["HLT_Mu50"]
    if(year == 2017 or year == 2018):
        triggers += ["HLT_TkMu100", "HLT_OldMu100"]
    else:
        triggers += ["HLT_TkMu50"]

    btag_cut = -1.

    #deepcsv medium WP's https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation
    if(year == 2018): btag_cut = 0.4168
    elif(year == 2017): btag_cut = 0.4506
    elif(year == 2016): btag_cut = 0.5847 #postVFP
    #elif(year == 2016APV): btag_cut = 0.6001 #preVFP



    nFiles = len(inputFileNames)
    print("Will run over %i files and output to %s with truth label %i" % (nFiles, outputFileName, process_flag))
    count = 0
    saved = 0

    n_JES_up = n_JES_down = n_JER_up = n_JER_down = 0

    out = Outputer_TTbar(outputFileName, truth_label =  process_flag, sample_type=sampleType, sort_pfcands=sort_pfcands, 
            include_systematics = include_systematics, year = year, do_top_ptrw = do_top_ptrw)

    if(include_systematics):
        out.get_weight_avgs(inTree)


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



                        
            
            AK8Jets = Collection(event, "FatJet")
            AK4Jets = Collection(event, "Jet")
            Mus = Collection(event, "Muon")
            nPVs = inTree.readBranch("PV_npvsGood")

            MET = inTree.readBranch('MET_pt')
            MET_phi = inTree.readBranch('MET_phi')

            if(len(Mus) == 0 or len(AK8Jets) == 0 or len(AK4Jets) == 0): continue

                


            ak4_min_pt = 30.
            ak8_min_pt = 200.

            pf_conts_start = 0 #keep track of indices for PF candidates
            jet_index = 0
            num_jets = 0

            sel_mu = None
            nMu = 0
            for mu in Mus:
                if(mu.tightId and mu.pt > 60. and abs(mu.eta) < 2.4):
                    nMu +=1
                    if(sel_mu is None): sel_mu = mu

            #cut on MET and muons
            if((sel_mu is None) or (nMu > 1) or (MET < 50.) or (MET + sel_mu.pt < 250.) or nPVs < 1 ): continue

            rel_pt_cut  = 25.
            ang_cut = 2. * ROOT.TMath.Pi() / 3.
            closest_jet = None
            min_jet_dR = 99999.
            nAK4s = 0
            j2_ak4 = None
            pass_btag = False
            for jet in AK4Jets:
                #jetId : bit1 = loose, bit2 = tight, bit3 = tightLepVeto
                #if(jet.jetId & 2 == 2 and jet.pt > ak4_min_pt and abs(jet.eta) < 2.4):
                if(jet.pt > ak4_min_pt and abs(jet.eta) < 2.4):
                    jet_dR = deltaR(jet, sel_mu)
                    nAK4s +=1
                    if jet_dR < min_jet_dR:
                        min_jet_dR = jet_dR
                        closest_jet = jet
                    #tightId and loose Pileup ID
                    if (jet.jetId & 2 == 2 and jet.puId % 2 == 1 and (ang_dist(sel_mu.phi, jet.phi)  < ang_cut) and jet.btagDeepB > btag_cut):
                        pass_btag = True
                        btag_jet = jet


            if(closest_jet is None): continue
            #print(min_jet_dR, rel_pt(sel_mu, closest_jet))
            mu_iso = (closest_jet is not None) and ((min_jet_dR > 0.4) or sel_mu.jetPtRelv2 > rel_pt_cut)

            ak4_cuts = nAK4s >= 2 and pass_btag





            j1_ak8 = None
            pf_cands_start = 0

            for i,jet in enumerate(AK8Jets):
                jet.pf_cands_start = pf_cands_start
                pf_cands_start += jet.nPFCand
                jet.idx = i
                #jetId : bit1 = loose, bit2 = tight, bit3 = tightLepVeto
                #want tight id
                if((jet.jetId & 2 == 2) and abs(jet.eta) < 2.4):
                    jet.PFConstituents_Start = pf_conts_start
                    if(jet.pt > 50): num_jets+=1
                    if((j1_ak8 is None or jet.pt > j1_ak8.pt) and jet.pt > 50. and ang_dist(jet.phi, sel_mu.phi) > ang_cut):
                        j1_ak8 = jet
                
                    jet.nPFConstituents = jet.nPFCand
                
                jet_index += 1
            
            
            

            ak8_cuts = (j1_ak8 is not None) and (j1_ak8.pt > ak8_min_pt) and not inHEMRegion(j1_ak8, year)
            #print(mu_iso, ak4_cuts, ak8_cuts)

            if(not mu_iso or not ak4_cuts or not ak8_cuts): continue


            saved+=1
            out.fill_event(inTree, event, j1_ak8, sel_mu, btag_jet)
            if(nEventsMax > 0 and saved >= nEventsMax): break
        print("Saved %i events" % saved)
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
            if(do_top_ptrw):
                print("Avg top pt rewight %.3f " % f['top_ptrw_avg'][0])
                

    print("Outputed to %s" % outputFileName)
    return saved


