# Read nanoAOD with PF constituents (aka pancakes), apply a pre-selection and output to an H5 file format
from H5_maker import *
from array import array
from correction_utils import *




def rel_pt(mu, jet):
    mu_vec = ROOT.TLorentzVector()
    jet_vec = ROOT.TLorentzVector()

    mu_vec.SetPtEtaPhiM(mu.pt, mu.eta, mu.phi, 0.)
    mu_vec.SetPz(0.) # pt vec only (?)
    jet_vec.SetPtEtaPhiM(jet.pt, jet.eta, jet.phi, jet.mass)
    return mu_vec.Perp(jet_vec.Vect())





class Outputer_TTbar(Outputer):
    def __init__(self, outputFileName="out.root", batch_size = 5000, truth_label = 0, sample_type="MC", 
            sort_pfcands = False, include_systematics = True, do_top_ptrw = False, year = "2018"):

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
        self.sys_weights = np.zeros((self.batch_size, 27), dtype=np.float32)
        self.jet1_JME_vars = np.zeros((self.batch_size, 12), dtype=np.float32)
        self.gen_parts = np.zeros((self.batch_size, 28), dtype=np.float32)


    
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
        PFCandsIdxs = list(Collection(event, "FatJetPFCands"))

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

        gen_parts = np.zeros(self.gen_parts.shape[1], dtype = np.float32)
        if(self.include_systematics):

            #JME corrections
            jet1.pt_corr = inTree.readBranch("FatJet_pt_nom")[jet1.idx]
            jet1.msoftdrop_corr = inTree.readBranch("FatJet_msoftdrop_nom")[jet1.idx]


            #JME systematics
            jet1_pt_JES_up = inTree.readBranch("FatJet_pt_jesTotalUp")[jet1.idx]
            jet1_msoftdrop_JES_up = inTree.readBranch("FatJet_msoftdrop_jesTotalUp")[jet1.idx]

            jet1_pt_JES_down = inTree.readBranch("FatJet_pt_jesTotalDown")[jet1.idx]
            jet1_msoftdrop_JES_down = inTree.readBranch("FatJet_msoftdrop_jesTotalDown")[jet1.idx]


            jet1_pt_JER_up = inTree.readBranch("FatJet_pt_jerUp")[jet1.idx]
            jet1_msoftdrop_JER_up = inTree.readBranch("FatJet_msoftdrop_jerUp")[jet1.idx]

            jet1_pt_JER_down = inTree.readBranch("FatJet_pt_jerDown")[jet1.idx]
            jet1_msoftdrop_JER_down = inTree.readBranch("FatJet_msoftdrop_jerDown")[jet1.idx]


            jet1_msoftdrop_JMS_up = inTree.readBranch("FatJet_msoftdrop_jmsUp")[jet1.idx]
            jet1_msoftdrop_JMS_down = inTree.readBranch("FatJet_msoftdrop_jmsDown")[jet1.idx]

            jet1_msoftdrop_JMR_up = inTree.readBranch("FatJet_msoftdrop_jmrUp")[jet1.idx]
            jet1_msoftdrop_JMR_down = inTree.readBranch("FatJet_msoftdrop_jmrDown")[jet1.idx]



            jet1.JME_vars = [jet1_pt_JES_up, jet1_msoftdrop_JES_up, jet1_pt_JES_down, jet1_msoftdrop_JES_down, 
                           jet1_pt_JER_up, jet1_msoftdrop_JER_up, jet1_pt_JER_down, jet1_msoftdrop_JER_down,
                           jet1_msoftdrop_JMS_up, jet1_msoftdrop_JMS_down, jet1_msoftdrop_JMR_up, jet1_msoftdrop_JMR_down]




            mu_weights =  get_lepton_weights(sel_mu, self.year)
            #print(sel_mu.pt, sel_mu.eta, mu_weights)


            #Renorm / Fac weights
            nScale = 0
            scale_fail = False
            if(not scale_fail):
                try:
                    nScale = inTree.readBranch("nLHEScaleWeight")
                except:
                    scale_fail = True

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







            #PDF's
            if(not scale_fail): pdf_up, pdf_down = get_pdf_weight(inTree)
            else: pdf_up = pdf_down = 1.0
            
            #Prefire
            if("2016" in self.year or "2017" in self.year):
                prefire_nom = inTree.readBranch("L1PreFiringWeight_Nom")
                prefire_up = inTree.readBranch("L1PreFiringWeight_Up")
                prefire_down = inTree.readBranch("L1PreFiringWeight_Dn")
            else:
                prefire_nom = prefire_up = prefire_down = 1.0

            #Pileup
            pileup_nom, pileup_up, pileup_down = get_pileup_weight(self.year, inTree.readBranch("Pileup_nTrueInt"))

            btag_nom = btag_up = btag_down = 1.0
            btag_nom = inTree.readBranch("Jet_btagSF_deepcsv_M") [jet1.idx]
            btag_up = inTree.readBranch("Jet_btagSF_deepcsv_M_up")[jet1.idx] / btag_nom
            btag_down = inTree.readBranch("Jet_btagSF_deepcsv_M_down")[jet1.idx] / btag_nom


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



            top_ptrw_nom = top_ptrw_up = top_ptrw_down = 1.0
            if(self.do_top_ptrw):
                
                #save gen particles
                top, anti_top, W, anti_W, fermion, anti_fermion, b_quark = get_ttbar_gen_parts(event, jet1)
                top_ptrw_nom, top_ptrw_up, top_ptrw_down = get_top_ptrw(event, top, anti_top)

                match = check_matching(jet1, fermion, anti_fermion, b_quark)
                #print(top, anti_top, W, anti_W, fermion, anti_fermion, b_quark)

                gen_parts = [match, top.pt, top.eta, top.phi, top.mass, 
                             anti_top.pt, anti_top.eta, anti_top.phi, anti_top.mass, 
                             W.pt, W.eta, W.phi, W.mass, 
                             anti_W.pt, anti_W.eta, anti_W.phi, anti_W.mass]

                #add quarks and b if they were found
                if(fermion is not None and anti_fermion is not None):
                    gen_parts += [fermion.pt, fermion.eta, fermion.phi, fermion.pdgId, anti_fermion.pt, anti_fermion.eta, anti_fermion.phi, anti_fermion.pdgId] 
                else: gen_parts += [0.]*6
                if(b_quark is not None): gen_parts += [b_quark.pt, b_quark.eta, b_quark.phi]
                else: gen_parts += [0.]*3



                gen_parts = np.array(gen_parts, dtype = np.float32)

                #print(gen_parts)

                






            gen_weight = prefire_nom * pileup_nom * btag_nom * top_ptrw_nom * mu_weights["nominal"] * np.sign(genWeight)
            sys_weights = [gen_weight, pdf_up, pdf_down, prefire_up, prefire_down, pileup_up, pileup_down, btag_up, btag_down, 
                            PS_ISR_up, PS_ISR_down, PS_FSR_up, PS_FSR_down, F_up, F_down, R_up, R_down, RF_up, RF_down, top_ptrw_up, top_ptrw_down,
                            mu_weights['trigger_up'], mu_weights['trigger_down'], mu_weights['id_up'], mu_weights['id_down'], mu_weights['iso_up'], mu_weights['iso_down']]
            #print(sys_weights)


            #clip extreme variations
            self.sys_weights[self.idx] = np.clip(np.array(sys_weights, dtype=np.float32), 1e-3, 1e3)

            #self.jet1_JME_vars[self.idx] = jet1.JME_vars


            

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
            idx = conv.pFCandsIdx
            if(i > j1_nPF): break
            cand = ROOT.Math.PtEtaPhiMVector(PFCands[idx].pt, PFCands[idx].eta, PFCands[idx].phi, PFCands[idx].mass)
            jet1_PFCands.append([cand.Px(), cand.Py(), cand.Pz(), cand.E()])


        #SV's

        self.event_info[self.idx] = np.array(event_info, dtype=np.float32)
        self.jet_kinematics[self.idx] = np.array(jet_kinematics, dtype = np.float32)
        self.jet1_extraInfo[self.idx] = np.array(jet1_extraInfo, dtype = np.float32)
        self.mu_info[self.idx] = np.array(mu_info, dtype = np.float32)
        self.btag_jet_info[self.idx] = np.array(btag_jet_info, dtype = np.float32)
        self.gen_parts[self.idx] = gen_parts
        
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
                    if(self.do_top_ptrw):
                        f.create_dataset("gen_parts", data=self.gen_parts, chunks = True, maxshape=(None, self.gen_parts.shape[1]), compression='gzip')

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
                    if(self.do_top_ptrw):
                        utils.append_h5(f, 'gen_parts', self.gen_parts)

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
                self.gen_parts = self.gen_parts[:self.idx]

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
                print("Weight avg:")
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

                #if(self.do_top_ptrw):
                    #avg_top_ptrw = np.mean(self.top_weights)
                    #f.create_dataset("top_ptrw_avg", data=np.array([avg_top_ptrw]))




def NanoReader_TTbar(process_flag, inputFileNames=["in.root"], outputFileName="out.root", json = '', year = "2018", nEventsMax = -1, sampleType = "MC", 
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
    if("2016" in year ): filters.append("Flag_CSCTightHaloFilter")
    if("2017" in year or "2018" in year): filters.append("Flag_ecalBadCalibFilter")

    triggers = ["HLT_Mu50"]
    if("2017" in year or "2018" in year):
        triggers += ["HLT_TkMu100", "HLT_OldMu100"]
    else:
        triggers += ["HLT_TkMu50"]

    btag_cut = -1.

    #deepcsv medium WP's https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation
    if("2018" in year): btag_cut = 0.4168
    elif("2017" in year): btag_cut = 0.4506
    elif("2016APV" in year): btag_cut = 0.6001 #preVFP
    elif("2016" in year): btag_cut = 0.5847 #postVFP
    else: 
        print("invalid year! %s")
        exit1(1)



    nFiles = len(inputFileNames)
    print("Will run over %i files and output to %s with truth label %i" % (nFiles, outputFileName, process_flag))
    count = 0
    saved = 0

    n_JES_up = n_JES_down = n_JER_up = n_JER_down = 0

    out = Outputer_TTbar(outputFileName, truth_label =  process_flag, sample_type=sampleType, sort_pfcands=sort_pfcands, 
            include_systematics = include_systematics, year = year, do_top_ptrw = do_top_ptrw)



#----------------- Begin loop over files ---------------------------------

    for fileName in inputFileNames:

        print("Opening file %s" % fileName)

        inputFile = TFile.Open(fileName)
        if(not inputFile): #check for null pointer
            print("Unable to open file %s, skipping \n" % fileName)
            continue

        #get input tree
        try:
            TTree = inputFile.Get("Events")
        except:
            print("Unable to get contents from file %s, skipping \n" % fileName)
            continue

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

        if(nTotal ==0): continue

        if(include_systematics):
            out.get_weight_avgs(inTree, ttbar = True)


        # Grab event tree from nanoAOD
        eventBranch = inTree.GetBranch('event')
        treeEntries = eventBranch.GetEntries()




# -------- Begin Loop over tree-------------------------------------


        #quick check
        event = Event(inTree, 0)
        AK8Jets = Collection(event, "FatJet")

        try:
            b = event.FatJet_nConstituents
            alt_lookup = False
        except:
            alt_lookup = True
            



        entries = inTree.entries
        for entry in range(entries):


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
            
            for trig in triggers:
                passTrigger = passTrigger or inTree.readBranch(trig)

            if(not passTrigger): continue



                        
            
            AK8Jets = Collection(event, "FatJet")
            AK4Jets = Collection(event, "Jet")
            Mus = Collection(event, "Muon")
            TrigObjs = Collection(event, "TrigObj")
            nPVs = inTree.readBranch("PV_npvsGood")

            MET = inTree.readBranch('MET_pt')
            MET_phi = inTree.readBranch('MET_phi')

            if(len(Mus) == 0 or len(AK8Jets) == 0 or len(AK4Jets) == 0): continue

                


            ak4_min_pt = 25.
            ak8_min_pt = 200.

            pf_conts_start = 0 #keep track of indices for PF candidates
            jet_index = 0
            num_jets = 0

            #trigger objs




            #select muon
            sel_mu = None
            nMu = 0
            for mu in Mus:
                if(mu.tightId and mu.pt > 60. and abs(mu.eta) < 2.4 and abs(mu.dxy) <= 0.2 and mu.miniPFRelIso_all < 0.10):
                    nMu +=1
                    if(sel_mu is None): sel_mu = mu

            if(sel_mu is None): continue


            #match muon to trig obj
            #Mu50 is 1024 and Mu100 is 2048 bit
            #https://cms-nanoaod-integration.web.cern.ch/integration/cms-swCMSSW_10_6_X/mc106Xul17v2_doc.html
            trig_obj = None
            min_trig_dR = 99999.
            for t_obj in TrigObjs:
                if(t_obj.id == 13 and t_obj.filterBits > 1024):

                    if(trig_obj is None or deltaR(t_obj, sel_mu) < min_trig_dR): 
                        trig_obj = t_obj
                        min_trig_dR = deltaR(trig_obj, sel_mu)

            if(trig_obj is None): 
                print("No trig obj?")
                continue

            if(deltaR(trig_obj, sel_mu) > 0.15):
                #print("trig reject dR %.3f" % deltaR(trig_obj, sel_mu))
                continue

            W_cand_px = sel_mu.pt * np.cos(sel_mu.phi) + MET * np.cos(MET_phi)
            W_cand_py = sel_mu.pt * np.sin(sel_mu.phi) + MET * np.sin(MET_phi)
            W_cand_pt = (W_cand_px**2 + W_cand_py**2)**(0.5)
            #print("W pt %.1f mu pt %.1f MET %.1f" % (W_cand_pt, sel_mu.pt, MET))



            #cut on MET and muons
            if((nMu > 1) or (MET < 50.) or (MET + sel_mu.pt < 100.) or nPVs < 1): continue

            ang_cut = 2.
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
                    #tightId and loose Pileup ID
                    if (jet.jetId & 2 == 2 and jet.puId % 2 == 1 and (abs(ang_dist(sel_mu.phi, jet.phi))  < ang_cut) and jet.btagDeepB > btag_cut):
                        pass_btag = True
                        btag_jet = jet


            ak4_cuts = nAK4s >= 2 and pass_btag





            j1_ak8 = None
            pf_cands_start = 0

            for i,jet in enumerate(AK8Jets):
                if(alt_lookup): jet.nConstituents = nPFCounter(jet_index, event)
                jet.pf_cands_start = pf_cands_start
                #pf_cands_start += jet.nPFCand
                pf_cands_start += jet.nConstituents
                jet.idx = i
                #jetId : bit1 = loose, bit2 = tight, bit3 = tightLepVeto
                #want tight id
                if((jet.jetId & 2 == 2) and abs(jet.eta) < 2.5):
                    jet.PFConstituents_Start = pf_conts_start
                    if(jet.pt > 50): num_jets+=1
                    if((j1_ak8 is None or jet.pt > j1_ak8.pt) and jet.pt > 50. and abs(ang_dist(jet.phi, sel_mu.phi)) > ang_cut):
                        j1_ak8 = jet
                
                    jet.nPFConstituents = jet.nConstituents
                
                jet_index += 1
            
            
            

            ak8_cuts = (j1_ak8 is not None) and (j1_ak8.pt > ak8_min_pt) and not inHEMRegion(j1_ak8, year)

            if(not ak4_cuts or not ak8_cuts): continue


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
            #if(do_top_ptrw):
                #print("Avg top pt rewight %.3f " % f['top_ptrw_avg'][0])
                

    print("Outputed to %s" % outputFileName)
    return saved


