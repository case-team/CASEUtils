from ttbar_h5_maker import *

def nGenCandCounter(index, event):
    count = 0
    jet_indices = event.GenFatJetCands_jetIdx
    length = len(list(event.GenFatJetCands_jetIdx))
    for i in range(length):
        if jet_indices[i] == index:
            count += 1
    
    return count 


class Outputer_Gen_TTbar(Outputer):
    def __init__(self, outputFileName="out.root", batch_size = 5000, truth_label = 0, sample_type="MC", 
            sort_pfcands = False, do_top_ptrw = False, year = "2018", herwig = False):

        self.batch_size = batch_size
        self.herwig = herwig
        self.output_name = outputFileName
        self.sample_type = sample_type
        self.first_write = False
        self.truth_label = np.array([[truth_label]]*batch_size, dtype=np.int8)
        self.idx = 0
        self.nBatch = 0
        self.n_pf_cands = 100 #how many PF candidates to save (max)
        self.do_top_ptrw = do_top_ptrw
        self.top_weights = []
        self.sort_pfcands = sort_pfcands
        self.year = year

        self.reset()

    def reset(self):
        self.idx = 0
        self.jet1_PFCands = np.zeros((self.batch_size, self.n_pf_cands,6), dtype=np.float16)
        self.jet_kinematics = np.zeros((self.batch_size, 4), dtype=np.float32)
        self.btag_jet_info = np.zeros((self.batch_size, 4), dtype=np.float32)
        self.mu_info = np.zeros((self.batch_size, 4), dtype=np.float32)
        self.event_info = np.zeros((self.batch_size, 6), dtype=np.float32)
        self.sys_weights = np.zeros((self.batch_size, 29), dtype=np.float32)
        self.gen_parts = np.zeros((self.batch_size, 28), dtype=np.float32)


    
    def fill_event(self, inTree, event, jet1, sel_mu, btag_jet, neutrino):
        #jet1 is ak8 jet

        if self.sample_type == "data":
            genWeight = 1
        else:
            genWeight = inTree.readBranch('genWeight')
        
        eventNum = inTree.readBranch('event')
        run = inTree.readBranch('run')
        PFCands = list(Collection(event, "GenCands"))
        PFCandsIdxs = list(Collection(event, "GenFatJetCands"))

        year_val = 2016.5 if 'APV' in self.year else int(self.year)

        event_info = [eventNum, neutrino.pt, neutrino.phi, genWeight, run, year_val]


        sys_weights = []
        jet1_JME_vars = []

        gen_parts = np.zeros(self.gen_parts.shape[1], dtype = np.float32)


        top_ptrw_nom = top_ptrw_up = top_ptrw_down = 1.0

            
        #save gen particles
        top, anti_top, W, anti_W, fermion, anti_fermion, b_quark, _,_,_ = get_ttbar_gen_parts(event, jet1, herwig = self.herwig)

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


        jet_kinematics = [jet1.pt, jet1.eta, jet1.phi, jet1.mass]
        btag_jet_info = [btag_jet.pt, btag_jet.eta, btag_jet.phi, btag_jet.mass]
        mu_charge = sel_mu.pdgId < 0
        mu_info = [sel_mu.pt, sel_mu.eta, sel_mu.phi, mu_charge]

        
        j1_nPF = min(self.n_pf_cands, jet1.nConstituents)
        range1 = PFCandsIdxs[jet1.pf_cands_start : jet1.pf_cands_start + j1_nPF] # indices of pf cands

        jet1_PFCands = []
        for i,conv in enumerate(range1):
            idx = conv.pFCandsIdx
            if(i > j1_nPF): break
            cand = ROOT.Math.PtEtaPhiMVector(PFCands[idx].pt, PFCands[idx].eta, PFCands[idx].phi, PFCands[idx].mass)
            jet1_PFCands.append([cand.Px(), cand.Py(), cand.Pz(), cand.E(), 1, 1])

        self.event_info[self.idx] = np.array(event_info, dtype=np.float32)
        self.jet_kinematics[self.idx] = np.array(jet_kinematics, dtype = np.float32)
        self.mu_info[self.idx] = np.array(mu_info, dtype = np.float32)
        self.btag_jet_info[self.idx] = np.array(btag_jet_info, dtype = np.float32)
        self.gen_parts[self.idx] = gen_parts
        
        # sort PFCands by pt
        if self.sort_pfcands:
            self.jet1_PFCands[self.idx,:jet1.nConstituents] = self.get_pfcands_sorted(np.array(jet1_PFCands, dtype = np.float32))
        else:
            self.jet1_PFCands[self.idx,:jet1.nConstituents] = np.array(jet1_PFCands, dtype = np.float32)

        nPS = inTree.readBranch("nPSWeight")
        if(nPS > 1):
            #order https://cms-nanoaod-integration.web.cern.ch/integration/cms-swCMSSW_10_6_X/mc106Xul17_doc.html
            PS_weights = inTree.readBranch("PSWeight")

            self.sys_weights[self.idx, 9] = PS_weights[0] #ISR_up
            self.sys_weights[self.idx, 10]  = PS_weights[2] #ISR_down
            self.sys_weights[self.idx, 11] = PS_weights[1] #FSR_up
            self.sys_weights[self.idx, 12]  = PS_weights[3] #FSR_down

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
                f.create_dataset("jet1_PFCands", data=self.jet1_PFCands, chunks = True, maxshape=(None, self.jet1_PFCands.shape[1], self.jet1_PFCands.shape[2]), compression='gzip')
                f.create_dataset("sys_weights", data=self.sys_weights, chunks = True, maxshape=(None, self.sys_weights.shape[1]))
                f.create_dataset("gen_parts", data=self.gen_parts, chunks = True, maxshape=(None, self.gen_parts.shape[1]), compression='gzip')

        else:
            with h5py.File(self.output_name, "a") as f:
                utils.append_h5(f,'truth_label',truth_label_write)
                utils.append_h5(f,'event_info',self.event_info)
                utils.append_h5(f,'jet_kinematics',self.jet_kinematics)
                utils.append_h5(f, "btag_jet_info", self.btag_jet_info)
                utils.append_h5(f,'jet1_PFCands',self.jet1_PFCands)
                utils.append_h5(f, 'mu_info', self.mu_info)
                utils.append_h5(f,'sys_weights',self.sys_weights)
                utils.append_h5(f, 'gen_parts', self.gen_parts)

        self.reset()

    def final_write_out(self, eff):
        if(self.idx < self.batch_size):
            print("Last batch only filled %i events, shortening arrays \n" % self.idx)
            self.jet1_PFCands = self.jet1_PFCands[:self.idx]
            self.jet_kinematics = self.jet_kinematics[:self.idx] 
            self.btag_jet_info = self.btag_jet_info[:self.idx]
            self.mu_info = self.mu_info[:self.idx]
            self.event_info = self.event_info[:self.idx]
            self.sys_weights = self.sys_weights[:self.idx]
            self.gen_parts = self.gen_parts[:self.idx]

        self.write_out()
        self.preselection_eff = eff
        with h5py.File(self.output_name, "a") as f:
            f.create_dataset("preselection_eff", data=np.array([eff]))

    def add_d_eta_eff(self, d_eta_cut = 1.3):
        with h5py.File(self.output_name, "a") as f:
            d_eta = f['jet_kinematics'][:, 1]
            d_eta_mask = d_eta < d_eta_cut
            d_eta_eff = np.mean(d_eta_mask)

            print("Delta eta cut (< %.2f) eff is %.3f " % (d_eta_cut, d_eta_eff))
            f.create_dataset("d_eta_eff", data=np.array([d_eta_eff]))



def NanoReader_Gen_TTbar(process_flag, inputFileNames=["in.root"], outputFileName="out.root", json = '', year = "2018", nEventsMax = -1, sampleType = "MC", 
        sort_pfcands=True,  do_top_ptrw = False, include_systematics = False, herwig = False):
    


    nFiles = len(inputFileNames)
    print("Will run over %i files and output to %s with truth label %i" % (nFiles, outputFileName, process_flag))
    count = 0
    saved = 0

    out = Outputer_Gen_TTbar(outputFileName, truth_label =  process_flag, sample_type=sampleType, sort_pfcands=sort_pfcands, 
            year = year, do_top_ptrw = do_top_ptrw, herwig = herwig)



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

        nTotal = TTree.GetEntries()
        inTree= InputTree(TTree) 
        print('Running over %i entries \n' % nTotal)

        if(nTotal ==0): continue

        # Grab event tree from nanoAOD
        eventBranch = inTree.GetBranch('event')
        treeEntries = eventBranch.GetEntries()

# -------- Begin Loop over tree-------------------------------------


        entries = inTree.entries
        printed = False
        for entry in range(entries):

            if count % 10000 == 0 :
                print('--------- Processing Event ' + str(count) +'   -- percent complete ' + str(100*count/nTotal/nFiles) + '% -- ')

            count +=1
            # Grab the event
            event = Event(inTree, entry)

                        
            
            AK8Jets = Collection(event, "GenJetAK8")
            AK4Jets = Collection(event, "GenJet")
            gen_parts = Collection(event, "GenPart")
            Mus = Collection(event, "Muon")


            top, anti_top, W, anti_W, fermion1, anti_fermion1, b_quark1, fermion2, anti_fermion2, b_quark2 = get_ttbar_gen_parts(event, None, verbose = False, herwig = herwig)
            if(top is None or anti_top is None): continue

            #PDGID
            MUON = 13 
            ELECTRON = 11

            sel_mu, neutrino = None, None
            if(fermion1 is not None and (abs(fermion1.pdgId) == MUON or  abs(fermion1.pdgId) == ELECTRON)): sel_mu, neutrino, b_quark = fermion1, anti_fermion1, b_quark1
            elif(anti_fermion1 is not None and (abs(anti_fermion1.pdgId) == MUON or abs(anti_fermion1.pdgId) == ELECTRON)): sel_mu, neutrino, b_quark = anti_fermion1, fermion1, b_quark1
            elif(fermion2 is not None and (abs(fermion2.pdgId) == MUON or abs(fermion2.pdgId) == ELECTRON)): sel_mu, neutrino, b_quark = fermion2, anti_fermion2, b_quark2
            elif(anti_fermion2 is not None and (abs(anti_fermion2.pdgId) == MUON or abs(anti_fermion2.pdgId) == ELECTRON)): sel_mu, neutrino, b_quark = anti_fermion2, fermion2, b_quark2

            if( len(AK8Jets) == 0 or len(AK4Jets) == 0): 
                #print("No jets")
                continue


            if(sel_mu is None): 
                #print("No Mu")
                continue

            ak4_min_pt = 25.
            ak8_min_pt = 200.

            pf_conts_start = 0 #keep track of indices for PF candidates
            jet_index = 0
            num_jets = 0
            muon_pt_cut = 60.
            MET_cut = 50.


            W_cand_px = sel_mu.pt * np.cos(sel_mu.phi) + neutrino.pt * np.cos(neutrino.phi)
            W_cand_py = sel_mu.pt * np.sin(sel_mu.phi) + neutrino.pt * np.sin(neutrino.phi)
            W_cand_pt = (W_cand_px**2 + W_cand_py**2)**(0.5)
            #print("W pt %.1f mu pt %.1f MET %.1f" % (W_cand_pt, sel_mu.pt, MET))



            #cut on MET and muons
            if(sel_mu.pt < muon_pt_cut or abs(sel_mu.eta) > 2.4 or neutrino.pt < MET_cut): 
                #print("Mu or MET")
                continue

            ang_cut = 2.
            min_jet_dR = 99999.
            nAK4s = 0
            btag_jet = None
            for jet in AK4Jets:
                if(jet.pt > ak4_min_pt and abs(jet.eta) < 2.4 and abs(ang_dist(sel_mu.phi, jet.phi))  < ang_cut and deltaR(jet, b_quark) < 0.1):
                    btag_jet = jet


            ak4_cuts = btag_jet is not None


            j1_ak8 = None
            pf_cands_start = 0

            for i,jet in enumerate(AK8Jets):
                jet.nConstituents = nGenCandCounter(jet_index, event)
                jet.pf_cands_start = pf_cands_start
                pf_cands_start += jet.nConstituents
                jet.idx = i
                #want tight id
                if(abs(jet.eta) < 2.5 and jet.pt > ak8_min_pt and abs(ang_dist(jet.phi, sel_mu.phi)) > ang_cut and deltaR(jet, W) < 0.4 or deltaR(jet, anti_W) < 0.4):
                    j1_ak8 = jet
                
                jet_index += 1
            
            
            

            ak8_cuts = (j1_ak8 is not None) and (j1_ak8.pt > ak8_min_pt)

            if(not ak4_cuts or not ak8_cuts): 
                #print("Jet cuts")
                continue


            saved+=1
            out.fill_event(inTree, event, j1_ak8, sel_mu, btag_jet, neutrino)
            if(nEventsMax > 0 and saved >= nEventsMax): break
        print("Saved %i events" % saved)
# -------- End Loop over tree-------------------------------------
# -------- End Loop over files-------------------------------------

    efficiency = float(saved)/count

    out.final_write_out(efficiency)
    print("Done. Selected %i events. Selection efficiency is %.3f \n" % (saved, out.preselection_eff))
                
    print("Outputed to %s" % outputFileName)
    return saved
