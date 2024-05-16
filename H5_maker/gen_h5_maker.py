# Read nanoAOD with PF constituents (aka pancakes), apply a pre-selection and output to an H5 file format

from ttbar_h5_maker import *

def nGenCandCounter(index, event):
    count = 0
    jet_indices = event.GenFatJetCands_jetIdx
    length = len(list(event.GenFatJetCands_jetIdx))
    for i in range(length):
        if jet_indices[i] == index:
            count += 1
    
    return count 


class Outputer_Gen(Outputer):
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
        self.n_pf_cands = 300 #how many PF candidates to save (max)
        self.do_top_ptrw = do_top_ptrw
        self.top_weights = []
        self.sort_pfcands = sort_pfcands
        self.year = year

        self.reset()

    def reset(self):
        self.idx = 0
        self.jet1_PFCands = np.zeros((self.batch_size, self.n_pf_cands,6), dtype=np.float16)
        self.jet2_PFCands = np.zeros((self.batch_size, self.n_pf_cands,6), dtype=np.float16)
        self.jet_kinematics = np.zeros((self.batch_size, 8), dtype=np.float32)
        self.event_info = np.zeros((self.batch_size, 6), dtype=np.float32)
        self.sys_weights = np.zeros((self.batch_size, 29), dtype=np.float32)
        #self.gen_info = np.zeros((self.batch_size, 14, 4), dtype=np.float32)
        self.gen_info = np.zeros((self.batch_size, 8, 4), dtype=np.float32)


    
    def fill_event(self, inTree, event, jet1, jet2):
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

        event_info = [eventNum, 0., 0., genWeight, run, year_val]


        sys_weights = []

        gen_info = np.zeros(self.gen_info.shape[1], dtype = np.float32)


        #save gen particles
        daughter1, daughter2, gen_qs = get_Wkk_gen_parts(event, herwig = self.herwig)
        #daughter1, daughter2, gen_qs = get_YtoHH_gen_parts(event, herwig = self.herwig)

        gen_info = [[daughter1.pt, daughter1.eta, daughter1.phi, daughter1.mass],
                    [daughter2.pt, daughter2.eta, daughter2.phi, daughter2.mass]] + gen_qs


        gen_info = np.array(gen_info, dtype = np.float32)


        jet_kinematics = [jet1.pt, jet1.eta, jet1.phi, jet1.mass, jet2.pt, jet2.eta, jet2.phi, jet2.mass]
        
        j1_nPF = min(self.n_pf_cands, jet1.nConstituents)
        j2_nPF = min(self.n_pf_cands, jet2.nConstituents)
        range1 = PFCandsIdxs[jet1.pf_cands_start : jet1.pf_cands_start + jet1.nConstituents] # indices of pf cands
        range2 = PFCandsIdxs[jet2.pf_cands_start : jet2.pf_cands_start + jet2.nConstituents] # indices of pf cands

        jet2_PFCands = []
        jet1_PFCands = []

        for i,conv in enumerate(range1):
            idx = conv.pFCandsIdx
            cand = ROOT.Math.PtEtaPhiMVector(PFCands[idx].pt, PFCands[idx].eta, PFCands[idx].phi, PFCands[idx].mass)
            jet1_PFCands.append([cand.Px(), cand.Py(), cand.Pz(), cand.E(), 1, 1])

        for i,conv in enumerate(range2):
            idx = conv.pFCandsIdx
            cand = ROOT.Math.PtEtaPhiMVector(PFCands[idx].pt, PFCands[idx].eta, PFCands[idx].phi, PFCands[idx].mass)
            jet2_PFCands.append([cand.Px(), cand.Py(), cand.Pz(), cand.E(), 1, 1])


        self.event_info[self.idx] = np.array(event_info, dtype=np.float32)
        self.jet_kinematics[self.idx] = np.array(jet_kinematics, dtype = np.float32)
        self.gen_info[self.idx] = gen_info

        
        # sort PFCands by pt
        if self.sort_pfcands:
            self.jet1_PFCands[self.idx,:j1_nPF] = self.get_pfcands_sorted(np.array(jet1_PFCands, dtype = np.float32))[:j1_nPF]
            self.jet2_PFCands[self.idx,:j2_nPF] = self.get_pfcands_sorted(np.array(jet2_PFCands, dtype = np.float32))[:j2_nPF]
        else:
            self.jet1_PFCands[self.idx,:j1_nPF] = np.array(jet1_PFCands, dtype = np.float32)[:j1_nPF]
            self.jet2_PFCands[self.idx,:j2_nPF] = np.array(jet2_PFCands, dtype = np.float32)[:j2_nPF]

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
                f.create_dataset("jet1_PFCands", data=self.jet1_PFCands, chunks = True, maxshape=(None, self.jet1_PFCands.shape[1], self.jet1_PFCands.shape[2]), compression='gzip')
                f.create_dataset("jet2_PFCands", data=self.jet2_PFCands, chunks = True, maxshape=(None, self.jet2_PFCands.shape[1], self.jet2_PFCands.shape[2]), compression='gzip')
                f.create_dataset("sys_weights", data=self.sys_weights, chunks = True, maxshape=(None, self.sys_weights.shape[1]))
                f.create_dataset("gen_info", data=self.gen_info, chunks = True, maxshape=(None, self.gen_info.shape[1], 4), compression='gzip')

        else:
            with h5py.File(self.output_name, "a") as f:
                utils.append_h5(f,'truth_label',truth_label_write)
                utils.append_h5(f,'event_info',self.event_info)
                utils.append_h5(f,'jet_kinematics',self.jet_kinematics)
                utils.append_h5(f,'jet1_PFCands',self.jet1_PFCands)
                utils.append_h5(f,'jet2_PFCands',self.jet2_PFCands)
                utils.append_h5(f,'sys_weights',self.sys_weights)
                utils.append_h5(f, 'gen_info', self.gen_info)

        self.reset()

    def final_write_out(self, eff):
        if(self.idx < self.batch_size):
            print("Last batch only filled %i events, shortening arrays \n" % self.idx)
            self.jet1_PFCands = self.jet1_PFCands[:self.idx]
            self.jet2_PFCands = self.jet2_PFCands[:self.idx]
            self.jet_kinematics = self.jet_kinematics[:self.idx] 
            self.event_info = self.event_info[:self.idx]
            self.sys_weights = self.sys_weights[:self.idx]
            self.gen_info = self.gen_info[:self.idx]

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



def NanoReader_Gen(process_flag, inputFileNames=["in.root"], outputFileName="out.root", json = '', year = "2018", nEventsMax = -1, sampleType = "MC", 
        sort_pfcands=True,  do_top_ptrw = False, include_systematics = False, herwig = False):
    


    nFiles = len(inputFileNames)
    print("Will run over %i files and output to %s with truth label %i" % (nFiles, outputFileName, process_flag))
    count = 0
    saved = 0

    out = Outputer_Gen(outputFileName, truth_label =  process_flag, sample_type=sampleType, sort_pfcands=sort_pfcands, 
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
            Mus = Collection(event, "Muon")


            #Swap to YtoHH too
            daughter1, daughter2, _ = get_Wkk_gen_parts(event, verbose = False, herwig = herwig)
            #daughter1, daughter2, _ = get_YtoHH_gen_parts(event, verbose = False, herwig = herwig)


            if(daughter1 is None or daughter2 is None): continue

            if( len(AK8Jets) == 0): continue


            ak8_min_pt = 300.

            pf_conts_start = 0 #keep track of indices for PF candidates
            jet_index = 0
            num_jets = 0

            j1_ak8 = j2_ak8 = None
            pf_cands_start = 0

            for i,jet in enumerate(AK8Jets):
                jet.nConstituents = nGenCandCounter(jet_index, event)
                jet.pf_cands_start = pf_cands_start
                pf_cands_start += jet.nConstituents
                jet.idx = i
                if(abs(jet.eta) < 2.5 and jet.pt > ak8_min_pt):
                    if(deltaR(jet, daughter1) < 0.4):
                        j1_ak8 = jet
                    elif(deltaR(jet, daughter2) < 0.4):
                        j2_ak8 = jet
                
                jet_index += 1

            if((j1_ak8 is None) or (j2_ak8 is None)): continue

            saved+=1
            out.fill_event(inTree, event, j1_ak8, j2_ak8)
            if(nEventsMax > 0 and saved >= nEventsMax): break
        print("Saved %i events" % saved)
# -------- End Loop over tree-------------------------------------
# -------- End Loop over files-------------------------------------

    efficiency = float(saved)/count

    out.final_write_out(efficiency)
    print("Done. Selected %i events. Selection efficiency is %.3f \n" % (saved, out.preselection_eff))
                
    print("Outputed to %s" % outputFileName)
    return saved
