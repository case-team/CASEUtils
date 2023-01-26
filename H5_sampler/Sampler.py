import h5py
import numpy as np 
from sys import exit

SCALAR_KEYS = ["preselection_eff", "d_eta_eff", "preselection_eff_JER_down", "preselection_eff_JER_up", "preselection_eff_JES_down", "preselection_eff_JES_up"]


class Sampler():
    def __init__(self, filename, xsec, lumi, holdout_frac = 0., isSignal = False):
        self.filename = filename
        with h5py.File(self.filename, "r") as h5_file:
            if isSignal:
                self.eff_xsec = 1.0
            else:
                self.eff_xsec = xsec * h5_file['preselection_eff'][0]
            self.keys = h5_file.keys()
            self.nEvents = h5_file['event_info'].shape[0]
        self.nSample = int(lumi * self.eff_xsec)
        
        if(isSignal): self.nSample = int(lumi * xsec)
        print("Will sample %i events from file %s" %(self.nSample, self.filename))
        self.n_to_sample_from = self.nEvents
        self.nHoldOut = 0
        self.with_replacement = False
        self.batch_setup = False
        if(holdout_frac > 0.):
            self.nHoldOut = int(self.nSample * holdout_frac)
            self.n_to_sample_from = self.nEvents  - self.nHoldOut
        if(self.nSample > self.n_to_sample_from):
            print("Warning not enough unique events in file %s to match cross section! Need to sample %i events but only have %i saved (%i saved for holdout)). \n" 
                    %(self.filename, self.nSample, self.n_to_sample_from, self.nHoldOut) + 
                    "Will sample with replacement. Each event will be reused an average of %.1f times" % (float(self.nSample) / self.n_to_sample_from)
                   )
            self.with_replacement = True

    def setup_batch(self, iBatch = 0, nBatches = 1):    
        #setup sampler so it can have the same random sample for all the keys in this batch
        
        n_batch_sample = self.nSample // nBatches #how many events to sample for this batch
        n_chunk_size = self.n_to_sample_from // nBatches #size of 'chunk' of events we are sampling from

        self.chunk_start = iBatch * n_chunk_size
        self.chunk_end = (iBatch+1) * n_chunk_size

        self.sample_idxs = np.random.choice(n_chunk_size, n_batch_sample, replace = self.with_replacement)
        self.batch_setup = True

    def sample(self, key):
        if(key not in self.keys):
            print("Error! Key %s not in list of keys for file %s. Keys are: \n" % (key, self.filename))
            print(self.keys)
            exit(1)
        if(not self.batch_setup):
            print("sampler not setup yet! Please call batch_setup before sample")
            exit(1)

        if(self.nSample < 1): return []
        with h5py.File(self.filename, "r") as h5_file:
            out = h5_file[key][self.chunk_start:self.chunk_end][self.sample_idxs]
            return out

    def sample_scalar(self, key):
        if(key not in self.keys):
            print("Error! Key %s not in list of keys for file %s. Keys are: \n" % (key, self.filename))
            print(self.keys)
            exit(1)

        with h5py.File(self.filename, "r") as h5_file:
            out = h5_file[key][:][0]
            return out

    def holdout(self, key):
        if(key not in self.keys):
            print("Error! Key %s not in list of keys for file %s. Keys are: \n" % (key, self.filename))
            print(self.keys)
            exit(1)
        if(self.nHoldOut == 0): return []
        with h5py.File(self.filename, "r") as h5_file:
            out = h5_file[key][self.n_to_sample_from:self.nEvents]
            return out


class MaxMultiSampler():
    """
    This class is essentially a container for multiple Sampler instances.
    It puts together multiple files with their own individual Sampler each.
    The Sampler's have a relative event count given by the given luminosities.
    And Their absolute count is determined by the maximum number of events
    that the input files have while respecting the relative composition.
    """
    def __init__(self, filename_list, lumi_list):
        self.filename_list = filename_list
        self.lumi_list = lumi_list
        assert len(filename_list) == len(lumi_list)

        # get the number of available evs and the normalized target composition
        event_counts = []
        composition = []
        for filename, lumi in zip(filename_list, lumi_list):
            with h5py.File(filename, "r") as h5_file:
                event_counts.append(h5_file['event_info'].shape[0])
                composition.append(h5_file['preselection_eff'][0] * lumi)

        norm_factor = sum(composition)
        for i in range(len(composition)):
            composition[i] /= norm_factor

        self.event_counts = event_counts
        self.composition = composition

        # compute largest possible event number from all samples
        evt_comp_ratios = [event_counts[i]/composition[i]
                           for i in range(len(event_counts))]
        min_ratio = min(evt_comp_ratios)
        target_event_counts = [weight*min_ratio for weight in composition]
        self.target_event_counts = target_event_counts

        # collect the samplers to be put into a BlackBox
        samplers = []
        for filename, lumi, n_events in zip(filename_list, lumi_list,
                                            target_event_counts):
            samplers.append(Sampler(filename,
                                    # Sampler will multiply the lumi again
                                    n_events/lumi,
                                    lumi,
                                    # Needs to be set to avoid re-multiplying
                                    #   the preselection efficiency
                                    isSignal=True,
                                    holdout_frac=-1))
        self.samplers = samplers

    def get_samplers(self):
        return self.samplers


class BlackBox():
    def __init__(self, samplers, keys = [], nBatches = 1):
        self.data = dict()
        self.holdout = dict()
        self.samplers = samplers
        if(len(keys) == 0):
            #empty list means keep everything
            self.keys = (self.samplers[0]).keys
            # for key in SCALAR_KEYS:
            #     if key in self.keys:
            #         self.keys.remove(key)
        else:
            self.keys = keys

        self.nEvents = 0
        for sample in self.samplers:
            self.nEvents += (sample.nSample // nBatches) * nBatches #avoid rounding issues 

        self.nHoldOut = 0
        for sample in self.samplers:
            self.nHoldOut += sample.nHoldOut

        self.nBatches = nBatches

    def writeOut(self, filename, batch_start =0):

        print("Creating a blackbox with %i events in %i batches" % (self.nEvents, self.nBatches))
        evts_per_batch = self.nEvents // self.nBatches
        for i in range(batch_start, self.nBatches):
            print("Starting batch %i \n" % i)
            shuffle_order = np.arange(evts_per_batch)

            if self.nBatches == 1:
                f_batch = filename
                if not f_batch.endswith(".h5"):
                    f_batch =  filename + ".h5"
            else:
                f_batch =  filename + "_batch%i.h5" % i
            np.random.shuffle(shuffle_order)

            for sam in self.samplers: #setup the sampler for this batch
                sam.setup_batch(iBatch = i, nBatches = self.nBatches)

            #Fill the data from the various samplers
            for key in self.keys:
                print("Getting data for key %s " % key)
                
                for j,sam in enumerate(self.samplers):
                    if key in SCALAR_KEYS:
                        if j == 0:
                            self.data[key] = sam.nSample * sam.sample_scalar(key)
                        else:
                            self.data[key] += sam.nSample * sam.sample_scalar(key)
                    else:
                        if j == 0: 
                            self.data[key] = sam.sample(key)

                        else:
                            self.data[key] = np.append(self.data[key], sam.sample(key), axis = 0)
                if key in SCALAR_KEYS:
                    # weighted average and transforming to np array scalar
                    self.data[key] = np.array([self.data[key] / self.nEvents])
                else:
                    #Shuffle the order
                    self.data[key] = self.data[key][shuffle_order]

            print("Writing out to file %s" % f_batch)
            with h5py.File(f_batch, "w") as f:
                for key in self.keys:
                    shape = list(self.data[key].shape)
                    shape[0] = None
                    f.create_dataset(key, data = self.data[key], chunks = True, maxshape = shape,   compression = 'gzip')

    def writeHoldOut(self, filename):
        if(self.nHoldOut <= 0):
            print("No holdout events specified!")
            return
        holdout_order = np.arange(self.nHoldOut)
        np.random.shuffle(holdout_order)
        for key in self.keys:
            for i,sam in enumerate(self.samplers):
                if(i==0):
                    self.holdout[key] = sam.holdout(key)
                else:
                    self.holdout[key] = np.append(self.holdout[key], sam.holdout(key), axis = 0)

            #shuffle
            self.holdout[key] = self.holdout[key][holdout_order]

        print("Writing out file with holdouts to %s" % filename)
        with h5py.File(filename, "w") as f:
            for key in self.keys:
                shape = list(self.holdout[key].shape)
                shape[0] = None
                f.create_dataset(key, data = self.holdout[key], chunks = True, maxshape = shape,   compression = 'gzip')
