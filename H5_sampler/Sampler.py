import h5py
import numpy as np 

class Sampler():
    def __init__(self, filename, xsec, lumi, isSignal = False):
        self.filename = filename
        with h5py.File(self.filename, "r") as h5_file:
            self.eff_xsec = xsec * h5_file['preselection_eff'][0]
            self.keys = h5_file.keys()
            self.nEvents = h5_file[self.keys[0]].shape[0]
        self.nSample = int(lumi * self.eff_xsec)
        if(isSignal): self.nSample = int(lumi * xsec)
        print("Will sample %i events from file %s" %(self.nSample, self.filename))
        if(self.nSample > self.nEvents):
            print("Warning not enough unique events in file %s to match cross section! Need to sample %i events but only have %i saved \n" 
                    %(self.filename, self.nSample, self.nEvents))
        #sample with replacement
        self.sample_idxs = np.random.choice(self.nEvents, self.nSample)

    def sample(self, key):
        if(key not in self.keys):
            print("Error! Key %s not in list of keys for file %s. Keys are: \n" % (key, self.filename))
            print(self.keys)
            exit(1)
        if(self.nSample < 1): return []
        with h5py.File(self.filename, "r") as h5_file:
            out = h5_file[key][()][self.sample_idxs]
            return out

class BlackBox():
    def __init__(self, samplers, keys = []):
        self.data = dict()
        if(len(keys) == 0):
            #empty list means keep everything
            keys = (samplers[0]).keys
            keys.remove('preselection_eff')
        self.keys = keys

        self.nEvents = 0
        for sample in samplers:
            self.nEvents += sample.nSample

        print("Creating a blackbox with %i events" % self.nEvents)
        shuffle_order = np.arange(self.nEvents)
        np.random.shuffle(shuffle_order)
        #Fill the data from the various samplers
        for key in keys:
            print("Getting data for key %s " % key)
            for i,sam in enumerate(samplers):
                if(i==0): self.data[key] = sam.sample(key)
                else: self.data[key] = np.append(self.data[key], sam.sample(key), axis = 0)

        #Shuffle the order
            self.data[key] = self.data[key][shuffle_order]

    def __getitem__(self, key):
       if(key not in self.keys):
           print("Error! Key %s not in list of keys for this dataset. Keys are: \n" % (key))
           print(self.keys)
           exit(1)
       return self.data[key]

    def writeOut(self, filename):
        print("Writing out to file %s" % filename)
        with h5py.File(filename, "w") as f:
            for key in self.keys:
                shape = list(self.data[key].shape)
                shape[0] = None
                f.create_dataset(key, data = self.data[key], chunks = True, maxshape = shape,   compression = 'gzip')




