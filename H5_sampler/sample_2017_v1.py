from Sampler import * 
import os

np.random.seed(123)

pbTofb = 1000.
lumi = 41.5
n_sig = 100000.
n_batches = 40
keys= [] #everything

base_dir = "/eos/cms/store/group/phys_b2g/CASE/h5_files/2017/"
qcd1_name = base_dir + "QCD_HT1000to15000_merge.h5"
qcd2_name = base_dir + "QCD_HT1500to2000_merge.h5"
qcd3_name = base_dir + "QCD_HT2000toInf_merge.h5"

sig1_name = base_dir +  "BulkGravToZZToZhadZhad_narrow_M-2500.h5"
sig2_name = base_dir +  "WprimeToWZToWhadZhad_narrow_M-3500.h5"
sig3_name = base_dir +  "WkkToWRadionToWWW_M2500-R0-08.h5"


for i in range(n_batches):
    print("starting batch %i" % i)
    qcd1 = Sampler( qcd1_name, pbTofb * 1088., lumi/n_batches)
    qcd2 = Sampler( qcd2_name, pbTofb * 99.11, lumi/n_batches)
    qcd3 = Sampler( qcd3_name, pbTofb * 20.23, lumi/n_batches)

    sig1 = Sampler(sig1_name, n_sig/n_batches, 1., isSignal = True)
    sig2 = Sampler(sig2_name, n_sig/n_batches, 1., isSignal = True)
    sig3 = Sampler(sig3_name, n_sig/n_batches, 1., isSignal = True)

    ws = [qcd1, qcd2, qcd3, sig1, sig2, sig3]
    BB = BlackBox(ws, keys)
    BB.writeOut(base_dir + 'BB_batch%i.h5' % i)
