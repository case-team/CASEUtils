from Sampler import * 
from optparse import OptionParser
import os

parser = OptionParser()
parser.add_option("-i", "--idx", type = int, default = 0, help="Batch index")
parser.add_option("-n", "--nBatch", type = int, default = 40, help="Number of Batches")
options, args = parser.parse_args()

np.random.seed(options.idx)

pbTofb = 1000.
lumi = 41.5
n_sig = 20000.
n_batches = options.nBatch
keys= [] #everything
holdout_frac = 0.4

base_dir = "/eos/cms/store/group/phys_b2g/CASE/h5_files/2017/"
qcd1_name = base_dir + "QCD_HT1000to15000_merge.h5"
qcd2_name = base_dir + "QCD_HT1500to2000_merge.h5"
qcd3_name = base_dir + "QCD_HT2000toInf_merge.h5"

sig1_name = base_dir +  "BulkGravToZZToZhadZhad_narrow_M-2500.h5"
sig2_name = base_dir +  "WprimeToWZToWhadZhad_narrow_M-2500.h5"
sig3_name = base_dir +  "WkkToWRadionToWWW_M2500-R0-08.h5"


print("starting batch %i" % i)
qcd1 = Sampler( qcd1_name, pbTofb * 1088., lumi/n_batches, holdout_frac = holdout_frac)
qcd2 = Sampler( qcd2_name, pbTofb * 99.11, lumi/n_batches, holdout_frac = holdout_frac)
qcd3 = Sampler( qcd3_name, pbTofb * 20.23, lumi/n_batches, holdout_frac = holdout_frac)

sig1 = Sampler(sig1_name, n_sig/n_batches, 1., isSignal = True, holdout_frac = holdout_frac)
sig2 = Sampler(sig2_name, n_sig/n_batches, 1., isSignal = True, holdout_frac = holdout_frac)
sig3 = Sampler(sig3_name, n_sig/n_batches, 1., isSignal = True, holdout_frac = holdout_frac)

ws = [qcd1, qcd2, qcd3, sig1, sig2, sig3]
BB = BlackBox(ws, keys)
BB.writeOut(base_dir + 'BB_batch%i.h5' % i)
if(options.idx == 0):
    BB.writeHoldOut("BB_testset.h5")
