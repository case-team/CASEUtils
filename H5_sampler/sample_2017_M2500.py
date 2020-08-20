from Sampler import * 
from optparse import OptionParser
import os

parser = OptionParser()
parser.add_option("-n", "--nBatch", type = int, default = 40, help="Number of Batches")
parser.add_option("-s", "--seed", type = int, default = 123, help="Seed")
options, args = parser.parse_args()

np.random.seed(options.seed)

pbTofb = 1000.
lumi = 41.5
n_sig = 10000.
bkg_holdout_frac = 0.01
sig_holdout_frac = 0.5
keys = ['event_info', 'jet1_PFCands', 'jet1_extraInfo', 'jet2_PFCands', 'jet2_extraInfo', 'jet_kinematics',  'truth_label']


base_dir = "/eos/cms/store/group/phys_b2g/CASE/h5_files/2017/"
out_dir = "/eos/cms/store/group/phys_b2g/CASE/h5_files/2017/BB_v2_2500/"
qcd1_name = base_dir + "QCD_HT1000to15000_merge.h5"
qcd2_name = base_dir + "QCD_HT1500to2000_merge.h5"
qcd3_name = base_dir + "QCD_HT2000toInf_merge.h5"

sig1_name = base_dir +  "BulkGravToZZToZhadZhad_narrow_M-2500.h5"
sig2_name = base_dir +  "WprimeToWZToWhadZhad_narrow_M-2500.h5"
sig3_name = base_dir +  "WkkToWRadionToWWW_M2500-R0-08.h5"
sig4_name = base_dir +  "BstarToTW_M-2600_LH.h5"


qcd1 = Sampler( qcd1_name, pbTofb * 1088., lumi, holdout_frac = bkg_holdout_frac)
qcd2 = Sampler( qcd2_name, pbTofb * 99.11, lumi, holdout_frac = bkg_holdout_frac)
qcd3 = Sampler( qcd3_name, pbTofb * 20.23, lumi, holdout_frac = bkg_holdout_frac)

sig1 = Sampler(sig1_name, n_sig, 1., isSignal = True, holdout_frac = sig_holdout_frac)
sig2 = Sampler(sig2_name, n_sig, 1., isSignal = True, holdout_frac = sig_holdout_frac)
sig3 = Sampler(sig3_name, n_sig, 1., isSignal = True, holdout_frac = sig_holdout_frac)
sig4 = Sampler(sig4_name, n_sig, 1., isSignal = True, holdout_frac = sig_holdout_frac)

ws = [qcd1, qcd2, qcd3, sig1, sig2, sig3]
BB = BlackBox(ws, keys, nBatches = options.nBatch)
os.system("mkdir %s" % out_dir)
f_out_name = out_dir + 'BB'
BB.writeOut(f_out_name)
h_name = out_dir + "BB_testset"
BB.writeHoldOut(h_name + ".h5")
#for i in range(options.nBatch):
#    os.system("python ../jet_images/make_jet_images.py -i %s_batch%i.h5 -o %s_batch%i_images.h5 " % (f_out_name, i,  f_out_name, i))
#
#os.system("python ../jet_images/make_jet_images.py -i %s.h5 -o %s_images.h5 " % (h_name, h_name))
