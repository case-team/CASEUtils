from Sampler import * 
from optparse import OptionParser
import os

parser = OptionParser()
parser.add_option("-n", "--nBatch", type = int, default = 20, help="Number of Batches")
parser.add_option("-s", "--seed", type = int, default = 123, help="Seed")
options, args = parser.parse_args()

np.random.seed(options.seed)

pbTofb = 1000.
lumi = 9.
n_sig = 10000.
#bkg_holdout_frac = 0.01
#sig_holdout_frac = 0.5
bkg_holdout_frac = -1.
sig_holdout_frac =-1.
keys = ['event_info', 'jet1_PFCands', 'jet1_extraInfo', 'jet2_PFCands', 'jet2_extraInfo', 'jet_kinematics',  'truth_label']


base_dir = "/eos/user/t/tloesche/CASE_data/CASE_pancakes/"
out_dir = "/eos/cms/store/group/phys_b2g/CASE/h5_files/2017/BB_norepeats_v3_2500/"
qcd1_name = base_dir + "QCD/qcd_1000to1500_merge.h5"
qcd2_name = base_dir + "QCD/qcd_1500to2000_merge.h5"
qcd3_name = base_dir + "QCD/qcd_2000toInf_merge.h5"

sig1_name = base_dir +  "signal/grav_merge.h5"
sig2_name = base_dir +  "signal/wprime_0.h5"
sig3_name = base_dir +  "signal/wkk_0.h5"
sig4_name = base_dir +  "signal/bstar_merge.h5"


qcd1 = Sampler( qcd1_name, pbTofb * 1088., lumi, holdout_frac = bkg_holdout_frac)
qcd2 = Sampler( qcd2_name, pbTofb * 99.11, lumi, holdout_frac = bkg_holdout_frac)
qcd3 = Sampler( qcd3_name, pbTofb * 20.23, lumi, holdout_frac = bkg_holdout_frac)

sig1 = Sampler(sig1_name, n_sig, 1., isSignal = True, holdout_frac = sig_holdout_frac)
sig2 = Sampler(sig2_name, n_sig, 1., isSignal = True, holdout_frac = sig_holdout_frac)
sig3 = Sampler(sig3_name, n_sig, 1., isSignal = True, holdout_frac = sig_holdout_frac)
sig4 = Sampler(sig4_name, n_sig, 1., isSignal = True, holdout_frac = sig_holdout_frac)

ws = [qcd1, qcd2, qcd3, sig1, sig2, sig3, sig4]
BB = BlackBox(ws, keys, nBatches = options.nBatch)
os.system("mkdir %s" % out_dir)
f_out_name = out_dir + 'BB'
BB.writeOut(f_out_name)
#h_name = out_dir + "BB_testset"
#BB.writeHoldOut(h_name + ".h5")
