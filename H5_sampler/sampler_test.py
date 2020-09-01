from Sampler import * 

np.random.seed(123)

pbTofb = 1000.
bkg_holdout_frac = 0.05
sig_holdout_frac = 0.2
nBatches = 2
lumi = 1.
qcd1 = Sampler("../H5_maker/test_files/QCD_HT1000to1500_test.h5", pbTofb * 1088., lumi, holdout_frac = bkg_holdout_frac)
qcd2 = Sampler("../H5_maker/test_files/QCD_HT1500to2000_test.h5", pbTofb * 99.11, lumi, holdout_frac = bkg_holdout_frac)
qcd3 = Sampler("../H5_maker/test_files/QCD_HT2000toInf_test.h5", pbTofb * 20.23, lumi, holdout_frac = bkg_holdout_frac)
sig = Sampler("../H5_maker/test_files/WprimeToWZToWhadZhad_narrow_M-3500_TuneCP5_13TeV-madgraph_test.h5", 20000.,lumi, isSignal = True, holdout_frac = sig_holdout_frac)

ws = [qcd1, qcd2, qcd3, sig]
#ws = [qcd1]
keys= ['event_info', 'jet_kinematics', 'truth_label', 'jet1_extraInfo', 'jet1_PFCands']
#keys= []
BB = BlackBox(ws, keys, nBatches = nBatches)
#print(BB['truth_label'].shape)
#print(BB['jet_kinematics'].shape)

#BB.writeOut('BB_test_Wprime')
BB.writeHoldOut('BB_test_Wprime_holdout.h5')
exit(1)

f1 = h5py.File("BB_test_Wprime_batch0.h5", "r")
f2 = h5py.File("BB_test_Wprime_batch1.h5", "r")
f_hout = h5py.File("BB_test_Wprime_holdout.h5", "r")
f_orig = h5py.File("../H5_maker/test_files/QCD_HT1000to1500_test.h5", "r")

e1 = f1['event_info'][:,0]
e2 = f2['event_info'][:,0]
e_hout = f_hout['event_info'][:,0]
e_orig = f_orig['event_info'][:,0]

print(e_orig[99000],  e_orig[75000], e_orig[0])
print(np.where(e_hout == e_orig[99000]))
print(np.where(e_hout == e_orig[0]))
print(np.where(e1 == e_orig[99000]))
print(np.where(e1 == e_orig[75000]))
print(np.where(e1 == e_orig[0]))
print(np.where(e2 == e_orig[99000]))
print(np.where(e2 == e_orig[75000]))
print(np.where(e2 == e_orig[0]))
print(np.where(e2 == e2[0]))
print(np.where(e1 == e2[0]))
print(np.where(e1 == e1[0]))
print(np.where(e2 == e1[0]))

