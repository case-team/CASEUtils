from Sampler import * 

np.random.seed(123)

pbTofb = 1000.
qcd1 = Sampler("../H5_maker/test_files/QCD_HT1000to1500_test.h5", pbTofb * 1088., 1.)
qcd2 = Sampler("../H5_maker/test_files/QCD_HT1500to2000_test.h5", pbTofb * 99.11, 1.)
qcd3 = Sampler("../H5_maker/test_files/QCD_HT2000toInf_test.h5", pbTofb * 20.23, 1.)
sig = Sampler("../H5_maker/test_files/WprimeToWZToWhadZhad_narrow_M-3500_TuneCP5_13TeV-madgraph_test.h5", 30000., 1.)

ws = [qcd1, qcd2, qcd3, sig]
#keys= ['event_info', 'jet_kinematics', 'truth_label', 'jet1_extraInfo', 'jet1_PFCands']
keys= []
BB = BlackBox(ws, keys)
#print(BB['truth_label'].shape)
#print(BB['jet_kinematics'].shape)

BB.writeOut('BB_test_Wprime.h5')
