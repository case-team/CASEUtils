from sampler import * 

np.random.seed(123)
w1 = Sampler("QCD_HT1000to1500_test.h5", 1088., 1.)
w2 = Sampler("QCD_HT1500to2000_test.h5", 99.11, 1.)
w3 = Sampler("QCD_HT2000toInf_test.h5", 20.23, 1.)

ws = [w1,w2,w3]
keys= ['event_info', 'jet_kinematics', 'truth_label', 'jet1_extraInfo', 'jet1_PFCands']
BB = BlackBox(ws, keys)
#print(BB['truth_label'].shape)
#print(BB['jet_kinematics'].shape)

BB.writeOut('BB_test.h5')
