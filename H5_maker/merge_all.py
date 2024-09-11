import os


categories = {

        "TT.h5" : ["TTToSemiLeptonic.h5", "TTToHadronic.h5", "TTTo2L2Nu.h5"],
        "diboson.h5" : ["WW.h5", "WZ.h5", "ZZ.h5"],
        "QCD_WJets.h5": ['QCD_MuEnriched_Pt120to170.h5', 
                         'QCD_MuEnriched_Pt170to300.h5', 
                         'QCD_MuEnriched_Pt300to470.h5', 
                         'QCD_MuEnriched_Pt470to600.h5', 
                         'QCD_MuEnriched_Pt600to800.h5', 
                          'QCD_MuEnriched_Pt800to1000.h5',
                          'QCD_MuEnriched_Pt1000toInf.h5',
                                'WJets_HT100to200.h5',
                                'WJets_HT200to400.h5',
                                'WJets_HT400to600.h5',
                                 'WJets_HT600to800.h5',
                                 'WJets_HT800to1200.h5', 
                                'WJets_HT1200to2500.h5',
                                 'WJets_HT2500toInf.h5' ],
        'SingleTop_merge.h5': ["SingleTop.h5", "SingleAntiTop.h5", "SingleTop_SChan.h5"],
        'TW.h5' : ['SingleTop_tW.h5', 'SingleTop_antitW.h5'],
        }

base_cmd = "python3 H5_merge.py "
base_dir = odir = "Lund_output_files17_may10/"



for key in categories.keys():
    cmd = base_cmd + odir + key + " " 
    for o in categories[key]:
        cmd += " " + base_dir + o + " " 
    print(cmd)
    os.system(cmd)

