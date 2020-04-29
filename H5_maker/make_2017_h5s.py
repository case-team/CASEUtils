from H5_maker import *



wjets = ("WJetsToQQ_HT-800toInf", [
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017//hadd/WJetsToQQ_HT-800toInf_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_0.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017//hadd/WJetsToQQ_HT-800toInf_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_1.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017//hadd/WJetsToQQ_HT-800toInf_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_10.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017//hadd/WJetsToQQ_HT-800toInf_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_11.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017//hadd/WJetsToQQ_HT-800toInf_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_12.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017//hadd/WJetsToQQ_HT-800toInf_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_13.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017//hadd/WJetsToQQ_HT-800toInf_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_14.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017//hadd/WJetsToQQ_HT-800toInf_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_15.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017//hadd/WJetsToQQ_HT-800toInf_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_16.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017//hadd/WJetsToQQ_HT-800toInf_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_17.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017//hadd/WJetsToQQ_HT-800toInf_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_18.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017//hadd/WJetsToQQ_HT-800toInf_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_19.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017//hadd/WJetsToQQ_HT-800toInf_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_2.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017//hadd/WJetsToQQ_HT-800toInf_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_3.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017//hadd/WJetsToQQ_HT-800toInf_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_4.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017//hadd/WJetsToQQ_HT-800toInf_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_5.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017//hadd/WJetsToQQ_HT-800toInf_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_6.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017//hadd/WJetsToQQ_HT-800toInf_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_7.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017//hadd/WJetsToQQ_HT-800toInf_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_8.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017//hadd/WJetsToQQ_HT-800toInf_qc19_3j_TuneCP5_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_9.root"
            ])

zjets = ("ZJetsToQQ_HT-800toInf", [
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017//hadd/ZJetsToQQ_HT-800toInf_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_0.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017//hadd/ZJetsToQQ_HT-800toInf_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_1.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017//hadd/ZJetsToQQ_HT-800toInf_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_10.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017//hadd/ZJetsToQQ_HT-800toInf_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_11.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017//hadd/ZJetsToQQ_HT-800toInf_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_12.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017//hadd/ZJetsToQQ_HT-800toInf_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_13.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017//hadd/ZJetsToQQ_HT-800toInf_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_14.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017//hadd/ZJetsToQQ_HT-800toInf_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_15.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017//hadd/ZJetsToQQ_HT-800toInf_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_16.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017//hadd/ZJetsToQQ_HT-800toInf_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_17.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017//hadd/ZJetsToQQ_HT-800toInf_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_18.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017//hadd/ZJetsToQQ_HT-800toInf_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_19.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017//hadd/ZJetsToQQ_HT-800toInf_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_2.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017//hadd/ZJetsToQQ_HT-800toInf_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_3.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017//hadd/ZJetsToQQ_HT-800toInf_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_4.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017//hadd/ZJetsToQQ_HT-800toInf_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_5.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017//hadd/ZJetsToQQ_HT-800toInf_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_6.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017//hadd/ZJetsToQQ_HT-800toInf_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_7.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017//hadd/ZJetsToQQ_HT-800toInf_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_8.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017//hadd/ZJetsToQQ_HT-800toInf_qc19_4j_TuneCP5_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_9.root"
            ])

qcd1000 = ("QCD_HT1000to1500", [
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_0.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_1.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_10.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_11.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_12.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_13.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_14.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_15.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_16.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_17.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_18.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_19.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_2.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_3.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_4.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_5.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_6.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_7.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_8.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_9.root"
            ])

qcd1500 = ("QCD_HT1500to2000", [
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_0.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_1.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_10.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_11.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_12.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_13.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_14.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_15.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_16.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_17.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_18.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_19.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_2.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_3.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_4.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_5.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_6.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_7.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_8.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_9.root"
            ])


qcd2000 = ("QCD_HT2000toInf", [
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_0.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_1.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_10.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_11.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_12.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_13.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_14.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_15.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_16.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_17.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_18.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_19.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_2.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_3.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_4.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_5.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_6.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_7.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_8.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8-hadd/nano_mc_2017_9.root"
            ])

ttbar = ("TTToHadronic", [
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/TTToHadronic_TuneCP5_13TeV-powheg-pythia8-hadd/nano_mc_2017_0.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/TTToHadronic_TuneCP5_13TeV-powheg-pythia8-hadd/nano_mc_2017_1.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/TTToHadronic_TuneCP5_13TeV-powheg-pythia8-hadd/nano_mc_2017_10.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/TTToHadronic_TuneCP5_13TeV-powheg-pythia8-hadd/nano_mc_2017_11.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/TTToHadronic_TuneCP5_13TeV-powheg-pythia8-hadd/nano_mc_2017_12.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/TTToHadronic_TuneCP5_13TeV-powheg-pythia8-hadd/nano_mc_2017_13.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/TTToHadronic_TuneCP5_13TeV-powheg-pythia8-hadd/nano_mc_2017_14.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/TTToHadronic_TuneCP5_13TeV-powheg-pythia8-hadd/nano_mc_2017_15.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/TTToHadronic_TuneCP5_13TeV-powheg-pythia8-hadd/nano_mc_2017_16.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/TTToHadronic_TuneCP5_13TeV-powheg-pythia8-hadd/nano_mc_2017_17.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/TTToHadronic_TuneCP5_13TeV-powheg-pythia8-hadd/nano_mc_2017_18.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/TTToHadronic_TuneCP5_13TeV-powheg-pythia8-hadd/nano_mc_2017_19.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/TTToHadronic_TuneCP5_13TeV-powheg-pythia8-hadd/nano_mc_2017_2.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/TTToHadronic_TuneCP5_13TeV-powheg-pythia8-hadd/nano_mc_2017_3.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/TTToHadronic_TuneCP5_13TeV-powheg-pythia8-hadd/nano_mc_2017_4.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/TTToHadronic_TuneCP5_13TeV-powheg-pythia8-hadd/nano_mc_2017_5.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/TTToHadronic_TuneCP5_13TeV-powheg-pythia8-hadd/nano_mc_2017_6.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/TTToHadronic_TuneCP5_13TeV-powheg-pythia8-hadd/nano_mc_2017_7.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/TTToHadronic_TuneCP5_13TeV-powheg-pythia8-hadd/nano_mc_2017_8.root",
                "root://cmseos.fnal.gov//store/user/lpcbacon/pancakes/02/2017/UL//hadd/TTToHadronic_TuneCP5_13TeV-powheg-pythia8-hadd/nano_mc_2017_9.root"
            ])

year = 2017
outdir = ""

NanoReader(0, inputFileNames = qcd1000[1], outputFileName = outdir + qcd1000[0] + ".h5", year = year)
