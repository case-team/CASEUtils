import h5py, math, commands, random
from array import array
import numpy as np
import time, sys, os, optparse, json

import ROOT
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
import CMS_lumi, tdrstyle
tdrstyle.setTDRStyle()
ROOT.gROOT.SetBatch(True)
ROOT.RooRandom.randomGenerator().SetSeed(random.randint(0, 1e+6))

from Fitter import Fitter
from DataCardMaker import DataCardMaker
from Utils import *
    
        
if __name__ == "__main__":
    parser = optparse.OptionParser()
    # 0.0003
    parser.add_option("--xsec","--xsec",dest="xsec",type=float,default=1,help="Injected signal cross section (suggested range 0-0.03)")
    parser.add_option("-M","-M",dest="mass",type=float,default=2500.,help="Injected signal mass")
    parser.add_option("-i","--inputFile",dest="inputFile",default='fit_inputs/no_selection_03p.h5',help="input h5 file")
    parser.add_option("-p","--plotDir",dest="plotDir",default='plots/',help="Where to put the plots")
    parser.add_option("--nToys", type=int, default =5, help='how many toys to generate')
    (options,args) = parser.parse_args()

    label = 'test'

    plot_dir = options.plotDir
    os.system("mkdir %s" % plot_dir)

    fine_bin_size = 4
     
    xsec = options.xsec
    mass = options.mass 
    nParsToTry = [2,3,4]

    binsx = [1530,1607,1687,1770,1856,1945,2037,2132,2231,2332,2438,2546,2659,2775,2895,3019,3147,3279,3416,3558,3704,3854,4010,4171,4337,4509,4686,4869,5058,5253,5500,5663,5877,6099,6328,6564,6808]
    #round to smallest precision we are storing mass values with, otherwise get weird effects related to bin size
    roundTo(binsx, fine_bin_size)

    roobins = ROOT.RooBinning(len(binsx)-1, array('d',binsx), "mjjbins")
    bins_fine = int(binsx[-1]-binsx[0])/fine_bin_size
    bins_sig_fit = array('f',truncate([binsx[0]+ib*fine_bin_size for ib in range(bins_fine+1)],0.8*mass,1.2*mass))
    large_bins_sig_fit = array('f',truncate(binsx,0.8*mass,1.2*mass))
    roobins_sig_fit = ROOT.RooBinning(len(large_bins_sig_fit)-1, array('d',large_bins_sig_fit), "mjjbins_sig")



     

     
    #Fit to signal events only (Can we make this optional?)

    #Signal data preparation 
    histos_sig = ROOT.TH1F("mjj_sig","mjj_sig",len(bins_sig_fit)-1,bins_sig_fit)
    load_h5_sig(options.inputFile,histos_sig, mass)
    print "************ Found",histos_sig.GetEntries(),"signal events"
    print


    print "########## FIT SIGNAL AND SAVE PARAMETERS ############"
    sig_outfile = ROOT.TFile("sig_fit.root","RECREATE")
    
    fitter=Fitter(['mjj_fine'])
    fitter.signalResonance('model_s',"mjj_fine",mass)
    fitter.w.var("MH").setVal(mass)
    fitter.importBinnedData(histos_sig,['mjj_fine'],'data')
    fres = fitter.fit('model_s','data',[ROOT.RooFit.Save(1)])
    #fres.Print()
    
    mjj_fine = fitter.getVar('mjj_fine')
    mjj_fine.setBins(len(bins_sig_fit))
    chi2_fine = fitter.projection("model_s","data","mjj_fine",plot_dir + "signal_fit.png")
    fitter.projection("model_s","data","mjj_fine",plot_dir + "signal_fit_log.png",0,True)
    chi2 = fitter.projection("model_s","data","mjj_fine",plot_dir + "signal_fit_binned.png",roobins_sig_fit)
    fitter.projection("model_s","data","mjj_fine",plot_dir + "signal_fit_log_binned.png",roobins_sig_fit,True)
    
    sig_outfile.cd()
    histos_sig.Write()

    graphs={'mean':ROOT.TGraphErrors(),'sigma':ROOT.TGraphErrors(),'alpha':ROOT.TGraphErrors(),'sign':ROOT.TGraphErrors(),'scalesigma':ROOT.TGraphErrors(),'sigfrac':ROOT.TGraphErrors()}
    for var,graph in graphs.iteritems():
        value,error=fitter.fetch(var)
        graph.SetPoint(0,mass,value)
        graph.SetPointError(0,0.0,error)

    sig_outfile.cd()
    for name,graph in graphs.iteritems(): graph.Write(name)
     
    sig_outfile.Close() 
    
    print "#############################"
    print "signal fit chi2 (fine binning)",chi2_fine
    print "signal fit chi2 (large binning)",chi2
    print "#############################"

    #Fit to background events only (Can we make this optional?)

    #Background data preparation
    histos_qcd = ROOT.TH1F("mjj_qcd","mjj_qcd",bins_fine,binsx[0],binsx[-1])
    load_h5_bkg(options.inputFile, histos_qcd)
    print "************ Found",histos_qcd.GetEntries(),"background events"
    print


    print
    print
    print "############# FIT BACKGROUND AND SAVE PARAMETERS ###########"
    chi2s = [0]*len(nParsToTry)
    ndofs = [0]*len(nParsToTry)
    probs = [0]*len(nParsToTry)
    fitters_QCD = [None] *len(nParsToTry)
    qcd_fnames = [""]*len(nParsToTry)

    for i,nPars in enumerate(nParsToTry):
        print("Trying %i parameter background fit" % nPars)
        qcd_fnames[i] = str(nPars) + 'par_qcd_fit.root'
        qcd_outfile = ROOT.TFile(qcd_fnames[i],'RECREATE')
         
        fitter_QCD=Fitter(['mjj_fine'])
        fitter_QCD.qcdShape('model_b','mjj_fine',nPars)
        fitter_QCD.importBinnedData(histos_qcd,['mjj_fine'],'data_qcd')
        fres = fitter_QCD.fit('model_b','data_qcd',[ROOT.RooFit.Save(1), ROOT.RooFit.Verbose(0)])
        #fres.Print()
        
        chi2_fine = fitter_QCD.projection("model_b","data_qcd","mjj_fine",plot_dir + str(nPars) + "par_qcd_fit.png",0,True)
        #chi2_binned = fitter_QCD.projection("model_b","data_qcd","mjj_fine",plot_dir + str(nPars) + "par_qcd_fit_binned.png",roobins,True)
        
        qcd_outfile.cd()
        histos_qcd.Write()

        mjj = fitter_QCD.getVar('mjj_fine')
        mjj.setBins(bins_fine)
        model = fitter_QCD.getFunc('model_b')
        dataset = fitter_QCD.getData('data_qcd')
        
        frame = mjj.frame()
        dataset.plotOn(frame,ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson), ROOT.RooFit.Name("data_qcd"),ROOT.RooFit.Invisible(),ROOT.RooFit.Binning(roobins))
        model.plotOn(frame,ROOT.RooFit.VisualizeError(fres,1),ROOT.RooFit.FillColor(ROOT.kRed-7),ROOT.RooFit.LineColor(ROOT.kRed-7),ROOT.RooFit.Name(fres.GetName()))
        model.plotOn(frame,ROOT.RooFit.LineColor(ROOT.kRed+1),ROOT.RooFit.Name("model_b"))
        
        framePulls = mjj.frame()
        #average bin edges instead of bin center
        useBinAverage = True
        hpull = frame.pullHist("data_qcd","model_b",useBinAverage)
        framePulls.addPlotable(hpull,"X0 P E1")
        
        dataset.plotOn(frame,ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson),ROOT.RooFit.Name("data_qcd"),ROOT.RooFit.XErrorSize(0),ROOT.RooFit.Binning(roobins))
        my_chi2,my_ndof = calculateChi2(histos_qcd,nPars,hpull)
        my_prob = ROOT.TMath.Prob(my_chi2,my_ndof)
        PlotFitResults(frame,fres.GetName(),nPars,framePulls,"data_qcd","model_b",my_chi2,my_ndof, str(nPars) + "par_qcd_fit_binned" , plot_dir)
        
        graphs = {}
        for p in range(nPars): graphs['p%i'%(p+1)] = ROOT.TGraphErrors()
        for var,graph in graphs.iteritems():
            print var
            value,error=fitter_QCD.fetch(var)
            graph.SetPoint(0,mass,value)
            graph.SetPointError(0,0.0,error)
             
        qcd_outfile.cd()              
        for name,graph in graphs.iteritems(): graph.Write(name)
            
        qcd_outfile.Close()

        print "#############################"
        print "% i Parameter results: " % nPars 
        print "bkg fit chi2/nbins (fine binning)",chi2_fine
        print "My chi2, ndof, prob", my_chi2, my_ndof, my_prob
        print "My chi/ndof, chi2/nbins", my_chi2/my_ndof, my_chi2/(my_ndof + nPars)
        print "#############################"

        chi2s[i] = my_chi2
        ndofs[i] = my_ndof
        probs[i] = my_prob
        fitters_QCD[i] = fitter_QCD

    best_i = f_test(nParsToTry, ndofs, chi2s)

    nPars = nParsToTry[best_i]
    fitter_QCD = fitters_QCD[best_i]
    qcd_fname = qcd_fnames[best_i]

    print("\n Chose %i parameters based on F-test ! \n" % nPars )


    #Fit to total data
    #if(options.no_toys):
    #    histos_sb = ROOT.TH1F("mjj_sb","mjj_sb",bins_fine,binsx[0],binsx[-1])
    #    load_h5_sb(options.inputFile, histos_sb)
    #    print "************ Found",histos_sb.GetEntries(),"total events"


    #    sb_outfile = ROOT.TFile('sb_fit.root','RECREATE')
    #    sb_outfile.cd()
    #    histos_sb.Write("mjj_sb")
    #    sb_outfile.Close()

    #    sig_data_name = 'mjj_sb'
    #else:

    print
    print 
    print "############# GENERATE SIGNAL+BACKGROUND DATA FROM PDFs ###########"

    signifs = []
    for i in range(options.nToys):
        toy_label = "toy%i" % i

        f_name = "cache%i.root" % random.randint(0, 1e+6)
        f = ROOT.TFile(f_name,"RECREATE")
        f.cd()
        w=ROOT.RooWorkspace("w","w")

        model_b = fitter_QCD.getFunc('model_b')
        model_s = fitter.getFunc('model_s')

        #model_b.Print("v")
        #model_s.Print("v")

        print "Generate",histos_qcd.Integral(),"background events from model_b"
        dataqcd = model_b.generateBinned(ROOT.RooArgSet(mjj),histos_qcd.Integral())
        hdataqcd = dataqcd.createHistogram("mjj_fine")
        hdataqcd.SetName("mjj_generate_qcd")

        if xsec != 0:
            print "Generate",int(histos_sig.Integral()*xsec),"signal events from model_s"
            datasig = model_s.generateBinned(ROOT.RooArgSet(mjj),int(histos_sig.Integral()*xsec))
            hdatasig = datasig.createHistogram("mjj_fine")
        else:
            hdatasig = ROOT.TH1F("mjj_generate_sig","mjj_generate_sig",histos_qcd.GetNbinsX(),histos_qcd.GetXaxis().GetXmin(),histos_qcd.GetXaxis().GetXmax())

        hdatasig.SetName("mjj_generate_sig")


        f_toy_name = 'sb_toy%i_fit.root' %i
        sb_outfile = ROOT.TFile(f_toy_name ,'RECREATE')
        sb_outfile.cd()
        htot = ROOT.TH1F()
        htot = hdataqcd.Clone("mjj_generate_tot")
        htot.Add(hdatasig)
        hdatasig.Write("mjj_generate_sig")
        hdataqcd.Write("mjj_generate_qcd")
        htot.Write("mjj_generate_tot")

        w.Delete()
        f.Close()
        f.Delete()
        os.system("rm %s" % f_name)
        sb_outfile.Close()
        sig_data_name = 'mjj_generate_tot'



        print
        print
        print "############ MAKE PER CATEGORY DATACARD AND WORKSPACE AND RUN COMBINE #############"
        
        card=DataCardMaker(toy_label)
        
        card.addSignalShape('model_signal_mjj','mjj','sig_fit.root',{'CMS_scale_j':1.0},{'CMS_res_j':1.0})
        constant = xsec
        card.addFixedYieldFromFile('model_signal_mjj',0,'sig_fit.root',histos_sig.GetName(),constant=constant)
        card.addSystematic("CMS_scale_j","param",[0.0,0.012])
        card.addSystematic("CMS_res_j","param",[0.0,0.08]) 
        
        card.addQCDShapeNoTag('model_qcd_mjj','mjj',qcd_fname,nPars)
        card.addFloatingYield('model_qcd_mjj',1,qcd_fname, histos_qcd.GetName())
        for i in range(1,nPars+1): card.addSystematic("CMS_JJ_p%i"%i,"flatParam",[])
        card.addSystematic("model_qcd_mjj_JJ_norm","flatParam",[])
        
        card.importBinnedData( f_toy_name, sig_data_name,["mjj"],'data_obs',1.0)
        card.makeCard()
        card.delete()
        
        cmd = 'text2workspace.py datacard_JJ_{l2}.txt -o workspace_JJ_{l1}_{l2}.root && combine -M Significance workspace_JJ_{l1}_{l2}.root -m {mass} -n significance_{l1}_{l2} && combine -M Significance workspace_JJ_{l1}_{l2}.root -m {mass} --pvalue -n pvalue_{l1}_{l2}'.format(mass=mass,l1=label,l2=toy_label)
        print cmd
        os.system(cmd)

        #run and visualize s+b fit as sanity check
        checkSBFit('workspace_JJ_{l1}_{l2}.root'.format(l1=label,l2=toy_label),toy_label,roobins,label + "_" + toy_label)
        f_signif_name = 'higgsCombinesignificance_{l1}_{l2}.Significance.mH{mass:.0f}.root'.format(mass = mass, l1=label, l2 = toy_label)
        f_pval_name = 'higgsCombinepvalue_{l1}_{l2}.Significance.mH{mass:.0f}.root'.format(mass = mass, l1=label, l2 = toy_label)

        f_signif = ROOT.TFile.Open(f_signif_name, "READ")
        lim = f_signif.Get("limit")
        lim.GetEntry(0)
        signifs.append(lim.limit)
        f_signif.Close()

    fitter.delete()

    for fitt in fitters_QCD:
        fitt.delete()


    print "Toy significances were: " 
    sig_sum = 0.
    for s in signifs:
        print s, " ", 
        sig_sum += s

    print "Avg is %.2f" % (sig_sum / options.nToys)

    check_rough_sig(options.inputFile, options.mass * 0.9, options.mass * 1.1)
    print("Done!")

