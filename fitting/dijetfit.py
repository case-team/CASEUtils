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

    sb1_edge = 2200
    sb2_edge = 2800

    roobins = ROOT.RooBinning(len(binsx)-1, array('d',binsx), "mjjbins")
    bins_fine = int(binsx[-1]-binsx[0])/fine_bin_size
     

    histos_sb = ROOT.TH1F("mjj_sb","mjj_sb",bins_fine,binsx[0],binsx[-1])
    #useful to check if doing blinding correctly
    #load_h5_sb(options.inputFile, histos_sb, sb1_edge = sb1_edge + 40. , sb2_edge  = sb2_edge - 40.)
    load_h5_sb(options.inputFile, histos_sb)
    print "************ Found",histos_sb.GetEntries(),"total events"
    print


    print
    print
    print "############# FIT BACKGROUND AND SAVE PARAMETERS ###########"
    chi2s = [0]*len(nParsToTry)
    ndofs = [0]*len(nParsToTry)
    probs = [0]*len(nParsToTry)
    fitters_QCD = [None] *len(nParsToTry)
    qcd_fnames = [""]*len(nParsToTry)

    regions = [("SB1", binsx[0], sb1_edge),
               ("SB2", sb2_edge, binsx[-1]),
               ("FULL", binsx[0], binsx[-1])]

    for i,nPars in enumerate(nParsToTry):
        print("Trying %i parameter background fit" % nPars)
        qcd_fnames[i] = str(nPars) + 'par_qcd_fit.root'
        qcd_outfile = ROOT.TFile(qcd_fnames[i],'RECREATE')
         
        fitter_QCD=Fitter(['mjj_fine'])
        fitter_QCD.qcdShape('model_b','mjj_fine',nPars)
        fitter_QCD.importBinnedData(histos_sb,['mjj_fine'],'data_qcd', regions = regions)
        fres = fitter_QCD.fit('model_b','data_qcd',[ROOT.RooFit.Save(1), ROOT.RooFit.Verbose(0), ROOT.RooFit.Range("SB1,SB2")])
        #fres.Print()
        
        chi2_fine = fitter_QCD.projection("model_b","data_qcd","mjj_fine",plot_dir + str(nPars) + "par_qcd_fit.png",0,True)
        #chi2_binned = fitter_QCD.projection("model_b","data_qcd","mjj_fine",plot_dir + str(nPars) + "par_qcd_fit_binned.png",roobins,True)
        
        qcd_outfile.cd()
        histos_sb.Write()

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
        my_chi2,my_ndof = calculateChi2(histos_sb,nPars,hpull)
        my_prob = ROOT.TMath.Prob(my_chi2,my_ndof)
        PlotFitResults(frame,fres.GetName(),nPars,framePulls,"data_qcd","model_b",my_chi2,my_ndof, str(nPars) + "par_qcd_fit_binned_blinded" , plot_dir)
        
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

    nPars_sig = 4
    nPars_QCD = nParsToTry[best_i]
    fitter_QCD = fitters_QCD[best_i]
    qcd_fname = qcd_fnames[best_i]

    #print("\n Chose %i parameters based on F-test ! \n" % nPars )


    ##fit to total data
    #print "************ Found",histos_sb.GetEntries(),"total events"


    #sb_outfile = ROOT.TFile('sb_fit.root','RECREATE')
    #sb_outfile.cd()
    #histos_sb.Write("mjj_sb")
    #sb_outfile.Close()

    #sig_data_name = 'mjj_sb'
    #

    ##Fit to SB
    #sb_label = "sb_fit"

#    print
#    print
#    print "############ MAKE PER CATEGORY DATACARD AND WORKSPACE AND RUN COMBINE #############"
#    
#    card=DataCardMaker(sb_label)
#    
#    card.addSignalShape('model_signal_mjj','mjj','sig_fit.root',{'CMS_scale_j':1.0},{'CMS_res_j':1.0})
#    constant = xsec
#    #TODO Take sig pdf from SB fit
#    card.addFixedYieldFromFile('model_signal_mjj',0,'sig_fit.root',"mjj_sig",constant=constant)
#    card.addSystematic("CMS_scale_j","param",[0.0,0.012])
#    card.addSystematic("CMS_res_j","param",[0.0,0.08]) 
#    
#    card.addQCDShapeNoTag('model_qcd_mjj','mjj',qcd_fname,nPars_QCD)
#    card.addFloatingYield('model_qcd_mjj',1,qcd_fname, "mjj_sb")
#    for i in range(1,nPars_QCD+1): card.addSystematic("CMS_JJ_p%i"%i,"flatParam",[])
#    card.addSystematic("model_qcd_mjj_JJ_norm","flatParam",[])
#    
#    card.importBinnedData( "sb_fit.root", sig_data_name,["mjj"],'data_obs',1.0)
#    card.makeCard()
#    card.delete()
#    
#    cmd = 'text2workspace.py datacard_JJ_{l2}.txt -o workspace_JJ_{l1}_{l2}.root && combine -M Significance workspace_JJ_{l1}_{l2}.root -m {mass} -n significance_{l1}_{l2} && combine -M Significance workspace_JJ_{l1}_{l2}.root -m {mass} --pvalue -n pvalue_{l1}_{l2}'.format(mass=mass,l1=label,l2=sb_label)
#    print cmd
#    os.system(cmd)
#    checkSBFit('workspace_JJ_{l1}_{l2}.root'.format(l1=label,l2=sb_label),sb_label,roobins,label + "_" + sb_label, nPars_QCD + nPars_sig)
#
#    f_signif_name = 'higgsCombinesignificance_{l1}_{l2}.Significance.mH{mass:.0f}.root'.format(mass = mass, l1=label, l2 = sb_label)
#    f_pval_name = 'higgsCombinepvalue_{l1}_{l2}.Significance.mH{mass:.0f}.root'.format(mass = mass, l1=label, l2 = sb_label)
#
#    f_signif = ROOT.TFile.Open(f_signif_name, "READ")
#    lim = f_signif.Get("limit")
#    lim.GetEntry(0)
#    signifs.append(lim.limit)
#    f_signif.Close()
#
#    fitter.delete()
#
#    for fitt in fitters_QCD:
#        fitt.delete()
#
#
#    print "Toy significances were: " 
#    sig_sum = 0.
#    for s in signifs:
#        print s, " ", 
#        sig_sum += s
#
#    print "Avg is %.2f" % (sig_sum / options.nToys)
#
#    check_rough_sig(options.inputFile, options.mass * 0.9, options.mass * 1.1)
#    print("Done!")
#
