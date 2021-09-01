from Fitter import Fitter
from DataCardMaker_counting import DataCardMaker
from Utils import *
import ROOT
import h5py
import numpy as np
import numdifftools        

## to numerically estimate the Jacobian of the fit functions, we'll need to reinstantiate it many times
## so here we put that QCD shape creation into a function that only depends on the parameters
def QCD_model(par_list, name="model"):

    assert len(par_list)>1, "At least two parameters need to be given!"

    ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
    poi = "mjj_fine"
    w = ROOT.RooWorkspace("w","w")
    w.factory(poi+"[1530, 6808]")

    if len(par_list) == 2:
        w.factory("p1["+str(par_list[0])+"]")
        w.factory("p2["+str(par_list[1])+"]")
        model = ROOT.RooGenericPdf(name, "pow(1-@0/13000., @1)/pow(@0/13000., @2)", ROOT.RooArgList(w.var(poi), w.var("p1"), w.var("p2")))
    elif len(par_list) == 3:
        w.factory("p1["+str(par_list[0])+"]")
        w.factory("p2["+str(par_list[1])+"]")
        w.factory("p3["+str(par_list[2])+"]")
        model = ROOT.RooGenericPdf(name, "pow(1-@0/13000., @1)/pow(@0/13000., @2+@3*log(@0/13000.))", ROOT.RooArgList(w.var(poi), w.var("p1"), w.var("p2"), w.var("p3")))
    elif len(par_list) == 4:
        w.factory("p1["+str(par_list[0])+"]")
        w.factory("p2["+str(par_list[1])+"]")
        w.factory("p3["+str(par_list[2])+"]")
        w.factory("p4["+str(par_list[3])+"]")
        model = ROOT.RooGenericPdf(name, "pow(1-@0/13000., @1)/ ( pow(@0/13000., @2+@3*log(@0/13000.)+@4*pow(log(@0/13000.),2)) )", ROOT.RooArgList(w.var(poi), w.var("p1"), w.var("p2"), w.var("p3"), w.var("p4")))
    else:
        raise NotImplementedError("Only implemented parameters up to 4!")
    return model, w


def dijetfit(options):

    sig_norm = 1.0 ##TODO adjust via command line input

    significances = {}
    p_values = {}
    limits_obs = {}
    limits_exp_median = {}
    limits_exp_plus1sigma = {}
    limits_exp_plus2sigma = {}
    limits_exp_minus1sigma = {}
    limits_exp_minus2sigma = {}

    plot_dir = options.plotDir
    if(plot_dir[-1] != '/'): plot_dir+= '/'
    os.system("mkdir %s" % plot_dir)
    os.system("mkdir %s" % plot_dir)

    fine_bin_size = 4
     
    mass = options.mass 
    nParsToTry = [2,3,4]

    binsx = [1530,1607,1687,1770,1856,1945,2037,2132,2231,2332,2438,2546,2659,2775,2895,3019,3147,3279,3416,3558,3704,3854,4010,4171,4337,4509,4686,4869,5058,5253,5500,5663,5877,6099,6328,6564,6808]

    ## finding closest bins to a SR of +/-10% mjj
    SR_low = mass-0.1*mass
    SR_high = mass+0.1*mass
    sb1_edge = binsx[np.argmin(np.abs(np.array(binsx)-SR_low))]+1
    sb2_edge = binsx[np.argmin(np.abs(np.array(binsx)-SR_high))]+1

    print "sb1_edge =", sb1_edge
    print "sb2_edge =", sb2_edge

    #round to smallest precision we are storing mass values with, otherwise get weird effects related to bin size
    roundTo(binsx, fine_bin_size)

    roobins = ROOT.RooBinning(len(binsx)-1, array('d',binsx), "mjjbins")
    bins_fine = int(binsx[-1]-binsx[0])/fine_bin_size
     
    histos_sb = ROOT.TH1F("mjj_sb","mjj_sb",bins_fine,binsx[0],binsx[-1])

    #useful to check if doing blinding correctly
    #load_h5_sb(options.inputFile, histos_sb, sb1_edge = sb1_edge + 40. , sb2_edge  = sb2_edge - 40.)
    load_h5_sb(options.inputFile, histos_sb)

    ## derive the number of signal events in the SR
    current_file = h5py.File(options.inputFile, "r")
    m_jj = current_file["mjj"][()].flatten()
    targets = current_file["truth_label"][()].flatten()
    mjj_SR_mask = (m_jj >= sb1_edge) & (m_jj <= sb2_edge)
    SR_targets = targets[mjj_SR_mask]
    n_SR_sig_truth = sum(SR_targets==1)
    n_SR_bkg_truth = sum(SR_targets==0)

    print
    print "############# FIT BACKGROUND AND SAVE PARAMETERS ###########"

    chi2s = [0]*len(nParsToTry)
    ndofs = [0]*len(nParsToTry)
    probs = [0]*len(nParsToTry)
    qcd_fnames = [""]*len(nParsToTry)

    regions = [("SB1", binsx[0], sb1_edge),
               ("SB2", sb2_edge, binsx[-1]),
               ("SR", sb1_edge, sb2_edge),
               ("FULL", binsx[0], binsx[-1])]

    blind_range = ROOT.RooFit.Range("SB1,SB2")
    full_range = ROOT.RooFit.Range("FULL")
    fit_ranges = [(binsx[0], sb1_edge), (sb2_edge, binsx[-1])]

    histos_sb_blind = apply_blinding(histos_sb, ranges = fit_ranges)
    num_blind = histos_sb_blind.Integral()
    norm = ROOT.RooFit.Normalization(num_blind, ROOT.RooAbsReal.NumEvent)

    for i,nPars in enumerate(nParsToTry):
        print("Trying %i parameter background fit" % nPars)
        qcd_fnames[i] = str(nPars) + 'par_qcd_fit%i.root' %i
        qcd_outfile = ROOT.TFile(qcd_fnames[i],'RECREATE')

        model_name = "model_b" + str(i)
    
        fitter_QCD=Fitter(['mjj_fine'])
        fitter_QCD.qcdShape(model_name,'mjj_fine',nPars)
        fitter_QCD.importBinnedData(histos_sb_blind,['mjj_fine'],'data_qcd_blind', regions = regions)

        #For some reason (???) the first fit gives some weird results, but if you fit it twice it works ok??? This should be fixed...
        fres = fitter_QCD.fit(model_name,'data_qcd_blind', options = [ROOT.RooFit.Save(1), ROOT.RooFit.Verbose(0), blind_range])
        #fres = fitter_QCD.fit(model_name,'data_qcd_blind', options = [ROOT.RooFit.Save(1), ROOT.RooFit.Verbose(0), blind_range])
        
        chi2_fine = fitter_QCD.projection(model_name,"data_qcd_blind","mjj_fine",plot_dir + str(nPars) + "par_qcd_fit.png",0,True)
        #chi2_binned = fitter_QCD.projection(model_name,"data_qcd","mjj_fine",plot_dir + str(nPars) + "par_qcd_fit_binned.png",roobins,True)
        
        qcd_outfile.cd()
        histos_sb_blind.Write()

        mjj = fitter_QCD.getVar('mjj_fine')
        mjj.setBins(bins_fine)
        model = fitter_QCD.getFunc(model_name)
        dataset = fitter_QCD.getData('data_qcd_blind')
        
        frame = mjj.frame()
        dataset.plotOn(frame,ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson), ROOT.RooFit.Name("data_qcd_blind"),ROOT.RooFit.Invisible(),ROOT.RooFit.Binning(roobins))
        model.plotOn(frame,ROOT.RooFit.VisualizeError(fres,1),ROOT.RooFit.FillColor(ROOT.kRed-7),ROOT.RooFit.LineColor(ROOT.kRed-7),ROOT.RooFit.Name(fres.GetName()), 
                ROOT.RooFit.Range("SB1,SB2"), norm)
        model.plotOn(frame,ROOT.RooFit.LineColor(ROOT.kRed+1),ROOT.RooFit.Name(model_name), 
                ROOT.RooFit.Range("SB1,SB2"), norm)

        #Compute integrals over signal region
        mjj_set = ROOT.RooArgSet(mjj)
        sig_region_integral = model.createIntegral(mjj_set, ROOT.RooFit.NormSet(mjj_set), ROOT.RooFit.Range("SR")).getVal()
        sb1_region_integral = model.createIntegral(mjj_set, ROOT.RooFit.NormSet(mjj_set), ROOT.RooFit.Range("SB1")).getVal()
        sb2_region_integral = model.createIntegral(mjj_set, ROOT.RooFit.NormSet(mjj_set), ROOT.RooFit.Range("SB2")).getVal()
        full_integral = model.createIntegral(mjj_set, ROOT.RooFit.NormSet(mjj_set), full_range).getVal()

        ##histos_sb is the full histogram containg SR and SB
        sideband_integral = sb1_region_integral + sb2_region_integral
        full_norm = full_integral/sideband_integral
        lower_bin = histos_sb.FindBin(sb1_edge)
        upper_bin = histos_sb.FindBin(sb2_edge)
        SR_integral = histos_sb.Integral(lower_bin, upper_bin)
        SB_integral = histos_sb.Integral(1, lower_bin-1) + histos_sb.Integral(upper_bin+1, histos_sb.GetNbinsX())
        print "SR data integral = ", SR_integral
        print "SB data integral = ", SB_integral
        print "their sum = ",  SR_integral + SB_integral
        print "full data integral = ", histos_sb.Integral()       

        bg_expected = (full_norm-1.0)*SB_integral
        print "bg expected =", bg_expected
        
        # 1.) Simple
        signal = SR_integral - bg_expected
        if signal < 0.0:
            signal = 0.0
        simple_significance = signal/np.sqrt(bg_expected)
        true_significance = n_SR_sig_truth / np.sqrt(n_SR_bkg_truth)

        print "Significance based on truth level information: ", true_significance
        print "Oversimplified significance based on fit: ", simple_significance

        # 2.) more complicated version, taking into account fit systematics
        def signal_fit_func(par_list):
            current_model, current_workspace = QCD_model(par_list, name="current_model")
            current_mjj = current_workspace.var("mjj_fine")
            current_mjj.setRange("SB1", binsx[0], sb1_edge)
            current_mjj.setRange("SB2", sb2_edge, binsx[-1])
            current_mjj.setRange("SR", sb1_edge, sb2_edge)
            current_mjj.setBins(bins_fine)
            current_mjj_set = ROOT.RooArgSet(current_mjj)
            current_sig_region_integral = current_model.createIntegral(current_mjj_set, ROOT.RooFit.NormSet(current_mjj_set), ROOT.RooFit.Range("SR")).getVal()
            current_sb1_region_integral = current_model.createIntegral(current_mjj_set, ROOT.RooFit.NormSet(current_mjj_set), ROOT.RooFit.Range("SB1")).getVal()
            current_sb2_region_integral = current_model.createIntegral(current_mjj_set, ROOT.RooFit.NormSet(current_mjj_set), ROOT.RooFit.Range("SB2")).getVal()
            return float(current_sig_region_integral)/(current_sb1_region_integral+current_sb2_region_integral)

        obs = np.array([SR_integral])
        exp = bg_expected

        covma = fres.covarianceMatrix()
        np_covma = np.zeros((nPars, nPars))
        for j in range(nPars):
            for k in range(nPars):
                np_covma[j,k] = covma[j][k]

        print "covariance matrix =", np_covma
       
        jac = numdifftools.core.Jacobian(signal_fit_func)

        fit_params = fres.floatParsFinal()
        fit_param_list = []
        for j in range(nPars):
            fit_param_list.append(fit_params[j].getValV())
        print "fit parameters list = ", fit_param_list

        x_signal_cov = np.dot(np.dot(jac(fit_param_list), np_covma), jac(fit_param_list).T)
        print "x_signal_cov =", x_signal_cov
        sigma = np.sqrt(x_signal_cov)[0][0]

        ## asymptotic formula enhanced with doing the profile fit analytically
        def Zscorefunc(obs,pred,sigma):
            pred2=pred-sigma**2
            LLR = 0.25*( (4*pred2) - (2*obs) + ((pred2**2)/(sigma**2)) + (2*(sigma**2))
                        - ((pred2*np.sqrt( (pred2**2) + (4*obs*(sigma**2)) ))/(sigma**2))
                        - 4*obs*np.log((1/(2*obs))*(pred2 + np.sqrt((pred2**2) + (4*obs*(sigma**2))))) )
            return np.sqrt(2*LLR)

        sigma_tot = np.sqrt(sigma**2 + obs)
        print "Observed: ", obs
        print "Expected: ", exp
        print "Difference in statistical sigma ", (obs-exp)/np.sqrt(obs)
        print "Systemtic bkg sigma: ", sigma

        if obs<exp+1:
            Zvaltryme = 0
            pvaltryme = 0.5
        else:
            Zvaltryme = Zscorefunc(obs,exp,sigma)
            pvaltryme = 1-ROOT.Math.normal_cdf(Zvaltryme)
        try:
            iter(Zvaltryme)
            Zvaltryme = Zvaltryme[0]
        except:
            Zvaltryme = Zvaltryme

        print "Z SCORE: ", Zvaltryme
        print "P VALUE: ", pvaltryme
        significances[i] = Zvaltryme
        p_values[i] = pvaltryme


        # compute observed CLs limits:

        ## analytically fitted background
        def theta_hat(obs, pred, sigma):
            pred2 = pred-sigma**2
            return 0.5*(pred2 + np.sqrt(pred2**2 + 4*obs*(sigma**2)))

        ## asymptotic formula for significance to use in limits setting
        def Zscore_limit(signal, bkg_pred, sigma, obs):
            if obs-bkg_pred > signal:
                return 0.
            theta_hat_pred = theta_hat(obs, signal+bkg_pred, sigma)
            theta_hat_obs = theta_hat(obs, obs, sigma)
            LLR = obs*np.log(theta_hat_obs/theta_hat_pred) - (theta_hat_obs-theta_hat_pred) \
                  - 0.5*(sigma**(-2))*(((theta_hat_obs - obs)**2) - (theta_hat_pred - signal - bkg_pred)**2)
            return np.sqrt(2*LLR)

        ## compute CLs value, assuming we have computed the background-only p-value already
        def CLs(signal, bkg_pred, sigma, obs, p_b):
            return (1-ROOT.Math.normal_cdf(Zscore_limit(signal, bkg_pred, sigma, obs)))/(1-p_b)

        ## scanning over different signals until we reach a CLs value closest to 5%
        CLs_vals = np.zeros(int(obs)*100)
        signal_pred_vals = np.linspace(obs-exp, obs, int(obs)*100)
        for j, signal_pred in enumerate(signal_pred_vals):
            CLs_vals[j] = CLs(signal_pred, exp, sigma, obs, pvaltryme)
            if CLs_vals[j] < 0.05:
                break
        CLs_idx = (np.abs(CLs_vals - 0.05)).argmin()
        
        print "upper limit on signal events (obs): {} at {:.2f}% CL".format(signal_pred_vals[CLs_idx].item(), (1-CLs_vals[CLs_idx])*100)
        limits_obs[i] = signal_pred_vals[CLs_idx].item()

        
        # compute expected CLs limits including uncertainty bands

        ## asymptotiv significance function for median expected upper limit
        def Zscore_limit_exp(signal, bkg_pred, sigma):
            pred2 = bkg_pred-sigma**2
            theta_hat_asimov = theta_hat(bkg_pred, signal+bkg_pred, sigma)
            LLR = bkg_pred*(np.log(pred2) - np.log(theta_hat_asimov) - 1) + theta_hat_asimov \
                  - 0.5*(sigma**(-2))*(theta_hat_asimov - signal- bkg_pred)**2 + 0.5*(sigma**2)
            return np.sqrt(2*LLR)

        ## scanning over expected median CLs values
        def expected_limit(bkg_pred, sigma):
            CLs_median_exp_vals = np.zeros(int(100*obs))
            signal_pred_median_exp_vals = np.linspace(0, obs, int(obs)*100)
            Zscore_asimov_vals = np.zeros(int(100*obs))
            for j, signal_pred in enumerate(signal_pred_median_exp_vals):
                Zscore_asimov_vals[j] = Zscore_limit_exp(signal_pred, bkg_pred, sigma)
                CLs_median_exp_vals[j] = 2*(1-ROOT.Math.normal_cdf(Zscore_asimov_vals[j]))
                if CLs_median_exp_vals[j] < 0.05:
                    break
            CLs_median_idx = np.nanargmin(np.abs(CLs_median_exp_vals - 0.05))
            return float(signal_pred_median_exp_vals[CLs_median_idx]), CLs_median_exp_vals[CLs_median_idx], Zscore_asimov_vals[CLs_median_idx]

        ## returning n sigma errorbands
        def expected_limit_errors(S_upper_median, Zscore_median, n_sigma):
            return float(S_upper_median/Zscore_median*(ROOT.Math.gaussian_quantile(1 - 0.05*ROOT.Math.normal_cdf(n_sigma), 1) + n_sigma))
           
        limit_exp_median, CL_val_exp_median, Zscore_asimov_median = expected_limit(exp, sigma)
        print "upper limit on signal events (exp): {} at {:.2f}% CL".format(limit_exp_median, (1-CL_val_exp_median)*100)

        limit_exp_plus1sigma = expected_limit_errors(limit_exp_median, Zscore_asimov_median, 1)
        print "upper limit on signal events (exp+1sigma): {} at {:.2f}% CL".format(limit_exp_plus1sigma, (1-CL_val_exp_median)*100)
        limit_exp_minus1sigma = expected_limit_errors(limit_exp_median, Zscore_asimov_median, -1)
        print "upper limit on signal events (exp-1sigma): {} at {:.2f}% CL".format(limit_exp_minus1sigma, (1-CL_val_exp_median)*100)
        limit_exp_plus2sigma = expected_limit_errors(limit_exp_median, Zscore_asimov_median, 2)
        print "upper limit on signal events (exp+2sigma): {} at {:.2f}% CL".format(limit_exp_plus2sigma, (1-CL_val_exp_median)*100)
        limit_exp_minus2sigma = expected_limit_errors(limit_exp_median, Zscore_asimov_median, -2)
        print "upper limit on signal events (exp-2sigma): {} at {:.2f}% CL".format(limit_exp_minus2sigma, (1-CL_val_exp_median)*100)

        limits_exp_median[i] = limit_exp_median
        limits_exp_plus1sigma[i] = limit_exp_plus1sigma
        limits_exp_plus2sigma[i] = limit_exp_plus2sigma
        limits_exp_minus1sigma[i] = limit_exp_minus1sigma
        limits_exp_minus2sigma[i] = limit_exp_minus2sigma

        dataset.plotOn(frame,ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson),ROOT.RooFit.Name("data_qcd_blind"),ROOT.RooFit.XErrorSize(0),ROOT.RooFit.Binning(roobins))
        model.plotOn(frame,ROOT.RooFit.Invisible(),ROOT.RooFit.Name(model_name), full_range, norm) #for correct pulls
        #average bin edges instead of bin center
        framePulls = mjj.frame()
        useBinAverage = True
        hpull = frame.pullHist("data_qcd_blind",model_name,useBinAverage)
        framePulls.addPlotable(hpull,"X0 P E1")
        
        my_chi2,my_ndof = calculateChi2(hpull,nPars, ranges = fit_ranges)
        my_prob = ROOT.TMath.Prob(my_chi2,my_ndof)
        PlotFitResults(frame,fres.GetName(),nPars,framePulls,"data_qcd_blind",model_name,my_chi2,my_ndof, str(nPars) + "par_qcd_fit_binned_blinded", plot_dir)
        
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
        fitter_QCD.delete()

    best_i = f_test(nParsToTry, ndofs, chi2s)


    nPars_QCD = nParsToTry[best_i]
    qcd_fname = qcd_fnames[best_i]


    print("\n Chose %i parameters based on F-test ! \n" % nPars_QCD )

    print "best fitting significance:", significances[best_i]
    print "best fitting p value:", p_values[best_i]
    print "best fitting 95% CL upper limit (obs):", limits_obs[best_i]
    print "best fitting 95% CL upper limit (exp):", limits_exp_median[best_i]
    print "best fitting 95% CL upper limit exp 1sigma bands:", limits_exp_minus1sigma[best_i], "to", limits_exp_plus1sigma[best_i]
    print "best fitting 95% CL upper limit exp 2sigma bands:", limits_exp_minus2sigma[best_i], "to", limits_exp_plus2sigma[best_i]

    results = dict()

    #QCD fit results
    results['chi2']  = chi2s[best_i]
    results['ndof']  = ndofs[best_i]
    results['fit_prob']  = probs[best_i]
    results['nPars_QCD']  = nPars_QCD

    results['signif'] = significances[best_i]
    results['p_value'] = p_values[best_i]
    results['obs_lim_events'] = limits_obs[best_i]*sig_norm
    results['exp_lim_events'] = limits_exp_median[best_i]*sig_norm
    results['exp_lim_plus1sigma_events'] = limits_exp_plus1sigma[best_i]*sig_norm
    results['exp_lim_minus1sigma_events'] = limits_exp_minus1sigma[best_i]*sig_norm
    results['exp_lim_plus2sigma_events'] = limits_exp_plus2sigma[best_i]*sig_norm
    results['exp_lim_minus2sigma_events'] = limits_exp_minus2sigma[best_i]*sig_norm

    return results
        

def fitting_options():
    parser = optparse.OptionParser()
    #parser.add_option("--sig_norm",type=float, default=1.0, help="Scale signal pdf normalization by this amount")
    parser.add_option("-M","-M",dest="mass",type=float,default=2500.,help="Injected signal mass")
    parser.add_option("-i","--inputFile",dest="inputFile",default='fit_inputs/no_selection_03p.h5',help="input h5 file")
    parser.add_option("-p","--plotDir",dest="plotDir",default='plots/',help="Where to put the plots")
    #arser.add_option("-l","--label",dest="label",default='test',help="Label for file names")
    return parser

if __name__ == "__main__":
    parser = fitting_options()
    (options,args) = parser.parse_args()
    dijetfit(options)


