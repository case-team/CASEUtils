import os
from Fitter import Fitter
from DataCardMaker import DataCardMaker
import Utils as utils
import numpy as np
from array import array
import pickle
import json
import random
import optparse

import tdrstyle
import ROOT
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
tdrstyle.setTDRStyle()
ROOT.gROOT.SetBatch(True)
ROOT.RooRandom.randomGenerator().SetSeed(random.randint(0, 1e+6))





def bias_test(options):

    label = options.label
    plot_dir = options.plotDir
    if(plot_dir[-1] != '/'):
        plot_dir += '/'

    if(not os.path.exists(plot_dir)):
        os.system("mkdir %s" % plot_dir)

    fout_name = plot_dir + "bias_test_results_%s_%s" % (options.mass, label)
    if(os.path.exists(fout_name + ".json")):
        #remove old results
        os.system("rm %s" % fout_name + ".json")

    fine_bin_size = 4
    mass = options.mass 

    binsx = [1460, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438,
             2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854,
             4010, 4171, 4337, 4509, 4700, 4900,  5100, 5300, 5500, 5800,
             6100, 6400, 6800]

    if(options.mjj_max < 0. and options.rebin): 
        options.mjj_max = utils.get_mjj_max(options.inputFile) + 5.0
        options.mjj_max = max(1.2*options.mass, options.mjj_max)
    #print("MJJ MAX %.2f" % options.mjj_max)
    
    if(options.mjj_min > 0 and options.mjj_min < binsx[-1]):
        start_idx = 0
        while(binsx[start_idx] < options.mjj_min):
            start_idx +=1
        binsx = binsx[start_idx:]

        if(abs(options.mjj_min - binsx[0]) < 50.):
            binsx[0] = options.mjj_min
        else:
            binsx.insert(0, options.mjj_min)
        print("Will start fit from %.0f GeV" % binsx[0])

    if(options.mjj_max > 0 and options.mjj_max < binsx[-1]):
        print("rebinning with max mjj %.2f" % options.mjj_max)
        end_idx = len(binsx)-1
        while(binsx[end_idx]   > options.mjj_max and end_idx > 0): 
            end_idx -=1
        binsx = binsx[:end_idx]

        if(abs(options.mjj_max - binsx[-1]) < 50.):
            binsx[-1] = options.mjj_max
        else:
            binsx.append(options.mjj_max)
        print("Will end fit at %.0f GeV" % binsx[-1])
        print(binsx)


    # round to smallest precision we are storing mass values with, otherwise
    # get weird effects related to bin size
    utils.roundTo(binsx, fine_bin_size)



    nbins_fine = int(binsx[-1] - binsx[0])/fine_bin_size

    histos_sb = ROOT.TH1F("mjj_sb", "mjj_sb" ,nbins_fine, binsx[0], binsx[-1])
    
    
    utils.load_h5_sb(options.inputFile, histos_sb)
    print("************ Found %i total events \n" % histos_sb.GetEntries())
    num_evts = histos_sb.GetEntries()

    

    if(options.rebin):
        bins_nonzero = utils.get_rebinning(binsx, histos_sb)
        print("Rebinning to avoid zero bins!")
        print("old", binsx)
        print("new", bins_nonzero)
        bins = bins_nonzero
    else:
        bins = binsx
    roobins = ROOT.RooBinning(len(bins)-1, array('d', bins), "mjjbins")


    if(not os.path.exists(options.sig_shape)):
        print("Sig file %s doesn't exist" % options.sig_shape)
        exit(1)
    sig_file_name = options.sig_shape







    fitting_histogram = histos_sb
    data_name = "data_qcd"

    #Do fit with alt bkg shape
    nPars = 4
    qcd_alt_fname = "altBkg_fit.root"


    qcd_outfile = ROOT.TFile(qcd_alt_fname, 'RECREATE')

    model_name = "model_alt"
    fitter_QCD = Fitter(['mjj_fine'], debug = False)
    if(options.alt_shape_ver != 1 and options.alt_shape_ver != 2):
        print("Unsupported alt_shape_ver %i " % options.alt_shape_ver)
        exit(1)
    alt_bkg_pdf = fitter_QCD.altBkgShape(model_name, 'mjj_fine', options.alt_shape_ver)
    fitter_QCD.importBinnedData(fitting_histogram, ['mjj_fine'], data_name)
    
    #Running fit two times seems to improve things (better initial guesses for params?)
    
    fres = fitter_QCD.fit(model_name, data_name, options=[ROOT.RooFit.Save(1), ROOT.RooFit.Verbose(0),  ROOT.RooFit.Minos(1)])
    fres = fitter_QCD.fit(model_name, data_name, options=[ROOT.RooFit.Save(1), ROOT.RooFit.Verbose(0),  ROOT.RooFit.Minos(1)])

    chi2_fine = fitter_QCD.projection(
        model_name, data_name, "mjj_fine",
        plot_dir + "alt_qcd_fit_%s.png" % label, 0, True)


    qcd_outfile.cd()

    mjj = fitter_QCD.getVar('mjj_fine')
    mjj.setBins(nbins_fine)
    model = fitter_QCD.getFunc(model_name)
    dataset = fitter_QCD.getData(data_name)

    #rescale so pdfs are in evts per 100 GeV
    low = roobins.lowBound()
    high = roobins.highBound()
    n = roobins.numBoundaries() - 1
    #RootFit default normalization is full range divided by number of bins
    default_norm = (high - low)/ n
    rescale = 100./ default_norm
    fit_norm = ROOT.RooFit.Normalization(rescale,ROOT.RooAbsReal.Relative)



    frame = mjj.frame()
    dataset.plotOn(frame, ROOT.RooFit.Name(data_name), ROOT.RooFit.Invisible(), ROOT.RooFit.Binning(roobins), ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2), 
            ROOT.RooFit.Rescale(rescale))

    model.plotOn(frame, ROOT.RooFit.VisualizeError(fres, 1), ROOT.RooFit.FillColor(ROOT.kRed - 7), ROOT.RooFit.LineColor(ROOT.kRed - 7), ROOT.RooFit.Name(fres.GetName()), 
                   fit_norm)

    model.plotOn(frame, ROOT.RooFit.LineColor(ROOT.kRed + 1), ROOT.RooFit.Name(model_name),  fit_norm)
    

    useBinAverage = True
    hpull = frame.pullHist(data_name, model_name, useBinAverage)
    hresid = frame.residHist(data_name, model_name, False, useBinAverage)
    dhist = ROOT.RooHist(frame.findObject(data_name, ROOT.RooHist.Class()))

    #print_roohist(dhist)

    #redraw data (so on top of model curves)
    if(options.rebin):
        dataset.plotOn(frame, ROOT.RooFit.Name(data_name),   ROOT.RooFit.Binning(roobins), ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2),
                   ROOT.RooFit.Rescale(rescale))
    else:
        dataset.plotOn(frame, ROOT.RooFit.Name(data_name),  ROOT.RooFit.XErrorSize(0), ROOT.RooFit.Binning(roobins), ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2),
                   ROOT.RooFit.Rescale(rescale))


    framePulls = mjj.frame()
    framePulls.addPlotable(hpull, "X0 P E1")


    #get fractional error on fit, evaluated at signal mass
    central = frame.getCurve(model_name);
    curve =  frame.getCurve("fitresults");
    upBound = ROOT.TGraph(central.GetN());
    loBound = ROOT.TGraph(central.GetN());

    for j in range(curve.GetN()):
        if( j < central.GetN() ): upBound.SetPoint(j, curve.GetX()[j], curve.GetY()[j]);
        else: loBound.SetPoint(j, curve.GetX()[j], curve.GetY()[j]);

    err_on_sig = (upBound.Eval(options.mass) - loBound.Eval(options.mass))/2.
    frac_err_on_sig = err_on_sig / central.Eval(options.mass)
    bkg_fit_frac_err = frac_err_on_sig



    my_chi2, my_ndof = utils.calculateChi2(hpull, nPars, excludeZeros = True, dataHist = dhist)
    my_prob = ROOT.TMath.Prob(my_chi2, my_ndof)
    utils.PlotFitResults(frame, fres.GetName(), nPars, framePulls, data_name,
                   [model_name], my_chi2, my_ndof,
                   "altBkg_qcd_fit_binned",
                   plot_dir, plot_label = label)

    graphs = {}
    for p in range(nPars):
        graphs['ap%i' % (p + 1)] = ROOT.TGraphErrors()


    #largest_frac_err = 0.
    for var, graph in graphs.iteritems():
        print(var)
        value, error = fitter_QCD.fetch(var)
        graph.SetPoint(0, mass, value)
        graph.SetPointError(0, 0.0, error)
        #frac_err = abs(error/value)
        #largest_frac_err = max(frac_err, largest_frac_err)

    qcd_outfile.cd()
    for name, graph in graphs.iteritems():
        graph.Write(name)
    qcd_outfile.Close()

    print("#############################")
    print("Alt bkg shape results: " )
    print("bkg fit chi2/nbins (fine binning) ", chi2_fine)
    print("My chi2, ndof, prob", my_chi2, my_ndof, my_prob)
    print("My chi/ndof, chi2/nbins", my_chi2/my_ndof, my_chi2/(my_ndof + nPars))
    print("Fit func fractional unc at sig mass ", bkg_fit_frac_err)
    print("#############################")

    
    # #################
    #Define PDF's for bias test

    reg_model_name = "model_reg_bkg"
    reg_bkg_pdf = fitter_QCD.qcdShape(reg_model_name, 'mjj_fine', nPars)



    sb_fname = "sb_altBkg_fit.root"
    sb_outfile = ROOT.TFile(sb_fname, 'RECREATE')
    sb_outfile.cd()
    histos_sb.Write("mjj_sb")
    sb_outfile.Close()
    sig_data_name = 'mjj_sb'
    sb_label = "raw"

    workspace = fitter_QCD.w

    dcb_gen_pdf, pars1 = fitter_QCD.addDCBSignalShape('model_signal_mjj_gen', 'mjj_fine', sig_file_name, {'CMS_scale_j': 1.0}, {'CMS_res_j': 1.0})
    dcb_fit_pdf, pars2 = fitter_QCD.addDCBSignalShape('model_signal_mjj_fit', 'mjj_fine', sig_file_name, {'CMS_scale_j': 1.0}, {'CMS_res_j': 1.0})

    constant = options.sig_norm

    #gen_mu = 0.1
    gen_mu = options.num_sig / (num_evts + options.num_sig)
    num_gen = options.num_sig + num_evts



    sig_strength_gen =  ROOT.RooRealVar("sig_strength_gen", "sig_strength_gen", gen_mu, -100., 100.)
    getattr(workspace, 'import')(sig_strength_gen, ROOT.RooFit.Rename("sig_strength_gen"))
    sig_strength_fit =  ROOT.RooRealVar("sig_strength_fit", "sig_strength_fit", gen_mu, -100., 100.)
    getattr(workspace, 'import')(sig_strength_fit, ROOT.RooFit.Rename("sig_strength_fit"))

    

    SB_alt_pdf = ROOT.RooAddPdf("SB_alt", "Sig + Alt bkg", ROOT.RooArgList(dcb_gen_pdf, alt_bkg_pdf), ROOT.RooArgList(sig_strength_gen))
    getattr(workspace, 'import')(SB_alt_pdf, ROOT.RooFit.Rename("SB_alt_pdf"))
    SB_fit_pdf = ROOT.RooAddPdf("SB_fit", "Sig + 4-par bkg", ROOT.RooArgList(dcb_fit_pdf, reg_bkg_pdf), ROOT.RooArgList(sig_strength_fit))
    getattr(workspace, 'import')(SB_fit_pdf, ROOT.RooFit.Rename("SB_fit_pdf"))

    # #################

    #Use RooMCStudy class to generate / fit many toys


    mcs = ROOT.RooMCStudy(SB_alt_pdf, ROOT.RooArgSet(fitter_QCD.w.var("mjj_fine")), ROOT.RooFit.FitModel(SB_fit_pdf),  ROOT.RooFit.Silence(), ROOT.RooFit.Binned(), 
        ROOT.RooFit.FitOptions(ROOT.RooFit.Save(True), ROOT.RooFit.PrintLevel(-1)))

    keepGenData = True

    mcs.generateAndFit(options.num_samples, int(round(num_gen)), keepGenData)

    genDataset = mcs.genData(0)
    fit_res = mcs.fitResult(0)
    toy_data_name = genDataset.GetName()
    toy_fit_name = fit_res.GetName()

    #set pdf pars equal to fit values from this toy for plotting
    pars = ROOT.RooArgSet()
    SB_fit_pdf.getParameters(pars)
    pars = fit_res.floatParsFinal()


    frame_toy = mjj.frame()
    genDataset.plotOn(frame_toy, ROOT.RooFit.Invisible(), ROOT.RooFit.Binning(roobins), ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2), ROOT.RooFit.Name(toy_data_name),
            ROOT.RooFit.Rescale(rescale))

    SB_fit_pdf.plotOn(frame_toy, ROOT.RooFit.VisualizeError(fit_res, 1), ROOT.RooFit.FillColor(ROOT.kRed - 7), ROOT.RooFit.LineColor(ROOT.kRed - 7), ROOT.RooFit.Name(fit_res.GetName()), 
                   fit_norm)

    SB_fit_pdf.plotOn(frame_toy,  ROOT.RooFit.Components(reg_model_name), ROOT.RooFit.LineColor(ROOT.kMagenta + 3),ROOT.RooFit.Name("Background"), fit_norm)
    SB_fit_pdf.plotOn(frame_toy, ROOT.RooFit.LineColor(ROOT.kRed + 1), ROOT.RooFit.Name("SB_fit"),  fit_norm)


    hpull_toy = frame_toy.pullHist()
    chi2 = frame_toy.chiSquare()

    frame_toyPulls = mjj.frame()
    frame_toyPulls.addPlotable(hpull_toy, "X0 P E1")


    SB_fit_pdf.plotOn(frame_toy, ROOT.RooFit.Components("model_signal_mjj_fit"), ROOT.RooFit.LineColor(ROOT.kBlue), fit_norm, ROOT.RooFit.Name("Signal"))

    if(options.rebin):
        genDataset.plotOn(frame_toy,ROOT.RooFit.Binning(roobins), ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2),ROOT.RooFit.Name(toy_data_name),
                   ROOT.RooFit.Rescale(rescale))
    else:
        genDataset.plotOn(frame_toy,  ROOT.RooFit.XErrorSize(0), ROOT.RooFit.Binning(roobins), ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2),ROOT.RooFit.Name(toy_data_name),
                   ROOT.RooFit.Rescale(rescale))

    #c_toy =ROOT.TCanvas("c_toy","",800,800)
    #frame_toy.Draw()
    #c_toy.Print( plot_dir + "toy_fit_test.png")





    utils.PlotFitResults(frame_toy, fit_res.GetName(), nPars, frame_toyPulls, toy_data_name,
                   ["SB_fit", "Signal", "Background"] , chi2, -1, "fit_toy1_%s" % label,
                   plot_dir, plot_label = label, has_sig = True, draw_sig = True)


    fitted_mu = []
    uncs = []
    pulls = []
    signifs = []

    for i in range(options.num_samples):

        mu = mcs.fitResult(i).floatParsFinal().find("sig_strength_fit")
        #mcs.fitResult(i).floatParsFinal().Print()
        #mu.Print()
        fitted_mu.append(mu.getVal())
        uncs.append(mu.getError())
        signifs.append(fitted_mu[i] / uncs[i])
        pulls.append((fitted_mu[i] - gen_mu) /uncs[i])


    resids = [(fit_mu - gen_mu)  for fit_mu in fitted_mu]

    mean_mu = np.mean(fitted_mu)
    mean_unc = np.mean(uncs)
    mean_resid = np.mean(resids)
    mean_signif = np.mean(signifs)
    mean_pull, std_pull = np.mean(pulls), np.std(pulls)
    print("Gen Mu %.3e " % gen_mu)
    print("Mean fitted Mu %.3e " % (mean_mu))
    print("Mean fitted unc %.3e " % mean_unc)
    print("Mean signif ~ %.2f " % (mean_signif))
    print("Pull Mean %.2f Std %.2f " % (mean_pull, std_pull))
    #check_rough_sig(options.inputFile, options.mass*0.9, options.mass*1.1)


    make_histo(resids, xmin = -5 * mean_unc, xmax = 5 * mean_unc, fout = plot_dir + "mu_residuals_%s.png" % label, xlabel = "(#mu_{fit} - #mu_{gen})", fit_gaus = True)
    if(gen_mu > 0.):
        resids_norm = [(fit_mu - gen_mu)/gen_mu  for fit_mu in fitted_mu]
        make_histo(resids, xmin = -5 * mean_unc/gen_mu, xmax = 5 * mean_unc/gen_mu, fout = plot_dir + "mu_residuals_norm_%s.png" % label, xlabel = "(#mu_{fit} - #mu_{gen})/#mu_{gen}", fit_gaus = True)
    make_histo(pulls, xmin = -5, xmax = 5, fout = plot_dir + "mu_pulls_%s.png" %label, xlabel = "(#mu_{fit} - #mu_{gen})/#sigma_{#mu}", fit_gaus = True)


    results = dict()
    results['script_options'] = vars(options)
    results['altbkg_chi2'] = my_chi2
    results['altbkg_ndof'] = my_ndof
    results['altbkg_prob'] = my_prob
    results['resids'] = resids
    results['uncs'] = uncs
    results['gen_mu'] = gen_mu
    results['mean_pull'] = mean_pull
    results['stddev_pull'] = std_pull
    results['mean_resid'] = mean_resid
    results['mean_signif'] = mean_signif
    results['mean_unc'] = mean_unc


    print("Saving fit results to %s" % ( fout_name + ".json"))
    with open(fout_name + ".json", "w") as jsonfile:
        json.dump(results, jsonfile, indent=4)

    return results

def make_histo(data, xmin = -99999., xmax = -99999., nbins = 20, fout = "", xlabel = "", fit_gaus = True):

    ROOT.gStyle.SetOptStat(2210)
    ROOT.gStyle.SetOptFit(11)

    if(xmin == -99999.): xmin = np.amin(data)
    if(xmax == -99999.): xmax = np.amax(data)

    h = ROOT.TH1F(fout + "_h", fout + "_h", nbins, xmin, xmax)
    for d in data:
        h.Fill(d)

    c = ROOT.TCanvas(fout + "_c", fout + "_c", 1000, 1000)
    h.GetXaxis().SetTitle(xlabel)
    h.GetXaxis().SetTitleSize(0.06)
    h.GetXaxis().SetLabelSize(0.04)
    h.SetLineColor(ROOT.kBlack)
    h.Draw("hist")
    if(fit_gaus):
        h.Fit("gaus")
        fit = h.GetFunction("gaus")
        fit.SetParName(1, "Gaussian : mean")
        fit.SetParName(2, "Gaussian : sigma")
        mean = fit.GetParameter(1)
        mean_err = fit.GetParError(1)
        sigma = fit.GetParameter(2)
        sigma_err = fit.GetParError(2)
        
        fit.SetLineColor(ROOT.kBlue)
        fit.SetLineWidth(3)
        fit.Draw("same")

        #tag1 = ROOT.TLatex(0.70,0.60,"#mu_{r}: %.2e +/- %.2e "%(mean, mean_err))
        #tag1.SetNDC(); tag1.SetTextFont(42); tag1.SetTextSize(0.04);
        #tag2 = ROOT.TLatex(0.70,0.57,"#sigma_{r}: %.2e +/- %.2e "%(sigma, sigma_err))
        #tag2.SetNDC(); tag1.SetTextFont(42); tag1.SetTextSize(0.04);
        #tag1.Draw()
        #tag2.Draw()


    if(len(fout) != 0): c.Print(fout)





def fitting_options():
    parser = optparse.OptionParser()
    parser.add_option("--mjj_min", type=float, default=-1.0,
                      help="Minimum mjj for the fit")
    parser.add_option("--num_sig", type=int, default=100, help="How many signal events to generate for the toys")
    parser.add_option("--num_samples", type=int, default=10, help="How many toys to generate + fit for the bias test" )
    parser.add_option("--alt_shape_ver", type=int, default=2, help="Which version of alt bkg shape pdf (1 or 2)" )
    parser.add_option("--mjj_max", type=float, default=-1.0,
                      help="Maximum mjj for the fit")
    parser.add_option("--sig_norm", type=float, default=1.0,
                      help="Scale signal pdf normalization by this amount")
    parser.add_option("--ftest_thresh", type=float, default=0.05,
                      help="Threshold to prefer a function in the f-test")
    parser.add_option("--err_thresh", type=float, default=0.5,
                      help="Threshold on fit unc to be included in f-test")
    parser.add_option("-s", "--sig_shape", default="signal_shape_m2500.root",
                      help="Pre-saved signal shape")
    parser.add_option("--rebin", default=False, action="store_true",
                      help="""Rebin dijet bins to make sure no bins less than 5 evts""")
    parser.add_option("-M", "-M", dest="mass", type=float, default=2500.,
                      help="Injected signal mass")
    parser.add_option("-i", "--inputFile", dest="inputFile",
                      default='fit_inputs/no_selection_03p.h5',
                      help="input h5 file")
    parser.add_option("-p", "--plotDir", dest="plotDir", default='plots/',
                      help="Where to put the plots")
    parser.add_option("-l", "--label", dest="label", default='test',
                      help="Label for file names")
    parser.add_option("--dcb-model", dest="dcb_model", action="store_true",
                      default=False,
                      help="""Whether to use double crystal ball model instead
                      of default model (gaussian core with single crystal ball)""")
    parser.add_option("--sig_norm_unc", dest="sig_norm_unc", type=float, default= -1.0, help="Fractional uncertainty on signal normalization")
    return parser


if __name__ == "__main__":
    parser = fitting_options()
    (options, args) = parser.parse_args()
    bias_test(options)

