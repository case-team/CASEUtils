import h5py, math, commands, random
from array import array
import numpy as np
import time, sys, os, optparse, json, copy
import copy

import ROOT
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
import CMS_lumi, tdrstyle
tdrstyle.setTDRStyle()
ROOT.gROOT.SetBatch(True)
ROOT.RooRandom.randomGenerator().SetSeed(random.randint(0, 1e+6))

# function to extract string from function -> in order write a proper json file
def returnString(func,ftype):
    if func.GetName().find("corr")!=-1:
        st = "("+str(func.GetParameter(0))+" + ("+str(func.GetParameter(1))+")*MJ1 + ("+str(func.GetParameter(2))+")*MJ2  + ("+str(func.GetParameter(3))+")*MJ1*MJ2)"
        if func.GetName().find("sigma")!=-1:
            st = "("+str(func.GetParameter(0))+" + ("+str(func.GetParameter(1))+")*MJ1 + ("+str(func.GetParameter(2))+")*MJ2 )"
        return st
    else:
        if ftype.find("pol")!=-1:
            st='(0'
            if func.GetName().find("corr")!=-1: 
                n = 1. #func.Integral(55,215)
                st = "(0"
                for i in range(0,func.GetNpar()):
                    st = st+"+("+str(func.GetParameter(i))+")"+("*(MJ1+MJ2)/2."*i)
                st+=")/"+str(n)
            else:
                for i in range(0,func.GetNpar()):
                    st=st+"+("+str(func.GetParameter(i))+")"+("*MH"*i)
                st+=")"
            return st
        if ftype.find("1/sqrt")!=-1:
            st='(0'
            if func.GetName().find("corr")!=-1:
                n = 1. # func.Integral(55,215)
                st = str(func.GetParameter(0))+"+("+str(func.GetParameter(1))+")*1/sqrt((MJ1+MJ2)/2.)/"+str(n)
            else:
                st = str(func.GetParameter(0))+"+("+str(func.GetParameter(1))+")"+")*1/sqrt(MH)"
                st+=")"
            return st
        if ftype.find("sqrt")!=-1 and ftype.find("1/")==-1:
            n =1.
            st='(0'
            if func.GetName().find("corr")!=-1: st = str(func.GetParameter(0))+"+("+str(func.GetParameter(1))+")"+"*sqrt((MJ1+MJ2)/2.))/"+str(n)
            else:
                st = str(func.GetParameter(0))+"+("+str(func.GetParameter(1))+")"+"*sqrt(MH)"
                st+=")"
            return st    
        if ftype.find("llog")!=-1:
            return str(func.GetParameter(0))+"+"+str(func.GetParameter(1))+"*log(MH)"
        if ftype.find("laur")!=-1:
            st='(0'
            for i in range(0,func.GetNpar()):
                st=st+"+("+str(func.GetParameter(i))+")"+"/MH^"+str(i)
            st+=")"
            return st    
        if ftype.find("spline")!=-1:
            print "write json for spline function: a list and not a string will be returned in this case"
            st=[]
            nnknots = func.GetNp()
            for i in range(0,nnknots):
                x = ROOT.Double(0) 
                y = ROOT.Double(0) 
                func.GetKnot(i,x,y)
                st.append([x,y])
            return st
        else:
            return ""

def get_palette(mode):

 palette = {}
 palette['gv'] = [] 
 
 colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
              '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
              '#bcbd22', '#17becf']

 for c in colors:
  palette['gv'].append(c)
 
 return palette[mode]

def convert_matrix(mat):
    mat_arr = mat.GetMatrixArray()
    return [[mat_arr[i + j*mat.GetNrows()] for j in range(mat.GetNcols())] for i in range(mat.GetNrows())]

def getBinning(binsMVV,minx,maxx,bins):
    l=[]
    if binsMVV=="":
        print " do this"
        print binsMVV
        for i in range(0,bins+1):
            l.append(minx + i* (maxx - minx)/bins)
    else:
        print "dot that"
        print binsMVV
        s = binsMVV.split(",")
        for w in s:
            l.append(int(w))
    return l

def truncate(binning,mmin,mmax):
    res=[]
    for b in binning:
        if b >= mmin and b <= mmax:
            res.append(b)
    return res


def PlotFitResults(frame,fitErrs,nPars,pulls,data_name,pdf_names,chi2,ndof,canvname, plot_dir, has_sig = False, plot_label = ""):

    c1 =ROOT.TCanvas("c1","",800,800)
    c1.SetLogy()
    c1.Divide(1,2,0,0,0)
    c1.SetLogy()
    c1.cd(1)
    p11_1 = c1.GetPad(1)
    p11_1.SetPad(0.01,0.26,0.99,0.98)
    p11_1.SetLogy()
    p11_1.SetRightMargin(0.05)



    p11_1.SetTopMargin(0.1)
    p11_1.SetBottomMargin(0.02)
    p11_1.SetFillColor(0)
    p11_1.SetBorderMode(0)
    p11_1.SetFrameFillStyle(0)
    p11_1.SetFrameBorderMode(0)
    frame.GetYaxis().SetTitleSize(0.06)
    frame.GetYaxis().SetTitleOffset(0.98)
    frame.SetMinimum(0.2)
    frame.SetMaximum(1E7)
    frame.SetName("mjjFit")
    frame.GetYaxis().SetTitle("Events / 100 GeV")
    frame.SetTitle("")
    frame.Draw()
        
    legend = ROOT.TLegend(0.45097293,0.54183362,0.6681766,0.779833)
    legend2 = ROOT.TLegend(0.45097293,0.54183362,0.6681766,0.779833)
    legend.SetTextSize(0.046)
    legend.SetLineColor(0)
    legend.SetShadowColor(0)
    legend.SetLineStyle(1)
    legend.SetLineWidth(1)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)
    legend.SetMargin(0.35)
    legend2.SetTextSize(0.038)
    legend2.SetLineColor(0)
    legend2.SetShadowColor(0)
    legend2.SetLineStyle(1)
    legend2.SetLineWidth(1)
    legend2.SetFillColor(0)
    legend2.SetFillStyle(0)
    legend2.SetMargin(0.35)
    legend.AddEntry(frame.findObject(data_name),"Data","lpe")
    if(not has_sig): 
        legend.AddEntry(frame.findObject(pdf_names[0]),"%i par. background fit"%nPars,"l")

    else: 
        legend.AddEntry(frame.findObject(pdf_names[0]),"Signal + Background Fit ","l")
    legend2.AddEntry("","","")
    legend2.AddEntry("","","")
    legend2.AddEntry("","","")
    legend2.AddEntry("","","")
    legend2.AddEntry(frame.findObject(fitErrs),"","f")
    legend2.AddEntry("","","")

    legend2.Draw("same")
    legend.Draw("same")

    pt = ROOT.TPaveText(0.18,0.06,0.54,0.17,"NDC")
    pt.SetTextFont(42)
    pt.SetTextAlign(12)
    pt.SetFillColor(0)
    pt.SetBorderSize(0)
    pt.SetFillStyle(0)
    pt.AddText("Chi2/ndf = %.2f/%i = %.2f"%(chi2,ndof,chi2/ndof))
    pt.AddText("Prob = %.3f"%ROOT.TMath.Prob(chi2,ndof))
    pt.Draw()

    pt2 = ROOT.TPaveText(0.5,0.8,0.6,0.9,"NDC")
    pt2.SetTextFont(42)
    pt2.SetTextAlign(22)
    pt2.SetFillColor(0)
    pt2.SetBorderSize(0)
    pt2.SetFillStyle(0)
    pt2.AddText(plot_label)
    pt2.Draw()

    
    c1.Update()

    c1.cd(2)
    p11_2 = c1.GetPad(2)
    p11_2.SetPad(0.01,0.02,0.99,0.27)
    p11_2.SetBottomMargin(0.35)
    p11_2.SetRightMargin(0.05)
    p11_2.SetGridx(0)
    p11_2.SetGridy(0)
    pulls.SetMinimum(-5)
    pulls.SetMaximum(5)
    pulls.SetTitle("")
    pulls.SetXTitle("Dijet invariant mass (GeV)")
    pulls.GetXaxis().SetTitleSize(0.06)
    pulls.SetYTitle("#frac{Data-Fit}{#sigma_{data}}")
    pulls.GetYaxis().SetTitleSize(0.15)
    pulls.GetYaxis().CenterTitle()
    pulls.GetYaxis().SetTitleOffset(0.30)
    pulls.GetYaxis().SetLabelSize(0.15)
    pulls.GetXaxis().SetTitleSize(0.17)
    pulls.GetXaxis().SetTitleOffset(0.91)
    pulls.GetXaxis().SetLabelSize(0.12)
    pulls.GetXaxis().SetNdivisions(906)
    pulls.GetYaxis().SetNdivisions(305)
    pulls.Draw("same")
    line = ROOT.TLine(frame.GetXaxis().GetXmin() , 0 , frame.GetXaxis().GetXmax(),0)
    line1  = ROOT.TLine(frame.GetXaxis().GetXmin(), 1 ,frame.GetXaxis().GetXmax(),1)
    line2  = ROOT.TLine(frame.GetXaxis().GetXmin(), -1 ,frame.GetXaxis().GetXmax(),-1)
    line1.SetLineStyle(2)
    line1.SetLineWidth(2)
    line2.SetLineStyle(2)
    line2.SetLineWidth(2)
    line.Draw("same")
    line1.Draw("same")
    line2.Draw("same")   
    c1.Update()

    canvname+='.png'
    c1.SaveAs(plot_dir + canvname)
    #c1.SaveAs(canvname.replace("png","C"),"C")

def print_roohist(h):
    a_x = array('d', [0.])
    a_val = array('d', [0.])
    a_data = array('d', [0.])

    for p in range (0,h.GetN()):
        h.GetPoint(p, a_x, a_val)
        e = h.GetErrorY(p)
        print("%i %.3f %.3f %.3f" % (p, a_x[0], a_val[0], e))



def calculateChi2(g_pulls, nPars, ranges = None, excludeZeros = True, dataHist = None):
     
    NumberOfVarBins = 0
    NumberOfObservations_VarBin = 0
    chi2_VarBin = 0.
     

    a_x = array('d', [0.])
    a_val = array('d', [0.])
    a_data = array('d', [0.])

    last_nonzero = 0

    if(excludeZeros):
        #find max non-zero bin
        for p in range (0,g_pulls.GetN()):
            g_pulls.GetPoint(p, a_x, a_val)
            dataHist.GetPoint(p, a_x, a_data)
            if(a_data[0] > 0.): last_nonzero = p

    for p in range (0,g_pulls.GetN()):
    
        g_pulls.GetPoint(p, a_x, a_val)
        x = a_x[0]
        pull = a_val[0]


        add = True
        if(ranges is not None and len(ranges) > 0):
            add = False
            for range_ in ranges:
                if(x >= range_[0] and x<= range_[1]):
                    add = True
         
        if(excludeZeros):
            #include 1 zero after last bin
            if(p > last_nonzero +1):
                add = False

        if(add):
            if (dataHist is not None ):
                dataHist.GetPoint(p, a_x, a_data)
                #print(x, a_data[0], pull)
            NumberOfObservations_VarBin+=1
            chi2_VarBin += pow(pull,2)
            
    ndf_VarBin = NumberOfObservations_VarBin - nPars
    return chi2_VarBin,ndf_VarBin

def apply_blinding(h, ranges = None):
    if(ranges is None or len(ranges) == 0):
        print("Must supply list of tuples specifying ranges to include")

    h_clone = h.Clone(h.GetName() + "_blinded")
    axis = h.GetXaxis()
    for i in range(axis.GetNbins()):
        low_edge = axis.GetBinLowEdge(i)
        high_edge = axis.GetBinLowEdge(i)

        inRange = False
        for interval in ranges:
            if(low_edge >= interval[0] and high_edge <= interval[1]):
                inRange = True

        if(not inRange):
            h_clone.SetBinContent(i, 0.)
            h_clone.SetBinError(i, 0.)
    return h_clone

def load(iFile,iVar,iName,iHist,iCut):
    lFile = ROOT.TFile(iFile)
    lTree = lFile.Get("Events")
    pTH  = ROOT.TH1F(iName+"tmp",  iName+"tmp", iHist.GetNbinsX(),iHist.GetXaxis().GetXmin(),iHist.GetXaxis().GetXmax())
    lTree.Draw(iVar+">>"+iName+"tmp",iCut)
    iHist.Add(pTH)
    lFile.Close()


def roundTo(arr, base):
    for i in range(len(arr)):
        x = arr[i]
        if(x % base == base/2): new_x = x
        else: new_x = int(base * round(float(x)/base)) + base/2
        arr[i] = new_x



def fill_hist(v, h, event_num = None):
    h.Sumw2()
    #print("%i events in the dataset " % len(v))
    if(event_num is None):
        for x in v: h.Fill(x)
    else:
        e_dict = dict()
        for idx,x in enumerate(v):
            e_num = event_num[idx]
            e_dict[e_num] = (x, e_dict.get(e_num, (0,0))[1] + 1)

        print("%i unique events"% len(e_dict.keys()))
        for k in e_dict.keys():
            x = e_dict[k][0]
            n = e_dict[k][1]
            h.Fill(x,n)

        print("%.0f filled events" % h.Integral())
    #h.Print("range")


def get_mjj_max(h_file):
    with h5py.File(h_file, "r") as f:
        mjj = np.array(f['mjj'][()])
        return np.amax(mjj)


def load_h5_sb(h_file, hist, correctStats=False, sb1_edge = -1., sb2_edge = -1.):
    event_num = None
    with h5py.File(h_file, "r") as f:
        mjj = np.array(f['mjj'][()])
        if(correctStats):
            event_num = f['event_num'][()]

    fill_hist(mjj, hist, event_num)

def load_h5_bkg(h_file, hist, correctStats = False):
    event_num = None
    with h5py.File(h_file, "r") as f:
        mjj = f['mjj'][()]
        is_sig = f['truth_label'][()]
        if(correctStats):
            event_num = f['event_num'][()]

    mask = (is_sig < 0.1)
    if(correctStats): event_num = event_num[mask]
    fill_hist(mjj[mask], hist, event_num)


def load_h5_sig(h_file, hist, sig_mjj, requireWindow = False, correctStats =False, mixed = False):
    event_num = None
    with h5py.File(h_file, "r") as f:
        try:
            mjj = f['jet_kinematics'][:, 0]
        except:
            mjj = f['mjj'][()]

        num_evts = mjj.shape[0]
        if(mixed):
            is_sig = f['truth_label'][()].flatten()
        else: 
            is_sig = np.ones_like(mjj)


        if(is_sig.shape[0] != mjj.shape[0]):
            #fix bug in old h5 maker where is_sig array would be too long
            is_sig = is_sig[:num_evts]
        
        if(correctStats):
            event_num = f['event_num'][()]


    if(requireWindow): mask = (mjj > 0.8*sig_mjj) & (mjj < 1.2*sig_mjj) & (is_sig > 0.9)
    else: mask = mjj > 0.
    if(correctStats): event_num = event_num[mask]
    fill_hist(mjj[mask], hist, event_num)

def check_rough_sig(h_file, m_low, m_high):
    with h5py.File(h_file, "r") as f:
        mjj = f['mjj'][()]
        is_sig = f['truth_label'][()]

    eps = 1e-6
    in_window = (mjj > m_low) & (mjj < m_high)
    sig_events = is_sig > 0.9
    bkg_events = is_sig < 0.1
    S = max(mjj[sig_events & in_window].shape[0], eps)
    B = max(mjj[bkg_events & in_window].shape[0], eps)
    print("Mjj window %f to %f " % (m_low, m_high))
    print("S = %i, B = %i, S/B %f, sigificance ~ %.1f " % (S, B, float(S)/B, S/np.sqrt(B)))

def get_below_bins(h, min_count = 5):
    out = []
    #skip first bin (can't merge left)
    for i in range(2, h.GetNbinsX()+1):
        c = h.GetBinContent(i)

        #width per 100 gev
        #width = h.GetBinWidth(i) / 100.
        #density = c / width
        #print(i, c, density) 

        #remove left edge of bin if below thresh
        if( c < min_count ): out.append(i-1)
    
    return out

def get_rebinning(binsx, histos_sb, min_count = 5):
    h_rebin = histos_sb.Clone("h_rebin_temp")
    rebins = copy.deepcopy(binsx)
    below_min = True

    while(below_min):
        h_rebin = h_rebin.Rebin(len(rebins)-1, "", array('d', rebins))
        below_bins = get_below_bins(h_rebin, min_count = min_count)
        if(len(below_bins) < 1): below_min = False
        else:
            rebins.pop(below_bins[-1])

    #print("new bins:", rebins)
    #print("Done with rebinning!")
    #h_rebin.Print("all")

    return rebins





    
def checkSBFit(filename,label,roobins,plotname, nPars, plot_dir, plot_label = ""):
    
    fin = ROOT.TFile.Open(filename,'READ')
    workspace = fin.w

    sig_name = 'shapeSig_model_signal_mjj_JJ_%s' % label
    model_tot = workspace.pdf('model_s')
    model_qcd = workspace.pdf('model_b')
    model_sig = workspace.pdf(sig_name)
    var = workspace.var('mjj')
    data = workspace.data('data_obs')
    data.Print("V")



    #sig_norm_var_total = workspace.function('n_exp_binJJ_%s_proc_model_signal_mjj' % label).getVal()
    #var.setRange("int_range", 1500., 6500)
    #set_ = ROOT.RooArgSet(var)
    #sig_norm_var_denom = model_sig.createIntegral(set_, set_, "int_range").getVal()

    #sig_norm_var_coeff = sig_norm_var_total/ sig_norm_var_denom
    #print(sig_norm_var_coeff, sig_norm_var_denom )
    #sig_norm_var_coeff = 1e-4

    #sig_norm_var = ROOT.RooRealVar("sig_norm_var", "Signal Normalization", sig_norm_var_coeff, sig_norm_var_coeff/100., sig_norm_var_coeff*100.)
    #sig_norm_var.Print("v")

    #sig_norm_pdf = ROOT.RooAddPdf("sig_norm", "Signal", ROOT.RooArgList(model_sig), ROOT.RooArgList(sig_norm_var))
    #sig_norm_pdf.Print("v")


    model = model_tot

    model_tot.Print("v")

    #default roofit normalization is total range divided by number of bins, we want per 100 GeV


    #rescale so pdfs are in evts per 100 GeV
    low = roobins.lowBound()
    high = roobins.highBound()
    n = roobins.numBoundaries() - 1

    #RootFit default normalization is full range divided by number of bins
    default_norm = (high - low)/ n

    rescale = 100./ default_norm

    fit_norm = ROOT.RooFit.Normalization(rescale,ROOT.RooAbsReal.Relative)


    
    fres = model.fitTo(data,ROOT.RooFit.SumW2Error(1),ROOT.RooFit.Minos(0),ROOT.RooFit.Verbose(0),ROOT.RooFit.Save(1),ROOT.RooFit.NumCPU(8)) 
    #fres.Print()
    
    frame = var.frame()
    
    data.plotOn(frame, ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson), ROOT.RooFit.Binning(roobins),ROOT.RooFit.Name("data_obs"),ROOT.RooFit.Invisible(), 
            ROOT.RooFit.Rescale(rescale))
    model.getPdf('JJ_%s'%label).plotOn(frame,ROOT.RooFit.VisualizeError(fres,1),ROOT.RooFit.FillColor(ROOT.kRed-7),ROOT.RooFit.LineColor(ROOT.kRed-7),ROOT.RooFit.Name(fres.GetName()),
            fit_norm)
    model.getPdf('JJ_%s'%label).plotOn(frame,ROOT.RooFit.LineColor(ROOT.kRed+1),ROOT.RooFit.Name("model_s"), fit_norm)

    #model_qcd.plotOn(frame,ROOT.RooFit.VisualizeError(fres,1),ROOT.RooFit.FillColor(ROOT.kGreen-7),ROOT.RooFit.LineColor(ROOT.kGreen-7), ROOT.RooFit.Name("Background"))
    #model_qcd.plotOn(frame,ROOT.RooFit.LineColor(ROOT.kRed+1),ROOT.RooFit.Name("Background"))
    #sig_norm_pdf.plotOn(frame,ROOT.RooFit.VisualizeError(fres,1),ROOT.RooFit.FillColor(ROOT.kBlue-7),ROOT.RooFit.LineColor(ROOT.kBlue-7), ROOT.RooFit.Name("Signal"))
    #sig_norm_pdf.plotOn(frame,ROOT.RooFit.LineColor(ROOT.kBlue+1), ROOT.RooFit.Name("Signal"))


    useBinAverage = True
    hpull = frame.pullHist("data_obs","model_s",useBinAverage)
    hresid = frame.residHist("data_obs", "model_s", False, useBinAverage)
    dhist = ROOT.RooHist(frame.findObject('data_obs', ROOT.RooHist.Class()))

    #print_roohist(dhist)


    chi2, ndof = calculateChi2(hpull, nPars + 1, excludeZeros = True, dataHist = dhist)

    #replot data on top 
    data.plotOn(frame,ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2), ROOT.RooFit.Binning(roobins),ROOT.RooFit.Name("data_obs"),ROOT.RooFit.XErrorSize(0), 
            ROOT.RooFit.Rescale(rescale))

    frame3 = var.frame()
    #average bin edges instead of bin center
    frame3.addPlotable(hpull,"X0 P E1")
    

    #chi2,ndof = calculateChi2(hpull, nPars +1)

    pdf_names = ["model_s"] 
    PlotFitResults(frame,fres.GetName(),nPars+1,frame3,"data_obs", pdf_names,chi2,ndof,'sbFit_'+plotname, plot_dir, has_sig = True, plot_label = plot_label)

    print "chi2,ndof are", chi2, ndof
    return chi2, ndof


def f_test(nParams, nDof, chi2, fit_errs, thresh = 0.05, err_thresh = 0.5):
    #assumes arrays are in increasing number of params order (ie nParams[0] is minimum number of params)
    print  "\n\n #################### STARTING F TEST #######################" 
    best_i = 0
    for i in range(1, len(nParams)):
        print("F test comparing %i to %i params" % (nParams[best_i], nParams[i]))
        nDof_base = nDof[best_i]
        chi2_base = chi2[best_i]

        nDof_new = nDof[i]
        chi2_new = chi2[i]

        F_num =   max((chi2_base - chi2_new), 0)/(abs(nDof_new - nDof_base))
        F_denom = chi2_new/nDof_new 
        F = F_num / F_denom

        prob = 1. - ROOT.TMath.FDistI(F, abs(nDof_new - nDof_base), nDof_new)

        print("Base chi2 was %.1f, new is %.1f" % (chi2_base, chi2_new))
        print("F is %.2f, prob is %.3f" % (F, prob))

        if(prob < thresh ):
            if(fit_errs[i] <  err_thresh or fit_errs[i] < fit_errs[best_i]):
                print("Prob below threshold, switching to %i parameters" % nParams[i])
                best_i = i
            else:
                print("Prob below threshold, but largest param error is too large(%.2f) so NOT adding parameters" % fit_errs[i])

        elif(fit_errs[best_i]  > err_thresh and fit_errs[i] < err_thresh):
                print("Prob not below threshold but previous best was above error threshold, so switch to %i params" % nParams[i])
                best_i = i



    return best_i
