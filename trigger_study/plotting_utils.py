from ROOT import TLatex, TCanvas, TLine, TFile, RooBinning, TLegend

def setTopPad(TopPad, r=4):
    TopPad.SetPad("TopPad", "", 0., 1./r, 1.0, 1.0, 0, -1, 0)
    TopPad.SetTopMargin(0.28/r)
    TopPad.SetBottomMargin(2*0.04/r)
    TopPad.SetRightMargin(0.05)
    TopPad.SetTicks(1, 1)

def setBotPad(BotPad, r=4):
    BotPad.SetPad("BotPad", "", 0., 0., 1.0, 1./r, 0, -1, 0)
    BotPad.SetTopMargin(0.5*r/100.)
    BotPad.SetBottomMargin(r/10.)
    BotPad.SetRightMargin(0.05)
    BotPad.SetTicks(1, 1)

def setPadStyle(h, r=1.2, isTop=False):
    h.GetXaxis().SetTitleSize(h.GetXaxis().GetTitleSize()*r*r)
    h.GetYaxis().SetTitleSize(h.GetYaxis().GetTitleSize()*r)
    h.GetXaxis().SetLabelSize(h.GetXaxis().GetLabelSize()*r)
    h.GetYaxis().SetLabelSize(h.GetYaxis().GetLabelSize()*r)
    if isTop: h.GetXaxis().SetLabelOffset(0.04)

def setBotStyle(h, r=4, fixRange=True):
    h.GetXaxis().SetLabelSize(h.GetXaxis().GetLabelSize()*(r-1));
    h.GetXaxis().SetLabelOffset(h.GetXaxis().GetLabelOffset()*(r-1));
    h.GetXaxis().SetTitleSize(h.GetXaxis().GetTitleSize()*(r-1));
    h.GetYaxis().SetLabelSize(h.GetYaxis().GetLabelSize()*(r-1));
    h.GetYaxis().SetNdivisions(505);
    h.GetYaxis().SetTitleSize(h.GetYaxis().GetTitleSize()*(r-1));
    h.GetYaxis().SetTitleOffset(h.GetYaxis().GetTitleOffset()/(r-1));
    if fixRange:
        h.GetYaxis().SetRangeUser(0., 2.)
        for i in range(1, h.GetNbinsX()+1):
            if h.GetBinContent(i)<1.e-6:
                h.SetBinContent(i, -1.e-6)

def drawCMS(lumi, text, onTop=False, year='', suppressCMS=False, suppress_year=False, large=False):
    latex = TLatex()
    latex.SetNDC()
    if large:
        latex.SetTextSize(0.06)
        upper_margin = 0.98
    else:
        latex.SetTextSize(0.045)
        upper_margin = 0.99
    latex.SetTextColor(1)
    latex.SetTextFont(42)
    latex.SetTextAlign(33)
    if (type(lumi) is float or type(lumi) is int):
        if float(lumi) > 0:
            if float(lumi)>100000.:
                latex.DrawLatex(0.9, upper_margin+0.012, "%.0f fb^{-1}  (13 TeV)" % (float(lumi)/1000.))
            else:
                latex.DrawLatex(0.9, upper_margin+0.012, "%.1f fb^{-1}  (13 TeV)" % (float(lumi)/1000.))
        if year!='':
            if year=="run2": year=""
            latex.DrawLatex(0.24, upper_margin, year)
        elif float(lumi) > 0:
            if lumi==35920.:
                year = '2016'
            elif lumi==41530.:
                year = '2017'
            elif lumi==59740.:
                year = '2018'
            elif lumi==137190.:
                year = ''
            if not suppress_year:
                latex.DrawLatex(0.24, upper_margin, year)
        else:
            latex.DrawLatex(0.9, upper_margin, "(13 TeV)")
    elif type(lumi) is str: latex.DrawLatex(0.95, 0.985, "%s  (13 TeV)" % lumi)
    if not onTop: latex.SetTextAlign(11)
    latex.SetTextFont(62)
    if large:
        latex.SetTextSize(0.07 if len(text)>0 else 0.08)
    else:
        latex.SetTextSize(0.05 if len(text)>0 else 0.06)
    if not suppressCMS:
        if not onTop:
            if large:
                latex.DrawLatex(0.15, 0.87 if len(text)>0 else 0.84, "CMS")
            else:
                latex.DrawLatex(0.15, 0.88 if len(text)>0 else 0.85, "CMS")
        else: latex.DrawLatex(0.24, 0.9925, "CMS")
    if large:
        latex.SetTextSize(0.05)
    else:
        latex.SetTextSize(0.04)
    latex.SetTextFont(52)
    if not onTop: latex.DrawLatex(0.15, 0.84, text)
    else: latex.DrawLatex(0.45, 0.98, text)
