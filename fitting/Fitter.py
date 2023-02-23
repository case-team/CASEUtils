import ROOT
import json
from numpy import random
from array import array
import os,sys,commands
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)

class Fitter(object):
    def __init__(self,poi = ['x'], debug = False):
        self.cache_name = "cache%i.root"%(random.randint(0, 1e+6))
        print("Making cache %s "  % self.cache_name)
        self.cache=ROOT.TFile(self.cache_name,"RECREATE")
        self.cache.cd()
        self.debug = debug
        self.cleanedup = False

        self.w=ROOT.RooWorkspace("w","w")
        self.dimensions = len(poi)
        self.poi=poi
        for v in poi:
            self.w.factory(v+"[1,161]")


    def __del__(self):
        if(not self.cleanedup): self.delete()

    def delete(self):
        if self.w:
            self.w.Delete()
        self.cache.Close()
        self.cache.Delete()
        os.system("rm %s" % self.cache_name)
        self.cleanedup = True

    def importBinnedData(self,histogram,poi = ["x"],name = "data", regions=[]):
        cList = ROOT.RooArgList()
        for i,p in enumerate(poi):
            var = self.w.var(p)
            cList.add(var)
            if i==0:
                axis=histogram.GetXaxis()
            elif i==1:
                axis=histogram.GetYaxis()
            elif i==2:
                axis=histogram.GetZaxis()
            else:
                print ('Asking for more than 3 D . ROOT doesnt support that, use unbinned data instead')
                return
            mini=axis.GetXmin()
            maxi=axis.GetXmax()
            bins=axis.GetNbins()
            binningx =[]
            for i in range(1,bins+2):
                #v = mmin + i * (mmax-mmin)/float(N)
                binningx.append(axis.GetBinLowEdge(i))
            if(len(regions) == 0):
                var.setMax(maxi)
                var.setMin(mini)
            else:
                for reg_name,reg_low,reg_high in regions:
                    var.setRange(reg_name, reg_low, reg_high)
            if(self.debug): 
                print (" set binning "+str(binningx)) 
            var.setBinning(ROOT.RooBinning(len(binningx)-1,array("d",binningx)))
            #a = self.w.var(p).getBinning()
            #for b in range(0,a.numBins()+1):
                #print a.binLow(b)
        dataHist=ROOT.RooDataHist(name,name,cList,histogram)
        getattr(self.w,'import')(dataHist,ROOT.RooFit.Rename(name))

    def fetch(self,var):
        #import pdb; pdb.set_trace()
        self.w.var(var).Print()
        print("Fetching value " ,self.w.var(var).getVal()  )
        print("Fetching error " ,self.w.var(var).getError())
        return (self.w.var(var).getVal(), self.w.var(var).getError())

    def getFunc(self,model = "model"):
        return self.w.pdf(model)

    def getData(self,data = "data"):
        return self.w.data(data)

    def getVar(self,var = "mjj"):
        return self.w.var(var)

    def getW(self):
        return self.w

    def fit(self,model = "model",data="data",options=[]):
        fit_data = self.w.data(data)

        if len(options)==0:
            fitresults = self.w.pdf(model).fitTo(fit_data)
        if len(options)==1:
            fitresults = self.w.pdf(model).fitTo(fit_data,options[0])	    
        if len(options)==2:
            fitresults = self.w.pdf(model).fitTo(fit_data,options[0],options[1])
        if len(options)==3:
            fitresults = self.w.pdf(model).fitTo(fit_data,options[0],options[1],options[2])
        if len(options)==4:
            fitresults = self.w.pdf(model).fitTo(fit_data,options[0],options[1],options[2],options[3])

        if fitresults:
            fitresults.Print() 
            f = ROOT.TFile.Open('fitresults.root','RECREATE')
            fitresults.SetName("fitresults")
            fitresults.Write()
            f.Close()	 
        return fitresults 

    def getLegend(self):
        self.legend = ROOT.TLegend(0.7510112,0.7183362,0.8502143,0.919833)
        self.legend.SetTextSize(0.032)
        self.legend.SetLineColor(0)
        self.legend.SetShadowColor(0)
        self.legend.SetLineStyle(1)
        self.legend.SetLineWidth(1)
        self.legend.SetFillColor(0)
        self.legend.SetFillStyle(0)
        self.legend.SetMargin(0.35)
        return self.legend

    def projection(self,model = "model",data="data",poi="x",filename="fit.root",binning=0,logy=False,xtitle='x',mass=1000):

        self.frame=self.w.var(poi).frame()

        f = ROOT.TFile.Open("fitresults.root",'READ')
        if f: 
            fr = f.Get('fitresults')
        else:
            fr = 0
            print("No fit result found (fitresults.root), plotting model only")

        if binning:
            self.w.data(data).plotOn(self.frame,ROOT.RooFit.Binning(binning),ROOT.RooFit.Invisible())
            if fr: 
                self.w.pdf(model).plotOn(self.frame,ROOT.RooFit.VisualizeError(fr,1),ROOT.RooFit.FillColor(ROOT.kRed-7),ROOT.RooFit.LineColor(ROOT.kRed-7),ROOT.RooFit.Name(fr.GetName()))
                self.w.pdf(model).plotOn(self.frame,ROOT.RooFit.LineColor(ROOT.kRed+1))	 
        else: 
            self.w.data(data).plotOn(self.frame,ROOT.RooFit.Invisible())
            if fr: 
                self.w.pdf(model).plotOn(self.frame,ROOT.RooFit.VisualizeError(fr,1),ROOT.RooFit.FillColor(ROOT.kRed-7),ROOT.RooFit.LineColor(ROOT.kRed-7),ROOT.RooFit.Name(fr.GetName()))
                self.w.pdf(model).plotOn(self.frame,ROOT.RooFit.LineColor(ROOT.kRed+1))

        if binning: self.w.data(data).plotOn(self.frame,ROOT.RooFit.Binning(binning))
        else: self.w.data(data).plotOn(self.frame)

        self.legend = self.getLegend()
        self.legend.AddEntry( self.w.pdf(model)," Full PDF","l")

        self.c=ROOT.TCanvas("c","c")
        if logy:
            self.frame.SetMinimum(1)
            self.frame.SetMaximum(1e+7)
            self.c.SetLogy()
        self.c.cd()
        self.frame.Draw()
        self.frame.GetYaxis().SetTitle('events')
        self.frame.GetXaxis().SetTitle(xtitle)
        self.frame.SetTitle('')
        self.c.Draw()

        self.legend.Draw("same")	    
        self.c.SaveAs(filename)
        pullDist = self.frame.pullHist()
        return self.frame.chiSquare()

    def signalResonance(self,name = 'model',poi="MVV",mass=0):

        ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")

        self.w.factory("MH[1000]")
        self.w.factory("mean[%.1f,%.1f,%.1f]"%(mass,0.8*mass,1.2*mass))
        self.w.factory("sigma[%.1f,%.1f,%.1f]"%(mass*0.05,mass*0.02,mass*0.10))
        self.w.factory("alpha[0.85,0.60,1.20]")
        self.w.factory("sign[6,0.1,150]")
        self.w.factory("scalesigma[2.0,1.2,3.6]")
        gsigma = ROOT.RooFormulaVar("gsigma","@0*@1", ROOT.RooArgList(self.w.var("sigma"),self.w.var("scalesigma")))
        getattr(self.w,'import')(gsigma,ROOT.RooFit.Rename('gsigma'))
        self.w.factory("Gaussian::gauss(%s,mean,gsigma)"%poi)
        self.w.factory("CBShape::cb(%s,mean,sigma,alpha,sign)"%poi)
        self.w.factory('SUM::'+name+'(sigfrac[0.0,0.0,0.850]*gauss,cb)')
    
    def signalDCB(self, name='model', poi="MVV", mass=0):
    
        ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
        
        self.w.factory("MH[1000]")
        self.w.factory("mean[%.1f,%.1f,%.1f]"%(mass,0.8*mass,1.2*mass))
        self.w.factory("sigma[%.1f,%.1f,%.1f]"%(mass*0.05,mass*0.02,mass*0.10))
        self.w.factory("alpha[1.2,0.0,18]")
        self.w.factory("alpha2[1.2,0.0,10]")
        self.w.factory("sign[5,0,600]")
        self.w.factory("sign2[5,0,50]")
        self.w.factory("DoubleCB::"+name+"(%s,mean,sigma,alpha,sign,alpha2,sign2)"%poi)

    def qcdShape(self,name = 'model',poi="MVV",nPars=2):

        if(type(poi) == str): poi = self.w.var(poi)

        ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")

        if nPars==2:
            self.w.factory("p1[9.28433e+00, -100. , 100.]")
            self.w.factory("p2[1.03641e+01, -200, 200]")
            model = ROOT.RooGenericPdf(name, "pow(1-@0/13000., @1)/pow(@0/13000., @2)", ROOT.RooArgList(poi, self.w.var("p1"), self.w.var("p2")))	 
        elif nPars==3: 
            self.w.factory("p1[9.28433e+00, -100. , 100.]")
            self.w.factory("p2[1.03641e+01, -200, 200]")
            self.w.factory("p3[2.35256e+00, -100., 100.]")
            model = ROOT.RooGenericPdf(name, "pow(1-@0/13000., @1)/pow(@0/13000., @2+@3*log(@0/13000.))", ROOT.RooArgList(poi, self.w.var("p1"), self.w.var("p2"), self.w.var("p3")))
        elif nPars==4:
            #default set
            self.w.factory("p1[9.28433e+00, -100. , 100.]")
            self.w.factory("p2[1.03641e+01, -200, 200]")	 
            self.w.factory("p3[2.35256e+00, -100., 100.]")
            self.w.factory("p4[4.17695e-01, -100., 100.]")
            #set-1
            #self.w.factory("p1[3.64825e-05, -1000. , 1000.]")
            #self.w.factory("p2[2.78348e+00, -100, 100]")	 
            #self.w.factory("p3[7.17321e+00, 0., 50.]")
            #self.w.factory("p4[-9.83737e-01, -1000., 1000.]")
            #set-2
            #self.w.factory("p1[3.45939e+02, -100. , 200.]")
            #self.w.factory("p2[-4.45967e+00, -100, 100]")	 
            #self.w.factory("p3[1.31413e+00, 0., 50.]")
            #self.w.factory("p4[1.91782e+02, -1000., 1000.]")
            #set-3
            #self.w.factory("p1[4.27857e+01, -100. , 1000.]")
            #self.w.factory("p2[-1.37580e-02, -100, 100]")	 
            #self.w.factory("p3[7.46099e+00, 0., 50.]")
            #self.w.factory("p4[1.78503e+00, -1000., 100.]")

            model = ROOT.RooGenericPdf(name, "pow(1-@0/13000., @1)/ ( pow(@0/13000., @2+@3*log(@0/13000.)+@4*pow(log(@0/13000.),2)) )", ROOT.RooArgList(poi, self.w.var("p1"), self.w.var("p2"), self.w.var("p3"), self.w.var("p4")))
            #alt func
            #model = ROOT.RooGenericPdf(name, "( @1*pow(1-@0/13000 + @4*pow(@0/13000,2),@2) ) / ( pow(@0/13000,@3) )", ROOT.RooArgList(self.w.var(poi), self.w.var("p1"), self.w.var("p2"), self.w.var("p3"), self.w.var("p4")))

        elif nPars ==5:
            self.w.factory("p1[9.28433e+00, -100. , 100.]")
            self.w.factory("p2[1.03641e+01, -200, 200]")	 
            self.w.factory("p3[2.35256e+00, -100., 100.]")
            self.w.factory("p4[4.17695e-01, -100., 100.]")
            self.w.factory("p5[1.00000e+01, -100., 100.]")

            model = ROOT.RooGenericPdf(name, "pow(exp(-@0/13000.),@5) *pow(1-@0/13000., @1)/ ( pow(@0/13000., @2+@3*log(@0/13000.)+@4*pow(log(@0/13000.),2)) )", 
                    ROOT.RooArgList(poi, self.w.var("p1"), self.w.var("p2"), self.w.var("p3"), self.w.var("p4"), self.w.var("p5")))

        elif nPars==6:

            self.w.factory("p1[6.02992e+00, -100. , 100]")
            self.w.factory("p2[6.28634e+00, -200., 200]")
            self.w.factory("p3[6.32552e-01, -100., 100.]")
            self.w.factory("p4[-8.53977e-02, -100., 100.]")
            self.w.factory("p5[1000, 1., 5000.]")
            self.w.factory("p6[1600, 0., 2500.]")
            model = ROOT.RooGenericPdf(name, "(0.5*tanh((@0-@6)/@5) + .5)*pow(1-@0/13000., @1) / ( pow(@0/13000., @2+@3*log(@0/13000.)+@4*pow(log(@0/13000.),2)) )", ROOT.RooArgList(poi, self.w.var("p1"), self.w.var("p2"), self.w.var("p3"), self.w.var("p4"), self.w.var("p5"), self.w.var("p6"))) 

        getattr(self.w,'import')(model,ROOT.RooFit.Rename(name))
        return model

    def altBkgShape(self,name = 'model',poi="MVV",ver= 1):

        if(type(poi) == str): poi = self.w.var(poi)


        if( ver == 1):
            self.w.factory("ap1[7.7, -100. , 100.]")
            self.w.factory("ap2[10, -200, 200]")	 
            self.w.factory("ap3[2.3, -100., 100.]")
            self.w.factory("ap4[0.4, -100., 100.]")

            model = ROOT.RooGenericPdf(name, "exp( -@1 * @0/13000.) / ( pow(@0/13000., @2+@3*log(@0/13000.)+@4*pow(log(@0/13000.),2)) )", 
                ROOT.RooArgList(poi, self.w.var("ap1"), self.w.var("ap2"), self.w.var("ap3"), self.w.var("ap4")))

        elif( ver == 2):
            self.w.factory("ap1[9.28433e+00, -100. , 100.]")
            self.w.factory("ap2[1.03641e+01, -200, 200]")	 
            self.w.factory("ap3[2.35256e+00, -100., 100.]")
            self.w.factory("ap4[4.17695e-01, -100., 100.]")

            model = ROOT.RooGenericPdf(name, " pow(@0/13000.,-@1) * exp( -@2 * @0/13000. - @3 * pow( @0/13000., 2) - @4 * pow(@0 / 13000., 3))", 
                ROOT.RooArgList(poi, self.w.var("ap1"), self.w.var("ap2"), self.w.var("ap3"), self.w.var("ap4")))

        getattr(self.w,'import')(model,ROOT.RooFit.Rename(name))
        return model

    def addDCBSignalShape(self, name, poi, jsonFile, scale={},
                          resolution={}):

        
        pdfName="_".join([name,])
        
        scaleStr='0'
        resolutionStr='0'

        scaleSysts=[]
        resolutionSysts=[]
        for syst,factor in scale.iteritems():
            self.w.factory(syst+"[0,-0.1,0.1]")
            scaleStr=scaleStr+"+{factor}*{syst}".format(factor=factor,syst=syst)
            scaleSysts.append(syst)
        for syst,factor in resolution.iteritems():
            self.w.factory(syst+"[0,-0.5,0.5]")
            resolutionStr=resolutionStr+"+{factor}*{syst}".format(factor=factor,syst=syst)
            resolutionSysts.append(syst)
            
        if(type(poi) == str): 
            poi = self.w.var(poi)
    
        f = ROOT.TFile(jsonFile,'READ')
        meanG = f.Get('mean')
        sigmaG = f.Get('sigma')
        alpha_oneG = f.Get('alpha')
        sign_oneG = f.Get('sign')
        alpha_twoG = f.Get('alpha2')
        sign_twoG = f.Get('sign2')

        x = ROOT.Double(0.)
        mean = ROOT.Double(0.)
        meanG.GetPoint(0,x,mean)
        sigma = ROOT.Double(0.)
        sigmaG.GetPoint(0,x,sigma)
        alpha_one = ROOT.Double(0.)
        alpha_oneG.GetPoint(0, x, alpha_one)
        alpha_two = ROOT.Double(0.)
        alpha_twoG.GetPoint(0, x, alpha_two)
        sign_one = ROOT.Double(0.)
        sign_oneG.GetPoint(0,x,sign_one)
        sign_two = ROOT.Double(0.)
        sign_twoG.GetPoint(0,x,sign_two)

        meanVar = "_".join(["MEAN", name, ])
        self.w.factory(
            "expr::{name}('{param}*(1+{vv_syst})',{vv_systs},{param})".format(
                name=meanVar, param=mean, vv_syst=scaleStr,
                vv_systs=','.join(scaleSysts)))

        sigmaVar = "_".join(["SIGMA", name, ])
        self.w.factory(
            "expr::{name}('{param}*(1+{vv_syst})',{vv_systs},{param})".format(
                name=sigmaVar, param=sigma, vv_syst=resolutionStr,
                vv_systs=','.join(resolutionSysts)))

        alphaOneVar = "_".join(["ALPHAONE", name, ])        
        alpha_one = ROOT.RooRealVar(alphaOneVar,alphaOneVar,alpha_one)
        getattr(self.w, 'import')(alpha_one, ROOT.RooFit.Rename(alphaOneVar))

        alphaTwoVar = "_".join(["ALPHATWO", name, ])        
        alpha_two = ROOT.RooRealVar(alphaTwoVar, alphaTwoVar, alpha_two)
        getattr(self.w, 'import')(alpha_two, ROOT.RooFit.Rename(alphaTwoVar))
        
        signOneVar = "_".join(["SIGNONE", name, ])
        sign_one = ROOT.RooRealVar(signOneVar, signOneVar, sign_one)    
        getattr(self.w, 'import')(sign_one, ROOT.RooFit.Rename(signOneVar))
        
        signTwoVar = "_".join(["SIGNTWO", name, ])
        sign_two = ROOT.RooRealVar(signTwoVar, signTwoVar, sign_two)    
        getattr(self.w, 'import')(sign_two, ROOT.RooFit.Rename(signTwoVar))
        
        #dcbFunc = "_".join(["dcb", name, ])
        #import pdb; pdb.set_trace()
        dcb = ROOT.RooDoubleCB(pdfName, pdfName, poi,
                                      self.w.function(meanVar),
                                      self.w.function(sigmaVar), alpha_one,
                                      sign_one, alpha_two, sign_two)

        getattr(self.w, 'import')(dcb, ROOT.RooFit.Rename(pdfName))

        return dcb, [alpha_one, alpha_two, sign_one, sign_two]


