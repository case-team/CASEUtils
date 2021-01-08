import ROOT
ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
import json
import sys
from array import array
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

class DataCardMaker:
    def __init__(self,tag):
        self.systematics=[]
	self.tag="JJ"+"_"+tag
        self.rootFile = ROOT.TFile("datacardInputs_%s.root"%self.tag,"RECREATE")
        self.rootFile.cd()
        self.w=ROOT.RooWorkspace("w","w")
	self.luminosity = 1.0
	self.contributions=[]
	self.systematics=[]

    def delete(self):
     if self.w:
      self.w.Delete()
      self.rootFile.Close()
      self.rootFile.Delete()
      
    def makeCard(self):

        f = open("datacard_"+self.tag+'.txt','w')
        f.write('imax 1\n')
        f.write('jmax {n}\n'.format(n=len(self.contributions)-1))
        f.write('kmax *\n')
        f.write('-------------------------\n')
        for c in self.contributions:
            f.write('shapes {name} {channel} {file}.root w:{pdf}\n'.format(name=c['name'],channel=self.tag,file="datacardInputs_"+self.tag,pdf=c['pdf']))
        f.write('shapes {name} {channel} {file}.root w:{name}\n'.format(name="data_obs",channel=self.tag,file="datacardInputs_"+self.tag))
        f.write('-------------------------\n')
        f.write('bin '+self.tag+'\n')
        f.write('observation  -1\n')
        f.write('-------------------------\n')
        f.write('bin\t') 

        for shape in self.contributions: f.write(self.tag+'\t')
        f.write('\n')

        #Sort the shapes by ID 
 
        shapes = sorted(self.contributions,key=lambda x: x['ID'])
        #print names
        f.write('process\t')
        for shape in shapes:
            f.write(shape['name']+'\t')
        f.write('\n')

        #Print ID
        f.write('process\t')
        for shape in shapes:
            f.write(str(shape['ID'])+'\t')
        f.write('\n')

        #print rates
        f.write('rate\t')
        for shape in shapes:
            f.write(str(shape['yield'])+'\t')
        f.write('\n')


        #Now systematics
        for syst in self.systematics:
            if syst['kind'] == 'param':
                f.write(syst['name']+'\t'+'param\t' +str(syst['values'][0])+'\t'+str(syst['values'][1])+'\n')
            elif syst['kind'] == 'flatParam':
                f.write(syst['name']+'\t'+'flatParam\n')
            elif 'rateParam' in syst['kind']:
                line = syst['name']+'\t'+str(syst['kind'])+"\t"+str(syst['bin'])+"\t"+str(syst['process'])+"\t"+str(syst['values'])+"\t"+str(syst['variables'])
                line+='\n' 
                f.write(line)
                
            elif syst['kind'] == 'discrete':
                f.write(syst['name']+'\t'+'discrete\n')

            elif syst['kind'] == 'lnN' :
                f.write(syst['name']+'\t'+ 'lnN\t' )
                for shape in shapes:
                    has=False
                    for name,v in syst['values'].iteritems():
                        if shape['name']==name:
                            f.write(str(v)+'\t' )
                            has=True
                            break;
                    if not has:
                            f.write('-\t' )
                f.write('\n' )
            elif syst['kind'] == 'lnU': 
                f.write(syst['name']+'\t'+ 'lnU\t' )
                for shape in shapes:
                    has=False
                    for name,v in syst['values'].iteritems():
                        if shape['name']==name:
                            f.write(str(v)+'\t' )
                            has=True
                            break;
                    if not has:
                            f.write('-\t' )
                f.write('\n' )
                            
                        
        f.close()


        self.rootFile.cd()
        self.w.Write("w",0,2000000000)
        self.rootFile.Close()
	
    def importBinnedData(self,filename,histoname,poi,name = "data_obs",scale=1):
        f=ROOT.TFile(filename)
        histogram=f.Get(histoname)
        histogram.Scale(scale)
        cList = ROOT.RooArgList()
        for i,p in enumerate(poi):
            cList.add(self.w.var(p))
            if i==0:
                axis=histogram.GetXaxis()
            elif i==1:
                axis=histogram.GetYaxis()
            elif i==2:
                axis=histogram.GetZaxis()
            else:
                print 'Asking for more than 3 D . ROOT doesnt support that, use unbinned data instead'
                return
            mini=axis.GetXmin()
            maxi=axis.GetXmax()
            bins=axis.GetNbins()
            self.w.var(p).setMin(mini)
            self.w.var(p).setMax(maxi)
            self.w.var(p).setBins(bins)
            self.w.var(p).setBins(bins,"cache")
	    mjj=self.w.var(p)
        dataHist=ROOT.RooDataHist(name,name,cList,histogram)
	#dataHist=ROOT.RooDataHist(name, name, ROOT.RooArgList(mjj), ROOT.RooFit.Import(histogram)) 
        getattr(self.w,'import')(dataHist,ROOT.RooFit.RenameVariable(name,name))
	
    def addSystematic(self,name,kind,values,bin="",process="",variables="",addPar = ""):
        if kind != 'rateParam': self.systematics.append({'name':name,'kind':kind,'values':values })
        else: self.systematics.append({'name':name,'kind':kind,'bin':bin,'process':process,'values':values,'variables':variables})
	
    def addFixedYieldFromFile(self,name,ID,filename,histoName,constant=1.0):
        pdfName="_".join([name,self.tag])
        f=ROOT.TFile(filename)
        histogram=f.Get(histoName)
        events=histogram.Integral()*self.luminosity*constant
        self.contributions.append({'name':name,'pdf':pdfName,'ID':ID,'yield':events})

    def addFloatingYield(self,name,ID,filename,histoName,mini=0,maxi=1e+9,constant=False):
        pdfName="_".join([name,self.tag])
        pdfNorm="_".join([name,self.tag,"norm"])
        f=ROOT.TFile(filename)
        histogram=f.Get(histoName)
        events=histogram.Integral()
        self.w.factory("{name}[{val},{mini},{maxi}]".format(name=pdfNorm,val=events,mini=mini,maxi=maxi))       
        if constant:
            self.w.var(pdfNorm).setConstant(1)
        self.contributions.append({'name':name,'pdf':pdfName,'ID':ID,'yield':1.0})
		
    def addSignalShape(self,name,variable,jsonFile,scale ={},resolution={}):
    
        pdfName="_".join([name,self.tag])
	
        #self.w.factory("MH[3000]")
        #self.w.var("MH").setConstant(1)
       
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
       
        self.w.factory(variable+"[0,13000]")
	
	f = ROOT.TFile(jsonFile,'READ')
	meanG = f.Get('mean')
	sigmaG = f.Get('sigma')
	alphaG = f.Get('alpha')
	scalesigmaG = f.Get('scalesigma')
	sigfracG = f.Get('sigfrac')
	signG = f.Get('sign')
	
	x = ROOT.Double(0.)
	mean = ROOT.Double(0.)
	meanG.GetPoint(0,x,mean)
	sigma = ROOT.Double(0.)
	sigmaG.GetPoint(0,x,sigma)	
	alpha = ROOT.Double(0.)
	alphaG.GetPoint(0,x,alpha)
	scalesigma = ROOT.Double(0.)
	scalesigmaG.GetPoint(0,x,scalesigma)
	sigfrac = ROOT.Double(0.)
	sigfracG.GetPoint(0,x,sigfrac)
	sign = ROOT.Double(0.)
	signG.GetPoint(0,x,sign)
		
	
	meanVar = "_".join(["MEAN",name,self.tag])
        self.w.factory("expr::{name}('{param}*(1+{vv_syst})',{vv_systs},{param})".format(name=meanVar,param=mean,vv_syst=scaleStr,vv_systs=','.join(scaleSysts)))

        sigmaVar = "_".join(["SIGMA",name,self.tag])
	self.w.factory("expr::{name}('{param}*(1+{vv_syst})',{vv_systs},{param})".format(name=sigmaVar,param=sigma,vv_syst=resolutionStr,vv_systs=','.join(resolutionSysts)))
		
	alphaVar = "_".join(["ALPHA",name,self.tag])		
        alpha = ROOT.RooRealVar(alphaVar,alphaVar,alpha)
	getattr(self.w,'import')(alpha,ROOT.RooFit.Rename(alphaVar))
	
	sigfracVar = "_".join(["SIGFRAC",name,self.tag])
	sigfrac = ROOT.RooRealVar(sigfracVar,sigfracVar,sigfrac)
	getattr(self.w,'import')(sigfrac,ROOT.RooFit.Rename(sigfracVar))
	
	scalesigmaVar = "_".join(["SCALESIGMA",name,self.tag])
	scalesigma = ROOT.RooRealVar(scalesigmaVar,scalesigmaVar,scalesigma)
	getattr(self.w,'import')(scalesigma,ROOT.RooFit.Rename(scalesigmaVar))
	
	signVar = "_".join(["SIGN",name,self.tag])
	sign = ROOT.RooRealVar(signVar,signVar,sign)	
	getattr(self.w,'import')(sign,ROOT.RooFit.Rename(signVar))

        gsigmaVar = "_".join(["GSIGMA",name,self.tag])		
        gsigma = ROOT.RooFormulaVar(gsigmaVar,"@0*@1", ROOT.RooArgList(self.w.function(sigmaVar),scalesigma))
        #getattr(self.w,'import')(gsigma,ROOT.RooFit.Rename(gsigmaVar))      

        gaussFunc = "_".join(["gauss",name,self.tag])	
        gauss = ROOT.RooGaussian(gaussFunc, gaussFunc, self.w.var(variable), self.w.function(meanVar), gsigma)
	cbFunc = "_".join(["cb",name,self.tag])
        cb    = ROOT.RooCBShape(cbFunc, cbFunc,self.w.var(variable), self.w.function(meanVar), self.w.function(sigmaVar), alpha, sign)
        model = ROOT.RooAddPdf(pdfName, pdfName, gauss, cb, self.w.var(sigfracVar))	
        getattr(self.w,'import')(model,ROOT.RooFit.Rename(pdfName))

    def addQCDShape(self,name,variable,preconstrains,nPars=4):

        pdfName=name+"_"+self.tag
	
        MVV=variable
        if self.w.var(MVV) == None: self.w.factory(MVV+"[0,10000]")
	
	f = ROOT.TFile.Open(preconstrains,'READ')
	parsG = [f.Get('p%i'%i) for i in range(1,nPars+1)]
	pars_val = [ROOT.Double(0.) for i in range(0,nPars)]       
        for i in range(1,nPars+1):
	 x = ROOT.Double(0.)
	 parsG[i-1].GetPoint(0,x,pars_val[i-1])
	 pName="_".join(["CMS_JJ_p%i"%i,self.tag])
	 errUp=pars_val[i-1]+parsG[i-1].GetErrorYhigh(0)*100.
	 errDown=pars_val[i-1]-parsG[i-1].GetErrorYlow(0)*100.
	 print i,pName,pars_val[i-1],parsG[i-1].GetErrorYhigh(0),parsG[i-1].GetErrorYlow(0),errUp,errDown
         self.w.factory("{name}[{val},{errDown},{errUp}]".format(name=pName,val=pars_val[i-1],errUp=errUp,errDown=errDown))
	
	if nPars==2: model = ROOT.RooGenericPdf(pdfName, "pow(1-@0/13000., @1)/pow(@0/13000., @2)", ROOT.RooArgList(self.w.var(MVV), self.w.var("CMS_JJ_p1_%s"%self.tag), self.w.var("CMS_JJ_p2_%s"%self.tag)))	 
	elif nPars==3: model = ROOT.RooGenericPdf(pdfName, "pow(1-@0/13000., @1)/pow(@0/13000., @2+@3*log(@0/13000.))", ROOT.RooArgList(self.w.var(MVV), self.w.var("CMS_JJ_p1_%s"%self.tag), self.w.var("CMS_JJ_p2_%s"%self.tag), self.w.var("CMS_JJ_p3_%s"%self.tag)))
	elif nPars==4: model = ROOT.RooGenericPdf(pdfName, "pow(1-@0/13000., @1)/ ( pow(@0/13000., @2+@3*log(@0/13000.)+@4*pow(log(@0/13000.),2)) )", ROOT.RooArgList(self.w.var(MVV), self.w.var("CMS_JJ_p1_%s"%self.tag), self.w.var("CMS_JJ_p2_%s"%self.tag), self.w.var("CMS_JJ_p3_%s"%self.tag), self.w.var("CMS_JJ_p4_%s"%self.tag)))

        getattr(self.w,'import')(model,ROOT.RooFit.Rename(pdfName))

    def addQCDShapeNoTag(self,name,variable,preconstrains,nPars=4):

        pdfName=name+"_"+self.tag
	
        MVV=variable
        if self.w.var(MVV) == None: self.w.factory(MVV+"[0,10000]")
	
	f = ROOT.TFile.Open(preconstrains,'READ')
	parsG = [f.Get('p%i'%i) for i in range(1,nPars+1)]
	pars_val = [ROOT.Double(0.) for i in range(0,nPars)]       
        for i in range(1,nPars+1):
	 x = ROOT.Double(0.)
	 parsG[i-1].GetPoint(0,x,pars_val[i-1])
	 pName="CMS_JJ_p%i"%i
	 errUp=pars_val[i-1]+parsG[i-1].GetErrorYhigh(0)*100.
	 errDown=pars_val[i-1]-parsG[i-1].GetErrorYlow(0)*100.
	 print i,pName,pars_val[i-1],parsG[i-1].GetErrorYhigh(0),parsG[i-1].GetErrorYlow(0),errUp,errDown
         self.w.factory("{name}[{val},{errDown},{errUp}]".format(name=pName,val=pars_val[i-1],errUp=errUp,errDown=errDown))
	
	if nPars==2: model = ROOT.RooGenericPdf(pdfName, "pow(1-@0/13000., @1)/pow(@0/13000., @2)", ROOT.RooArgList(self.w.var(MVV), self.w.var("CMS_JJ_p1"), self.w.var("CMS_JJ_p2")))	 
	elif nPars==3: model = ROOT.RooGenericPdf(pdfName, "pow(1-@0/13000., @1)/pow(@0/13000., @2+@3*log(@0/13000.))", ROOT.RooArgList(self.w.var(MVV), self.w.var("CMS_JJ_p1"), self.w.var("CMS_JJ_p2"), self.w.var("CMS_JJ_p3")))
	elif nPars==4: model = ROOT.RooGenericPdf(pdfName, "pow(1-@0/13000., @1)/ ( pow(@0/13000., @2+@3*log(@0/13000.)+@4*pow(log(@0/13000.),2)) )", ROOT.RooArgList(self.w.var(MVV), self.w.var("CMS_JJ_p1"), self.w.var("CMS_JJ_p2"), self.w.var("CMS_JJ_p3"), self.w.var("CMS_JJ_p4")))

        getattr(self.w,'import')(model,ROOT.RooFit.Rename(pdfName))
