from HiggsAnalysis.CombinedLimit.PhysicsModel import *
import ROOT
from array import array

class Function:
    def __init__(self, x, y):
        f = ROOT.TSpline3("spline", x, y, len(x))
        self.f = f

    def __call__(self, x, par=None):
        x = x[0]
        #y = self.f.Eval(x)
        y = 0.5
        return y




class VariableSigEff(PhysicsModel):
    def __init__(self):
        super(VariableSigEff, self).__init__()
        return

    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""
        maxR = 100

        self.modelBuilder.doVar("r[1,-%i, %i]" % (maxR, maxR));
        self.modelBuilder.doSet("POI","r")

        r_var = self.modelBuilder.out.var('r')

        #TODO Dynamically load these or something
        #assume flat below zero and after max range
        #xs = array('d', [-maxR, 0, 0.5, 1.0, 1.5, 2.0, maxR])
        #effs = array('d', [0.01, 0.01, 0.01, 0.05, 0.2, 0.25, 0.3, 0.3, 0.3])

        #xs = array('d', [-maxR, maxR])
        #effs = array('d', [0.5, 0.5])

        #RooSpline only in very recent root versions so do this work around with TF1

        #f = ROOT.TF1("tf1_sig_eff_spline", Function(xs, effs), min(xs), max(xs), 0)

        #Using binded TF1
        f = ROOT.TF1("tf1_sig_eff", "0.5 * x", -2*maxR, 2*maxR)
        func = ROOT.RooFit.bindFunction("norm", f, r_var)

        #Using native RooFit Function
        #c = ROOT.RooRealVar("c", "c", 0.5)
        #func = ROOT.RooPolyVar("norm", "norm",  r_var, ROOT.RooArgList(c), 1)

        getattr(self.modelBuilder.out, 'import')(func)
        self.modelBuilder.out.ls()

    def getYieldScale(self,bin,process):
        if 'signal' in process: return "norm"
        #if 'signal' in process: return "halfR"
        else: return 1




sig_eff_model = VariableSigEff()
