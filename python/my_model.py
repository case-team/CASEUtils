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


def vectorize(x):
    v = ROOT.vector('float')()
    for o in x:
        v.push_back(o)
    return v


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

        #TODO Dynamically load these values or something

        #assume effiency flat below zero and after max injection
        #xs = array('d', [-maxR, 0, 0.5, 1.0, 1.5, 2.0, maxR])
        #effs = array('d', [0.01, 0.01, 0.01, 0.05, 0.2, 0.25, 0.3, 0.3, 0.3])

        xs = [-maxR, maxR]
        effs = [0.5, 0.5]

        xs = vectorize(xs)
        
        #normalization is effs times cross section times 'norm factor' (TODO)
        norm_factor = 1.0
        ys = vectorize([norm_factor * effs[i] * xs[i] for i in range(len(effs))])


        #RooSpline between different injection points
        func = ROOT.RooNCSpline_1D_fast("norm", "norm", r_var, xs, ys)

        getattr(self.modelBuilder.out, 'import')(func)
        self.modelBuilder.out.ls()

    def getYieldScale(self,bin,process):
        if 'signal' in process: return "norm"
        else: return 1




sig_eff_model = VariableSigEff()
