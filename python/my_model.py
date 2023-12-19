from HiggsAnalysis.CombinedLimit.PhysicsModel import *
import matplotlib.pyplot as plt
import ROOT
from array import array
import h5py
import os,sys
import numpy as np

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
    def __init__(self, map_file = ""):
        super(VariableSigEff, self).__init__()
        self.map_file = map_file
        return

    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""
        #maximum cross section
        maxR = 100

        #self.modelBuilder.doVar("xs[1,-%i, %i]" % (maxR, maxR));
        self.modelBuilder.doVar("r[1,-%i, %i]" % (-maxR, maxR));
        self.modelBuilder.doSet("POI","r")

        r_var = self.modelBuilder.out.var('r')

        if(not os.path.exists(self.map_file)):
            print("Signal Efficiency map %s is missing! Exiting" % self.map_file) 

        f = h5py.File(self.map_file)
        xs = f['injected_xsecs'][()]
        effs = f['effs'][()]


        const = 100
        norm_factor = f['norm_factor'][0] #lumi * preselection eff
        norm_factor /= const #nominal yield is 100, so scale by this amount

        floor_effs = False
        do_plot = True

        if(floor_effs):
            #add artificial points to the eff map to keep interpolation flat in between different values
            effs_new = []
            xs_new = []
            for i in range(len(effs) - 1 ):
                shift = 0.1
                xs_new.append(xs[i])
                effs_new.append(effs[i])

                #add ficticious point right before next point
                xs_new.append(xs[i+1]* (1. - shift))
                effs_new.append(effs[i])
            xs = xs_new
            effs = effs_new
            

        if(do_plot):
            plt.figure(figsize=(15,10))
            plt.plot(xs, effs,  markersize = 15, c = 'blue', marker='o', linestyle = 'solid')

            plt.xlim([0, None])
            plt.ylim([0, None])
            plt.xlabel(r"Cross Section (fb)")
            plt.ylabel("Signal Efficiency")
            plt.savefig("sig_eff_map.png" , bbox_inches="tight")


        #assume efficiencies flat before first and after last injection
        xs = np.insert(xs, 0, -maxR)
        xs = np.append(xs, maxR)
        effs = np.insert(effs, 0, effs[0])
        effs = np.append(effs, effs[-1])

        #xs = np.arange(-maxR, maxR, 1)
        #effs = [0.2]*len(xs)



        print('LOADING cross sections: ')
        print(xs)
        print('LOADING effs: ')
        print(effs)

        yields = [norm_factor*xs[i]*effs[i] for i in range(len(xs))]
        xs_full = np.arange(-maxR, maxR, 1)
        interp_yields = np.interp(xs_full, xs, yields)

        #assume effiency flat below zero and after max injection
        #xs = array('d', [-maxR, 0, 0.5, 1.0, 1.5, 2.0, maxR])
        #effs = array('d', [0.01, 0.01, 0.01, 0.05, 0.2, 0.25, 0.3, 0.3, 0.3])

        xs_final = vectorize(xs_full)
        #yield is effs times cross section times 'norm factor'
        #yields_final = vectorize([norm_factor*xs[i]*effs[i] for i in range(len(xs))])
        yields_final = vectorize(interp_yields)

        BC = ROOT.RooNCSplineCore.bcClamped

        #RooSpline between different injection points
        #func = ROOT.RooNCSpline_1D_fast("norm", "norm", r_var, xs, yields, BC,BC)
        func = ROOT.RooNCSpline_1D_fast("norm", "norm", r_var, xs_final, yields_final)

        getattr(self.modelBuilder.out, 'import')(func)
        self.modelBuilder.out.ls()

    def getYieldScale(self,bin,process):
        if 'signal' in process: return "norm"
        else: return 1




sig_eff_model = VariableSigEff("/uscms_data/d3/oamram/CASE_analysis/src/CASE/CASEUtils/python/sig_eff_map.h5")
