import argparse
import json
import os
import ROOT
from Utils import roundTo
from dijetfit import fit_signalmodel


def fit_signals(options):

    out_dir = os.path.abspath(options.outDir)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    fine_bin_size = 4
    masses = options.masses

    binsx = [1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438,
             2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854,
             4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5500, 5663, 5877,
             6099, 6328, 6564, 6808]

    # round to smallest precision we are storing mass values with, otherwise
    # get weird effects related to bin size
    roundTo(binsx, fine_bin_size)
    bins_fine = int(binsx[-1]-binsx[0])/fine_bin_size

    if options.dcbModel:
        full_graphs = {'mean': ROOT.TGraphErrors(),
                       'sigma': ROOT.TGraphErrors(),
                       'alpha': ROOT.TGraphErrors(),
                       'alpha2': ROOT.TGraphErrors(),
                       'sign': ROOT.TGraphErrors(),
                       'sign2': ROOT.TGraphErrors()}
    else:
        full_graphs = {'mean': ROOT.TGraphErrors(),
                       'sigma': ROOT.TGraphErrors(),
                       'alpha': ROOT.TGraphErrors(),
                       'sign': ROOT.TGraphErrors(),
                       'scalesigma': ROOT.TGraphErrors(),
                       'sigfrac': ROOT.TGraphErrors()}

    for i, mass in enumerate(masses):
        print("########## FIT SIGNAL AND SAVE PARAMETERS ############")
        sig_file_name = os.path.join(out_dir, "sig_fit_{}.root".format(mass))
        current_fit = fit_signalmodel(options.inputFiles[i], sig_file_name,
                                      mass, binsx, bins_fine, out_dir + "/",
                                      return_fit=True,
                                      dcb_model=options.dcbModel, 
                                      fit_range = options.fitRange)

        parameters = dict()
        for var, graph in full_graphs.iteritems():
            val, err = current_fit.fetch(var)
            parameters[var] = val
            parameters['%s-err' % var] = err
            graph.SetPoint(i, mass, val)
            graph.SetPointError(i, 0.0, err)

        sig_file_name = os.path.join(out_dir, "sig_fit_{}.json".format(mass))
        with open(sig_file_name, 'w') as f:
            json.dump(parameters, f, indent=4)

    full_file_name = os.path.join(out_dir, "full_fit.root")
    full_outfile = ROOT.TFile(full_file_name, "RECREATE")
    full_outfile.cd()

    for name, graph in full_graphs.iteritems():
        graph.Write(name)

    full_outfile.Close()


def fitting_options():
    parser = argparse.ArgumentParser()
    parser.add_argument("-M", "-M", dest="masses", type=int, nargs="+",
                        help="List of injected signal masses")
    parser.add_argument("-i", "--inputFiles", dest="inputFiles", nargs="+",
                        help="""List of input h5 files.
                        Must have same length as mass list""")
    parser.add_argument("-o", "--outDir", dest="outDir", default='plots/',
                        help="Where to store output files")
    parser.add_argument("--dcb-model", "--dcbModel", dest="dcbModel", action="store_true",
                        default=False,
                        help="Whether or not to use double crystal ball model")
    parser.add_argument("--fitRange", dest="fitRange", type = float, default=0.2,
                        help="What mass range to perform fit to signal shape over (in terms of frac. of signal mass)")
    return parser


if __name__ == "__main__":
    parser = fitting_options()
    options = parser.parse_args()
    fit_signals(options)
