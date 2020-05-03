from H5_maker import *
import os
import sys

# merge multiple h5's together (like hadd)
# syntax is python H5_merge.py output_name.h5 input1.h5 input2.h5 input3.h5 ...

def merge(fin_name, fout_name):
    fin = h5py.File(fin_name, "r")
    fout = h5py.File(fout_name, "r+")

    fin_keys = fin.keys()
    fout_keys = fout.keys()

    if(fin_keys != fout_keys):
        print("Input and output files have different datasets!")
        print("fin %s: " % fin_name, fin_keys)
        print("fout %s: " %fout_name, fout_keys)
        print("skipping this dataset")
        return

    fin_keys.remove("preselection_eff")
    n_fin = float(fin[fin_keys[0]].shape[0])
    n_fout = float(fout[fin_keys[0]].shape[0])
    #preselection efficiency is weighted average of two files
    fout['preselection_eff'][0] = (n_fin * fin['preselection_eff'][0] + n_fout * fout['preselection_eff'][0]) / (n_fin + n_fout)
    for key in fin_keys:
        append_h5(fout, key, fin[key])




def merge_multiple(fout_name, fs):
    print("Merging H5 files: ", fs)
    print("Dest %s" % fout_name)
    os.system("cp %s %s" % (fs[0], fout_name))
    for fin_name in fs[1:]:
        merge(fin_name, fout_name)



#merge_multiple("test.h5", ["output_files/QCD_HT1000to1500_0.h5", "output_files/QCD_HT1000to1500_1.h5", "output_files/QCD_HT1000to1500_2.h5"])
if __name__ == "__main__":
    #print(sys.argv[1], sys.argv[2:])
    merge_multiple(sys.argv[1], sys.argv[2:])


