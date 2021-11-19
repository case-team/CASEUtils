from __future__ import absolute_import
import utils
import h5py
import os
import sys
import copy
from os import listdir
from os.path import isfile, join

# merge multiple h5's together (like hadd)
# syntax is python H5_merge.py output_name.h5 input1.h5 input2.h5 input3.h5 ...

def merge(fin_name, fout_name):
    fin = h5py.File(fin_name, "r")
    fout = h5py.File(fout_name, "r+")

    fin_keys = list(fin.keys())
    fout_keys = list(fout.keys())

    if(fin_keys != fout_keys):
        print("Input and output files have different datasets!")
        print("fin %s: " % fin_name, fin_keys)
        print("fout %s: " %fout_name, fout_keys)
        print("skipping this dataset")
        return

    for key in copy.copy(fin_keys):
        if('_eff' in key):
            fin_keys.remove(key)
            n_fin = float(fin[fin_keys[0]].shape[0])
            n_fout = float(fout[fin_keys[0]].shape[0])
            #efficiency is weighted average of two files
            fout[key][0] = (n_fin * fin[key][0] + n_fout * fout[key][0]) / (n_fin + n_fout)


    for key in fin_keys:
        utils.append_h5(fout, key, fin[key])
        


def my_copy(fin_name, fout_name):
    fin = h5py.File(fin_name, "r")
    fout = h5py.File(fout_name, "w")

    fin_keys = fin.keys()

    for key in fin_keys:
        shape = list(fin[key].shape)
        shape[0] = None
        fout.create_dataset(key, data = fin[key], chunks = True, maxshape = shape, compression = 'gzip')



def merge_multiple(fout_name, fs):
    print("Merging H5 files: ", fs)
    print("Dest %s" % fout_name)
    os.system("cp %s %s" % (fs[0], fout_name))
    #my_copy(fs[0], fout_name)
    for fin_name in fs[1:]:
        print("Merging %s" % fin_name)
        merge(fin_name, fout_name)



#merge_multiple("test.h5", ["output_files/QCD_HT1000to1500_0.h5", "output_files/QCD_HT1000to1500_1.h5", "output_files/QCD_HT1000to1500_2.h5"])
if __name__ == "__main__":
    #print(sys.argv[1], sys.argv[2:])
    
    pt_bin = sys.argv[3]
    input_dir = sys.argv[2]
    
    onlyfiles = [join(input_dir, f) for f in listdir(input_dir) if (isfile(join(input_dir, f)) & (pt_bin in f))]
    merge_multiple(sys.argv[1], onlyfiles)
    
    print("Done!")

