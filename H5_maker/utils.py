import h5py
import numpy as np
import os


def get_file_list(dname):
    if('store' in dname and 'eos' not in dname):
        cmd = "xrdfs root://cmseos.fnal.gov/ ls -u %s > temp_flist.txt" % dname
    else:
        cmd = "find %s/*.root > temp_flist.txt" % dname
    os.system(cmd)
    print(cmd)
    f = open("temp_flist.txt")
    f_list = f.read().splitlines()
    f_list = list(filter(lambda x: ".root" in x, f_list))
    f.close()
    os.system("rm temp_flist.txt")
    return f_list


def append_h5(f, name, data):
    prev_size = f[name].shape[0]
    f[name].resize(( prev_size + data.shape[0]), axis=0)
    f[name][prev_size:] = data
