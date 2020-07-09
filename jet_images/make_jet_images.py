from __future__ import print_function, division
from optparse import OptionParser

import h5py
from ImageUtils import *
import numpy as np
import math

parser = OptionParser()
parser.add_option("--npix", default =32, type = int, help="Number of pixels")
parser.add_option("-i", "--input", dest = 'fin_name',  default = 'test.h5', help="Input file name")
parser.add_option("-o", "--output", dest='fout_name',  default = '', help="Output file name (leave blank for adding images to input file")
options, args = parser.parse_args()

# Load files
batch_size = 5000
fin_name = options.fin_name
fout_name = options.fout_name
if(fout_name == ""): fout_name = fin_name
excludes = []
npix = options.npix
img_width = 1.2
rotate = False
overwrite = True


if(fin_name != fout_name):
    #output to different file
    fin = h5py.File(fin_name, 'r')
    fout = h5py.File(fout_name, 'w')
    print("Going to copy all the data from %s to %s, will add jet images after " % (fin_name, fout_name))
    for key in fin.keys():
        if key in excludes:
            continue
        print("Copying key %s" % key)
        fin.copy(key, fout)

    print(fout.keys())

else:
    #add jet images to existing file
    print("going to add jet images to file %s" %fin_name)
    fin = h5py.File(fin_name, 'r+')
    if(overwrite):
        if('j1_images' in fin.keys()):
            print("deleting existing j1 images")
            del fin['j1_images']
        if('j2_images' in fin.keys()):
            print("deleting existing j2 images")
            del fin['j2_images']
    fout = fin

total_size = fin[fin.keys()[0]].shape[0]
iters = int(math.ceil(float(total_size)/batch_size))

print("going to make jet images for %i events with batch size %i (%i batches) \n \n" % (total_size, batch_size, iters))
for i in range(iters):

    print("batch %i \n" %i)
    start_idx = i*batch_size
    end_idx = min(total_size, (i+1)*batch_size)

    jet1_PFCands = fin['jet1_PFCands'][start_idx:end_idx]
    jet2_PFCands = fin['jet2_PFCands'][start_idx:end_idx]

    j1_4vec = fin['jet_kinematics'][start_idx:end_idx,2:6]
    j2_4vec = fin['jet_kinematics'][start_idx:end_idx,6:10]

    j1_images = np.zeros((end_idx - start_idx, npix, npix), dtype = np.float16)
    j2_images = np.zeros((end_idx - start_idx, npix, npix), dtype = np.float16)
    for j in range(end_idx - start_idx):

        j1_image = make_image(j1_4vec[j], jet1_PFCands[j], npix = npix, img_width = img_width, norm = True, rotate = rotate)
        j2_image = make_image(j2_4vec[j], jet2_PFCands[j], npix = npix, img_width = img_width, norm = True, rotate = rotate)


        j1_images[j] = j1_image
        j2_images[j] = j2_image


    if(i == 0):
        fout.create_dataset("j1_images", data=j1_images, chunks = True, maxshape=(None,npix,npix), compression = 'gzip')
        fout.create_dataset("j2_images", data=j2_images, chunks = True, maxshape=(None,npix,npix), compression = 'gzip')
    else:
        fout['j1_images'].resize((fout['j1_images'].shape[0] + j1_images.shape[0]), axis=0)
        fout['j1_images'][-j1_images.shape[0]:] = j1_images
        fout['j2_images'].resize((fout['j2_images'].shape[0] + j2_images.shape[0]), axis=0)
        fout['j2_images'][-j2_images.shape[0]:] = j2_images



fin.close()
if(fin_name != fout_name):
    fout.close()
print("Finished all batches! Output file saved to %s" %(fout_name))


            
                              

