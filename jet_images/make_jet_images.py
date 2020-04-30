from __future__ import print_function, division

import h5py
from ImageUtils import *
import numpy as np


# Load files
batch_size = 5000
fin_name = "../data/QCD_HT1000to1500_test.h5"
fout_name = fin_name
#fout_name = "../data/images_test.h5"
excludes = ['jet1_PFCands', 'jet2_PFCands']
npix = 40
img_width = 1.2


do_plots = False
has_sig_bit = False

if(fin_name != fout_name):
    #output to different file
    fin = h5py.File(fin_name, 'r')
    fout = h5py.File(fout_name, 'w')
    print("Going to copy all the data from %s to %s, will add jet images after ")
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
    fout = fin

total_size = fin[fin.keys()[0]].shape[0]
iters = total_size//batch_size




print("going to make jet images for %i events with batch size %i \n \n" % (total_size, batch_size))
for i in range(iters):

    print("batch %i \n" %i)
    start_idx = i*batch_size
    end_idx = min(total_size, (i+1)*batch_size)

    jet1_PFCands = fin['jet1_PFCands'][start_idx:end_idx]
    jet2_PFCands = fin['jet2_PFCands'][start_idx:end_idx]

    j1_4vec = fin['jet_kinematics'][start_idx:end_idx,2:6]
    j2_4vec = fin['jet_kinematics'][start_idx:end_idx,6:10]

    j1_images = np.zeros((batch_size, npix, npix), dtype = np.float16)
    j2_images = np.zeros((batch_size, npix, npix), dtype = np.float16)
    for j in range(end_idx - start_idx):

        j1_image = make_image(j1_4vec[j], jet1_PFCands[j], npix = npix, img_width = img_width, norm = True)
        j2_image = make_image(j2_4vec[j], jet2_PFCands[j], npix = npix, img_width = img_width, norm = True)


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



print("Finished all batches! Output file saved to %s" %(fout_name))


            
                              

