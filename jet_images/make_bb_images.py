import os

n = 40

for i in range(n):
    bb_file = "/eos/cms/store/group/phys_b2g/CASE/h5_files/2017/BB_files/BB_batch%i.h5" % i
    bb_file_out = "/eos/cms/store/group/phys_b2g/CASE/h5_files/2017/BB_files_images/BB_batch%i.h5" % i
    print("Start file %s " % bb_file)
    os.system("cp %s temp.h5" % bb_file)
    os.system("python make_jet_images.py -i temp.h5")
    os.system("cp temp.h5 %s" % bb_file_out)
    os.system("rm -f temp.h5")
    

final_output_name = "/eos/cms/store/group/phys_b2g/CASE/h5_files/2017/BB_v1_with_images.h5"
output_dir = "/eos/cms/store/group/phys_b2g/CASE/h5_files/2017/BB_files_images/*.h5" 
os.system("python ../H5_maker/H5_merge.py %s %s" %(final_output_name, output_dir))
