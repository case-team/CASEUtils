from H5_merge import *


if __name__ == "__main__":
    #print(sys.argv[1], sys.argv[2:])
    
    pt_bin = sys.argv[3]
    input_dir = sys.argv[2]
    
    onlyfiles = [join(input_dir, f) for f in listdir(input_dir) if (isfile(join(input_dir, f)) & (pt_bin in f))]
    merge_multiple(sys.argv[1], onlyfiles)
    
    print("Done!")

