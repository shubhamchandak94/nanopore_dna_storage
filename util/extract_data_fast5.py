import json
import os
import h5py
import sys
import numpy as np
if len(sys.argv) != 4:
    print('Incorrect usage')
    print('python extract_data_fast5.py SAMFILENAME FAST5DIR OUTFILE.hdf5')
    sys.exit(1)

infile_sam = sys.argv[1]
fast5_dir = sys.argv[2]
outfile = sys.argv[3]

print('SAM file:',infile_sam)
print('FAST5 dir:',fast5_dir)
print('HDF5 outfile:',outfile)

# first read sam file and make a dictionary mapping the read id to the reference sequence
sam_dict = {}
with open(infile_sam,'r') as f:
    for line in f:
        if line[0] == '@':
            continue
        arr = line.split()
        read_id = 'read_'+(arr[0].split(' '))[0]
        ref = arr[2]
        sam_dict[read_id] = ref

print('number of read ids in sam', len(sam_dict))

# now go through each of the fast5 files and if the read_id is in sam_dict then put the reference and the raw signals (after conversion to float) into an hdf5 file
fout = h5py.File(outfile,'w')
num_done = 0
for filename in os.listdir(fast5_dir):
    with h5py.File(fast5_dir+'/'+filename) as f:
        for read_id in f.keys():
            if read_id in sam_dict:
                fout.create_group(read_id)
                fout[read_id].create_dataset("raw_signal",data=f[read_id]['Raw']['Signal'])
                fout[read_id].attrs.create("ref",data = np.string_(sam_dict[read_id]))
                num_done += 1
                if num_done%1000 == 0:
                    print(str(num_done),'done')

print('number of matching read ids found in fast5', num_done)

fout.close()
