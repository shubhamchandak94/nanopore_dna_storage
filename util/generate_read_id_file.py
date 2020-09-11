import os
import subprocess
import binascii
import struct
import numpy as np
import h5py
import warnings
import random
warnings.filterwarnings("ignore", category=DeprecationWarning)
EXP_ID = 17
NUM_READS = 10000
READ_IDS_FILE = '/raid/nanopore/shubham/20200804_nanopore_pool_data/data/read_ids/exp_'+str(EXP_ID)+'_read_ids.'+str(NUM_READS)+'.txt'
HDF5_FILE = '/raid/nanopore/shubham/20200804_nanopore_pool_data/data/raw_signal/raw_signal_'+str(EXP_ID)+'.hdf5'

f_raw = h5py.File(HDF5_FILE)
readid_list = list(f_raw.keys())
print('len(readid_list):',len(readid_list))
readid_list_small = random.sample(readid_list,NUM_READS)
f_read_id = open(READ_IDS_FILE,'w')
for rid in readid_list_small:
    f_read_id.write(rid+'\n')
f_read_id.close()
