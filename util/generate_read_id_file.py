import os
import subprocess
import binascii
import struct
import numpy as np
import h5py
import warnings
import random
import sys
warnings.filterwarnings("ignore", category=DeprecationWarning)
EXP_ID = sys.argv[1]
NUM_READS = int(sys.argv[2])
DATA_PATH = sys.argv[3]
if len(sys.argv) == 5:
    REPLICATE = sys.argv[4]
    READ_IDS_FILE = DATA_PATH+'/read_ids/'+REPLICATE+'/exp_'+str(EXP_ID)+'_read_ids.'+str(NUM_READS)+'.txt'
    HDF5_FILE = DATA_PATH+'/raw_signal/'+REPLICATE+'/raw_signal_'+str(EXP_ID)+'.hdf5'
else:
    READ_IDS_FILE = DATA_PATH+'/read_ids/exp_'+str(EXP_ID)+'_read_ids.'+str(NUM_READS)+'.txt'
    HDF5_FILE = DATA_PATH+'/raw_signal/raw_signal_'+str(EXP_ID)+'.hdf5'

f_raw = h5py.File(HDF5_FILE)
readid_list = list(f_raw.keys())
print('len(readid_list):',len(readid_list))
readid_list_small = random.sample(readid_list,NUM_READS)
f_read_id = open(READ_IDS_FILE,'w')
for rid in readid_list_small:
    f_read_id.write(rid+'\n')
f_read_id.close()
