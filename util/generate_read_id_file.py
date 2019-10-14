import os
import subprocess
import binascii
import struct
import numpy as np
import h5py
import warnings
import random
warnings.filterwarnings("ignore", category=DeprecationWarning)
OUT_DIR = '../../nanopore_dna_storage_data/decoded_lists/exp_7/'
HDF5_FILE = '../../nanopore_dna_storage_data/raw_signal/raw_signal_7.hdf5'
NUM_READS = 20000

f_raw = h5py.File(HDF5_FILE)
readid_list = list(f_raw.keys())
readid_list_small = random.sample(readid_list,NUM_READS)
f_read_id = open(OUT_DIR+'read_ids.txt','w')
for rid in readid_list_small:
    f_read_id.write(rid+'\n')
f_read_id.close()
