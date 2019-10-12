import os
import subprocess
import util
import binascii
import struct
import numpy as np
import h5py
import warnings
import random
warnings.filterwarnings("ignore", category=DeprecationWarning)
OUT_DIR = '/raid/nanopore/shubham/20190629_nanopore_data/decoded_lists/11_5/'
HDF5_FILE = '/raid/nanopore/shubham/20190629_nanopore_data/raw_signal/raw_signal_7.hdf5'
NUM_READS = 20000
NUM_READS_PER_FILE = 2000

f_raw = h5py.File(HDF5_FILE)
readid_list = list(f_raw.keys())
readid_list_small = random.sample(readid_list,NUM_READS)

start_id_num = 0
i = 0
done = False
while True:
    num_read_ids_file = NUM_READS_PER_FILE
    if start_id_num + num_read_ids_file >= NUM_READS:
        done = True
        num_read_ids_file = NUM_READS - start_id_num
    if num_read_ids_file == 0:
        break
    f_read_id = open(OUT_DIR+'read_ids_'+str(i)+'.txt','w')
    for j in range(start_id_num, start_id_num+num_read_ids_file):
        f_read_id.write(readid_list_small[j]+'\n')
    f_read_id.close()
    i += 1
    start_id_num += num_read_ids_file
    if done:
        break
