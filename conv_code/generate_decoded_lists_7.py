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
LIST_SIZE = 64
NUM_THREADS = 10
START_BARCODE = 'GCTAGTACGCGAACAGAGTGCAGTA'
END_BARCODE = 'ACAGATGCAGTAATTCTCACGAACT'
START_BARCODE_RC = util.reverse_complement(END_BARCODE)
END_BARCODE_RC = util.reverse_complement(START_BARCODE)
OUT_DIR = '/raid/nanopore/shubham/20190629_nanopore_data/decoded_lists/8_5/'
HDF5_FILE = '/raid/nanopore/shubham/20190629_nanopore_data/raw_signal/raw_signal_6.hdf5'
PATH_TO_CPP_EXEC = "/raid/nanopore/shubham/nanopore_dna_storage/conv_code/viterbi_nanopore.out"
PATH_TO_FLAPPIE = "/raid/nanopore/shubham/nanopore_dna_storage/flappie/flappie"
NUM_READS = 20000
MSG_LEN = 180
MEM_CONV = 8
RATE_CONV = 5

f_info = open(OUT_DIR+'info.txt','w') # store the ref and read id for future reference
f_raw = h5py.File(HDF5_FILE)
readid_list = list(f_raw.keys())
readid_list_small = random.sample(readid_list,NUM_READS)
for i,readid in enumerate(readid_list_small):
    print("i:",i)
    rnd = str(np.random.randint(10000000))
    raw_data = f_raw[readid]['raw_signal']
    print(readid)
    print(f_raw[readid].attrs['ref'])
    f_info.write(readid+'\t'+f_raw[readid].attrs['ref'].decode("utf-8")+'\n')
    # create fast5 from raw data
    fast5_filename='tmp.'+rnd+'.fast5'
    util.create_fast5(raw_data,fast5_filename)

    # call flappie to generate transition posterior table
    post_filename = 'tmp.'+rnd+'.post'
    fastq_filename = 'tmp.'+rnd+'.fastq'
    trans_filename = 'tmp.'+rnd+'.trans'
    subprocess.run([PATH_TO_FLAPPIE, fast5_filename, '--post-output-file', post_filename, '--trans-output-file', trans_filename, '-o',fastq_filename])

    # truncate post according to barcode
    (start_pos, end_pos, dist_start, dist_end) = util.find_barcode_pos_in_post(trans_filename,fastq_filename,START_BARCODE,END_BARCODE)
    (start_pos_RC, end_pos_RC, dist_start_RC, dist_end_RC) = util.find_barcode_pos_in_post(trans_filename,fastq_filename,START_BARCODE_RC,END_BARCODE_RC)
    rc = False
    if dist_start + dist_end > dist_start_RC + dist_end_RC:
        rc = True
        start_pos = start_pos_RC
        end_pos = end_pos_RC

    if start_pos == -1 or end_pos - start_pos + 1 < MEM_CONV+MSG_LEN+1:
        print ("Failure in barcode removing.")
        os.remove(fast5_filename)
        os.remove(post_filename)
        os.remove(trans_filename)
        os.remove(fastq_filename)
        continue
    new_post_filename = 'tmp.'+rnd+'.post.new'
    util.truncate_post_file(post_filename, new_post_filename, start_pos, end_pos)
    decoded_filename = OUT_DIR + 'list_' + str(i)
    rc_flag = ''
    if rc:
        rc_flag = '--rc'
    print(rc_flag)
    subprocess.run([PATH_TO_CPP_EXEC,'-m','decode','-i',new_post_filename,'-o',decoded_filename,'--msg-len',str(MSG_LEN),'-l',str(LIST_SIZE),'-t',str(NUM_THREADS),'--mem-conv',str(MEM_CONV),rc_flag,'--max-deviation','20','-r',str(RATE_CONV)])

    os.remove(fast5_filename)
    os.remove(post_filename)
    os.remove(new_post_filename)
    os.remove(trans_filename)
    os.remove(fastq_filename)

f_info.close()
