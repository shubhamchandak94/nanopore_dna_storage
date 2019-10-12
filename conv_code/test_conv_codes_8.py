import os
import subprocess
import util
import binascii
import struct
import numpy as np
import h5py
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
LIST_SIZE = 1
NUM_THREADS = 16
START_BARCODE = 'CACCTGTGCTGCGTCAGGCTGTGTC'
END_BARCODE = 'GCTGTCCGTTCCGCATTGACACGGC'
START_BARCODE_RC = util.reverse_complement(END_BARCODE)
END_BARCODE_RC = util.reverse_complement(START_BARCODE)
HDF5_FILE = '/raid/nanopore/shubham/20190629_nanopore_data/raw_signal/raw_signal_7.hdf5'
CONV_INPUT_FILE = '/raid/nanopore/shubham/20190629_nanopore_data/oligo_files/conv_input_7.txt'
PATH_TO_CPP_EXEC = "/raid/nanopore/shubham/nanopore_dna_storage/conv_code/viterbi_nanopore.out"
PATH_TO_FLAPPIE = "/raid/nanopore/shubham/nanopore_dna_storage/flappie/flappie"
NUM_ATTEMPTS = 100 
MSG_LEN = 180
MEM_CONV = 11
RATE_CONV = 5

f_raw = h5py.File(HDF5_FILE)
readid_list = list(f_raw.keys())
with open(CONV_INPUT_FILE) as f:
    conv_input_list = [s.rstrip('\n') for s in f.readlines()] 
num_success = 0
num_attempted = 0
for _ in range(NUM_ATTEMPTS):
    num_attempted += 1
    print("num_attempted:",num_attempted)
    rnd = str(np.random.randint(10000000))

    readid = np.random.choice(readid_list)
    raw_data = f_raw[readid]['raw_signal']
    print(f_raw[readid].attrs['ref'])
    bit_string = conv_input_list[int(((f_raw[readid].attrs['ref']).decode('utf-8')).split('_')[-1])]
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

    if start_pos == -1 or end_pos - start_pos + 1 < MEM_CONV+MSG_LEN+1:
        print ("Failure in barcode removing.")
        os.remove(fast5_filename)
        os.remove(post_filename)
        os.remove(trans_filename)
        os.remove(fastq_filename)
        continue
    new_post_filename = 'tmp.'+rnd+'.post.new'
    if not rc:
        util.truncate_post_file(post_filename, new_post_filename, start_pos, end_pos)
    else:
        util.truncate_post_file(post_filename, new_post_filename, start_pos_RC, end_pos_RC)
    decoded_filename = 'tmp.'+rnd+'.dec'
    rc_flag = ''
    if rc:
        rc_flag = '--rc'
    print(rc_flag)
    subprocess.run([PATH_TO_CPP_EXEC,'-m','decode','-i',new_post_filename,'-o',decoded_filename,'--msg-len',str(MSG_LEN),'-l',str(LIST_SIZE),'-t',str(NUM_THREADS),'--mem-conv',str(MEM_CONV),rc_flag,'--max-deviation','20','-r',str(RATE_CONV)])

    with open(decoded_filename) as f:
        decoded_msg_list = [l.rstrip('\n') for l in f.readlines()]
    print(decoded_msg_list[0])
    print(bit_string)
    if decoded_msg_list[0] == bit_string:
        num_success += 1
        print('success')
    os.remove(fast5_filename)
    os.remove(post_filename)
    os.remove(decoded_filename)
    os.remove(new_post_filename)
    os.remove(trans_filename)
    os.remove(fastq_filename)

print('num success:',num_success)
