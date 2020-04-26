import os
import subprocess
import helper
import binascii
import struct
import numpy as np
import h5py
import warnings
import random
import argparse
import shutil 
import h5py

parser = argparse.ArgumentParser(description='generate decoded lists from raw signal')
parser.add_argument('--hdf_file',type=str,required=True)
parser.add_argument('--out_prefix',type=str,required=True)
parser.add_argument('--read_id_file',type=str,required=True)
parser.add_argument('--info_file',type=str,required=True)
parser.add_argument('--mem_conv',type=int,required=True)
parser.add_argument('--msg_len',type=int,required=True)
parser.add_argument('--rate_conv',type=int,required=True)
parser.add_argument('--list_size',type=int,required=True)
parser.add_argument('--start_barcode',type=str,required=True)
parser.add_argument('--end_barcode',type=str,required=True)
parser.add_argument('--barcode_search_extend_len',type=int,help='default: set to min of start and end barcode len')
parser.add_argument('--barcode_extend_penalty',type=float,default=0.6)
parser.add_argument('--num_threads',type=int,default=1)
parser.add_argument('--bonito_model_path',type=str)

args = parser.parse_args()
print(args)

# warnings.filterwarnings("ignore", category=DeprecationWarning)

LIST_SIZE = args.list_size
NUM_THREADS = args.num_threads
START_BARCODE = args.start_barcode
END_BARCODE = args.end_barcode
START_BARCODE_RC = helper.reverse_complement(END_BARCODE)
END_BARCODE_RC = helper.reverse_complement(START_BARCODE)
OUT_PREFIX = args.out_prefix
HDF5_FILE = args.hdf_file
PATH_TO_CPP_EXEC = "viterbi/viterbi_nanopore.out"
MSG_LEN = args.msg_len
MEM_CONV = args.mem_conv
RATE_CONV = args.rate_conv
info_file = args.info_file
read_id_file = args.read_id_file
if args.bonito_model_path is None:
    bonito_model_path = 'dna_r9.4.1'
else:
    bonito_model_path = args.bonito_model_path
if args.barcode_search_extend_len is None:
    barcode_search_extend_len = min(len(START_BARCODE),len(END_BARCODE))
else:
    barcode_search_extend_len = args.barcode_search_extend_len

with open(read_id_file) as f:
    readid_list = [l.rstrip('\n') for l in f.readlines()]
f_info = open(info_file,'w') # store the ref and read id for future reference
f_raw = h5py.File(HDF5_FILE,'r')

for i,readid in enumerate(readid_list):
    print("i:",i)
    rnd = str(np.random.randint(10000000))
    raw_data = f_raw[readid]['raw_signal']
    print(readid)
    print(f_raw[readid].attrs['ref'])
    f_info.write(readid+'\t'+f_raw[readid].attrs['ref'].decode("utf-8")+'\n')
    # create fast5 from raw data
    fast5_dir = 'tmp_input_' + rnd + '_' + str(readid)
    os.mkdir(fast5_dir)

    fast5_filename = os.path.join(fast5_dir, 'tmp.'+rnd+'.fast5')
    helper.create_fast5(raw_data,fast5_filename)

    # call bonito to generate CTC posterior table
    post_filename = 'tmp.'+rnd+'_'+readid+'.post'
    trans_filename = 'tmp.'+rnd+'_'+readid+'.trans'
    fastq_filename = 'tmp.'+rnd+'_'+readid+'.fastq'

    subprocess.run(['bonito','basecaller', bonito_model_path, fast5_dir, '--post_file', post_filename])

    # convert bonito post file to fastq and move (trans) files
    helper.bonito_basecall_generate_move(post_filename,fastq_filename,trans_filename)

    # truncate post according to barcode
    (start_pos, end_pos, dist_start, dist_end) = helper.find_barcode_pos_in_post(trans_filename,fastq_filename,START_BARCODE,END_BARCODE,barcode_search_extend_len,args.barcode_extend_penalty)
    (start_pos_RC, end_pos_RC, dist_start_RC, dist_end_RC) = helper.find_barcode_pos_in_post(trans_filename,fastq_filename,START_BARCODE_RC,END_BARCODE_RC,barcode_search_extend_len,args.barcode_extend_penalty)
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
    helper.truncate_post_file(post_filename, new_post_filename, start_pos, end_pos)
    
    decoded_filename = OUT_PREFIX + '_' + str(i)
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
    shutil.rmtree(fast5_dir)
f_info.close()
