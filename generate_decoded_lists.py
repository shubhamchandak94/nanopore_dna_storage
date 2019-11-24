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
parser.add_argument('--num_threads',type=int,default=1)
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
PATH_TO_FLAPPIE = "flappie/flappie"
MSG_LEN = args.msg_len
MEM_CONV = args.mem_conv
RATE_CONV = args.rate_conv
info_file = args.info_file
read_id_file = args.read_id_file

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
    fast5_dir = 'tmp_input_' + str(readid)
    os.mkdir(fast5_dir)

    fast5_filename = os.path.join(fast5_dir, 'tmp.'+rnd+'.fast5')
    helper.create_fast5(raw_data,fast5_filename)

    # call flappie to generate transition posterior table
    post_filename = 'tmp.'+rnd+'.post'
    trans_filename = 'tmp.'+rnd+'.trans'
    #subprocess.run([PATH_TO_FLAPPIE, fast5_filename, '--post-output-file', post_filename, '--trans-output-file', trans_filename, '-o',fastq_filename])

    tmp_output_dir = 'tmp_output_' + str(readid) 
    fastq_filename = os.path.join(tmp_output_dir, "fastq_runid_test_0_0.fastq") 
   
    GUPPY_MODEL_PATH='/home/kedart/code/taiyaki/taiyaki_guppy_default/training_1/model_checkpoint_00010.checkpoint.json'
    subprocess.run([helper.PATH_TO_GUPPY, '--input_path', fast5_dir, '--save_path', tmp_output_dir, '--flowcell', 'FLO-MIN106', '--kit', 'SQK-LSK109', '--device', 'auto', '--post_out', '--fast5_out', '--num_callers', '30'])
    #subprocess.run([helper.PATH_TO_GUPPY, '--trim_strategy', 'none', '--model_file', GUPPY_MODEL_PATH, '--input_path', fast5_dir, '--save_path', tmp_output_dir, '--flowcell', 'FLO-MIN106', '--kit', 'SQK-LSK109', '--device', 'auto', '--post_out', '--fast5_out', '--num_callers', '30'])
    
    # Convert guppy output to post_file, trans_file
    guppy_output_fast5 = tmp_output_dir + "/workspace/tmp." + rnd + ".fast5"
    print(guppy_output_fast5)
    #helper.guppy_output_transitions(guppy_output_fast5, trans_filename)
    helper.guppy_output_state_data(guppy_output_fast5, post_filename) 

    with open(trans_filename,'w') as outfile:
        subprocess.run(['./viterbi_basecall.out','decode',post_filename,'tmp.basecall'], stdout=outfile)

    # truncate post according to barcode
    (start_pos, end_pos, dist_start, dist_end) = helper.find_barcode_pos_in_post(trans_filename,fastq_filename,START_BARCODE,END_BARCODE)
    (start_pos_RC, end_pos_RC, dist_start_RC, dist_end_RC) = helper.find_barcode_pos_in_post(trans_filename,fastq_filename,START_BARCODE_RC,END_BARCODE_RC)
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
    
    subprocess.run(['./a.out','decode',new_post_filename,'tmp.basecall'])
    subprocess.run(['cat','tmp.basecall'])

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
    shutil.rmtree(tmp_output_dir)
f_info.close()
