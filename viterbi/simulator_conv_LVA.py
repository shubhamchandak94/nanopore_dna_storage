import sys
import numpy as np
import scrappy
from uuid import uuid4
import scipy.stats as st
import subprocess
from fast5_research import Fast5
import os
import distance
import util
import argparse

parser = argparse.ArgumentParser(description='Simulation for convolutional code.')
parser.add_argument('--num_trials',type=int,default=100)
parser.add_argument('--list_size',type=int,default=1)
parser.add_argument('--num_thr',type=int,default=8)
parser.add_argument('--mem_conv',type=int,default=6)
parser.add_argument('--rate',type=int,default=1)
parser.add_argument('--msg_len',type=int,default=100)
parser.add_argument('--deepsimdwell',type=str,default="False")
parser.add_argument('--reversecomp',type=str,default="False")
args = parser.parse_args()
print(args)

PATH_TO_CPP_EXEC = "/raid/nanopore/shubham/flappie_new/viterbi/viterbi_nanopore.out"
PATH_TO_FLAPPIE = "/raid/nanopore/shubham/flappie_new/flappie"
NUM_TRIALS = args.num_trials
LIST_SIZE = args.list_size
NUM_THR = args.num_thr
MEM_CONV = args.mem_conv
RATE = args.rate
msg_len = args.msg_len # after attaching sync markers
# Note: markers not attached for the padding.

SYN_SUB_PROB = 0.005
SYN_DEL_PROB = 0.005
SYN_INS_PROB = 0.0005
if args.deepsimdwell == "False":
    deepsimdwell = False
else:
    deepsimdwell = True

if args.reversecomp == "False":
    revcomp = False
else:
    revcomp = True

sync_marker = ''
sync_marker_period = 67 # first marker at 0, second at sync_marker_period
msg_len_before_sync_markers = msg_len - (msg_len//sync_marker_period)*len(sync_marker) - min([msg_len%sync_marker_period,len(sync_marker)]) # subtract marker lengths from complete periods and last incomplete period
print('msg_len_before_sync_markers',msg_len_before_sync_markers)

correct_list_top = []
correct_list_list = []

for _ in range(NUM_TRIALS):
    msg_before_markers = np.random.choice(['0','1'], msg_len_before_sync_markers)
    msg_after_markers = []
    pos_in_msg_vec = 0
    for i in range(msg_len):
        if i%sync_marker_period < len(sync_marker):
            msg_after_markers.append(sync_marker[i%sync_marker_period])
        else:
            msg_after_markers.append(msg_before_markers[pos_in_msg_vec])
            pos_in_msg_vec += 1
    assert pos_in_msg_vec == msg_len_before_sync_markers
    msg = ''.join(msg_after_markers)
    print(msg)
    assert len(msg) == msg_len
    rnd = str(np.random.randint(10000000))
    file_msg = 'tmp.'+rnd
    file_seq = 'tmp.enc.'+rnd
    with open(file_msg,'w') as f:
        f.write(msg)
    subprocess.run([PATH_TO_CPP_EXEC,'-m','encode','-i',file_msg,'-o',file_seq,'--mem-conv',str(MEM_CONV),'-r',str(RATE),'--msg-len',str(msg_len)])
    
    seq = util.read_seq(file_seq)
    print('len(seq):',len(seq))

    if revcomp:
        seq = util.reverse_complement(seq)

    fast5_filename='tmp.'+rnd+'.fast5'
    util.simulate_read(seq, SYN_SUB_PROB, SYN_DEL_PROB, SYN_INS_PROB, fast5_filename, deepSimDwellFlag = deepsimdwell)

    # call flappie to generate transition posterior table
    post_filename = 'tmp.'+rnd+'.post'
    decoded_filename = 'tmp.'+rnd+'.dec'
    subprocess.run([PATH_TO_FLAPPIE, fast5_filename, '--post-output-file', post_filename])
    rc_flag = ''
    if revcomp:
        rc_flag = '--rc'
    subprocess.run([PATH_TO_CPP_EXEC,'-m', 'decode','-i',post_filename,'-o',decoded_filename,'--mem-conv',str(MEM_CONV),'--msg-len',str(msg_len),'-l',str(LIST_SIZE),'-t',str(NUM_THR),'-r',str(RATE),rc_flag])
    with open(decoded_filename) as f:
        decoded_msg_list = [l.rstrip('\n') for l in f.readlines()]

    print('Top message:')
    print(decoded_msg_list[0])
    print('List size:',len(decoded_msg_list))
    top_correct = (decoded_msg_list[0] == msg)
    list_correct = (msg in decoded_msg_list)
    print('Top correct:', top_correct)
    print('List correct:', list_correct)
    correct_list_top.append(top_correct)
    correct_list_list.append(list_correct)

    # remove temporary files
    os.remove(file_msg)
    os.remove(file_seq)
    os.remove(fast5_filename)
    os.remove(post_filename)
    os.remove(decoded_filename)


print('Summary statistics:')
print('Number total:', NUM_TRIALS)
print('Number top correct:', sum(correct_list_top))
print('Number list correct:', sum(correct_list_list))
