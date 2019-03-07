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
PATH_TO_CPP_EXEC = "/raid/nanopore/shubham/flappie/shubham/viterbi_nanopore.out"
PATH_TO_FLAPPIE = "/raid/nanopore/shubham/flappie/flappie"
NUM_TRIALS = 10
LIST_DECODING_MODE = 1
LIST_SIZE = 1
NUM_THR = 30

SYN_SUB_PROB = 0.0
SYN_DEL_PROB = 0.0
SYN_INS_PROB = 0.0
msg_len = 142 # after attaching sync markers
# Note: markers not attached for the padding.
sync_marker = '110'
sync_marker_period = 9 # first marker at 0, second at sync_marker_period
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
    subprocess.run([PATH_TO_CPP_EXEC,'encode',file_msg,file_seq])
    
    seq = util.read_seq(file_seq)
    os.remove(file_msg)
    os.remove(file_seq)

    fast5_filename='tmp.'+rnd+'.fast5'
    util.simulate_read(seq, SYN_SUB_PROB, SYN_DEL_PROB, SYN_INS_PROB, fast5_filename, deepSimDwellFlag = True)

    # call flappie to generate transition posterior table
    post_filename = 'tmp.'+rnd+'.post'
    decoded_filename = 'tmp.'+rnd+'.dec'
    subprocess.run([PATH_TO_FLAPPIE, fast5_filename, '--post-output-file', post_filename])
    subprocess.run([PATH_TO_CPP_EXEC,'decode',post_filename,decoded_filename,str(msg_len),'-l',str(LIST_DECODING_MODE),str(LIST_SIZE),str(NUM_THR)])
    with open(decoded_filename) as f:
        decoded_msg_list = [l.rstrip('\n') for l in f.readlines()]

    # remove temporary files
    os.remove(fast5_filename)
    os.remove(post_filename)
    os.remove(decoded_filename)

    print('Top message:')
    print(decoded_msg_list[0])
    print('List size:',len(decoded_msg_list))
    top_correct = (decoded_msg_list[0] == msg)
    list_correct = (msg in decoded_msg_list)
    print('Top correct:', top_correct)
    print('List correct:', list_correct)
    correct_list_top.append(top_correct)
    correct_list_list.append(list_correct)
print('Summary statistics:')
print('Number total:', NUM_TRIALS)
print('Number top correct:', sum(correct_list_top))
print('Number list correct:', sum(correct_list_list))