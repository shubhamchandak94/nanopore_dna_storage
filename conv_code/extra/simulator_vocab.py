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
infile_vocab = "/raid/nanopore/shubham/flappie/shubham/vocab_files/vocab_10mer_256"
with open(infile_vocab,'r') as f:
    vocab = [l.rstrip('\n') for l in f.readlines()]
vocab_size = len(vocab)
NUM_TRIALS = 100
SYN_SUB_PROB = 0.02
SYN_DEL_PROB = 0.02
SYN_INS_PROB = 0.005
msg_len = 12 # number of words in oligo 
print('msg_len:',msg_len)

hamming_list = []
correct_list = []

for _ in range(NUM_TRIALS):
    msg = np.random.randint(vocab_size,size=msg_len).tolist()
    print(msg)
    assert len(msg) == msg_len
    rnd = str(np.random.randint(10000000))
    seq = ''.join([vocab[i] for i in msg])

    fast5_filename='tmp.'+rnd+'.fast5'
    util.simulate_read(seq, SYN_SUB_PROB, SYN_DEL_PROB, SYN_INS_PROB, fast5_filename, deepSimDwellFlag = True)

    # call flappie to generate transition posterior table
    post_filename = 'tmp.'+rnd+'.post'
    decoded_filename = 'tmp.'+rnd+'.dec'
    subprocess.run([PATH_TO_FLAPPIE, fast5_filename, '--post-output-file', post_filename])
    subprocess.run([PATH_TO_CPP_EXEC,'decode',post_filename,decoded_filename,str(msg_len),infile_vocab])
    with open(decoded_filename) as f:
        decoded_msg = [int(l.rstrip('\n')) for l in f.readlines()]
    assert len(decoded_msg) == msg_len

    # remove temporary files
    os.remove(fast5_filename)
    os.remove(post_filename)
    os.remove(decoded_filename)

    print(decoded_msg)
    hamming_list.append(sum([msg[i] != decoded_msg[i] for i in range(msg_len)]))
    print('Hamming distance',hamming_list[-1])
    correct_list.append((hamming_list[-1] == 0))
print('Summary statistics:')
print('Number total:', NUM_TRIALS)
print('Number correct:', sum(correct_list))
print('Average symbol error rate:', sum(hamming_list)/(msg_len*NUM_TRIALS))
