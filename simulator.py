import sys
import numpy as np
import scrappy
import subprocess
from fast5_research import Fast5
import os
import distance
import helper
import argparse
import math
import warnings
import time
warnings.filterwarnings("ignore", category=DeprecationWarning) 


parser = argparse.ArgumentParser(description='Simulation for convolutional code.')
parser.add_argument('--num_trials',type=int,default=100)
parser.add_argument('--list_size',type=int,default=1)
parser.add_argument('--num_thr',type=int,default=8)
parser.add_argument('--mem_conv',type=int,default=6)
parser.add_argument('--rate',type=int,default=1)
parser.add_argument('--msg_len',type=int,default=100)
parser.add_argument('--deepsimdwell',type=str,default="False")
parser.add_argument('--reversecomp',type=str,default="False")
parser.add_argument('--syn_sub_prob',type=float,default=0.002)
parser.add_argument('--syn_del_prob',type=float,default=0.0085)
parser.add_argument('--syn_ins_prob',type=float,default=0.0005)
args = parser.parse_args()
print(args)

PATH_TO_CPP_EXEC = "viterbi/viterbi_nanopore.out"
NUM_TRIALS = args.num_trials
LIST_SIZE = args.list_size
NUM_THR = args.num_thr
MEM_CONV = args.mem_conv
RATE = args.rate
msg_len = args.msg_len 

SYN_SUB_PROB = args.syn_sub_prob
SYN_DEL_PROB = args.syn_del_prob
SYN_INS_PROB = args.syn_ins_prob
if args.deepsimdwell == "False":
    deepsimdwell = False
else:
    deepsimdwell = True

if args.reversecomp == "False":
    revcomp = False
else:
    revcomp = True

correct_list_top = []
correct_list_list = []
hamming_list = []
hamming_list_8 = []
hamming_list_16 = []
edit_list = []

for _ in range(NUM_TRIALS):
    msg = ''.join(np.random.choice(['0','1'], msg_len))
    print(msg)
    rnd = str(np.random.randint(10000000))
    file_msg = 'tmp.'+rnd
    file_seq = 'tmp.enc.'+rnd
    with open(file_msg,'w') as f:
        f.write(msg+'\n')
    subprocess.run([PATH_TO_CPP_EXEC,'-m','encode','-i',file_msg,'-o',file_seq,'--mem-conv',str(MEM_CONV),'-r',str(RATE),'--msg-len',str(msg_len)])
    
    seq = helper.read_seq(file_seq)
    print('len(seq):',len(seq))

    if revcomp:
        seq = helper.reverse_complement(seq)

    tmp_output_dir = 'tmp_output_' +rnd + '/'
    os.mkdir(tmp_output_dir)
    fast5_filename= tmp_output_dir+'tmp.'+rnd+'.fast5'
    helper.simulate_read(seq, SYN_SUB_PROB, SYN_DEL_PROB, SYN_INS_PROB, fast5_filename, deepSimDwellFlag = deepsimdwell)

    # call bonito to generate posterior table
    post_filename = 'tmp.'+rnd+'.post'
    decoded_filename = 'tmp.'+rnd+'.dec'

    subprocess.run(['bonito','basecaller', 'dna_r9.4.1', tmp_output_dir, '--post_file', post_filename])
    rc_flag = ''
    if revcomp:
        rc_flag = '--rc'
    start = time.time()
    subprocess.run([PATH_TO_CPP_EXEC,'-m', 'decode','-i',post_filename,'-o',decoded_filename,'--mem-conv',str(MEM_CONV),'--msg-len',str(msg_len),'-l',str(LIST_SIZE),'-t',str(NUM_THR),'-r',str(RATE),rc_flag,'--max-deviation','20'])
    stop = time.time()
    print("Decoding time:",stop-start,"s")
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
    hamming_list.append(distance.hamming(msg,decoded_msg_list[0]))
    hamming_list_8.append(np.sum([(decoded_msg_list[0][i*8:(i+1)*8] != msg[i*8:(i+1)*8]) for i in range(math.ceil(len(msg)/8))]))
    hamming_list_16.append(np.sum([(decoded_msg_list[0][i*16:(i+1)*16] != msg[i*16:(i+1)*16]) for i in range(math.ceil(len(msg)/16))]))
    print('Hamming distance of top:',hamming_list[-1])
    print('Hamming distance of top (8 blocks):',hamming_list_8[-1])
    print('Hamming distance of top (16 blocks):',hamming_list_16[-1])
    edit_list.append(distance.levenshtein(msg,decoded_msg_list[0]))
    print('Edit distance',edit_list[-1])
    if not top_correct:
        print('Error pattern (original, errors):')
        print(msg)
        print(''.join([msg[i] if (msg[i] == decoded_msg_list[0][i]) else '*' for i in range(len(msg))]))

    # remove temporary files
    os.remove(file_msg)
    os.remove(file_seq)
    os.remove(fast5_filename)
    os.rmdir(tmp_output_dir)
    os.remove(post_filename)
    os.remove(decoded_filename)


print('Summary statistics:')
print('Number total:', NUM_TRIALS)
print('Number top correct:', sum(correct_list_top))
print('Number list correct:', sum(correct_list_list))
print('Average bit error rate of top:', sum(hamming_list)/(msg_len*NUM_TRIALS))
print('Average 8 block error rate of top:', sum(hamming_list_8)/(math.ceil(msg_len/8)*NUM_TRIALS))
print('Average 16 block error rate of top:', sum(hamming_list_16)/(math.ceil(msg_len/16)*NUM_TRIALS))
print('Average edit distance rate of top:', sum(edit_list)/(msg_len*NUM_TRIALS))
