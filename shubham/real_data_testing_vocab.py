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
import json

START_BARCODE = 'CTGGCGGCCTTGGCCGACTATCTGC'
END_BARCODE = 'ACCATGTCGTACAGTCGTTGTAACA'
JSON_FILE_WITH_RAW_DATA = '/raid/nanopore/pass_rl_250_sorted_coded_04_no_reverse_raw.json'
OLIGO_FILE = '/raid/nanopore/shubham/Nanopore-DNA-storeage/final_oligos/coded_04_uncoded_ccsds_trivialmodulation_rotation_305.txt'
PATH_TO_CPP_EXEC = "/raid/nanopore/shubham/flappie/shubham/viterbi_nanopore.out"
PATH_TO_FLAPPIE = "/raid/nanopore/shubham/flappie/flappie"
NUM_TRIALS = 100 # in each trial pick one read at random from real data
KMER_LEN = 10 # kmer length for when we artificially create a vocabulary
ADD_VOCAB = 80 # add additional random kmers to vocab to make things harder

with open(OLIGO_FILE,'r') as f:
    oligos = [l.rstrip('\n') for l in f.readlines()]
# ASSUME ALL OLIGOS HAVE SAME LENGTH
oligo_len = len(oligos[0])
for oligo in oligos:
    assert len(oligo) == oligo_len

assert oligo_len%KMER_LEN == 0
msg_len = oligo_len//KMER_LEN
print('msg_len',msg_len)
msg = list(range(msg_len)) # message is always trivial, we just change the vocabulary
print('msg',msg)

# read raw data JSON
with open(JSON_FILE_WITH_RAW_DATA,'r') as f:
    raw_data = json.loads(f.read())

num_reads = len(raw_data)
print('num_reads',num_reads)

hamming_list = []
correct_list = []

for _ in range(NUM_TRIALS):
    rnd = str(np.random.randint(10000000))
    read_idx = np.random.randint(num_reads)
    print('read_idx',read_idx)
    print('read_id',raw_data[read_idx][0])
    print('oligo_id',raw_data[read_idx][1])
    read_idx_raw_data = raw_data[read_idx][2]
    read_idx_oligo_num = int(((raw_data[read_idx][1]).split('_'))[-1]) # oligo num is given by part after last _
    correct_oligo = oligos[read_idx_oligo_num]
    # create vocab from correct oligo and write to file
    vocab = [correct_oligo[i*KMER_LEN:(i+1)*KMER_LEN] for i in range(msg_len)]
    vocab = vocab + [''.join(np.random.choice(['A','C','G','T'],KMER_LEN)) for _ in range(ADD_VOCAB)]
    file_vocab = "tmp."+rnd+".vocab"
    with open(file_vocab,'w') as f:
        f.write(''.join([word+'\n' for word in vocab]))
    vocab_size = len(vocab)
    print('vocab_size',vocab_size)

    # create fast5 from raw data
    fast5_filename='tmp.'+rnd+'.fast5'
    util.create_fast5(read_idx_raw_data,fast5_filename)

    # call flappie to generate transition posterior table
    post_filename = 'tmp.'+rnd+'.post'
    fastq_filename = 'tmp.'+rnd+'.fastq'
    trans_filename = 'tmp.'+rnd+'.trans'
    subprocess.run([PATH_TO_FLAPPIE, fast5_filename, '--post-output-file', post_filename, '--trans-output-file', trans_filename, '-o',fastq_filename])
    
    # truncate post according to barcode
    (start_pos, end_pos) = util.find_barcode_pos_in_post(trans_filename,fastq_filename,START_BARCODE,END_BARCODE)
    if start_pos == -1:
        raise NameError('error in find_barcode_pos_in_post')
    if end_pos - start_pos + 1 < msg_len:
        raise NameError('Too short start_pos - end_pos')
    new_post_filename = 'tmp.'+rnd+'.post.new'
    util.truncate_post_file(post_filename, new_post_filename, start_pos, end_pos)

    decoded_filename = 'tmp.'+rnd+'.dec'
    subprocess.run([PATH_TO_CPP_EXEC,'decode',new_post_filename,decoded_filename,str(msg_len),file_vocab])
    with open(decoded_filename) as f:
        decoded_msg = [int(l.rstrip('\n')) for l in f.readlines()]
    assert len(decoded_msg) == msg_len

    # remove temporary files
    os.remove(fast5_filename)
    os.remove(post_filename)
    os.remove(decoded_filename)
    os.remove(new_post_filename)
    os.remove(trans_filename)
    os.remove(fastq_filename)
    os.remove(file_vocab)

    print(decoded_msg)
    hamming_list.append(sum([msg[i] != decoded_msg[i] for i in range(msg_len)]))
    print('Hamming distance',hamming_list[-1])
    correct_list.append((hamming_list[-1] == 0))
print('Summary statistics:')
print('Number total:', NUM_TRIALS)
print('Number correct:', sum(correct_list))
print('Average symbol error rate:', sum(hamming_list)/(msg_len*NUM_TRIALS))
