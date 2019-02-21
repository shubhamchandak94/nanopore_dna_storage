import sys
import numpy as np
import scrappy
import scipy.stats as st
import subprocess
import os
import distance
import util
PATH_TO_CPP_EXEC = "/raid/nanopore/shubham/flappie/shubham/viterbi_nanopore.out"
PATH_TO_FLAPPIE = "/raid/nanopore/shubham/flappie/flappie"
NUM_TRIALS = 100
SYN_SUB_PROB = 0.0
SYN_DEL_PROB = 0.0
SYN_INS_PROB = 0.0
msg_len = 139 # after attaching sync markers
# Note: markers not attached for the padding.
sync_marker = '110'
sync_marker_period = 9 # first marker at 0, second at sync_marker_period
msg_len_before_sync_markers = msg_len - (msg_len//sync_marker_period)*len(sync_marker) - min([msg_len%sync_marker_period,len(sync_marker)]) # subtract marker lengths from complete periods and last incomplete period
print('msg_len_before_sync_markers',msg_len_before_sync_markers)

# from deepsimulator (probably more realistic dwell times)
def rep_rvs(size,a):
    # array_1 = np.ones(int(size*0.075)).astype(int)
    # samples = st.alpha.rvs(3.3928495261646932, 
    #     -7.6451557771999035, 50.873948369526737, 
    #     size=(size-int(size*0.075))).astype(int)
    # samples = np.concatenate((samples, array_1), 0)
    # samples[samples<1] = 2
    # print(a)
    a = a*5
    array_1 = np.ones(int(size*(0.075-0.015*a))).astype(int)
    samples = st.alpha.rvs(3.3928495261646932+a, 
        -7.6451557771999035+(2*a), 50.873948369526737, 
        size=(size-int(size*(0.075-0.015*a)))).astype(int)
    samples = np.concatenate((samples, array_1), 0)
    samples[samples<1] = 2
    np.random.shuffle(samples)
    return samples

hamming_list = []
edit_list = []
correct_list = []

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
    with open('tmp.'+rnd,'w') as f:
        f.write(msg)
    subprocess.run([PATH_TO_CPP_EXEC,'encode','tmp.'+rnd,'tmp.enc.'+rnd])

    infile_seq = 'tmp.enc.'+rnd
    # can handle fasta header, but also works without it
    f = open(infile_seq)
    seq = f.readline().rstrip('\n')
    if seq[0] == '>':
            seq = f.readline().rstrip('\n')
    f.close()
    len_seq = len(seq)
    print('Length of seq: ', len_seq)
    print(seq)
    os.remove('tmp.'+rnd)
    os.remove('tmp.enc.'+rnd)
    syn_seq = util.simulate_indelsubs(seq, sub_prob = SYN_SUB_PROB, del_prob = SYN_DEL_PROB, ins_prob = SYN_INS_PROB)
    print('Length of synthesized sequence', len(syn_seq))
    print(syn_seq)
    squiggle_array = scrappy.sequence_to_squiggle(syn_seq,rescale=True).data(as_numpy=True)
    raw_data = np.array([])

    # for dwell time use deepsimulator since the one provided by scrappie is way too good
    # scrappie dwell gives around 3-4% edit distance while ds dwell gives around 15%
    ds_alpha = 0.1 # 0.1 is default parameter in deepsim
    squiggle_array[:,0] = rep_rvs(squiggle_array.shape[0], ds_alpha)

    for squig in squiggle_array:
            mean = squig[1]
            stdv = squig[2]
            dwell = squig[0]
            raw_data = np.append(raw_data, np.random.laplace(mean, stdv/np.sqrt(2), int(round(dwell))))

    print('Length of raw signal: ', len(raw_data))

    # flappie


    decoded_msg = scrappy.basecall_raw_viterbi_conv(raw_data,PATH_TO_CPP_EXEC, msg_len)
    print(decoded_msg)
    hamming_list.append(distance.hamming(msg,decoded_msg))
    print('Hamming distance',hamming_list[-1])
    edit_list.append(distance.levenshtein(msg,decoded_msg))
    print('Edit distance',edit_list[-1])
    correct_list.append((hamming_list[-1] == 0))
print('Summary statistics:')
print('Number total:', NUM_TRIALS)
print('Number correct:', sum(correct_list))
print('Average bit error rate:', sum(hamming_list)/(msg_len_before_sync_markers*NUM_TRIALS))
