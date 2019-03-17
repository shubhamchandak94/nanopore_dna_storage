import sys
import numpy as np
import scrappy
from uuid import uuid4
import scipy.stats as st
import subprocess
from fast5_research import Fast5
import distance

def simulate_indelsubs(read, sub_prob = 0.0, del_prob = 0.0, ins_prob = 0.0):
    '''
    add iid indels and substitions to read
    '''
    char_list = [c for c in read]
    pos_in_char_list = 0
    new_char_list = []
    alphabet = {}
    alphabet['all'] = ['A','C','G','T']
    alphabet['A'] = ['C','G','T']
    alphabet['C'] = ['A','G','T']
    alphabet['G'] = ['C','A','T']
    alphabet['T'] = ['C','G','A']
    while True:
        ins = (np.random.random_sample()<ins_prob)
        if ins:
            new_char_list.append(np.random.choice(alphabet['all']))
        else:
            if pos_in_char_list == len(char_list):# end of original read and not inserting
                break
            _del = (np.random.random_sample()<del_prob) 
            if _del:
                pos_in_char_list += 1
            else:
                sub = (np.random.random_sample()<sub_prob)
                if sub:
                    new_char_list.append(np.random.choice(alphabet[char_list[pos_in_char_list]]))
                else:
                    new_char_list.append(char_list[pos_in_char_list])
                pos_in_char_list += 1
    return ''.join(new_char_list)

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

def create_fast5(raw_data, fast5_filename):
    raw_data = np.array(raw_data)
    # create fast5 (from https://nanoporetech.github.io/fast5_research/examples.html)
    # example of how to digitize data
    start, stop = int(min(raw_data - 1)), int(max(raw_data + 1))
    rng = stop - start
    digitisation = 8192.0
    bins = np.arange(start, stop, rng / digitisation)
    # np.int16 is required, the library will refuse to write anything other
    raw_data_binned = np.digitize(raw_data, bins).astype(np.int16)

    # The following are required meta data
    channel_id = {
        'digitisation': digitisation,
        'offset': 0,
        'range': rng,
        'sampling_rate': 4000,
        'channel_number': 1,
        }
    read_id = {
        'start_time': 0,
        'duration': len(raw_data),
        'read_number': 1,
        'start_mux': 1,
        'read_id': str(uuid4()),
        'scaling_used': 1,
        'median_before': 0,
    }
    tracking_id = {
        'exp_start_time': '1970-01-01T00:00:00Z',
        'run_id': str(uuid4()).replace('-',''),
        'flow_cell_id': 'FAH00000',
    }
    context_tags = {}

    with Fast5.New(fast5_filename, 'w', tracking_id=tracking_id, context_tags=context_tags, channel_id=channel_id) as h:
        h.set_raw(raw_data_binned, meta=read_id, read_number=1)

def simulate_read(seq, SYN_SUB_PROB, SYN_DEL_PROB, SYN_INS_PROB, fast5_filename, deepSimDwellFlag = True, deepSimAlpha = 0.1):
    syn_seq = simulate_indelsubs(seq, sub_prob = SYN_SUB_PROB, del_prob = SYN_DEL_PROB, ins_prob = SYN_INS_PROB)
    print('Length of synthesized sequence', len(syn_seq))
    print(syn_seq)
    squiggle_array = scrappy.sequence_to_squiggle(syn_seq,rescale=True).data(as_numpy=True)
    raw_data = np.array([])

    if deepSimDwellFlag:
        # for dwell time use deepsimulator since the one provided by scrappie is way too good
        # scrappie dwell gives around 3-4% edit distance while ds dwell gives around 15%
        ds_alpha = deepSimAlpha # 0.1 is default parameter in deepsim
        squiggle_array[:,0] = rep_rvs(squiggle_array.shape[0], ds_alpha)

    for squig in squiggle_array:
            mean = squig[1]
            stdv = squig[2]
            dwell = squig[0]
            raw_data = np.append(raw_data, np.random.laplace(mean, stdv/np.sqrt(2), int(round(dwell))))

    print('Length of raw signal: ', len(raw_data))
    create_fast5(raw_data, fast5_filename)

def read_seq(infile_seq):
    # can handle fasta header, but also works without it
    f = open(infile_seq)
    seq = f.readline().rstrip('\n')
    if seq[0] == '>':
        seq = f.readline().rstrip('\n') 
    f.close()
    len_seq = len(seq)
    print('Length of seq: ', len_seq)
    print(seq)
    return seq

def find_barcode_pos_in_post(trans_filename,fastq_filename,start_barcode,end_barcode):
    '''
    find position of best edit distance match for barcodes in the post matrix
    looks at fastq to find the best match for barcode_start and barcode_end and then finds
    corresponding entries in trans_filename. Returns a tuple (start_pos,end_pos) which represents
    start and end position of actual payload in the post matrix (both inclusive, zero-indexed). 
    One could then slightly extend these or not, depending on what works best. 
    If things fail, return (-1,-1)
    '''
    # load basecalled read from fastq
    with open(fastq_filename,'r') as f:
        _ = f.readline()
        basecall = f.readline().rstrip('\n')
    basecall_len = len(basecall)
    # load entries in trans_filename
    with open(trans_filename,'r') as f:
        trans_arr = [int(l.rstrip('\n')) for l in f.readlines()]
    
    start_barcode_len = len(start_barcode)
    end_barcode_len = len(end_barcode)
    if start_barcode_len + end_barcode_len > basecall_len:
        print('Too short read')
        return (-1,-1)

    start_bc_edit_distance = []
    for i in range(basecall_len//2+1-start_barcode_len): # only search in first half for start barcode
        start_bc_edit_distance.append(distance.levenshtein(start_barcode,basecall[i:i+start_barcode_len]))

    end_bc_edit_distance = []
    for i in range(basecall_len//2,basecall_len-end_barcode_len): # only search in first half for end barcode (so that things don't break if start and end barcodes are same)
        end_bc_edit_distance.append(distance.levenshtein(end_barcode,basecall[i:i+end_barcode_len]))

    # find best match positions
    start_bc_first_base = start_bc_edit_distance.index(min(start_bc_edit_distance))
    end_bc_first_base = basecall_len//2+end_bc_edit_distance.index(min(end_bc_edit_distance))
    start_bc_last_base = start_bc_first_base + start_barcode_len - 1
    start_pos = trans_arr[start_bc_last_base+1]-1
    end_pos = trans_arr[end_bc_first_base-1]-1
    print('basecall_len', basecall_len)
    print('start_bc_last_base', start_bc_last_base)
    print('end_bc_first_base', end_bc_first_base)
    print('start_pos_in_post',start_pos)
    print('end_pos_in_post',end_pos)
    print('basecall',basecall)
    print('start_barcode',start_barcode)
    print('start_bcmatch',basecall[start_bc_first_base:start_bc_first_base+start_barcode_len])
    print('end_barcode',end_barcode)
    print('end_bcmatch',basecall[end_bc_first_base:end_bc_first_base+end_barcode_len])

    assert end_pos >= start_pos
    return (start_pos,end_pos)

def truncate_post_file(old_post_filename, new_post_filename, start_pos, end_pos, bytes_per_blk = 160):
    '''
    Truncate post file to [start_pos,end_pos] and write to new_post_filename. 
    bytes_per_blk is 160 by default
    160 = sizeof(float)*40 (40 entries in transition matrix per timestep for flappie)
    '''
    with open(old_post_filename,'rb') as f:
        data = f.read()
    assert len(data)%bytes_per_blk == 0
    assert end_pos >= start_pos
    assert len(data) >= (end_pos+1)*bytes_per_blk
    with open(new_post_filename,'wb') as f:
        f.write(data[start_pos*bytes_per_blk:(end_pos+1)*bytes_per_blk])
    return
