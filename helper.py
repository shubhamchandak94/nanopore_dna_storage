import sys
import numpy as np
import scrappy
from uuid import uuid4
import scipy.stats as st
import subprocess
from fast5_research import Fast5
import distance
import math
import json
import binascii 
import struct
import os
import crc8
import filecmp
import h5py

REPO_PATH = os.path.dirname(os.path.realpath(__file__))+'/'
PATH_TO_RS_CODE = REPO_PATH+'RSCode_schifra/'
PATH_TO_VITERBI_NANOPORE = REPO_PATH+'viterbi/viterbi_nanopore.out'
sys.path.insert(0, PATH_TO_RS_CODE)
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

import RSCode_16bit_fileio as RSCode

# parameters for PRP x -> ax+b mod 2^12 to randomize index (12 = index_len)
prp_a = 1751
prp_b = 2532
prp_a_inv = 3303 # calculated using util.modinv in dna_storage
index_len = 12
crc_len = 8

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

# from deepsimulator (different dwell time distribution - higher variance)
def rep_rvs(size,a):
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
    if raw_data.dtype == np.int16:
        raw_data_binned = raw_data
    else:
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
        'read_id': 'test',
        'scaling_used': 1,
        'median_before': 0,
    }
    tracking_id = {
        'exp_start_time': '1970-01-01T00:00:00Z',
        'run_id': 'test',
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

def find_barcode_pos_in_post(trans_filename,fastq_filename,start_barcode,end_barcode,extend_len=0, extend_penalty=1.0):
    '''
    find position of best edit distance match for barcodes in the post matrix
    looks at fastq to find the best match for barcode_start and barcode_end and then finds
    corresponding entries in trans_filename. Returns a tuple (start_pos,end_pos) which represents
    start and end position of actual payload in the post matrix (both inclusive, zero-indexed). 
    One could then slightly extend these or not, depending on what works best. 
    If things fail, return (-1,-1)
    extend_len: extra length at start and end to search for barcode match (useful if basecaller cuts
    a bit of the barcode, e.g., with guppy)
    extend_penalty: float between 0.0 and 1.0 telling the penalty we impose per base on extend operation (i.e., barcode hanging off the sides). Setting to 1 just means we compute edit distance as it is (this can penalize a bit much and miss perfect match of say 10 bases out of 25 base barcode). Setting to 0 is bad because then there is chance that the empty match (completely hanging off) is selected.
    '''
    assert extend_len >= 0
    assert 0.0 <= extend_penalty <= 1.0
    one_minus_extend_penalty = 1.0-extend_penalty
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
        return (-1,-1,np.inf,np.inf)

    start_bc_edit_distance = []
    for i in range(-extend_len,0):
        start_bc_edit_distance.append(distance.levenshtein(start_barcode,basecall[:start_barcode_len+i])+i*one_minus_extend_penalty)
    for i in range(basecall_len-start_barcode_len):
        start_bc_edit_distance.append(distance.levenshtein(start_barcode,basecall[i:i+start_barcode_len]))

    # find best match positions
    start_bc_first_base = start_bc_edit_distance.index(min(start_bc_edit_distance))-extend_len
    start_bc_last_base = start_bc_first_base + start_barcode_len - 1

    end_bc_edit_distance = []
    for i in range(start_bc_last_base,basecall_len-end_barcode_len+1): # end_bc must be after start_bc 
        end_bc_edit_distance.append(distance.levenshtein(end_barcode,basecall[i:i+end_barcode_len]))
    for i in range(1,extend_len+1):
        end_bc_edit_distance.append(distance.levenshtein(end_barcode,basecall[basecall_len-end_barcode_len+i:basecall_len])-i*one_minus_extend_penalty)

    # find best match positions
    # need sanity check in case range(start_bc_last_base,basecall_len-end_barcode_len) is empty
    if start_bc_last_base >= basecall_len-end_barcode_len:
      end_bc_first_base = basecall_len-end_barcode_len + 1 + end_bc_edit_distance.index(min(end_bc_edit_distance))
    else:
      end_bc_first_base = start_bc_last_base+end_bc_edit_distance.index(min(end_bc_edit_distance))
    print('basecall_len', basecall_len)
    print('start_bc_last_base', start_bc_last_base)
    print('end_bc_first_base', end_bc_first_base)
    
    start_pos = trans_arr[start_bc_last_base+1]-1
    end_pos = trans_arr[end_bc_first_base-1]-1
    print('start_pos_in_post',start_pos)
    print('end_pos_in_post',end_pos)
    print('basecall',basecall)
    print('start_barcode',start_barcode)
    print('start_bcmatch',basecall[max(0,start_bc_first_base):start_bc_first_base+start_barcode_len])
    print('end_barcode',end_barcode)
    print('end_bcmatch',basecall[end_bc_first_base:end_bc_first_base+end_barcode_len])

    if end_pos < start_pos:
        print('Barcode removal failure')
        return (-1,-1,np.inf,np.inf)
    return (start_pos,end_pos,min(start_bc_edit_distance),min(end_bc_edit_distance))

def truncate_post_file(old_post_filename, new_post_filename, start_pos, end_pos, bytes_per_blk = 20):
    '''
    Truncate post file to [start_pos,end_pos] and write to new_post_filename. 
    bytes_per_blk is 20 by default
    20 = sizeof(float)*5 (5 entries in CTC posterior per timestep for bonito)
    '''
    with open(old_post_filename,'rb') as f:
        data = f.read()
    assert len(data)%bytes_per_blk == 0
    assert end_pos >= start_pos
    assert len(data) >= (end_pos+1)*bytes_per_blk
    with open(new_post_filename,'wb') as f:
        f.write(data[start_pos*bytes_per_blk:(end_pos+1)*bytes_per_blk])
    return

# below function from https://codereview.stackexchange.com/questions/151329/reverse-complement-of-a-dna-string
def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return ''.join([complement[base] for base in dna[::-1]])

def encode(data_file, oligo_file, bytes_per_oligo, RS_redundancy, conv_m, conv_r, pad=False):
    # single bit padding might be needed depending on puncturing pattern for successful encoding
    # RS_redundancy is a float representing the additional percent redundancy added (e.g. 0.1 for 10%)
    # conv_r can be 1, 2, 3, 4, 5, 7
    # conv_m can be 6, 8, 11, 14
    # bytes_per_oligo is the number of message bytes per oligo. It should be a multiple of 2.

    assert bytes_per_oligo%2 == 0
    assert conv_m in [6,8,11,14]
    assert conv_r in [1,2,3,4,5,7]

    with open(data_file,'rb') as f:
        data = f.read()
    
    # pad data to multiple of bytes_per_oligo
    data_size = len(data)
    data_size_padded = math.ceil(data_size/bytes_per_oligo)*bytes_per_oligo
    (msg_len, num_oligos_data, num_oligos_RS, num_oligos) = compute_parameters(bytes_per_oligo, RS_redundancy, data_size_padded, pad)
    data_padded = data.ljust(data_size_padded,b'0')
    segmented_data = [data_padded[i*bytes_per_oligo:(i+1)*bytes_per_oligo] for i in range(num_oligos_data)]
    segemented_data_with_RS = RSCode.MainEncoder(segmented_data, num_oligos_RS)

    conv_input_file = oligo_file+'.conv_input'
    with open(conv_input_file, 'w') as f:
        # attach index, CRC and pad to each oligo and write to f
        for index, oligo in enumerate(segemented_data_with_RS):
            index_prp = (prp_a*index+prp_b)%(2**index_len)
            bin_index_string = bin(index_prp)[2:].zfill(index_len)
            index_bytes = bitstring2bytestring(bin_index_string, 8*math.ceil(index_len/8))
            crc = crc8.crc8(index_bytes+oligo)
            bit_string_oligo = bin_index_string + bytestring2bitstring(oligo+crc.digest(),8*bytes_per_oligo+crc_len)
            if pad:
                bit_string_oligo = bit_string_oligo + '0'
            f.write(bit_string_oligo + '\n')
    
    # apply convolutional encoding to each oligo
    subprocess.run([PATH_TO_VITERBI_NANOPORE,'-m', 'encode','-i',conv_input_file,'-o',oligo_file,'--mem-conv',str(conv_m),'--msg-len',str(msg_len),'-r',str(conv_r)])
    
    with open(oligo_file) as f:
        oligo_len = len(f.readline().rstrip('\n'))
        print('oligo_len',oligo_len)
    print('writing rate (bits per base):', data_size*8/(oligo_len*num_oligos))
    return 


def rotate_left(input_str, rot_val):
    rot_val = rot_val % len(input_str)
    out1 = input_str[rot_val:] 
    out2 = input_str[:rot_val]
    return out1+out2

def rotate_right(input_str, rot_val):
    rot_val = rot_val % len(input_str)
    out1 = input_str[-rot_val:] 
    out2 = input_str[:-rot_val]
    return out1+out2


def encode_2crc(data_file, oligo_file, bytes_per_oligo, RS_redundancy, conv_m, conv_r, pad=False):
    # single bit padding might be needed depending on puncturing pattern for successful encoding
    # RS_redundancy is a float representing the additional percent redundancy added (e.g. 0.1 for 10%)
    # conv_r can be 1, 2, 3, 4, 5, 7
    # conv_m can be 6, 8, 11, 14
    # bytes_per_oligo is the number of message bytes per oligo. It should be a multiple of 2.

    assert bytes_per_oligo%2 == 0
    assert conv_m in [6,8,11,14]
    assert conv_r in [1,2,3,4,5,7]

    with open(data_file,'rb') as f:
        data = f.read()
    
    # pad data to multiple of bytes_per_oligo
    data_size = len(data)
    data_size_padded = math.ceil(data_size/bytes_per_oligo)*bytes_per_oligo
    (msg_len, num_oligos_data, num_oligos_RS, num_oligos) = compute_parameters(bytes_per_oligo, RS_redundancy, data_size_padded, pad)
    msg_len += crc_len
    data_padded = data.ljust(data_size_padded,b'0')
    segmented_data = [data_padded[i*bytes_per_oligo:(i+1)*bytes_per_oligo] for i in range(num_oligos_data)]
    segemented_data_with_RS = RSCode.MainEncoder(segmented_data, num_oligos_RS)

    conv_input_file = oligo_file+'.conv_input'
    with open(conv_input_file, 'w') as f:
        # attach index, CRC and pad to each oligo and write to f
        for index, oligo in enumerate(segemented_data_with_RS):
            index_prp = (prp_a*index+prp_b)%(2**index_len)
            bin_index_string = bin(index_prp)[2:].zfill(index_len)
            index_bytes = bitstring2bytestring(bin_index_string, 8*math.ceil(index_len/8))

            # Rotate the oligo
            # RS code operates on 2 bytes at a time, so the shift is proportional
            oligo = rotate_left(oligo, index_prp*2) 

            # Generate 2 CRCs with first half and second half of the oligo
            oligo_len = len(oligo)
            oligo1 = oligo[:2*(oligo_len//4)]
            oligo2 = oligo[2*(oligo_len//4):]
            crc1 = crc8.crc8(index_bytes+oligo1)
            crc2 = crc8.crc8(index_bytes+oligo2)
            bit_string_oligo = bin_index_string + bytestring2bitstring(crc1.digest()+oligo+crc2.digest(),8*bytes_per_oligo+crc_len+crc_len)
            if pad:
                bit_string_oligo = bit_string_oligo + '0'
            f.write(bit_string_oligo + '\n')
    
    # apply convolutional encoding to each oligo
    print('msg_len',msg_len)
    print('conv_m',conv_m)
    print('conv_r',conv_r)
    print('conv_input_file',conv_input_file)
    subprocess.run([PATH_TO_VITERBI_NANOPORE,'-m', 'encode','-i',conv_input_file,'-o',oligo_file,'--mem-conv',str(conv_m),'--msg-len',str(msg_len),'-r',str(conv_r)])
    
    with open(oligo_file) as f:
        oligo_len = len(f.readline().rstrip('\n'))
        print('oligo_len',oligo_len)
    print('writing rate (bits per base):', data_size*8/(oligo_len*num_oligos))
    return 

def compute_parameters(bytes_per_oligo, RS_redundancy, data_size_padded, pad):
    msg_len = index_len + crc_len + 8*bytes_per_oligo + pad
    assert data_size_padded%bytes_per_oligo == 0
    num_oligos_data = data_size_padded//bytes_per_oligo 
    num_oligos_RS = int(num_oligos_data*RS_redundancy) 
    num_oligos = num_oligos_data + num_oligos_RS
    print('msg_len', msg_len)
    print('num_oligos_data',num_oligos_data)
    print('num_oligos_RS',num_oligos_RS)
    print('num_oligos',num_oligos)
    return (msg_len, num_oligos_data, num_oligos_RS, num_oligos)

def bitstring2bytestring(bitstring, bitstring_len):
    return binascii.unhexlify(((hex(int(bitstring,2)))[2:]).zfill(bitstring_len//4))

def bytestring2bitstring(bytestring, bitstring_len):
    return bin(int(binascii.hexlify(bytestring), 16))[2:].zfill(bitstring_len)

def decode_list_CRC_index(decoded_msg_list, bytes_per_oligo, num_oligos, pad):
    for decoded_msg_ in decoded_msg_list:
        # remove padding, if any
        if pad:
            decoded_msg = decoded_msg_[:-1]
        else:
            decoded_msg = decoded_msg_
        length_with_crc = math.ceil(len(decoded_msg)/8)*8
        bytestring_with_crc = bitstring2bytestring(decoded_msg, length_with_crc)
        crc = crc8.crc8(bytestring_with_crc[:-crc_len//8])
        if crc.digest() == bytestring_with_crc[-crc_len//8:]:
            index_bit_string = bytestring2bitstring(bytestring_with_crc[:math.ceil(index_len/8)], 8*math.ceil(index_len/8))
            index_bit_string = index_bit_string[-index_len:]
            index = (prp_a_inv*((int(index_bit_string,2))-prp_b))%(2**index_len)
            payload_bytes = bitstring2bytestring(decoded_msg[index_len:-crc_len], bytes_per_oligo*8)
            if index < num_oligos:
                return (index, payload_bytes, decoded_msg_)
    return (None, None, None)


def get_output_list(decoded_oligo1, decoded_oligo2, decoded_index, oligo_len):
    
    if decoded_index is None:
        return []
    
    # Get index
    index_bit_string = bytestring2bitstring(decoded_index,8*math.ceil(index_len/8))
    index_bit_string = index_bit_string[-index_len:]
    index_prp = int(index_bit_string,2)
    index = (prp_a_inv*(index_prp-prp_b))%(2**index_len)
    oligo_RS_segments = []
    oligo1_len = (2*(oligo_len//4))
    oligo2_len = oligo_len - oligo1_len


    if decoded_oligo1 is None:
        oligo_RS_segments.extend([None for i in range(oligo1_len//2)])
    else:
        oligo_RS_segments.extend([decoded_oligo1[i:i+2] for i in range(0, len(decoded_oligo1), 2)])

    if decoded_oligo2 is None:
        oligo_RS_segments.extend([None for i in range(oligo2_len//2)])
    else:
        oligo_RS_segments.extend([decoded_oligo2[i:i+2] for i in range(0, len(decoded_oligo2), 2)])


    # rotate back the oligos
    oligo_RS_segments_rotated = rotate_right(oligo_RS_segments, index_prp)
    # generate output lists
    output_list = []
    for pos, segment in enumerate(oligo_RS_segments_rotated):
        if segment is not None:
            output_list.append((index, pos, segment))

    return output_list


def decode_list_2CRC_index(decoded_msg_list, bytes_per_oligo, num_oligos, pad):
    decoded_oligo1 = None
    decoded_oligo2 = None
    decoded_index = None
    for decoded_msg_ in decoded_msg_list:
        # remove padding, if any
        if pad:
            decoded_msg = decoded_msg_[:-1]
        else:
            decoded_msg = decoded_msg_

        length_with_crc = math.ceil(len(decoded_msg)/8)*8
        bytestring_with_crc = bitstring2bytestring(decoded_msg, length_with_crc)
        index_bytes = bytestring_with_crc[:math.ceil(index_len/8)]
        oligo_len = length_with_crc//8 - 2*crc_len//8 - math.ceil(index_len/8)
        oligo1_len = (2*(oligo_len//4))
        crc_len_in_bytes = crc_len//8
        index_crc_len = len(index_bytes) + crc_len_in_bytes
        crc1_bytes = bytestring_with_crc[len(index_bytes): index_crc_len]
        oligo1 = bytestring_with_crc[index_crc_len: oligo1_len+index_crc_len ]
        oligo2 = bytestring_with_crc[index_crc_len + oligo1_len: -crc_len_in_bytes]
        crc2_bytes = bytestring_with_crc[-crc_len_in_bytes:]

        index_oligo1 = index_bytes + oligo1
        index_oligo2 = index_bytes + oligo2
        crc1 = crc8.crc8(index_oligo1) 
        crc2 = crc8.crc8(index_oligo2)
       
        index_bit_string = bytestring2bitstring(index_bytes,8*math.ceil(index_len/8))
        index_bit_string = index_bit_string[-index_len:]
        index_prp = int(index_bit_string,2)
        index = (prp_a_inv*(index_prp-prp_b))%(2**index_len)
        # If both the crcs match, extract the index etc.
        if crc1.digest() == crc1_bytes and crc2.digest() == crc2_bytes:
            if index < num_oligos:
                decoded_oligo1 = oligo1
                decoded_oligo2 = oligo2
                decoded_index = index_bytes
                break
        
         
        # Try to extract the first oligo
        if crc1.digest() == crc1_bytes:
            if (index < num_oligos) and (decoded_oligo1 is None):
                decoded_oligo1 = oligo1
                decoded_index = index_bytes

        # Try to extract the first oligo
        if crc2.digest() == crc2_bytes:
            if (index < num_oligos) and (decoded_oligo2 is None):
                decoded_oligo2 = oligo2
                decoded_index = index_bytes

        
    # Create segments from the oligo_output
    output_list = get_output_list(decoded_oligo1, decoded_oligo2, decoded_index, oligo_len)

    return output_list


########
# Bonito helper functions
########

def bonito_basecall_generate_move(post_file, fastq_file, trans_file):
    '''
    perform CTC decoding (greedy) to generate fastq file and a file with 
    the time steps when each new base starts.
    The reason for working with greedy is that it allows us to make a trans 
    (move) matrix which is not possible with beam search.
    Note that the decoder still uses beam search later on, this is just for
    the barcode removal
    '''
    labels = [ "N", "A", "C", "G", "T",]
    posteriors = np.fromfile(post_file, dtype=np.float32)
    posteriors = np.reshape(posteriors, (-1,5))
    basecall = ""
    trans = [] # times when new base arrives
    prev_base = 'N'
    for i in range(posteriors.shape[0]):
        next_base = labels[np.argmax(posteriors[i,:])]
        if next_base != prev_base:
            prev_base = next_base
            if next_base != 'N':
                trans.append(i)
                basecall = basecall + next_base
    with open(fastq_file,'w') as f:
        f.write('@read\n')
        f.write(basecall+'\n')
        f.write('+\n')
        f.write('#'*len(basecall)+'\n')
    with open(trans_file,'w') as f:
        f.write(''.join([str(val)+'\n' for val in trans]))
    return

