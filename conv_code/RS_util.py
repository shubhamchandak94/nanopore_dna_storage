import math
import subprocess 
import json
import binascii 
import struct
import numpy as np
import os
import crc8
import filecmp
import util
import sys
PATH_TO_RS_CODE = '/raid/nanopore/shubham/nanopore_dna_storage/RSCode_schifra/'
sys.path.insert(0, PATH_TO_RS_CODE)
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

import RSCode_16bit_fileio as RSCode

# parameters for PRP x -> ax+b mod 2^12 to randomize index (12 = index_len)
prp_a = 1751
prp_b = 2532
prp_a_inv = 3303 # calculated using util.modinv in dna_storage

PATH_TO_VITERBI_NANOPORE = '/raid/nanopore/shubham/flappie_new/viterbi/viterbi_nanopore.out'
PATH_TO_FLAPPIE = "/raid/nanopore/shubham/flappie_new/flappie"

index_len = 12
crc_len = 8

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
    segemented_data_with_RS = RSCode.MainEncoder(segmented_data, num_oligos_RS, path=PATH_TO_RS_CODE)

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

def simulate_and_decode(oligo_file, decoded_data_file,  num_reads, data_file_size, bytes_per_oligo, RS_redundancy, conv_m, conv_r, pad=False, syn_sub_prob = 0.005, syn_del_prob = 0.005, syn_ins_prob = 0.0005, deepsimdwell = False, num_thr = 16, list_size = 1):
    # data_file_size in bytes
    data_size_padded = math.ceil(data_file_size/bytes_per_oligo)*bytes_per_oligo
    (msg_len, num_oligos_data, num_oligos_RS, num_oligos) = compute_parameters(bytes_per_oligo, RS_redundancy, data_size_padded, pad)
    with open(oligo_file) as f:
        oligo_list = [l.rstrip('\n') for l in f.readlines()]
    print('oligo_len',len(oligo_list[0]))
    decoded_dict = {}
    num_success = 0
    num_attempted = 0
    for _ in range(num_reads):
        num_attempted += 1
        print('num_attempted:',num_attempted)
        rnd = str(np.random.randint(10000000))
        oligo = np.random.choice(oligo_list)
        rc = np.random.choice([True, False])
        if rc:
            oligo = util.reverse_complement(oligo)
        fast5_filename='tmp.'+rnd+'.fast5'
        util.simulate_read(oligo, syn_sub_prob, syn_del_prob, syn_ins_prob, fast5_filename, deepSimDwellFlag = deepsimdwell)
        # call flappie to generate transition posterior table
        post_filename = 'tmp.'+rnd+'.post'
        decoded_filename = 'tmp.'+rnd+'.dec'
        subprocess.run([PATH_TO_FLAPPIE, fast5_filename, '--post-output-file', post_filename])
        post_filename = 'tmp.'+rnd+'.post'
        decoded_filename = 'tmp.'+rnd+'.dec'
        subprocess.run([PATH_TO_FLAPPIE, fast5_filename, '--post-output-file', post_filename])
        rc_flag = ''
        if rc:
            rc_flag = '--rc'
        subprocess.run([PATH_TO_VITERBI_NANOPORE,'-m', 'decode','-i',post_filename,'-o',decoded_filename,'--mem-conv',str(conv_m),'--msg-len',str(msg_len),'-l',str(list_size),'-t',str(num_thr),'-r',str(conv_r),rc_flag,'--max-deviation','20'])
        with open(decoded_filename) as f:
            decoded_msg_list = [l.rstrip('\n') for l in f.readlines()]
        for decoded_msg in decoded_msg_list:
            # remove padding, if any
            if pad:
                decoded_msg = decoded_msg[:-1]
            length_with_crc = math.ceil(len(decoded_msg)/8)*8
            bytestring_with_crc = bitstring2bytestring(decoded_msg, length_with_crc)
            crc = crc8.crc8(bytestring_with_crc[:-crc_len//8])
            if crc.digest() == bytestring_with_crc[-crc_len//8:]:
                index_bit_string = bytestring2bitstring(bytestring_with_crc[:math.ceil(index_len/8)], 8*math.ceil(index_len/8))
                index_bit_string = index_bit_string[-index_len:]
                index = (prp_a_inv*((int(index_bit_string,2))-prp_b))%(2**index_len)
                payload_bytes = bitstring2bytestring(decoded_msg[index_len:-crc_len], bytes_per_oligo*8)
                if index < num_oligos:
                    num_success += 1
                    print('Success')
                    print('num success:',num_success)
                    if index not in decoded_dict:
                        decoded_dict[index] = payload_bytes
                        print('New index!',index)
                        print('num_unique',len(decoded_dict))
                    else:
                        print('Already seen index',index)
                    break
                else:
                    print('Index out of range')
            else:
                print('CRC failed')
        os.remove(fast5_filename)
        os.remove(post_filename)
        os.remove(decoded_filename)
    # do RS decoding
    decoded_list = [[k, decoded_dict[k]] for k in decoded_dict]
    print(decoded_list)
    print(num_oligos_RS)
    print(num_oligos)
    RS_decoded_list = RSCode.MainDecoder(decoded_list, num_oligos_RS, num_oligos, path=PATH_TO_RS_CODE)
    assert len(RS_decoded_list) == num_oligos_data
    decoded_data = b''.join(RS_decoded_list)
    decoded_data = decoded_data[:data_file_size]
    with open(decoded_data_file,'wb') as f:
        f.write(decoded_data)
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
# testing encode + simulate_and_decode
# infile = 'myfile_1K'
# infile_size = 1000
# RS_redundancy = 1
# encode(data_file = infile, oligo_file = infile+'.oligos', bytes_per_oligo = 12, RS_redundancy = RS_redundancy, conv_m = 6, conv_r = 1, pad=False)
# simulate_and_decode(oligo_file = infile+'.oligos', decoded_data_file = infile+'.decoded', num_reads = 500, data_file_size = infile_size, bytes_per_oligo = 12, RS_redundancy = RS_redundancy, conv_m = 6, conv_r = 1, pad=False)
# print('filecmp',filecmp.cmp(infile, infile+'.decoded'))
