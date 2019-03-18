import json
import os
import subprocess
import util
import binascii
import struct
import numpy as np

LIST_SIZE = 1
NUM_THREADS = 32
START_BARCODE = 'ATGTTCGGAACGTCAAGACCGAGGA'
END_BARCODE = 'TGGCTCCATTATGCTACAATCACTA'
JSON_FILE_WITH_RAW_DATA = '/raid/nanopore/shubham/20190306_2044_MN19956_FAK54380_79e29d9c/conv_m6_9sync110/merged_rl_250_no_reverse.ids.shuf.20000.raw.json'
RAPTOR_HEADER_FILE = '/raid/nanopore/shubham/Nanopore-DNA-storeage/oligos_2_7_19/code/conv_m6_9sync110_init_32payload_init100101_oligos_raptor_header'
DECODED_FILE = 'conv_m6_9sync110_decoded_l1.bz2'
PATH_TO_CPP_EXEC = "/raid/nanopore/shubham/flappie_new/viterbi/viterbi_nanopore.out"
PATH_TO_FLAPPIE = "/raid/nanopore/shubham/flappie_new/flappie"
NUM_READS_BEFORE_RAPTOR = 1100 # number of unique reads to obtain before calling raptor decoding (should be slightly larger than the optimum to guarantee recovery
MEM_CONV = 6
SYNC_MARKER = '110'
SYNC_MARKER_PERIOD = 9
index_len = 16
crc_len = 16
payload_len = 32
MSG_LEN_AFTER_MARKERS = 97
bin_block_len = index_len + crc_len + payload_len

# parameters for PRP x -> ax+b mod 2^16 to randomize index (16 = index_len)
prp_a = 30009
prp_b = 9901
prp_a_inv = 59657 # calculated using util.modinv in dna_storage

# read raw data JSON
with open(JSON_FILE_WITH_RAW_DATA,'r') as f:
    raw_data_list = json.loads(f.read())

decode_raptor_script = "/raid/nanopore/shubham/dna_storage/rq --debug decode "

f_in = open(RAPTOR_HEADER_FILE,'r')
data = json.loads(f_in.read()) 
f_in.close()

data['symbols'] = [] # for raptor

index_set = set([])
num_attempted = 0 
num_success = 0
for raw_data in raw_data_list:
    num_attempted += 1
    print("num_attempted:",num_attempted)
    rnd = str(np.random.randint(10000000))
    # create fast5 from raw data
    fast5_filename='tmp.'+rnd+'.fast5'
    util.create_fast5(raw_data[1],fast5_filename)

    # call flappie to generate transition posterior table
    post_filename = 'tmp.'+rnd+'.post'
    fastq_filename = 'tmp.'+rnd+'.fastq'
    trans_filename = 'tmp.'+rnd+'.trans'
    subprocess.run([PATH_TO_FLAPPIE, fast5_filename, '--post-output-file', post_filename, '--trans-output-file', trans_filename, '-o',fastq_filename])

    # truncate post according to barcode
    (start_pos, end_pos) = util.find_barcode_pos_in_post(trans_filename,fastq_filename,START_BARCODE,END_BARCODE)
    if start_pos == -1 or end_pos - start_pos + 1 < MEM_CONV+MSG_LEN_AFTER_MARKERS+1:
        print ("Failure in barcode removing.")
        os.remove(fast5_filename)
        os.remove(post_filename)
        os.remove(trans_filename)
        os.remove(fastq_filename)
        continue
    new_post_filename = 'tmp.'+rnd+'.post.new'
    util.truncate_post_file(post_filename, new_post_filename, start_pos, end_pos)
    decoded_filename = 'tmp.'+rnd+'.dec'
    if SYNC_MARKER == '':
        subprocess.run([PATH_TO_CPP_EXEC,'-m','decode_conv','-i',new_post_filename,'-o',decoded_filename,'--msg-len',str(MSG_LEN_AFTER_MARKERS),'-l',str(LIST_SIZE),'-t',str(NUM_THREADS),'--mem-conv',str(MEM_CONV)])
    else:
        subprocess.run([PATH_TO_CPP_EXEC,'-m','decode_conv','-i',new_post_filename,'-o',decoded_filename,'--msg-len',str(MSG_LEN_AFTER_MARKERS),'-l',str(LIST_SIZE),'-t',str(NUM_THREADS),'--mem-conv',str(MEM_CONV),'--sync-marker',str(SYNC_MARKER),'--sync-period',str(SYNC_MARKER_PERIOD)])

    with open(decoded_filename) as f:
        decoded_msg_list = [l.rstrip('\n') for l in f.readlines()]

    for decoded_msg in decoded_msg_list:
        # first remove markers
        if len(SYNC_MARKER) == 0:
            bit_string_with_crc = decoded_msg
        else:
            bit_string_with_crc = ''.join([decoded_msg[i] for i in range(len(decoded_msg)) if (i%SYNC_MARKER_PERIOD >= len(SYNC_MARKER))])
            assert len(bit_string_with_crc) == bin_block_len
        byte_with_crc = binascii.unhexlify(((hex(int(bit_string_with_crc,2)))[2:]).zfill(bin_block_len//4))
        int_crc = struct.unpack('H',byte_with_crc[-(crc_len//8):])[0]
        if binascii.crc_hqx(byte_with_crc[:-(crc_len//8)],0) == int_crc:
            num_success += 1
            print('Success')
            print('num success:',num_success)
            bit_string_with_index = bin(int(binascii.hexlify(byte_with_crc[:-(crc_len//8)]), 16))[2:].zfill(bin_block_len-crc_len)
            index = (prp_a_inv*((int(bit_string_with_index[:index_len],2))-prp_b))%(2**index_len)
            if index not in index_set:
                # new 
                index_set.add(index)
                data['symbols'].append([index,bit_string_with_index[index_len:]])
                print('New index!',index)
                print('num_unique',len(index_set))
            else:
                # already seen
                print('Already seen index',index)
            break
    os.remove(fast5_filename)
    os.remove(post_filename)
    os.remove(decoded_filename)
    os.remove(new_post_filename)
    os.remove(trans_filename)
    os.remove(fastq_filename)
    if len(index_set) == NUM_READS_BEFORE_RAPTOR:
        break
#    if num_attempted == 100:
#        exit(0)

tmpfile = 'tmpfile'+str(np.random.randint(1000000))
f_tmp = open(tmpfile,'w')
f_tmp.write(json.dumps(data, indent=2, separators=(',',': ')))
f_tmp.close()
ret = subprocess.call([decode_raptor_script+" "+tmpfile +" "+DECODED_FILE], shell=True)
os.remove(tmpfile)
