import helper
import os 
import random
import math
import filecmp

# SET THESE PARAMETERS BEFORE RUNNING
NUM_TRIALS = 10
LIST_SIZE = 8
NUM_READS_TOTAL = 10000
NUM_READS_TO_USE_START = 2750 # starting point
NUM_READS_TO_USE_STEP = 250 # step by which we increase num reads at each step
bytes_per_oligo = 18
DECODED_LISTS_DIR = '/mnt/ix1/Projects_lite/20180528_HJL_DGI/riskyKmer/shubham/nanopore_dna_storage/data_20210304_MIN_0964/decoded_lists/barcode01/exp_0_bonito_default_bcp_0.6/'
RS_REDUNDANCY = 0.25
pad = True
ORIGINAL_FILE = '/mnt/ix1/Projects_lite/20180528_HJL_DGI/riskyKmer/shubham/nanopore_dna_storage/nanopore_dna_storage/oligos_8_4_20/data_files.tar.bz2.enc.1' # for checking if decoding was successful

data_file_size = os.path.getsize(ORIGINAL_FILE)
data_size_padded = math.ceil(data_file_size/bytes_per_oligo)*bytes_per_oligo

(msg_len, num_oligos_data, num_oligos_RS, num_oligos) = helper.compute_parameters(bytes_per_oligo, RS_REDUNDANCY, data_size_padded, pad)

print('list size:',LIST_SIZE)

NUM_READS_TO_USE = NUM_READS_TO_USE_START
while True:
    print('NUM_READS_TO_USE:',NUM_READS_TO_USE)
    num_successes = 0
    random.seed(999)
    for _ in range(NUM_TRIALS):
        list_read_ids = random.sample(list(range(NUM_READS_TOTAL)),NUM_READS_TO_USE)
        decoded_index_dict = {}
        # map from index to [[payload_bytes,count]]
        for i in list_read_ids:
            list_file = DECODED_LISTS_DIR+'list_'+str(i)
            if os.path.isfile(list_file):
                with open(list_file) as f:
                     decoded_msg_list = [l.rstrip('\n') for l in f.readlines()]
                decoded_msg_list = decoded_msg_list[:LIST_SIZE]
                (index, payload_bytes, decoded_msg) = helper.decode_list_CRC_index(decoded_msg_list,bytes_per_oligo,num_oligos,pad)
                if index is not None:
                    if index in decoded_index_dict:
                        found = False
                        for tup in decoded_index_dict[index]:
                            if tup[0] == payload_bytes:
                                tup[1] += 1
                                found = True
                                break
                        if not found:
                            decoded_index_dict[index].append([payload_bytes,1])
                        decoded_index_dict[index] = sorted(decoded_index_dict[index],key=lambda x: -x[1])
                    else:
                        decoded_index_dict[index] = [[payload_bytes,1]]
        decoded_list = [[k, decoded_index_dict[k][0][0]] for k in decoded_index_dict]
        RS_decoded_list = helper.RSCode.MainDecoder(decoded_list, num_oligos_RS, num_oligos)
        decoded_data = b''.join(RS_decoded_list)
        decoded_data = decoded_data[:data_file_size]
        decoded_data_file = 'tmpfile.'+str(random.randint(0,1000000))
        with open(decoded_data_file,'wb') as f:
            f.write(decoded_data)
        if filecmp.cmp(ORIGINAL_FILE, decoded_data_file):
            num_successes += 1
            print('Success')
        else:
            print('Failure')
            os.remove(decoded_data_file)
            break
        os.remove(decoded_data_file)

    if num_successes == NUM_TRIALS:
        print('min number of reads:',NUM_READS_TO_USE)
        break
    else:
        NUM_READS_TO_USE += NUM_READS_TO_USE_STEP
        if NUM_READS_TO_USE > NUM_READS_TOTAL:
            break
