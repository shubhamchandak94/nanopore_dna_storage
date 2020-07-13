import helper
import os 
import random
import math
import filecmp

# SET THESE PARAMETERS BEFORE RUNNING
NUM_TRIALS = 10
LIST_SIZE = 8
NUM_READS_TOTAL = 10000
NUM_READS_TO_USE = 1500
bytes_per_oligo  = 20
DECODED_LISTS_DIR = 'data_real_exp/11_5_model_train_on_oligo_1_2_4_5_8_9_10_11_on_pretrained_lr_e-6_test_on_oligo_3_global_accuracy_weights_100_0.5_full/'
RS_REDUNDANCY = 0.3
pad = False
ORIGINAL_FILE = '/raid/nanopore/shubham/nanopore_dna_storage_data/encoded_file/data_files.tar.bz2.enc' # for checking if decoding was successful


data_file_size = os.path.getsize(ORIGINAL_FILE)
data_size_padded = math.ceil(data_file_size/bytes_per_oligo)*bytes_per_oligo

(msg_len, num_oligos_data, num_oligos_RS, num_oligos) = helper.compute_parameters(bytes_per_oligo, RS_REDUNDANCY, data_size_padded, pad)

print('NUM_READS_TO_USE:',NUM_READS_TO_USE)
print('list size:',LIST_SIZE)

num_successes = 0

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

print('NUM_TRIALS',NUM_TRIALS)
print('num_successes',num_successes)
