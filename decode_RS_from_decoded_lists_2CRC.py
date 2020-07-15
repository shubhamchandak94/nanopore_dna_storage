import helper
import os 
import random
import math
import filecmp

# SET THESE PARAMETERS BEFORE RUNNING
NUM_TRIALS = 10
LIST_SIZE = 8
NUM_READS_TOTAL = 10000
NUM_READS_TO_USE = 5500
bytes_per_oligo  = 18
DECODED_LISTS_DIR = '../nanopore_dna_storage_data/decoded_lists/exp_7/'
RS_REDUNDANCY = 0.3
pad = False
ORIGINAL_FILE = '../nanopore_dna_storage_data/encoded_file/data_files.tar.bz2.enc' # for checking if decoding was successful


data_file_size = os.path.getsize(ORIGINAL_FILE)
data_size_padded = math.ceil(data_file_size/bytes_per_oligo)*bytes_per_oligo

(msg_len, num_oligos_data, num_oligos_RS, num_oligos) = helper.compute_parameters(bytes_per_oligo, RS_REDUNDANCY, data_size_padded, pad)

print('NUM_READS_TO_USE:',NUM_READS_TO_USE)
print('list size:',LIST_SIZE)

num_successes = 0

for _ in range(NUM_TRIALS):
    list_read_ids = random.sample(list(range(NUM_READS_TOTAL)),NUM_READS_TO_USE)
    decoded_index_dict = {}

    num_RS_segments = bytes_per_oligo//2
    decoded_index_RS_segments_dict = [{} for i in range(num_RS_segments)]
    # map from index to [[payload_bytes,count]]
    for i in list_read_ids:
        list_file = DECODED_LISTS_DIR+'list_'+str(i)
        if os.path.isfile(list_file):
            with open(list_file) as f:
                decoded_msg_list = [l.rstrip('\n') for l in f.readlines()]
                decoded_msg_list = decoded_msg_list[:LIST_SIZE]
            
            # Output a list of index, pos and payload_bytes
            output_list = helper.decode_list_2CRC_index(decoded_msg_list,bytes_per_oligo,num_oligos,pad)
            for (index, pos, payload_bytes) in output_list:
                if index in decoded_index_dict[pos]:
                    found = False
                    for tup in decoded_index_dict[pos][index]:
                        if tup[0] == payload_bytes:
                            tup[1] += 1
                            found = True
                            break
                        if not found:
                            decoded_index_dict[pos][index].append([payload_bytes,1])
                    decoded_index_dict[pos][index] = sorted(decoded_index_dict[pos][index],key=lambda x: -x[1])
                else:
                    decoded_index_dict[pos][index] = [[payload_bytes,1]]
    
    # Prepare data to be sent to RS decoder
    decoded_list = []
    for segment_id in range(num_RS_segments):
        decoded_list.append([[k, decoded_index_dict[segment_id][k][0][0]] for k in decoded_index_dict[segment_id]])
    
    # Decode each segment separately
    RS_decoded_list = []
    for segment_id in range(num_RS_segments):
        RS_decoded_list.append(helper.RSCode.MainDecoder(decoded_list[segment_id], num_oligos_RS, num_oligos))

    # Combine the segments into a list
    decoded_oligos = [''.join(list(i)) for i in zip(*l)]

    decoded_data = b''.join(decoded_oligos)
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
