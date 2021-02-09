import helper
import os 
import random
import math
import filecmp

# SET THESE PARAMETERS BEFORE RUNNING
NUM_TRIALS = 10
LIST_SIZE = 8
NUM_READS_TOTAL = 8949
NUM_READS_TO_USE_START = 5500 # starting point
NUM_READS_TO_USE_STEP = 250 # step by which we increase num reads at each step
bytes_per_oligo  = 20
DECODED_LISTS_DIR = '/mnt/ix1/Projects_lite/20180528_HJL_DGI/riskyKmer/shubham/nanopore_dna_storage/data_20201104_MIN_0923/decoded_lists/barcode22/exp_12_bonito_default_bcp_0.6/'
RS_REDUNDANCY = 0.25
pad = True
ORIGINAL_FILE = '/mnt/ix1/Projects_lite/20180528_HJL_DGI/riskyKmer/shubham/nanopore_dna_storage/nanopore_dna_storage/oligos_8_4_20/data_files.tar.bz2.enc.2' # for checking if decoding was successful

data_file_size = os.path.getsize(ORIGINAL_FILE)
data_size_padded = math.ceil(data_file_size/bytes_per_oligo)*bytes_per_oligo

(msg_len, num_oligos_data, num_oligos_RS, num_oligos) = helper.compute_parameters(bytes_per_oligo, RS_REDUNDANCY, data_size_padded, pad)
crc_len = 8
msg_len += crc_len # 2 CRC

print('list size:',LIST_SIZE)

NUM_READS_TO_USE = NUM_READS_TO_USE_START
while True:
    print('NUM_READS_TO_USE:',NUM_READS_TO_USE)
    num_successes = 0
    random.seed(999)
    for _ in range(NUM_TRIALS):
        list_read_ids = random.sample(list(range(NUM_READS_TOTAL)),NUM_READS_TO_USE)
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
                    if index in decoded_index_RS_segments_dict[pos]:
                        found = False
                        for tup in decoded_index_RS_segments_dict[pos][index]:
                            if tup[0] == payload_bytes:
                                tup[1] += 1
                                found = True
                                break
                            if not found:
                                decoded_index_RS_segments_dict[pos][index].append([payload_bytes,1])
                        decoded_index_RS_segments_dict[pos][index] = sorted(decoded_index_RS_segments_dict[pos][index],key=lambda x: -x[1])
                    else:
                        decoded_index_RS_segments_dict[pos][index] = [[payload_bytes,1]]
        
        # Prepare data to be sent to RS decoder
        decoded_list = []
        for segment_id in range(num_RS_segments):
            decoded_list.append([[k, decoded_index_RS_segments_dict[segment_id][k][0][0]] for k in decoded_index_RS_segments_dict[segment_id]])
        
        # Decode each segment separately
        RS_decoded_list = []
        for segment_id in range(num_RS_segments):
            RS_decoded_list.append(helper.RSCode.MainDecoder(decoded_list[segment_id], num_oligos_RS, num_oligos))
    
        # Combine the segments into a list
        decoded_oligos = [b''.join(list(i)) for i in zip(*RS_decoded_list)]
    
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
            break
        os.remove(decoded_data_file)

    if num_successes == NUM_TRIALS:
        print('min number of reads:',NUM_READS_TO_USE)
        break
    else:
        NUM_READS_TO_USE += NUM_READS_TO_USE_STEP
        if NUM_READS_TO_USE > NUM_READS_TOTAL:
            break
