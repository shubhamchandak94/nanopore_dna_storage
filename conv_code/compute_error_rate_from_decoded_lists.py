import RS_util
import os 

LIST_SIZE = 1
DECODED_LISTS_DIR = '/raid/nanopore/shubham/20190629_nanopore_data/decoded_lists/8_1/'
CONV_INPUT_FILE = '/raid/nanopore/shubham/20190629_nanopore_data/oligo_files/conv_input_0.txt'
NUM_READS = 6000
pad = False
bytes_per_oligo = 10
print('list size:',LIST_SIZE)
with open(CONV_INPUT_FILE) as f:
    conv_input_list = [s.rstrip('\n') for s in f.readlines()]
num_oligos = len(conv_input_list)
print('num_oligos',num_oligos)

num_correct = 0
num_erasure_CRC_index = 0
num_barcode_failure = 0
num_error_CRC_index = 0
decoded_index_dict = {}
# map from index to [[decoded_msg,count]]
for i in range(NUM_READS):
    list_file = DECODED_LISTS_DIR+'list_'+str(i)
    if not os.path.isfile(list_file):
        num_barcode_failure += 1
    else:
        with open(list_file) as f:
            decoded_msg_list = [l.rstrip('\n') for l in f.readlines()]
            decoded_msg_list = decoded_msg_list[:LIST_SIZE]
        (index, payload_bytes, decoded_msg) = RS_util.decode_list_CRC_index(decoded_msg_list,bytes_per_oligo,num_oligos,pad)
        if index == None:
            num_erasure_CRC_index += 1
        else:
            if index in decoded_index_dict:
                found = False
                for tup in decoded_index_dict[index]:
                    if tup[0] == decoded_msg:
                        tup[1] += 1
                        found = True
                        break
                if not found:
                    decoded_index_dict[index].append([decoded_msg,1])
                decoded_index_dict[index] = sorted(decoded_index_dict[index],key=lambda x: -x[1])
            else:
                decoded_index_dict[index] = [[decoded_msg,1]]
            if decoded_msg == conv_input_list[index]:
                num_correct += 1
            else:
                num_error_CRC_index += 1

print('NUM_READS:',NUM_READS)
print('num_correct:',num_correct)
print('num_barcode_failure:',num_barcode_failure)
print('num_erasure_CRC_index:',num_erasure_CRC_index)
print('num_error_CRC_index:',num_error_CRC_index)

# print overall RS stats
num_erasures = num_oligos - len(decoded_index_dict)
num_errors = 0
for index in decoded_index_dict:
    if decoded_index_dict[index][0][0] != conv_input_list[index]:
        num_errors += 1

print('Stats for RS decoding')
print('num_oligos:',num_oligos)
print('num_errors:',num_errors)
print('num_erasures:',num_erasures)
print('num_erasures+2*num_errors:',num_erasures+2*num_errors)