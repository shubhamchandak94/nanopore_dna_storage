import RS_util
import os 

LIST_SIZE = 1
DECODED_LISTS_DIR = '/raid/nanopore/shubham/20190629_nanopore_data/decoded_lists/8_5/'
CONV_INPUT_FILE = '/raid/nanopore/shubham/20190629_nanopore_data/oligo_files/conv_input_6.txt'
NUM_READS = 1400
pad = False
bytes_per_oligo = 20

f_info = open(DECODED_LISTS_DIR+'info.txt')
with open(CONV_INPUT_FILE) as f:
    conv_input_list = [s.rstrip('\n') for s in f.readlines()]
num_oligos = len(conv_input_list)
print('num_oligos',num_oligos)

num_correct = 0
num_erasure_CRC_index = 0
num_barcode_failure = 0
num_error_CRC_index = 0

for i in range(NUM_READS):
    correct_bit_string = conv_input_list[int(((f_info.readline().split())[1]).split('_')[-1])]
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
            if decoded_msg == correct_bit_string:
                num_correct += 1
            else:
                num_error_CRC_index += 1

print('NUM_READS:',NUM_READS)
print('num_correct:',num_correct)
print('num_barcode_failure:',num_barcode_failure)
print('num_erasure_CRC_index:',num_erasure_CRC_index)
print('num_error_CRC_index:',num_error_CRC_index)
