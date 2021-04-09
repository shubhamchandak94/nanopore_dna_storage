import helper
import os 

# PARAMETERS TO BE SET BEFORE RUNNING THE CODE
LIST_SIZE = 8
DECODED_LISTS_DIR = '/mnt/ix1/Projects_lite/20180528_HJL_DGI/riskyKmer/shubham/nanopore_dna_storage/data_20210304_MIN_0964/decoded_lists/barcode01/exp_0_bonito_default_bcp_0.6/'
CONV_INPUT_FILE = '/mnt/ix1/Projects_lite/20180528_HJL_DGI/riskyKmer/shubham/nanopore_dna_storage/nanopore_dna_storage/oligos_8_4_20/reads.0.conv_input'
pad = True
bytes_per_oligo = 18

print('list size:',LIST_SIZE)
with open(CONV_INPUT_FILE) as f:
    conv_input_list = [s.rstrip('\n') for s in f.readlines()]
num_oligos = len(conv_input_list)
print('num_oligos',num_oligos)

num_reads = 0
num_correct = 0
num_erasure_CRC_index = 0
num_error_CRC_index = 0
decoded_index_dict = {}
# map from index to [[decoded_msg,count]]

for filename in os.listdir(DECODED_LISTS_DIR):
    if not filename.startswith("list_"):
        continue
    num_reads += 1
    with open(os.path.join(DECODED_LISTS_DIR,filename)) as f:
        decoded_msg_list = [l.rstrip('\n') for l in f.readlines()]
        decoded_msg_list = decoded_msg_list[:LIST_SIZE]
    (index, payload_bytes, decoded_msg) = helper.decode_list_CRC_index(decoded_msg_list,bytes_per_oligo,num_oligos,pad)
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

print('num_reads:',num_reads)
print('num_correct:',num_correct)
print('num_erasure_CRC_index:',num_erasure_CRC_index)
print('num_error_CRC_index:',num_error_CRC_index)
