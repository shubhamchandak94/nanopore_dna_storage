#!/bin/bash


HDF5_INPUT="/raid/nanopore/shubham/20190629_nanopore_data/raw_signal/raw_signal_3.hdf5"
python3 -u generate_decoded_lists.py \
--hdf_file $HDF5_INPUT \
--out_prefix data_real_exp/8_3_default_100_approach_1/list \
--read_id_file data_real_exp/read_id_lists/read_id_8_3.txt.100 \
--info_file data_real_exp/8_3_default_100_approach_1/info.txt \
--mem_conv 8 \
--msg_len 164 \
--rate_conv 3 \
--list_size 8 \
--start_barcode CTGGCGGCCTTGGCCGACTATCTGC \
--end_barcode ACCATGTCGTACAGTCGTTGTAACA \
--num_threads 10 \
--barcode_search_extend_len 10 \
