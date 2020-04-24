#!/bin/bash


HDF5_INPUT="/raid/nanopore/shubham/20190629_nanopore_data/raw_signal/raw_signal_6.hdf5"
python3 -u generate_decoded_lists.py \
--hdf_file $HDF5_INPUT \
--out_prefix data_real_exp/8_5_default_100_approach_1/list \
--read_id_file data_real_exp/read_id_lists/read_id_8_5.txt.100 \
--info_file data_real_exp/8_5_default_100_approach_1/info.txt \
--mem_conv 8 \
--msg_len 180 \
--rate_conv 5 \
--list_size 8 \
--start_barcode GCTAGTACGCGAACAGAGTGCAGTA \
--end_barcode ACAGATGCAGTAATTCTCACGAACT \
--num_threads 10 \
--barcode_search_extend_len 10 \
