#!/bin/bash


HDF5_INPUT="/raid/nanopore/shubham/20190629_nanopore_data/raw_signal/raw_signal_0.hdf5"
python3 -u generate_decoded_lists.py \
--hdf_file $HDF5_INPUT \
--out_prefix data/output_3_rescaled_no_trim_default_model/list \
--read_id_file data/read_id_tiny.txt \
--info_file data/output_1/info.txt \
--mem_conv 8 \
--msg_len 100 \
--rate_conv 1 \
--list_size 1 \
--start_barcode CTGGCTCCTCTGTATGTTGGAGAAT \
--end_barcode CCTATATGTACCTCTATCGTAAGTC \
--num_threads 10
