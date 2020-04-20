#!/bin/bash


HDF5_INPUT="/raid/nanopore/shubham/20190629_nanopore_data/raw_signal/raw_signal_3.hdf5"
python3 -u generate_decoded_lists.py \
--hdf_file $HDF5_INPUT \
--out_prefix data_guppy/8_3_checkpoint_0/list \
--read_id_file read_id_files_1000/read_id_8_3.txt.100 \
--info_file data_guppy/8_3_checkpoint_0/info.txt \
--mem_conv 8 \
--msg_len 164 \
--rate_conv 3 \
--list_size 8 \
--start_barcode CTGGCGGCCTTGGCCGACTATCTGC \
--end_barcode ACCATGTCGTACAGTCGTTGTAACA \
--num_threads 10 \
--barcode_search_extend_len 10 \
--guppy_model_path /home/kedart/code/taiyaki/taiyaki_exp2/training_over_pretrained/model_checkpoint_00000.checkpoint.json
#--guppy_model_path /home/kedart/code/taiyaki/taiyaki_guppy_default/training_1/model_checkpoint_00010.checkpoint.json
# comment guppy_model_path line to use default model