#!/bin/bash

ROOT_PATH=/raid/nanopore/shubham/20200804_nanopore_pool_data/

penalty=0.5
model_folder=/raid/nanopore/shubham/20200214_nanopore_pool_data/nanopore_dna_storage/bonito_model/model_train_on_oligo_1_2_4_5_8_9_10_11_on_pretrained_lr_e-6_test_on_oligo_3_global_accuracy_weights_100/

HDF5_INPUT=$ROOT_PATH/data/raw_signal/raw_signal_4.hdf5
DIRNAME=$ROOT_PATH/data/decoded_lists/exp_4_bonito_trained_bcp_0.5/
READ_ID_FILE=$ROOT_PATH/data/read_ids/exp_4_read_ids.1000.txt
mkdir -p $DIRNAME
python3 -u generate_decoded_lists.py \
--hdf_file $HDF5_INPUT \
--out_prefix $DIRNAME/list \
--read_id_file $READ_ID_FILE \
--info_file $DIRNAME/info.txt \
--mem_conv 8 \
--msg_len 180 \
--rate_conv 5 \
--list_size 8 \
--start_barcode CATCAGCAGTAGAGAGTAGCGCGAT \
--end_barcode TATCATCGACGCTAGCAGTGTCTGC \
--num_threads 2 \
--barcode_extend_penalty $penalty \
--bonito_model_path $model_folder
# comment out bonito_model_path line for default model
