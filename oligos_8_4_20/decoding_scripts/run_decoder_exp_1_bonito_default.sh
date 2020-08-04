#!/bin/bash

ROOT_PATH=/raid/nanopore/shubham/20200214_nanopore_pool_data/

penalty=0.6
#model_folder=$ROOT_PATH/nanopore_dna_storage/bonito_model/model_train_on_oligo_1_2_4_5_8_9_10_11_on_pretrained_lr_e-6_test_on_oligo_3_global_accuracy_weights_100/

HDF5_INPUT=$ROOT_PATH/data/raw_signal/raw_signal_1.hdf5
DIRNAME=$ROOT_PATH/data/decoded_lists/exp_1_bonito_default_bcp_0.6/
READ_ID_FILE=$ROOT_PATH/data/read_ids/exp_1_read_ids.1000.txt
mkdir -p $DIRNAME
python3 -u generate_decoded_lists.py \
--hdf_file $HDF5_INPUT \
--out_prefix $DIRNAME/list \
--read_id_file $READ_ID_FILE \
--info_file $DIRNAME/info.txt \
--mem_conv 8 \
--msg_len 164 \
--rate_conv 3 \
--list_size 32 \
--start_barcode TCTATCTACTCGTGCTCGCTAGCTG \
--end_barcode TGTCTGCACTGCACTAGTCGCATGT \
--num_threads 4 \
--barcode_extend_penalty $penalty \
#--bonito_model_path $model_folder
# comment out bonito_model_path line for default model
