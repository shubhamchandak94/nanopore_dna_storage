#!/bin/bash

ROOT_PATH=/mnt/ix1/Projects_lite/20180528_HJL_DGI/riskyKmer/shubham/nanopore_dna_storage/data_20210304_MIN_0964/

penalty=0.6
#model_folder=$ROOT_PATH/nanopore_dna_storage/bonito_model/model_train_on_oligo_1_2_4_5_8_9_10_11_on_pretrained_lr_e-6_test_on_oligo_3_global_accuracy_weights_100/
REPLICATE=$1

HDF5_INPUT=$ROOT_PATH/raw_signal/$REPLICATE/raw_signal_1.hdf5
DIRNAME=$ROOT_PATH/decoded_lists/$REPLICATE/exp_1_bonito_default_bcp_0.6/
READ_ID_FILE=$ROOT_PATH/read_ids/$REPLICATE/exp_1_read_ids.10000.txt
mkdir -p $DIRNAME
python3 -u generate_decoded_lists.py \
--hdf_file $HDF5_INPUT \
--out_prefix $DIRNAME/list \
--read_id_file $READ_ID_FILE \
--info_file $DIRNAME/info.txt \
--mem_conv 8 \
--msg_len 164 \
--rate_conv 3 \
--list_size 8 \
--start_barcode TCTATCTACTCGTGCTCGCTAGCTG \
--end_barcode TGTCTGCACTGCACTAGTCGCATGT \
--num_threads 6 \
--barcode_extend_penalty $penalty \
#--bonito_model_path $model_folder
# comment out bonito_model_path line for default model
