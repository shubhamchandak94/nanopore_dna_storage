#!/bin/bash
for i in {0..9}
do
    python generate_decoded_lists.py \
        --hdf_file /raid/nanopore/shubham/20190629_nanopore_data/raw_signal/raw_signal_1.hdf5 \
        --out_prefix /raid/nanopore/shubham/20190629_nanopore_data/decoded_lists/11_1/list_$i \
        --read_id_file \
            /raid/nanopore/shubham/20190629_nanopore_data/decoded_lists/11_1/read_ids_$i.txt \
        --info_file /raid/nanopore/shubham/20190629_nanopore_data/decoded_lists/11_1/info_$i.txt \
        --mem_conv 11 \
        --msg_len 100 \
        --rate_conv 1 \
        --list_size 8 \
        --start_barcode TGCGGATGCGGAAGTATGGTCCTCG \
        --end_barcode CACTAGAAGCATGTCGCTATCGAGT \
        &
done
