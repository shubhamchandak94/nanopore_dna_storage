#!/bin/bash
for i in {0..9}
do
    python generate_decoded_lists.py \
        --hdf_file /raid/nanopore/shubham/20190629_nanopore_data/raw_signal/raw_signal_2.hdf5 \
        --out_prefix /raid/nanopore/shubham/20190629_nanopore_data/decoded_lists/14_1/list_$i \
        --read_id_file \
            /raid/nanopore/shubham/20190629_nanopore_data/decoded_lists/14_1/read_ids_$i.txt \
        --info_file /raid/nanopore/shubham/20190629_nanopore_data/decoded_lists/14_1/info_$i.txt \
        --mem_conv 14 \
        --msg_len 100 \
        --rate_conv 1 \
        --list_size 4 \
        --start_barcode AGTAACGCCTATTGATAACGAAGCA \
        --end_barcode TAACCTTCGCTGCTAGGAACTGTCT \
        &
done
