#!/bin/bash
for i in {0..14}
do
    python generate_decoded_lists.py \
        --hdf_file /raid/nanopore/shubham/20190629_nanopore_data/raw_signal/raw_signal_9.hdf5 \
        --out_prefix /raid/nanopore/shubham/20190629_nanopore_data/decoded_lists/11_3_RS_0.2/list_$i \
        --read_id_file \
            /raid/nanopore/shubham/20190629_nanopore_data/decoded_lists/11_3_RS_0.2/read_ids_$i.txt \
        --info_file /raid/nanopore/shubham/20190629_nanopore_data/decoded_lists/11_3_RS_0.2/info_$i.txt \
        --mem_conv 11 \
        --msg_len 164 \
        --rate_conv 3 \
        --list_size 8 \
        --start_barcode GTATACATTCCTTGCCAACATAGTA \
        --end_barcode TGAATCTGGATACGCGTTCCTCAAC \
        &
done
