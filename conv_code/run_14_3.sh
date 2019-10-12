#!/bin/bash
for i in {10..39}
do
    python generate_decoded_lists.py \
        --hdf_file /raid/nanopore/shubham/20190629_nanopore_data/raw_signal/raw_signal_5.hdf5 \
        --out_prefix /raid/nanopore/shubham/20190629_nanopore_data/decoded_lists/14_3/list_$i \
        --read_id_file \
            /raid/nanopore/shubham/20190629_nanopore_data/decoded_lists/14_3/read_ids_$i.txt \
        --info_file /raid/nanopore/shubham/20190629_nanopore_data/decoded_lists/14_3/info_$i.txt \
        --mem_conv 14 \
        --msg_len 164 \
        --rate_conv 3 \
        --list_size 4 \
        --start_barcode ATGTTCGGAACGTCAAGACCGAGGA \
        --end_barcode TGGCTCCATTATGCTACAATCACTA \
        &
done
