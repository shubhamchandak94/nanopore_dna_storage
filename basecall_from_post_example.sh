#!/bin/bash
set -e
./viterbi/viterbi_basecall.out decode test_files/test_file_post test_files/test_file_basecall.1 > test_files/test_file_move.1
cmp test_files/test_file_basecall.1 test_files/test_file_basecall
cmp test_files/test_file_move.1 test_files/test_file_move
rm test_files/test_file_basecall.1 test_files/test_file_move.1
