make clean
make
g++ viterbi/viterbi_convolutional_code.cpp -std=c++11 -o viterbi/viterbi_nanopore.out -Wall  -Iviterbi -fopenmp -O3 -march=native
