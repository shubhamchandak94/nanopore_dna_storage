#make clean
##make
g++ viterbi/viterbi_nanopore.cpp -std=c++11 -o viterbi/viterbi_nanopore.out -O3 -march=native -Wall -fopenmp -Iviterbi
