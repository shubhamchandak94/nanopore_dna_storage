cd RSCode_schifra
make clean
make
cd ../
g++ viterbi/viterbi_convolutional_code.cpp -std=c++11 -o viterbi/viterbi_nanopore.out -Wall  -Iviterbi -fopenmp -O3 -march=native
g++ viterbi/viterbi_basecall.cpp -std=c++11 -o viterbi/viterbi_basecall.out -Wall -fopenmp -O3 -march=native
g++ util/read_length_distribution.cpp -std=c++11 -o util/read_length_distribution.out -Wall -O3 -march=native
