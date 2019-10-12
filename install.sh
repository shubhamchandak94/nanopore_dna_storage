cd flappie/
make clean
make
cd ../
cd RSCode_schifra
make clean
make
cd ../
g++ conv_code/viterbi_convolutional_code.cpp -std=c++11 -o conv_code/viterbi_nanopore.out -Wall  -Iviterbi -fopenmp -O3 -march=native
g++ util/read_length_distribution.cpp -std=c++11 -o util/read_length_distribution.out -Wall -O3 -march=native
