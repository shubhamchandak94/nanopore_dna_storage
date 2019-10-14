# Nanopore DNA storage
DNA storage for nanopore sequencing using convolutional coding and flappie basecaller

### [Supplementary material](https://github.com/shubhamchandak94/nanopore_dna_storage/blob/master/supplementary_material.pdf)

### Data: https://github.com/shubhamchandak94/nanopore_dna_storage_data

Code tested on Ubuntu 18.04.1
## Download and install
Download:
```
git clone https://github.com/shubhamchandak94/nanopore_dna_storage/
```
Flappie models are stored using git lfs, see instructions in `flappie/` directory README in case of issues.

Install:
```
./install.sh
```
Other than this, you might need to install the following Python3 packages: crc8, distance, fast5_research, h5py, numpy, scipy, scrappy, struct. Also, flappie needs certain dependencies listed in the `flappie/` directory.
