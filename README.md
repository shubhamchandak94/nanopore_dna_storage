# Nanopore DNA storage
DNA storage for nanopore sequencing using convolutional coding and basecaller-decoder integration

### [Supplementary material](https://github.com/shubhamchandak94/nanopore_dna_storage/blob/master/supplementary_material.pdf)

### Data: https://github.com/shubhamchandak94/nanopore_dna_storage_data

Code tested on Ubuntu 18.04.1 and Python3.
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
Other than this, you might need to install the following Python3 packages: crc8, distance, fast5_research, h5py, numpy, scipy, scrappy, struct (see `install_python_packages.sh`). Also, flappie needs certain dependencies listed in the `flappie/` directory.

## General instructions
In many of the scripts, you need to set the path for the corresponding data directories as well as the encoding parameters. 
The current paths assume that the data is stored in `../nanopore_dna_storage_data/`

## Parameters for experiments
The file `encode_experiments.py` was used for generating the oligos and contains the parameters required for the decoding.

## Convolutional code decoding of raw signal data
First, generate a list of read ids to decode using the script `util/generate_read_ids.py` (changing the parameters as needed). Then, run `generate_decoded_lists.py` with the relevant parameters. Below is an example execution for m=11, r=5/6 and list size 8.
```
python3 generate_decoded_lists.py \
--hdf_file ../nanopore_dna_storage_data/raw_signal/raw_signal_7.hdf5 \
--out_prefix ../nanopore_dna_storage_data/decoded_lists/exp_7/list \
--read_id_file ../nanopore_dna_storage_data/decoded_lists/exp_7/read_ids.txt \
--info_file ../nanopore_dna_storage_data/decoded_lists/exp_7/info.txt \
--mem_conv 11 \
--msg_len 180 \
--rate_conv 5 \
--list_size 8 \
--start_barcode CACCTGTGCTGCGTCAGGCTGTGTC \
--end_barcode GCTGTCCGTTCCGCATTGACACGGC
```
Note that the `rate_conv` parameter is set to 1, 2, 3, 4, 5 and 7 for convolutional code rates of 1/2, 2/3, 3/4, 4/5, 5/6 and 7/8, respectively. The `msg_len` parameter is the length of the binary input to the convolutional code encoder. This writes the decoded lists to files named `list_1, list_2, ...` in the directory `../nanopore_dna_storage_data/decoded_lists/exp_7/` and also generates an `info.txt` file listing the decoded reads.
The parameters for the different experiments can be found in the files `encode_experiments.py` and `logs/encoding_log`.

## Computing error rates of convolutional code decoding
To compute the number of errors (due to no CRC match being found and due to incorrect CRC match) from the decoded lists, use `compute_error_rate_from_decoded_lists.py` after setting the parameters.

## Performing RS decoding of the lists
To perform RS decoding from the lists, use `decode_RS_from_decoded_lists.py` after setting the parameters. `NUM_TRIALS` denotes the number of decoding trials performed. Each trial involves a subsampling of `NUM_READS_TO_USE` reads from the set of decoded lists of size `NUM_READS_TOTAL`, an attempt at RS decoding and comparison with the original encoded file to ensure successful decoding.

## Running simulations
The script `simulator.py` can be used to perform simulations to test various parameters for the convolutional code. The simulation of the raw signal is performed using the [scrappie](https://github.com/nanoporetech/scrappie) simulator with an optional mode to use dwell time distribution from [DeepSimulator](https://github.com/lykaust15/DeepSimulator). An example execution is shown below:
```
python3 simulator.py \
--num_trials 100 \
--list_size 8 \
--mem_conv 11 \
--rate 5 \
--msg_len 180 \
--deepsimdwell False \
--reversecomp False \
--syn_sub_prob 0.004 \
--syn_del_prob 0.0085 \
--syn_ins_prob 0.0005
```
The `mem_conv` parameter can be set to 6, 8, 11 or 14. The `rate` parameter can be set to 1, 2, 3, 4, 5 and 7 for convolutional code rates of 1/2, 2/3, 3/4, 4/5, 5/6 and 7/8, respectively. The `msg_len` parameter decides the length of the binary input to the convolutional encoder, and should be set depending upon the desired oligo length. Setting `deepsimdwell` to `True` uses the dwell time distribution from DeepSim which has higher variance (based on our experience, keeping this `False` gives results closer to reality). `reversecomp` can be set to True to simulate reverse complemented read decoding. Finally, the parameters `syn_sub_prob`, `syn_del_prob` and `syn_ins_prob` decide the iid substitution, deletion and insertion error rates introduced during the synthesis.
