# Nanopore DNA storage
DNA storage for nanopore sequencing using convolutional coding and basecaller-decoder integration

### [BioRxiv](https://www.biorxiv.org/content/10.1101/2019.12.20.871939v2)

### [Supplementary material](https://github.com/shubhamchandak94/nanopore_dna_storage/blob/master/supplementary_material.pdf)

### Data: https://github.com/shubhamchandak94/nanopore_dna_storage_data

Code tested on Ubuntu 18.04.1 and Python3.
## Download and install
Download:
```
git clone --recursive https://github.com/shubhamchandak94/nanopore_dna_storage/
```

Install:
```
./install.sh
```
All the steps below should be run in a virtual environment created for bonito (see `bonito/` directory for details).
Other than this, you might need to install the following Python3 packages: crc8, distance, fast5_research, h5py, numpy, scipy, scrappy, struct (see `install_python_packages_bonito_venv.sh`).

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

## Notes on bonito integration
We have implemented two approaches for bonito integration. Note that bonito outputs probabilities instead of log-probabilities, which caused issues in the earlier implementations that assumed the wrong thing.

Both approaches are based on ideas from Parallel LVA algorithm as described in 
https://github.com/shubhamchandak94/kBestViterbi/blob/master/kBestViterbi.py
or in ieeexplore.ieee.org/iel1/26/12514/00577040.pdf and on beam search for CTC as described in https://distill.pub/2017/ctc/
and in https://gist.github.com/awni/56369a90d03953e370f3964c826ed4b0.
### Approach 1
Implemented in `viterbi_convolutional_code_approach_1.cpp`. We have states correspoding to conv code state and pos in msg. Each state stores a list of messages with their score of ending in blank and non-blank. At next step, we add one character, add (logsumexp) the scores for all the ways a new message can be obtained, and take the top ones. This approach is very closely tied to the way beam search usually works, just keeping a beam for each state.

One issue is that if we replace logsumexp with max and set list size to 1, this does not reduce to finding the best path satisfying the code constraint.

### Approach 2
Implemented in `viterbi_convolutional_code_approach_2.cpp`. We have states correspoding to conv code state, pos in msg and CTC state. Each state stores a list of messages with their score. At each time step, we consider the valid transitions (stay transitions to blank, stay transitions to same nonblank base, nonstay transitions to nonblank base). We also combine paths corresponding to same message using logsumexp. This approach is very similar to what we did for flappie integration (although the implementation is a bit different since we do logsumexp instead of max, so the heap strategy is no longer applicable).

This does reduce to the best path satisfying code constraint if logsumexp is changed to max and list size is set to 1.

### Improvements to barcode removal
Some analysis showed that the previous barcode removal strategy led to suboptimal performance. The following issues were identified and fixed:
- Earlier the start barcode was being searched only in first half of read, and end barcode only in second half. This is bad when read is long (adapters, chimeric, something else). So we allow search of start barcode throughout and end barcode anywhere after the start barcode.
- Earlier barcode was searched only within the read, i.e., the possibility that part of the barcode is cut off was not entertained. While working on guppy, a barcode extend len parameter was included that allowed this to happen. The parameter determined how much overhang of barcode was allowed and was set to 10 arbitrarily. Now this has been set to the barcode length by default, which improves things a bit when more than 10 bases of the barcode are cut off. Note that the overhang portion contributes to the edit distance and is not "free".
- As mentioned above, the extended barcode portion that is not part of read contributes to the edit distance. In some cases, this can be a bit bad, e.g., consider that last 15 bases of 25 length barcode perfectly match the first 25 bases of the read. In this case the edit distance for the appropriate shift is still 10 which can lead to a overall worse match being selected from within the read. To resolve this issue, we add a barcode extend penalty parameter. When the parameter is 1.0, this corresponds to the usual edit distance. When it is 0.0, that means we don't penalize the extension at all, which is bas because then the "match" with barcode completely outside read is very likely to get selected. We set it to 0.6 based on some experiments. In general, this can be chosen based on the likelihood of the barcode being not present in the read and amount of barcode generally present in the read.
