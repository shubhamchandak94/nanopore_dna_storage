# Nanopore DNA storage
DNA storage for nanopore sequencing using convolutional coding and basecaller-decoder integration. This branch improves upon the original results described in the papers below by switching to bonito basecalling and some other improvements. The updates are briefly described [later](#updates) on this page.

### [BioRxiv](https://www.biorxiv.org/content/10.1101/2019.12.20.871939v2)

### [ICASSP 2020](https://ieeexplore.ieee.org/document/9053441)

### Data for bonito branch: https://github.com/shubhamchandak94/nanopore_dna_storage_data/tree/bonito

Code tested on Ubuntu 18.04 with Python3. 

#### Table of Contents
- [Download and installation](#download-and-installation)
- [Usage of provided scripts](#usage-of-provided-scripts)
- [Analysis workflow](#analysis-workflow)
- [Updates](#updates)

## Download and installation
Download:
```
git clone --recursive -b bonito https://github.com/shubhamchandak94/nanopore_dna_storage/
cd nanopore_dna_storage/
```

Install (compile RS code and convolutional code):
```
./install.sh
```
All the commands should be run in a virtual environment created for bonito (see `bonito/` directory for details, follow the steps in the section Developer Quickstart). For convenience, we provide the steps below as well (you can use conda environments instead of venv if you wish):
```
cd bonito
python3 -m venv venv3
source venv3/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
python setup.py develop
```

Other than this, you might need to install the following Python3 packages: crc8, distance, fast5_research, h5py, numpy, scipy, scrappy, struct (just run `./install_python_packages_bonito_venv.sh`).

Finally, to verify the dependencies are installed correctly, run `python helper.py` which runs some basic encoding and decoding roundtrip tests.

## Usage of provided scripts
In many of the scripts, you need to set the path for the corresponding data directories as well as the encoding parameters. Details about the experiments and corresponding files and scripts can be found in `oligos_8_4_20/`. The file `helper.py` contains the major functions, and these are called from the different scripts. It might be helpful to study the [Analysis workflow](#analysis-workflow) section below which describes the decoding steps in more detail.

### Parameters for experiments
The file `oligos_8_4_20/encode_experiments.py` was used for generating the oligos and contains the parameters required for the decoding.

### Convolutional code decoding of raw signal data
First, generate a list of read ids to decode using the script `util/generate_read_ids.py` (changing the parameters as needed). Then, run `generate_decoded_lists.py` with the relevant parameters. See the directory `oligos_8_4_20/decoding_scripts/` for some examples.

Note that the `rate_conv` parameter is set to 1, 2, 3, 4, 5 and 7 for convolutional code rates of 1/2, 2/3, 3/4, 4/5, 5/6 and 7/8, respectively. The `msg_len` parameter is the length of the binary input to the convolutional code encoder. This writes the decoded lists to files named `list_0, list_1, ...` in the output directory and also generates an `info.txt` file listing the decoded reads.
The parameters for the different experiments can be found in the files `oligos_8_4_20/encode_experiments.py` and `oligos_8_4_20/encoding_log.txt`.

### Computing error rates of convolutional code decoding
To compute the number of errors (due to no CRC match being found and due to incorrect CRC match) from the decoded lists, use `compute_error_rate_from_decoded_lists.py` after setting the parameters. For experiments employing the two-CRC strategy, use `compute_error_rate_from_decoded_lists_2CRC.py`.

### Performing RS decoding of the lists
To perform RS decoding from the lists, use `decode_RS_from_decoded_lists.py` after setting the parameters. `NUM_TRIALS` denotes the number of decoding trials performed. Each trial involves a subsampling of `NUM_READS_TO_USE` reads from the set of decoded lists of size `NUM_READS_TOTAL`, an attempt at RS decoding and comparison with the original encoded file to ensure successful decoding. For experiments employing the two-CRC strategy, use `decode_RS_from_decoded_lists_2CRC.py`.

### Finding minimum number of reads required for RS decoding of the lists
To find the minimum number of reads required for successful RS decoding from the lists, use `compute_min_reads_for_decoding.py` after setting the parameters. `NUM_TRIALS` denotes the number of decoding trials performed (success declared when all trials succeed). Each trial involves a subsampling of `NUM_READS_TO_USE` reads from the set of decoded lists of size `NUM_READS_TOTAL`, an attempt at RS decoding and comparison with the original encoded file to ensure successful decoding. You should specify the `NUM_READS_TO_USE_START` (the initial value of `NUM_READS_TO_USE`) and `NUM_READS_TO_USE_STEP` (the value by which `NUM_READS_TO_USE` is incremented until we get success). For experiments employing the two-CRC strategy, use `compute_min_reads_for_decoding_2CRC.py`.

### Running simulations
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

In addition, we provide functions to simulate the entire pipeline including the Reed-Solomon outer code in `helper.py`. A small sample simulation is included in the `if __name__ == '__main__'` block and can be run by executing `python helper.py`.

### Utility scripts
The [`util/`](util/) directory contains several utility scripts used for preparing data in a format that `generate_decoded_lists.py` can interpret, and also for performing various forms of statistical analysis of the data based on alignment.

## Analysis workflow
In this section, we describe the decoding workflow (as mentioned above, the file [`oligos_8_4_20/encode_experiments.py`](oligos_8_4_20/encode_experiments.py) was used for generating the oligos and contains the parameters required for the decoding). We start from the fast5 files and discuss the steps involved in getting to the decoded file. We have provided the raw files as well as several intermediate analysis files at https://github.com/shubhamchandak94/nanopore_dna_storage_data/tree/bonito. _Note that the intermediate analysis files can be used to test parts of the pipeline without needing to go through the entire process._ 

For concreteness, we will focus on an example of a specific sequencing run with multiple replicates/barcodes (`20210304_MIN_0964` in the data repo), and specifically the subpool experiment 0 (convolutional code m=6, r=3/4) for `barcode01`. We also include the minor modifications (mostly in path names) required for sequencing run with no barcode multiplexing. 

For context, note that we included 18 experiments in a single synthesis pool with different PCR primers. The 18 experiments correspond to different convolutional code memory, rate and number of CRCs used, with all the parameters described in [`oligos_8_4_20/encode_experiments.py`](oligos_8_4_20/encode_experiments.py)

### Directory structure
We will assume the fast5 files are in `$FAST5_PATH` and the analysis data is in `$DATA`. The `$DATA` directory has further subdirectories, initially just the `$DATA/oligo_files/` that contains the files included [here](https://github.com/shubhamchandak94/nanopore_dna_storage_data/tree/bonito/oligo_files).

### Basecalling
First step is to basecall the fast5 so the reads can be separated according to the different experiments and barcodes. We generated the fastq files from the fast5 files using Guppy 4.0.14 basecaller. The results should be similar with more recent versions of the basecaller.

Downloading basecaller:
```
wget https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy_4.0.14_linux64.tar.gz
tar -xzvf ont-guppy_4.0.14_linux64.tar.gz
```
To basecall, run:
```
./ont-guppy/bin/guppy_basecaller -i $FAST5_PATH -s $DATA/fastq/ont_guppy -c dna_r9.4.1_450bps_hac.cfg -x cuda:all
```
For multiplexed experiments (pool `20210304_MIN_0964`), run the barcoder:
```
./ont-guppy/bin/guppy_barcoder --input_path $DATA/fastq/ont_guppy --save_path $DATA/fastq/ont_guppy_demultiplexed --config configuration.cfg --barcode_kits "EXP-NBD104 EXP-NBD114" -x cuda:0
```
Next we need to merge the fastq files produced by a basecaller into a single file using a command like:
```
cat $DATA/fastq/ont_guppy_demultiplexed/barcode01/*.fastq > $DATA/fastq/barcode01/merged.fastq
```
In case of no barcoding, just run
```
cat $DATA/fastq/ont_guppy/*.fastq > $DATA/fastq/merged.fastq
```

### Alignment, separating by experiment and preparing raw data for decoding
The next step involves aligning the reads to the original oligo files (only to separate the experiments, the original oligos are not used for decoding). In addition we generate some statistics and extract a random subset of raw signals for each experiment that will be used for the decoding in subsequent steps. We provide a script for this in [`util/align_compute_stats_replicates.sh`](util/align_compute_stats_replicates.sh). First step is to set the variables at the top of the script:
- `MINIMAP2`: path to minimap2 aligner executable
- `SAMTOOLS`: path to samtools executable
- `DATA_PATH`: the `$DATA` directory
- `FAST5_PATH`: the `$FAST5_PATH` directory
- `EXPERIMENTS`: which experiments are included in the sequencing run (e.g., `"0 1 2 5 8"` for the `20210304_MIN_0964` run)
- `NUM_READS_FOR_DECODING`: how many reads we decode with the Viterbi decoder. Since Viterbi decoding is a slow process, we keep this to a small number (10000) which is more than sufficient for successful decoding and subsampling experiments.

Once these parameters are set, you can just run
```
./util/align_compute_stats_replicates.sh barcode01
```
to process the `barcode01` replicate.

In case of no barcoding, use the `util/align_compute_stats.sh` script with same variables to be set, and run it as
```
./util/align_compute_stats.sh
```

After this step, the `$DATA` directory contains the following subdirectories (with further subsubdirectories for each barcode if applicable):
- `oligo_files`: the original oligos (made available [here](https://github.com/shubhamchandak94/nanopore_dna_storage_data/tree/bonito/oligo_files))
- `fastq`: basecalled and demultiplexed reads
- `aligned`: aligned reads, before and after separating by experiments
- `raw_signal`: Contains the raw signals in hdf5 format for the different experiments (made available [here](https://github.com/shubhamchandak94/nanopore_dna_storage_data/tree/bonito/raw_signal/))
- `stats`: Contains statistics for the basecalling errors and alignment (made available [here](https://github.com/shubhamchandak94/nanopore_dna_storage_data/tree/bonito/stats))

One further analysis step described in the paper is the coverage analysis. We describe an example execution below:
```
grep "^[^@]" $DATA/aligned/barcode01/exp_aligned_0.filtered.sam | cut -f 3 > $DATA/aligned/barcode01/exp_aligned_0.filtered.txt
python util/compute_coverage_stats.py $DATA/aligned/barcode01/exp_aligned_0.filtered.txt 880 5
```
The first step extracts the reference sequence of the different aligned oligos, and the second step uses this to determine the coverage variance and the fraction of oligos with zero coverage. For consistency across experiments with different overall coverage, we determine these metrics for a subsampling with 5x coverage (the third parameter). The second parameter `880` denotes the number of oligos for the specific experiment which can be seen in the supplementary tables provided, or in [`oligos_8_4_20/encoding_log.txt`](oligos_8_4_20/encoding_log.txt).

### Extracting read ids for decoding
In this step, we extract the read ids that will be used for the decoding. These are extracted from the `raw_signal` obtained in the previous step. This is a relic of an older analysis pipeline, but can still be useful for testing decoding on a small number of reads instead of all 10000 reads extracted in the previous step. The first step is to `mkdir -p $DATA/read_ids/barcode01/`. Then run
```
python util/generate_read_id_file.py 0 10000 $DATA barcode01
```
This will generate a file `$DATA/read_ids/barcode01/exp_0_read_ids.10000.txt` which contains the list of read ids that will be decoded in the next step (obtained from `$DATA/raw_signal/barcode01/raw_signal_0.hdf5` corresponding to experiment 0). 

In case of no barcoding, do the following:
```
mkdir -p $DATA/read_ids/
python util/generate_read_id_file.py 0 10000 $DATA
```
This will generate a file `$DATA/read_ids/exp_0_read_ids.10000.txt` which contains the list of read ids that will be decoded in the next step (obtained from `$DATA/raw_signal/raw_signal_0.hdf5` corresponding to experiment 0). 

We have made available all the read id files used in our analysis [here](https://github.com/shubhamchandak94/nanopore_dna_storage_data/tree/bonito/raw_signal).

### Convolutional code Viterbi decoding
This is the main decoding step in the pipeline, and the one that takes the most time. Note that all the inputs as well as the outputs are already available in the [data repo](https://github.com/shubhamchandak94/nanopore_dna_storage_data/). All the decoding scripts are available in this repository, the ones for the `20200814_MIN_0880` pool (all 18 experiments, no barcoding) are at [`oligos_8_4_20/decoding_scripts`](oligos_8_4_20/decoding_scripts), and the ones for the `20210304_MIN_0964` pool (5 experiments, 6 replicates) are at [`oligos_8_4_20/replicates`](oligos_8_4_20/replicates). The only difference is the path convention, otherwise the scripts are similar. Note that we have one script per experiment since they pass different convolutional code parameters and primer sequences to [`generate_decoded_lists.py`](generate_decoded_lists.py) which actually does the decoding.

Let's focus on an example for `barcode01` and experiment 0, with the decoding script at [`oligos_8_4_20/replicates/run_decoder_exp_0_bonito_default.sh`](oligos_8_4_20/replicates/run_decoder_exp_0_bonito_default.sh). We first make the following change to this file: on line 3, set `ROOT_PATH` to the `$DATA` directory. You might also want to change the `--num_threads` or the `--list_size` parameter to manage the speed. We also provide the possibility of passing a non-default model, but we weren't able to get any benefits from the finetuning, and we used just the default parameters throughout. Finally, you should verify that the `HDF5_INPUT` and the `READ_ID_FILE` files follow the correct path convention you used (and edit them if not). You can also modify the `DIRNAME` line to obtain the output in a different directory. Next we run,
```
./oligos_8_4_20/replicates/run_decoder_exp_0_bonito_default.sh barcode01/
```
This can take up to a few days depending on the parameters and the number of reads in question. Once this finishes, we obtain the decoded lists for all the reads in the `DIRNAME` directory, with one file per read, named as `list_0`, `list_1`, etc. Missing file denotes that the primer removal failed and the convolutional code decoding was not performed. Note that the list files correspond to the read ids in the file provided as input to the script. 

All the decoded lists from this step are made available [here](https://github.com/shubhamchandak94/nanopore_dna_storage_data/tree/bonito/decoded_lists/).

### Analyzing decoded lists
The final step involves obtaining the final results from the decoded lists from the previous step (made available [here](https://github.com/shubhamchandak94/nanopore_dna_storage_data/tree/bonito/decoded_lists/)). 

The script [`compute_min_reads_for_decoding.py`](compute_min_reads_for_decoding.py) can be used to find the minimum number of reads needed for successful decoding, which is used to calculate our chief metric (reading cost). The script attempts decoding with progressively larger number of reads until we have success. There are several parameters in the script that need to be set before running:

- `NUM_READS_TOTAL`: total number of reads used for decoding, set to 10000 for this example.
- `NUM_READS_TO_USE_START`: initial number of reads to start testing the decoding (can set this as small as desired, but it's better to start with a reasonable value so we get success sooner). We set `NUM_READS_TO_USE = NUM_READS_TO_USE_START` at the start.
- `NUM_READS_TO_USE_STEP`: the value by which `NUM_READS_TO_USE` is incremented until we get success
- `NUM_TRIALS`: number of decoding trials performed at each value of `NUM_READS_TO_USE` (success declared when all trials succeed) 
- `LIST_SIZE`: list size (should be less than or equal to the list size in the Viterbi decoding)
- `bytes_per_oligo`, `pad`: related to the convolutional code parameters, use the value for your experiment from [`oligos_8_4_20/encoding_log.txt`](oligos_8_4_20/encoding_log.txt).
- `DECODED_LISTS_DIR`: set to the `DIRNAME` parameter from the previous step.
- `ORIGINAL_FILE`: set to `oligos_8_4_20/data_files.tar.bz2.enc.1` (this is used to verify that the decoding succeeded)

Another useful script is [`compute_error_rate_from_decoded_lists.py`](compute_error_rate_from_decoded_lists.py) which can help calculate the number of reads that were decoded successfully. The parameters in this are a subset of the parameters in [`compute_min_reads_for_decoding.py`](compute_min_reads_for_decoding.py), and should be set in the same way.

Note that there is another set of scripts for the 2CRC case (experiments 9 to 17 in the oligo pool), with similar parameters.

## Updates 
### Notes on bonito integration
We have implemented two approaches for bonito integration. Note that bonito outputs probabilities instead of log-probabilities, which caused issues in the earlier implementations that assumed the wrong thing. 

Both approaches are based on ideas from Parallel LVA algorithm as described in 
https://github.com/shubhamchandak94/kBestViterbi/blob/master/kBestViterbi.py
or in ieeexplore.ieee.org/iel1/26/12514/00577040.pdf and on beam search for CTC as described in https://distill.pub/2017/ctc/
and in https://gist.github.com/awni/56369a90d03953e370f3964c826ed4b0.
#### Approach 1 (default)
Implemented in `viterbi_convolutional_code_approach_1.cpp`. We have states correspoding to conv code state and pos in msg. Each state stores a list of messages with their score of ending in blank and non-blank. At next step, we add one character, add (logsumexp) the scores for all the ways a new message can be obtained, and take the top ones. This approach is very closely tied to the way beam search usually works, just keeping a beam for each state.

One issue is that if we replace logsumexp with max and set list size to 1, this does not reduce to finding the best path satisfying the code constraint.

#### Approach 2
Implemented in `viterbi_convolutional_code_approach_2.cpp`. We have states correspoding to conv code state, pos in msg and CTC state. Each state stores a list of messages with their score. At each time step, we consider the valid transitions (stay transitions to blank, stay transitions to same nonblank base, nonstay transitions to nonblank base). We also combine paths corresponding to same message using logsumexp. This approach is very similar to what we did for flappie integration (although the implementation is a bit different since we do logsumexp instead of max, so the heap strategy is no longer applicable).

This does reduce to the best path satisfying code constraint if logsumexp is changed to max and list size is set to 1.

### Improvements to barcode removal
Some analysis showed that the previous barcode removal strategy led to suboptimal performance. The following issues were identified and fixed:
- Earlier the start barcode was being searched only in first half of read, and end barcode only in second half. This is bad when read is long (adapters, chimeric, something else). So we allow search of start barcode throughout and end barcode anywhere after the start barcode.
- Earlier barcode was searched only within the read, i.e., the possibility that part of the barcode is cut off was not entertained. While working on guppy, a barcode extend len parameter was included that allowed this to happen. The parameter determined how much overhang of barcode was allowed and was set to 10 arbitrarily. Now this has been set to the barcode length by default, which improves things a bit when more than 10 bases of the barcode are cut off. Note that the overhang portion contributes to the edit distance and is not "free".
- As mentioned above, the extended barcode portion that is not part of read contributes to the edit distance. In some cases, this can be a bit bad, e.g., consider that last 15 bases of 25 length barcode perfectly match the first 25 bases of the read. In this case the edit distance for the appropriate shift is still 10 which can lead to a overall worse match being selected from within the read. To resolve this issue, we add a barcode extend penalty parameter. When the parameter is 1.0, this corresponds to the usual edit distance. When it is 0.0, that means we don't penalize the extension at all, which is bas because then the "match" with barcode completely outside read is very likely to get selected. We set it to 0.6 based on some experiments. In general, this can be chosen based on the likelihood of the barcode being not present in the read and amount of barcode generally present in the read.

### Two-CRC stategy
Idea is to use 2 CRCs so that we can use part of the decoded message from convolutional coding if the whole thing doesn't succeed. This can help because the errors are often localized, and so as long as the index + first/second half of the payload is correct, we get useful information. Since the index is closer to the first half, we observed that the correct decoding in the index happened more frequently with the correct decoding of the first half. This is not ideal because the RS outer coding is applied per columnar segment, and we need a threshold number of correct segments in each column. To make sure each column has roughly equal number of correctly decoded segments, we cyclically rotate the segments in the oligos after RS encoding so that each RS code has equal number of segments in different parts of the oligo.  
