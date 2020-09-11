This directory contains some of the data associated with the decoding. The fast5 files will be made available at TODO and they can be used to generate the intermediate raw signal and aligned files using the script available at `../../util/align_compute_stats.sh`. The intermediate files can be made available separately due to their large size. 
- `read_ids/`: list of read ids that were decoded for each experiment. These were generated using `../../util/generate_read_id_file.py`.
- `stats/`: basecalling stats for the overall pool and the individual experiments.
- `decoded_lists/`: decoded lists from convolutional decoding. These were generated using the scripts at `../decoded_scripts/` and can be used to recover the data after doing RS decoding.
