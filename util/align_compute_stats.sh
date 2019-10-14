#!/bin/bash

# exit when any command fails
set -e

# script for aligning the fastq files to the oligos
# needs the paths below to minimap2, samtools and the directory containing the oligos and fastq files
# we used samtools 1.9 (https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2)
# and minimap2 version 2.13 (https://github.com/lh3/minimap2/releases/tag/v2.13)
MINIMAP2='/raid/shubham/minimap2/minimap2'
SAMTOOLS='/raid/shubham/samtools-1.9/samtools'
DATA_PATH='../../nanopore_dna_storage/'

# create .bed files from .fa files (will be used later for separating reads from different experiments)
for i in {0..12}; do
    # first create fai file
    $SAMTOOLS faidx $DATA_PATH/oligo_files/oligos_$i.fa
    awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' $DATA_PATH/oligo_files/oligos_$i.fa.fai > $DATA_PATH/oligo_files/oligos_$i.bed
done

# align fastq
mkdir -p $DATA_PATH/aligned
$MINIMAP2 -ax map-ont $DATA_PATH/oligo_files/merged.fa $DATA_PATH/fastq/merged.fastq > $DATA_PATH/aligned/merged.sam

 now separate aligned sam file 
for i in {0..12}; do
   $SAMTOOLS view -h -L $DATA_PATH/oligo_files/oligos_$i.bed $DATA_PATH/aligned/merged.sam > $DATA_PATH/aligned/exp_aligned_$i.sam
done

# filter by removing secondary alignments and reads with chimeric alignments (due to sequences getting ligated to each other)
for i in {0..12}; do
    $SAMTOOLS view -h -F 256 $DATA_PATH/aligned/exp_aligned_$i.sam | grep -v 'SA:Z:' > $DATA_PATH/aligned/exp_aligned_$i.filtered.sam
done

# extract raw signal by experiment 
mkdir -p $DATA_PATH/raw_signal
for i in {0..12}; do
    python3 extract_data_fast5.py $DATA_PATH/aligned/exp_aligned_$i.filtered.sam $DATA_PATH/fast5_pass $DATA_PATH/raw_signal/raw_signal_$i.hdf5
done

# convert to BAM and sort (for collecting stats)
$SAMTOOLS view -b -o $DATA_PATH/aligned/merged.bam $DATA_PATH/aligned/merged.sam
$SAMTOOLS sort -O BAM -o $DATA_PATH/aligned/merged.sorted.bam --reference $DATA_PATH/olig_files/merged.fa $DATA_PATH/aligned/merged.bam

# run samtools stats for computing error rates
mkdir -p $DATA_PATH/stats
$SAMTOOLS stats -r $DATA_PATH/oligo_files/merged.fa $DATA_PATH/aligned/merged.sorted.bam > $DATA_PATH/stats/merged.sorted.stats
grep ^MPC $DATA_PATH/stats/merged.sorted.stats | cut -f 2- > $DATA_PATH/stats/merged.sorted.substitution_stats
grep ^IC $DATA_PATH/stats/merged.sorted.stats | cut -f 2- > $DATA_PATH/stats/merged.sorted.insertion_stats
grep ^IC $DATA_PATH/stats/merged.sorted.stats | cut -f 2- > $DATA_PATH/stats/merged.sorted.deletion_stats
num_mapped=$(grep "reads mapped:" $DATA_PATH/stats/merged.sorted.stats | cut -f 3)
python3 compile_plot_stats.py $DATA_PATH/stats/merged.sorted $num_mapped
