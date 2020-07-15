#!/bin/bash

# exit when any command fails
set -e

# guppy commands used
# Download: wget https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy_4.0.11_linux64.tar.gz
# untar: tar -xzvf ont-guppy_4.0.11_linux64.tar.gz
# Basecall: ./ont-guppy/bin/guppy_basecaller --input_path 20200710_MIN_0866/20200708_dnastorage_convcode2_resynth/20200708_2225_MN19956_FAN43057_9c0c7c71/fast5 --device cuda:3 --config dna_r9.4.1_450bps_hac.cfg --save_path data/fastq/guppy_out
# cd data/fastq/
# cat guppy_out/*.fastq > merged.fastq


# script for aligning the fastq files to the oligos
# needs the paths below to minimap2, samtools and the directory containing the oligos and fastq files
# we used samtools 1.9 (https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2)
# and minimap2 version 2.13 (https://github.com/lh3/minimap2/releases/tag/v2.13)
MINIMAP2='/raid/shubham/minimap2/minimap2'
SAMTOOLS='/raid/shubham/samtools-1.9/samtools'
DATA_PATH='../../data/'
FAST5_PATH='/raid/nanopore/shubham/20200214_nanopore_pool_data/20200710_MIN_0866/20200708_dnastorage_convcode2_resynth/20200708_2225_MN19956_FAN43057_9c0c7c71/fast5/'
# create .bed files from .fa files (will be used later for separating reads from different experiments)
for i in {0..17}; do
    # first create fai file
    $SAMTOOLS faidx $DATA_PATH/oligo_files/oligos_$i.fa
    awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' $DATA_PATH/oligo_files/oligos_$i.fa.fai > $DATA_PATH/oligo_files/oligos_$i.bed
done

# align fastq
mkdir -p $DATA_PATH/aligned
$MINIMAP2 -ax map-ont $DATA_PATH/oligo_files/merged.fa $DATA_PATH/fastq/merged.fastq > $DATA_PATH/aligned/merged.sam

# now separate aligned sam file 
for i in {0..17}; do
   $SAMTOOLS view -h -L $DATA_PATH/oligo_files/oligos_$i.bed $DATA_PATH/aligned/merged.sam > $DATA_PATH/aligned/exp_aligned_$i.sam
done

# filter by removing secondary alignments and reads with chimeric alignments (due to sequences getting ligated to each other)
for i in {0..17}; do
    $SAMTOOLS view -h -F 256 $DATA_PATH/aligned/exp_aligned_$i.sam | grep -v 'SA:Z:' > $DATA_PATH/aligned/exp_aligned_$i.filtered.sam
done

# extract raw signal by experiment 
mkdir -p $DATA_PATH/raw_signal
for i in {0..17}; do
    python3 extract_data_fast5.py $DATA_PATH/aligned/exp_aligned_$i.filtered.sam $FAST5_PATH $DATA_PATH/raw_signal/raw_signal_$i.hdf5
done

# convert to BAM and sort (for collecting stats)
$SAMTOOLS view -b -o $DATA_PATH/aligned/merged.bam $DATA_PATH/aligned/merged.sam
$SAMTOOLS sort -O BAM -o $DATA_PATH/aligned/merged.sorted.bam --reference $DATA_PATH/oligo_files/merged.fa $DATA_PATH/aligned/merged.bam
for i in {0..17}; do
    $SAMTOOLS view -b -o $DATA_PATH/aligned/exp_aligned_$i.filtered.bam $DATA_PATH/aligned/exp_aligned_$i.filtered.sam
    $SAMTOOLS sort -O BAM -o $DATA_PATH/aligned/exp_aligned_$i.filtered.sorted.bam --reference $DATA_PATH/oligo_files/merged.fa $DATA_PATH/aligned/exp_aligned_$i.filtered.bam
done

# run samtools stats for computing error rates
mkdir -p $DATA_PATH/stats
$SAMTOOLS stats -r $DATA_PATH/oligo_files/merged.fa $DATA_PATH/aligned/merged.sorted.bam > $DATA_PATH/stats/merged.sorted.stats.txt
for i in {0..17}; do
    $SAMTOOLS stats -r $DATA_PATH/oligo_files/merged.fa $DATA_PATH/aligned/exp_aligned_$i.filtered.sorted.bam > $DATA_PATH/stats/exp_aligned_$i.filtered.sorted.stats.txt
done
