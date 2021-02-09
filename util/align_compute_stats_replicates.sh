#!/bin/bash

# exit when any command fails
set -e

# guppy commands used
# Download: wget https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy_4.0.14_linux64.tar.gz
# untar: tar -xzvf ont-guppy_4.0.14_linux64.tar.gz
# Basecall: ./ont-guppy/bin/guppy_basecaller --input_path 20200710_MIN_0866/20200708_dnastorage_convcode2_resynth/20200708_2225_MN19956_FAN43057_9c0c7c71/fast5 --device cuda:3 --config dna_r9.4.1_450bps_hac.cfg --save_path data/fastq/guppy_out
# cd data/fastq/
# cat guppy_out/*.fastq > merged.fastq


# script for aligning the fastq files to the oligos
# needs the paths below to minimap2, samtools and the directory containing the oligos and fastq files
# we used samtools 1.9 (https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2)
# and minimap2 version 2.13 (https://github.com/lh3/minimap2/releases/tag/v2.13)
MINIMAP2='/raid/shubham/minimap2/minimap2'
SAMTOOLS='/raid/shubham/samtools-1.9/samtools'
DATA_PATH='/raid/nanopore/shubham/20200804_nanopore_pool_data/data_20210205_MIN_0952/'
FAST5_PATH='//raid/nanopore/20210205_MIN_0952/dnastorage_hifi_apex/20210204_0014_MN19956_FAO47706_eb4e84b9/fast5'

REPLICATE=$1
EXPERIMENTS="0 1 2 5 8"
NUM_READS_FOR_DECODING=10000

# create .bed files from .fa files (will be used later for separating reads from different experiments)
for i in $EXPERIMENTS; do
    # first create fai file
    $SAMTOOLS faidx $DATA_PATH/oligo_files/oligos_$i.fa
    awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' $DATA_PATH/oligo_files/oligos_$i.fa.fai > $DATA_PATH/oligo_files/oligos_$i.bed
done

# align fastq
mkdir -p $DATA_PATH/aligned/$REPLICATE
$MINIMAP2 -ax map-ont $DATA_PATH/oligo_files/merged.fa $DATA_PATH/fastq/$REPLICATE/merged.fastq > $DATA_PATH/aligned/$REPLICATE/merged.sam

# now separate aligned sam file 
for i in $EXPERIMENTS; do
   $SAMTOOLS view -h -L $DATA_PATH/oligo_files/oligos_$i.bed $DATA_PATH/aligned/$REPLICATE/merged.sam > $DATA_PATH/aligned/$REPLICATE/exp_aligned_$i.sam
done

# filter by removing secondary alignments and reads with chimeric alignments (due to sequences getting ligated to each other)
for i in $EXPERIMENTS; do
    $SAMTOOLS view -h -F 256 $DATA_PATH/aligned/$REPLICATE/exp_aligned_$i.sam | grep -v 'SA:Z:' > $DATA_PATH/aligned/$REPLICATE/exp_aligned_$i.filtered.sam
done

# convert to BAM and sort (for collecting stats)
$SAMTOOLS view -b -o $DATA_PATH/aligned/$REPLICATE/merged.bam $DATA_PATH/aligned/$REPLICATE/merged.sam
$SAMTOOLS sort -O BAM -o $DATA_PATH/aligned/$REPLICATE/merged.sorted.bam --reference $DATA_PATH/oligo_files/merged.fa $DATA_PATH/aligned/$REPLICATE/merged.bam
for i in $EXPERIMENTS; do
    $SAMTOOLS view -b -o $DATA_PATH/aligned/$REPLICATE/exp_aligned_$i.filtered.bam $DATA_PATH/aligned/$REPLICATE/exp_aligned_$i.filtered.sam
    $SAMTOOLS sort -O BAM -o $DATA_PATH/aligned/$REPLICATE/exp_aligned_$i.filtered.sorted.bam --reference $DATA_PATH/oligo_files/merged.fa $DATA_PATH/aligned/$REPLICATE/exp_aligned_$i.filtered.bam
done

# run samtools stats for computing error rates
mkdir -p $DATA_PATH/stats/$REPLICATE
$SAMTOOLS stats -r $DATA_PATH/oligo_files/merged.fa $DATA_PATH/aligned/$REPLICATE/merged.sorted.bam > $DATA_PATH/stats/$REPLICATE/merged.sorted.stats.txt
for i in $EXPERIMENTS; do
    $SAMTOOLS stats -r $DATA_PATH/oligo_files/merged.fa $DATA_PATH/aligned/$REPLICATE/exp_aligned_$i.filtered.sorted.bam > $DATA_PATH/stats/$REPLICATE/exp_aligned_$i.filtered.sorted.stats.txt
done

# extract raw signal by experiment 
mkdir -p $DATA_PATH/raw_signal/$REPLICATE/
for i in $EXPERIMENTS; do
    python3 extract_data_fast5.py $DATA_PATH/aligned/$REPLICATE/exp_aligned_$i.filtered.sam $FAST5_PATH $DATA_PATH/raw_signal/$REPLICATE/raw_signal_$i.hdf5 $NUM_READS_FOR_DECODING
done
