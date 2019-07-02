#!/bin/bash

# exit when any command fails
set -e

# script for aligning the fastq files to the oligos
# needs the paths below to minimap2, samtools and the directory containing the oligos and fastq files
# we used samtools 1.9 (https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2)
# and minimap2 version 2.13 (https://github.com/lh3/minimap2/releases/tag/v2.13)
MINIMAP2='/raid/shubham/minimap2/minimap2'
SAMTOOLS='/raid/shubham/samtools-1.9/samtools'
DATA_PATH='/raid/nanopore/shubham/20190629_nanopore_data/'

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
    python extract_data_fast5.py $DATA_PATH/aligned/exp_aligned_$i.filtered.sam $DATA_PATH/aligned/fast5_pass $DATA_PATH/raw_signal/raw_signal_$i.hdf5
done

## convert to BAM and sort (for collecting stats)
#for i in {1..12}; do
#    $SAMTOOLS view -b -o $DATA_PATH/aligned/exp_aligned_$i.filtered.bam $DATA_PATH/aligned/exp_aligned_$i.filtered.sam
#    $SAMTOOLS sort -O BAM -o $DATA_PATH/aligned/exp_aligned_sorted_$i.filtered.bam --reference $DATA_PATH/oligo_files/oligos_$i.fa $DATA_PATH/aligned/exp_aligned_$i.filtered.bam
#done
#
## run samtools stats for computing error rates 
#mkdir -p $DATA_PATH/stats
#for i in {1..9}; do
#    $SAMTOOLS stats -r $DATA_PATH/oligo_files/oligos_$i.fa $DATA_PATH/aligned/exp_aligned_sorted_$i.bam > $DATA_PATH/stats/exp_aligned_sorted_$i.stats
#    grep ^MPC $DATA_PATH/stats/exp_aligned_sorted_$i.stats | cut -f 2- > $DATA_PATH/stats/exp_aligned_sorted_$i.substitution_stats
#    grep ^IC $DATA_PATH/stats/exp_aligned_sorted_$i.stats | cut -f 2,4 > $DATA_PATH/stats/exp_aligned_sorted_$i.insertion_stats 
#    grep ^IC $DATA_PATH/stats/exp_aligned_sorted_$i.stats | cut -f 2,6 > $DATA_PATH/stats/exp_aligned_sorted_$i.deletion_stats 
#    num_mapped=$(grep "reads mapped:" $DATA_PATH/stats/exp_aligned_sorted_$i.stats | cut -f 3)
#    python3 compile_plot_stats.py $DATA_PATH/stats/exp_aligned_sorted_$i $num_mapped
#done
#
## compute coverage stats (subsample to 5x and compute the coefficient of variation)
#subsampling_coverage=5
#for i in {1..9}; do
#    # first extract the ref name from sam
#    awk '$1 !~ "^@"' < $DATA_PATH/aligned/exp_aligned_$i.sam | cut -f 3 > $DATA_PATH/aligned/exp_aligned_$i.ref
#    num_oligos=$(wc -l $DATA_PATH/oligo_files/oligos_$i.bed | cut -f 1 -d ' ')
#    python3 compute_coverage_cv.py $DATA_PATH/aligned/exp_aligned_$i.ref $num_oligos $subsampling_coverage > $DATA_PATH/stats/exp_aligned_cvg_$subsampling_coverage.$i.cv.txt
#done
#
