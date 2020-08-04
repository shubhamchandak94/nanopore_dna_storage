#!/bin/bash

# $1 -> name of the basecalled output prefix
# $2 -> The Fasta file

if [[ $# != 3 ]]; then
    echo "There should be four arguments"
    echo $#
    exit 1
    
fi

OLIGO=reads

mkdir $1
# run minimap2 alignment
MINIMAP2=/raid/shubham/minimap2/minimap2
SAMTOOLS=/raid/shubham/samtools-1.9/samtools
FASTA=$2
BASECALL=$3
$MINIMAP2 -I 16G -x map-ont -t 32 -a --secondary=no $FASTA $BASECALL | $SAMTOOLS view -bST $FASTA - > $1/basecalls.bam
$SAMTOOLS sort -O BAM -o $1/basecalls_sorted.bam --reference $FASTA $1/basecalls.bam 
$SAMTOOLS index $1/basecalls_sorted.bam  
$SAMTOOLS stats -r $FASTA $1/basecalls_sorted.bam > $1/stats.txt
head -40  $1/stats.txt

