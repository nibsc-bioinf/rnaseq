#!/bin/bash

dir=$1
reference=$2
date=$3
sample=$4
threads=$5

cd $dir
mkdir -p processed/$sample/$reference

t=`grep -P "^$sample\t" input/samples.$date | cut -f2`
fq1=`echo $t | sed 's/;/_R1_001.fastq.gz,/g' | sed "s|$date|input/$date|g" | sed 's/,$//'`
fq2=`echo $t | sed 's/;/_R2_001.fastq.gz,/g' | sed "s|$date|input/$date|g" | sed 's/,$//'`
tophat2 -p $threads -o processed/$sample/$reference -G input/reference/$reference.gtf input/reference/$reference $fq1 $fq2
cufflinks -o processed/$sample/$reference -G input/reference/$reference.gtf processed/$sample/$reference/accepted_hits.bam
