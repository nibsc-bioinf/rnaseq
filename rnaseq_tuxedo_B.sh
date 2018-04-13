#!/bin/bash

dir=$1
reference=$2
date=$3
sample=$4
threads=$5

cd $dir
mkdir -p processed/$sample/$reference

cuffquant -o processed/$sample/$reference -p $threads processed/cuffmerge.$reference/merged.gtf processed/$sample/$reference/accepted_hits.bam
