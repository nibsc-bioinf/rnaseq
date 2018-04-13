#!/bin/bash

dir=$1
date=$2
sample2=$3

sample1=${date}_$sample2

cd $dir

mkdir -p qc/images
mkdir -p qc/$sample1

for extension in L001_R1 L001_R2 L002_R1 L002_R2 L003_R1 L003_R2 L004_R1 L004_R2
do
  fastqc input/$date/fastq/raw/${sample2}_${extension}_001.fastq.gz -o qc/$sample1/
  unzip -j qc/$sample1/${sample2}_${extension}_001_fastqc.zip ${sample2}_${extension}_001_fastqc/Images/per_base_quality.png
  mv per_base_quality.png qc/images/${sample1}_${extension}.png
done

