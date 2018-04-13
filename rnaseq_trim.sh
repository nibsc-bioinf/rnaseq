#!/bin/bash

function fastqStats () {
  stats=`echo $1 | sed 's/.fastq.gz/.stats/'`
  file=`echo $1 | rev | cut -d'/' -f1-2 | rev`
  echo $stats
  zcat $1 | sed -n '2~4p' | wc | awk -v fastq=$file -F ' ' '{ printf "%s\t%d\t%d\t%0.2f\n", fastq, $1, ($3-$1), ($3/$1-1); }' > $stats
}

dir=$1
date=$2
sample=$3

fastq="input/$date/fastq/raw/${sample}_L001_R1_001.fastq.gz"

adapter="TruSeq3-PE.fa"

cd $dir

for L00x in L001 L002 L003 L004
do
  fastq1=`echo $fastq  | sed "s/_L001_/_${L00x}_/"`
  fastq2=`echo $fastq1 | sed 's/_R1_/_R2_/'`
  trimp1=`echo $fastq1 | sed 's/raw/30/' | sed 's/_R1_.*/.p1.fastq.gz/'`
  trimp2=`echo $fastq1 | sed 's/raw/30/' | sed 's/_R1_.*/.p2.fastq.gz/'`
  trims1=`echo $fastq1 | sed 's/raw/30/' | sed 's/_R1_.*/.s1.fastq.gz/'`
  trims2=`echo $fastq1 | sed 's/raw/30/' | sed 's/_R1_.*/.s2.fastq.gz/'`

  java -jar /usr/local/bin/trimmomatic-0.33.jar PE -threads 16 -phred33 $fastq1 $fastq2 $trimp1 $trims1 $trimp2 $trims2 \
      ILLUMINACLIP:/usr/local/bin/adapters/$adapter:2:30:10 SLIDINGWINDOW:5:30 LEADING:10 TRAILING:10 MINLEN:50

  fastqStats $fastq1
  fastqStats $fastq2
  fastqStats $trimp1
  fastqStats $trimp2
  fastqStats $trims1
  fastqStats $trims2
done
