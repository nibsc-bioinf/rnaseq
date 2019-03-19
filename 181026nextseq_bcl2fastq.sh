#!/bin/bash

dirOutput=$1 #/sequencing/nextseq/output/180711_NB501506_0020_AHKGF7BGX7
dirProcessed=$2 #/sequencing/nextseq/processed/180711
dirZips=$3 #/sequencing/zips
date=$4 #180711

cd $dirProcessed

#
#  Process sequencing data into FASTQ.
#
bcl2fastq --barcode-mismatches 0 -R $dirOutput --interop-dir $dirProcessed/InterOp --sample-sheet SampleSheet.csv -o bcl2fastq

#
#  Link from processed output data into processed
#
cd fastq/raw
for s in `find ../../bcl2fastq -name "*.fastq.gz"` ; do ln -s $s ; done
cd ../.. 

projects=`find "$dirProcessed/fastq/raw/" -name "*.fastq.gz" | cut -f8 -d'/' | cut -b1-3 | sort -u`
for project in $projects; do
  echo $project
  ls $dirProcessed/fastq/raw/${project}*.fastq.gz
  zip -j $dirZips/$project-$date.zip $dirProcessed/fastq/raw/${project}*.fastq.gz
done

chown -R root.hpc_sequencing_writers "$dirProcessed"
chmod -R 755 $dirProcessed
