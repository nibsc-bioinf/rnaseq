#!/bin/bash

dir=$1
project=$2
reference=$3
threads=$4

cd $dir

experiment=`grep "$project" /sequencing/projects/rnaseq.experiments | tail -1 | cut -f4-`
design=`echo "processed/$experiment/$reference/abundances.cxb" | sed "s:\([ \t,]\):/$reference/abundances.cxb\1processed/:g"`

echo /usr/bin/time cuffdiff -o processed/cuffdiff.$reference -p $threads processed/cuffmerge.$reference/merged.gtf $design
echo /usr/bin/time cuffnorm -o processed/cuffnorm.$reference -p $threads processed/cuffmerge.$reference/merged.gtf $design
/usr/bin/time cuffdiff -o processed/cuffdiff.$reference -p $threads processed/cuffmerge.$reference/merged.gtf $design
/usr/bin/time cuffnorm -o processed/cuffnorm.$reference -p $threads processed/cuffmerge.$reference/merged.gtf $design

