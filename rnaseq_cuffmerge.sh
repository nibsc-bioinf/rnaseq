#!/bin/bash

dir=$1
reference=$2
threads=$3

cd $1

cuffmerge -p $threads -o processed/cuffmerge.$reference -s input/reference/$reference.fa -g input/reference/$reference.gtf manifest.$reference.gtf
