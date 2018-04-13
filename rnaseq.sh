#!/bin/bash

project=$1
date=$2
reference=$3
#reference="ucsc_hg38"

dir=/sequencing/projects/$project
threads=16


cd $dir
mkdir -p analysis input log processed

cd input
ln -sf /sequencing/references/rnaseq reference
ln -sf /sequencing/nextseq/processed/$date
rm -f samples.$date
find $date/fastq/raw/ -name $project*_R1_001.fastq.gz | sed 's/_R1_001.fastq.gz//' | sort > a
keys=`cut -f4 -d'/' a | sed 's/_L00.$//' | uniq`
for k in $keys
do
  s="${date}_${k}"
  t=`grep $k a | tr '\n' ';'`
  echo -e "$s\t$t" >> samples.$date
done
cd ..

rm -f qsub.rnaseq.$project.$date.sh $dir/manifest.$reference.gtf $dir/abundances.$reference.list

#
#  TRIM
#
for sample1 in `cut -f1 input/samples.$date`
do
  sample2=`echo $sample1 | cut -f2- -d'_'`
  echo qsub -N "rnaseq.trim.$project.$date"            \
            -pe parallel 32 -cwd -S /bin/bash          \
            -o $dir/log/qsub.log.trim.$date.$sample1 \
            -e $dir/log/qsub.err.trim.$date.$sample1 \
            /usr/local/mark/rnaseq_trim.sh $dir $date $sample2 >> qsub.rnaseq.$project.$date.sh
done

#
#  FASTQC
#
for sample1 in `cut -f1 input/samples.$date`
do
  sample2=`echo $sample1 | cut -f2- -d'_'`
  echo qsub -N "rnaseq.fastqc.$project.$date"          \
            -pe parallel 32 -cwd -S /bin/bash          \
            -o $dir/log/qsub.log.fastqc.$date.$sample1 \
            -e $dir/log/qsub.err.fastqc.$date.$sample1 \
            /usr/local/mark/rnaseq_fastqc.sh $dir $date $sample2 >> qsub.rnaseq.$project.$date.sh
done

#
#  TUXEDO (A)
#
for sample1 in `cut -f1 input/samples.$date`
do
  sample2=`echo $sample1 | cut -f2- -d'_'`
  echo qsub -N "rnaseq.tuxedoA.$project.$date"                       \
            -pe parallel $threads -cwd -S /bin/bash                  \
            -o $dir/log/qsub.log.tophat2.$date.$reference.$sample1.A \
            -e $dir/log/qsub.err.tophat2.$date.$reference.$sample1.A \
            /usr/local/mark/rnaseq_tuxedo_A.sh $dir $reference $date $sample1 $threads >> qsub.rnaseq.$project.$date.sh
done

#
#  CUFFMERGE
#
echo qsub -N "rnaseq.cuffmerge.$project.$date"            \
          -hold_jid "rnaseq.tuxedoA.$project.$date"       \
          -pe parallel 32 -cwd -S /bin/bash               \
          -o $dir/log/qsub.log.cuffmerge.$date.$reference \
          -e $dir/log/qsub.err.cuffmerge.$date.$reference \
          /usr/local/mark/rnaseq_cuffmerge.sh $dir $reference $threads >> qsub.rnaseq.$project.$date.sh

#
#  TUXEDO (B)
#
for sample1 in `cut -f1 input/samples.$date`
do
  sample2=`echo $sample1 | cut -f2- -d'_'`
  echo qsub -N "rnaseq.tuxedoB.$project.$date"                       \
            -hold_jid "rnaseq.cuffmerge.$project.$date"              \
            -pe parallel $threads -cwd -S /bin/bash                  \
            -o $dir/log/qsub.log.tophat2.$date.$reference.$sample1.B \
            -e $dir/log/qsub.err.tophat2.$date.$reference.$sample1.B \
            /usr/local/mark/rnaseq_tuxedo_B.sh $dir $reference $date $sample1 $threads >> qsub.rnaseq.$project.$date.sh
  echo "processed/$sample1/$reference/transcripts.gtf" >> $dir/manifest.$reference.gtf
  echo "processed/$sample1/$reference/abundances.cxb"  >> $dir/abundances.$reference.list
done

#
#  CUFFDIFF
#
echo qsub -N "rnaseq.cuffdiff.$project.$date"            \
          -hold_jid "rnaseq.tuxedoB.$project.$date"      \
          -pe parallel 32 -cwd -S /bin/bash              \
          -o $dir/log/qsub.log.cuffdiff.$date.$reference \
          -e $dir/log/qsub.err.cuffdiff.$date.$reference \
          /usr/local/mark/rnaseq_cuffdiff.sh $dir $project $reference $threads >> qsub.rnaseq.$project.$date.sh

#
#  REPORT
#
echo qsub -N "rnaseq.report.$project.$date"            \
          -hold_jid "rnaseq.fastqc.$project.$date,rnaseq.cuffdiff.$project.$date"   \
          -pe parallel 32 -cwd -S /usr/bin/Rscript     \
          -o $dir/log/qsub.log.report.$date.$reference \
          -e $dir/log/qsub.err.report.$date.$reference \
          /usr/local/mark/rnaseq_report.r $dir $project $reference $threads >> qsub.rnaseq.$project.$date.sh

chmod +x qsub.rnaseq.$project.$date.sh

