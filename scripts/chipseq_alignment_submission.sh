#!/bin/bash/
## Script for submitting chipseq_alignment.sh jobs to the hpc.
## Based on job submission template.

## Date: 26/09/2018
## Last updated: 07/09/2021

# Specify necessary parameters for scripts.
dirLoc="/data/RUNX_KO_ChIPseq"
geneBuild="hg38"

fq=()

for f in ${dirLoc}/fastq/*.fastq.gz
do
        fq+=( $(echo $f | cut -d '/' -f 6 | cut -d '.' -f 1) )
done

mkdir ${dirLoc}/aligned_bam
mkdir ${dirLoc}/fastQC

for f in ${fq[@]}
do
  lib=$f
  script2Run=$dirLoc"/scripts/chipseq_alignment.sh"
  echo "Submitting $lib for alignment"
  qsub $script2Run -N $lib"_bowtie2_align" -v lib="$lib",dirLoc="$dirLoc",geneBuild="$geneBuild"
done
