#!/bin/bash
## Script to perform alignment of ChIP-seq data using bowtie2.
##

## Date: 20/09/2018
## Last updated: 07/01/2020


#PBS -m ae
#PBS -l nodes=1:ppn=4
#PBS -S /bin/bash

###############################################################################
## Required qsub parameters:
## $lib = sample/library unique identifier (name) without filepaths
## $dirLoc = filepath of root directory in which the analysis is taking place
## $geneBuild = genome build to be used in experiment/analysis
###############################################################################

# Run fastQC to generate fastq quality report.
/usr/local/FastQC/fastqc ${dirLoc}/fastq/${lib}.fastq.gz --outdir=${dirLoc}/fastQC

# Run bowtie2 aligner on fastq files.
/usr/local/bowtie2-2.1.0/bowtie2 -x /data/Resources/${geneBuild}/${geneBuild} -U ${dirLoc}/fastq/${lib}.fastq.gz -k 1 -q -p 4 | /usr/local/samtools-1.3.1/bin/samtools view -bS - > ${dirLoc}/aligned_bam/${lib}.bam

# Filter aligned reads based on mappability. Use a threshold of 30 to filter reads.
/usr/local/samtools-1.3.1/bin/samtools view -b ${dirLoc}/aligned_bam/${lib}.bam > ${dirLoc}/aligned_bam/${lib}.nofilter.bam
/usr/local/samtools-1.3.1/bin/samtools view -b -q 30 ${dirLoc}/aligned_bam/${lib}.bam > ${dirLoc}/aligned_bam/${lib}.q30.bam

# Sort bams.
/usr/local/samtools-1.3.1/bin/samtools sort ${dirLoc}/aligned_bam/${lib}.nofilter.bam > ${dirLoc}/aligned_bam/${lib}.nofilter.sorted.bam
/usr/local/samtools-1.3.1/bin/samtools sort ${dirLoc}/aligned_bam/${lib}.q30.bam > ${dirLoc}/aligned_bam/${lib}.q30.sorted.bam

# Index bams.
/usr/local/samtools-1.3.1/bin/samtools index ${dirLoc}/aligned_bam/${lib}.nofilter.sorted.bam
/usr/local/samtools-1.3.1/bin/samtools index ${dirLoc}/aligned_bam/${lib}.q30.sorted.bam

# Remove intermediate files
rm ${dirLoc}/aligned_bam/${lib}.bam
rm ${dirLoc}/aligned_bam/${lib}.nofilter.bam
rm ${dirLoc}/aligned_bam/${lib}.q30.bam