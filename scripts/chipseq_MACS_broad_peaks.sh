#!/bin/bash
## Script to perform MACS modelling analysis of ChIP-seq data using bowtie2.
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
###############################################################################

cellLine=$(echo $lib | cut -d '_' -f 1)

rep=$(echo $lib | cut -d '_' -f 3)

module load macs2/2.1.1

# Run MACS2 on treatment/control pairs of BAM files. 
macs2 callpeak -t ${dirLoc}/aligned_bam/${lib}.nofilter.sorted.bam -c ${dirLoc}/aligned_bam/${cellLine}_input_${rep}.nofilter.sorted.bam -f BAM -g hs --broad --keep-dup auto --outdir ${dirLoc}/MACS2_broad_results -n ${lib}_nofilter -q 0.05 -B
# Run MACS2 on filtered treatment/control pairs of BAM files. 
macs2 callpeak -t ${dirLoc}/aligned_bam/${lib}.q30.sorted.bam -c ${dirLoc}/aligned_bam/${cellLine}_input_${rep}.q30.sorted.bam -f BAM -g hs --broad --keep-dup auto --outdir ${dirLoc}/MACS2_broad_results -n ${lib}_q30 -q 0.05 -B
# Create GFF files for ROSE. 
awk -v OFS="\t" 'NR > 29 {print $1, $9, ".", $2, $3, ".", ".", ".", $9}' ${dirLoc}/MACS2_broad_results/${lib}_nofilter_peaks.xls > ${dirLoc}/MACS2_broad_results/${lib}_nofilter_peaks.gff
awk -v OFS="\t" 'NR > 29 {print $1, $9, ".", $2, $3, ".", ".", ".", $9}' ${dirLoc}/MACS2_broad_results/${lib}_q30_peaks.xls > ${dirLoc}/MACS2_broad_results/${lib}_q30_peaks.gff
# Generate track files
macs2 bdgcmp -t ${dirLoc}/MACS2_broad_results/${lib}_nofilter_treat_pileup.bdg -c ${dirLoc}/MACS2_broad_results/${lib}_nofilter_control_lambda.bdg -o ${dirLoc}/MACS2_broad_results/${lib}_nofilter_FE_track.bdg -m FE
macs2 bdgcmp -t ${dirLoc}/MACS2_broad_results/${lib}_nofilter_treat_pileup.bdg -c ${dirLoc}/MACS2_broad_results/${lib}_nofilter_control_lambda.bdg -o ${dirLoc}/MACS2_broad_results/${lib}_nofilter_logLR_track.bdg -m logLR -p 0.00001
macs2 bdgcmp -t ${dirLoc}/MACS2_broad_results/${lib}_q30_treat_pileup.bdg -c ${dirLoc}/MACS2_broad_results/${lib}_q30_control_lambda.bdg -o ${dirLoc}/MACS2_broad_results/${lib}_q30_FE_track.bdg -m FE
macs2 bdgcmp -t ${dirLoc}/MACS2_broad_results/${lib}_q30_treat_pileup.bdg -c ${dirLoc}/MACS2_broad_results/${lib}_q30_control_lambda.bdg -o ${dirLoc}/MACS2_broad_results/${lib}_q30_logLR_track.bdg -m logLR -p 0.00001
