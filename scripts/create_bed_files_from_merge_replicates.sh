#!/bin/bash
## Script to convert merge_replicate files to bed files
##
## Date: 18/11/2021
## Last updated: 18/11/2021

# Specify necessary parameters.
dirLoc="/data/RUNX_KO_ChIPseq"

for f in ${dirLoc}/replicate_overlap_results/*.txt
do
	fq=( $(echo $f | cut -d '/' -f 6 | cut -d '.' -f 1) )
	awk -F"\t" '{OFS="\t"};NR > 1 {print $1,$2,$3,$4,".","."}' $f > ${dirLoc}/replicate_overlap_results/${fq}.bed
done