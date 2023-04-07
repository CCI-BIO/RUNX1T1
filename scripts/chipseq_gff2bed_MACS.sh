#!/bin/bash
## Script to convert MACS2 GFF files to bed files compatible with homer. 
##

## Date: 17/10/2018
## Last updated: 13/09/2021

# Specify necessary parameters.
dirLoc="/data/RUNX_KO_ChIPseq"

for f in ${dirLoc}/MACS2_broad_results/*.gff
do
	fq=( $(echo $f | cut -d '/' -f 6 | cut -d '.' -f 1) )
	awk -F"\t" '{OFS="\t"};{print $1,$4,$5,$2,".","."}' $f > ${dirLoc}/MACS2_broad_results/${fq}_macs.bed
done

for f in ${dirLoc}/MACS2_narrow_results/*.gff
do
	fq=( $(echo $f | cut -d '/' -f 6 | cut -d '.' -f 1) )
	awk -F"\t" '{OFS="\t"};{print $1,$4,$5,$2,".","."}' $f > ${dirLoc}/MACS2_narrow_results/${fq}_macs.bed
done