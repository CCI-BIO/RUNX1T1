#!/bin/bash
## intersect_macs_peaks.sh
## Script to obtain intersection of bed files from merge replicates output.
##
## Date: 06/11/2018
## Last updated: 19/11/2021

dirLoc="/data/RUNX_KO_ChIPseq"

mkdir ${dirLoc}/consensus_overlaps

### SORT ALL BED FILES ###
sort -k1,1 -k2,2n ${dirLoc}/replicate_overlap_results/CONT_H3K27ac_final_overlapped_peak_table.bed > ${dirLoc}/consensus_overlaps/CONT_H3K27ac_consensus_peaks.sorted.bed
sort -k1,1 -k2,2n ${dirLoc}/replicate_overlap_results/CONT_H3K27me3_final_overlapped_peak_table.bed > ${dirLoc}/consensus_overlaps/CONT_H3K27me3_consensus_peaks.sorted.bed
sort -k1,1 -k2,2n ${dirLoc}/replicate_overlap_results/CONT_H3K4me1_final_overlapped_peak_table.bed > ${dirLoc}/consensus_overlaps/CONT_H3K4me1_consensus_peaks.sorted.bed
sort -k1,1 -k2,2n ${dirLoc}/replicate_overlap_results/CONT_H3K4me2_final_overlapped_peak_table.bed > ${dirLoc}/consensus_overlaps/CONT_H3K4me2_consensus_peaks.sorted.bed
sort -k1,1 -k2,2n ${dirLoc}/replicate_overlap_results/CONT_H3K4me3_final_overlapped_peak_table.bed > ${dirLoc}/consensus_overlaps/CONT_H3K4me3_consensus_peaks.sorted.bed
sort -k1,1 -k2,2n ${dirLoc}/replicate_overlap_results/POS_H3K27ac_final_overlapped_peak_table.bed > ${dirLoc}/consensus_overlaps/POS_H3K27ac_consensus_peaks.sorted.bed
sort -k1,1 -k2,2n ${dirLoc}/replicate_overlap_results/POS_H3K27me3_final_overlapped_peak_table.bed > ${dirLoc}/consensus_overlaps/POS_H3K27me3_consensus_peaks.sorted.bed
sort -k1,1 -k2,2n ${dirLoc}/replicate_overlap_results/POS_H3K4me1_final_overlapped_peak_table.bed > ${dirLoc}/consensus_overlaps/POS_H3K4me1_consensus_peaks.sorted.bed
sort -k1,1 -k2,2n ${dirLoc}/replicate_overlap_results/POS_H3K4me2_final_overlapped_peak_table.bed > ${dirLoc}/consensus_overlaps/POS_H3K4me2_consensus_peaks.sorted.bed
sort -k1,1 -k2,2n ${dirLoc}/replicate_overlap_results/POS_H3K4me3_final_overlapped_peak_table.bed > ${dirLoc}/consensus_overlaps/POS_H3K4me3_consensus_peaks.sorted.bed

# POS H3K27ac 
bedtools intersect -a ${dirLoc}/consensus_overlaps/POS_H3K27ac_consensus_peaks.sorted.bed \
	-b ${dirLoc}/consensus_overlaps/POS_H3K27me3_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/POS_H3K4me1_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/POS_H3K4me2_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/POS_H3K4me3_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/CONT_H3K27ac_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/CONT_H3K27me3_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/CONT_H3K4me1_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/CONT_H3K4me2_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/CONT_H3K4me3_consensus_peaks.sorted.bed \
	-names POS_H3K27me3 POS_H3K4me1 POS_H3K4me2 POS_H3K4me3 CONT_H3K27ac CONT_H3K27me3 CONT_H3K4me1 CONT_H3K4me2 CONT_H3K4me3 \
	-sorted -wa -wb > ${dirLoc}/consensus_overlaps/POS_H3K27ac_consensus_peaks_overlaps.bed

# POS H3K27me3
bedtools intersect -a ${dirLoc}/consensus_overlaps/POS_H3K27me3_consensus_peaks.sorted.bed \
	-b ${dirLoc}/consensus_overlaps/POS_H3K27ac_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/POS_H3K4me1_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/POS_H3K4me2_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/POS_H3K4me3_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/CONT_H3K27ac_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/CONT_H3K27me3_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/CONT_H3K4me1_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/CONT_H3K4me2_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/CONT_H3K4me3_consensus_peaks.sorted.bed \
	-names POS_H3K27ac POS_H3K4me1 POS_H3K4me2 POS_H3K4me3 CONT_H3K27ac CONT_H3K27me3 CONT_H3K4me1 CONT_H3K4me2 CONT_H3K4me3 \
	-sorted -wa -wb > ${dirLoc}/consensus_overlaps/POS_H3K27me3_consensus_peaks_overlaps.bed

# POS H3K4me1
bedtools intersect -a ${dirLoc}/consensus_overlaps/POS_H3K4me1_consensus_peaks.sorted.bed \
	-b ${dirLoc}/consensus_overlaps/POS_H3K27ac_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/POS_H3K27me3_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/POS_H3K4me2_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/POS_H3K4me3_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/CONT_H3K27ac_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/CONT_H3K27me3_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/CONT_H3K4me1_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/CONT_H3K4me2_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/CONT_H3K4me3_consensus_peaks.sorted.bed \
	-names POS_H3K27ac POS_H3K27me3 POS_H3K4me2 POS_H3K4me3 CONT_H3K27ac CONT_H3K27me3 CONT_H3K4me1 CONT_H3K4me2 CONT_H3K4me3 \
	-sorted -wa -wb > ${dirLoc}/consensus_overlaps/POS_H3K4me1_consensus_peaks_overlaps.bed

# POS H3K4me2
bedtools intersect -a ${dirLoc}/consensus_overlaps/POS_H3K4me2_consensus_peaks.sorted.bed \
	-b ${dirLoc}/consensus_overlaps/POS_H3K27ac_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/POS_H3K27me3_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/POS_H3K4me1_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/POS_H3K4me3_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/CONT_H3K27ac_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/CONT_H3K27me3_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/CONT_H3K4me1_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/CONT_H3K4me2_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/CONT_H3K4me3_consensus_peaks.sorted.bed \
	-names POS_H3K27ac POS_H3K27me3 POS_H3K4me1 POS_H3K4me3 CONT_H3K27ac CONT_H3K27me3 CONT_H3K4me1 CONT_H3K4me2 CONT_H3K4me3 \
	-sorted -wa -wb > ${dirLoc}/consensus_overlaps/POS_H3K4me2_consensus_peaks_overlaps.bed


# POS H3K4me3
bedtools intersect -a ${dirLoc}/consensus_overlaps/POS_H3K4me3_consensus_peaks.sorted.bed \
	-b ${dirLoc}/consensus_overlaps/POS_H3K27ac_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/POS_H3K27me3_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/POS_H3K4me1_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/POS_H3K4me2_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/CONT_H3K27ac_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/CONT_H3K27me3_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/CONT_H3K4me1_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/CONT_H3K4me2_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/CONT_H3K4me3_consensus_peaks.sorted.bed \
	-names POS_H3K27ac POS_H3K27me3 POS_H3K4me1 POS_H3K4me2 CONT_H3K27ac CONT_H3K27me3 CONT_H3K4me1 CONT_H3K4me2 CONT_H3K4me3 \
	-sorted -wa -wb > ${dirLoc}/consensus_overlaps/POS_H3K4me3_consensus_peaks_overlaps.bed

# CONT H3K27ac
bedtools intersect -a ${dirLoc}/consensus_overlaps/CONT_H3K27ac_consensus_peaks.sorted.bed \
	-b ${dirLoc}/consensus_overlaps/POS_H3K27ac_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/POS_H3K27me3_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/POS_H3K4me1_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/POS_H3K4me2_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/POS_H3K4me3_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/CONT_H3K27me3_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/CONT_H3K4me1_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/CONT_H3K4me2_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/CONT_H3K4me3_consensus_peaks.sorted.bed \
	-names POS_H3K27ac POS_H3K27me3 POS_H3K4me1 POS_H3K4me2 POS_H3K4me3 CONT_H3K27me3 CONT_H3K4me1 CONT_H3K4me2 CONT_H3K4me3 \
	-sorted -wa -wb > ${dirLoc}/consensus_overlaps/CONT_H3K27ac_consensus_peaks_overlaps.bed

# CONT H3K27me3
bedtools intersect -a ${dirLoc}/consensus_overlaps/CONT_H3K27me3_consensus_peaks.sorted.bed \
	-b ${dirLoc}/consensus_overlaps/POS_H3K27ac_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/POS_H3K27me3_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/POS_H3K4me1_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/POS_H3K4me2_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/POS_H3K4me3_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/CONT_H3K27ac_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/CONT_H3K4me1_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/CONT_H3K4me2_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/CONT_H3K4me3_consensus_peaks.sorted.bed \
	-names POS_H3K27ac POS_H3K27me3 POS_H3K4me1 POS_H3K4me2 POS_H3K4me3 CONT_H3K27ac CONT_H3K4me1 CONT_H3K4me2 CONT_H3K4me3 \
	-sorted -wa -wb > ${dirLoc}/consensus_overlaps/CONT_H3K27me3_consensus_peaks_overlaps.bed

# CONT H3K4me1
bedtools intersect -a ${dirLoc}/consensus_overlaps/CONT_H3K4me1_consensus_peaks.sorted.bed \
	-b ${dirLoc}/consensus_overlaps/POS_H3K27ac_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/POS_H3K27me3_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/POS_H3K4me1_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/POS_H3K4me2_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/POS_H3K4me3_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/CONT_H3K27ac_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/CONT_H3K27me3_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/CONT_H3K4me2_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/CONT_H3K4me3_consensus_peaks.sorted.bed \
	-names POS_H3K27ac POS_H3K27me3 POS_H3K4me1 POS_H3K4me2 POS_H3K4me3 CONT_H3K27ac CONT_H3K27me3 CONT_H3K4me2 CONT_H3K4me3 \
	-sorted -wa -wb > ${dirLoc}/consensus_overlaps/CONT_H3K4me1_consensus_peaks_overlaps.bed

# CONT H3K4me2
bedtools intersect -a ${dirLoc}/consensus_overlaps/CONT_H3K4me2_consensus_peaks.sorted.bed \
	-b ${dirLoc}/consensus_overlaps/POS_H3K27ac_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/POS_H3K27me3_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/POS_H3K4me1_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/POS_H3K4me2_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/POS_H3K4me3_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/CONT_H3K27ac_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/CONT_H3K27me3_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/CONT_H3K4me1_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/CONT_H3K4me3_consensus_peaks.sorted.bed \
	-names POS_H3K27ac POS_H3K27me3 POS_H3K4me1 POS_H3K4me2 POS_H3K4me3 CONT_H3K27ac CONT_H3K27me3 CONT_H3K4me1 CONT_H3K4me3 \
	-sorted -wa -wb > ${dirLoc}/consensus_overlaps/CONT_H3K4me2_consensus_peaks_overlaps.bed


# CONT H3K4me3
bedtools intersect -a ${dirLoc}/consensus_overlaps/CONT_H3K4me3_consensus_peaks.sorted.bed \
	-b ${dirLoc}/consensus_overlaps/POS_H3K27ac_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/POS_H3K27me3_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/POS_H3K4me1_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/POS_H3K4me2_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/POS_H3K4me3_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/CONT_H3K27ac_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/CONT_H3K27me3_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/CONT_H3K4me1_consensus_peaks.sorted.bed \
	${dirLoc}/consensus_overlaps/CONT_H3K4me2_consensus_peaks.sorted.bed \
	-names POS_H3K27ac POS_H3K27me3 POS_H3K4me1 POS_H3K4me2 POS_H3K4me3 CONT_H3K27ac CONT_H3K27me3 CONT_H3K4me1 CONT_H3K4me2 \
	-sorted -wa -wb > ${dirLoc}/consensus_overlaps/CONT_H3K4me3_consensus_peaks_overlaps.bed
