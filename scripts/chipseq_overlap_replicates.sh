#!/bin/bash

dirLoc="/data/RUNX_KO_ChIPseq"

mkdir ${dirLoc}/loj_overlaps

for lib in H3K4me2 H3K4me3
do
	for treatment in POS CONT
	do
		bedtools intersect -a ${dirLoc}/MACS2_narrow_results/${treatment}_${lib}_rep1_q30_peaks_macs.bed \
		-b ${dirLoc}/MACS2_narrow_results/${treatment}_${lib}_rep2_q30_peaks_macs.bed \
		${dirLoc}/MACS2_narrow_results/${treatment}_${lib}_rep3_q30_peaks_macs.bed \
		-names rep2 rep3 -loj -sorted > ${dirLoc}/loj_overlaps/${treatment}_${lib}_rep1_loj_overlaps.bed

		bedtools intersect -a ${dirLoc}/MACS2_narrow_results/${treatment}_${lib}_rep2_q30_peaks_macs.bed \
		-b ${dirLoc}/MACS2_narrow_results/${treatment}_${lib}_rep1_q30_peaks_macs.bed \
		${dirLoc}/MACS2_narrow_results/${treatment}_${lib}_rep3_q30_peaks_macs.bed \
		-names rep1 rep3 -loj -sorted > ${dirLoc}/loj_overlaps/${treatment}_${lib}_rep2_loj_overlaps.bed

		bedtools intersect -a ${dirLoc}/MACS2_narrow_results/${treatment}_${lib}_rep3_q30_peaks_macs.bed \
		-b ${dirLoc}/MACS2_narrow_results/${treatment}_${lib}_rep1_q30_peaks_macs.bed \
		${dirLoc}/MACS2_narrow_results/${treatment}_${lib}_rep2_q30_peaks_macs.bed \
		-names rep1 rep2 -loj -sorted > ${dirLoc}/loj_overlaps/${treatment}_${lib}_rep3_loj_overlaps.bed
	done
done


for lib in H3K27ac H3K27me3 H3K4me1
do
	for treatment in POS CONT
	do
		bedtools intersect -a ${dirLoc}/MACS2_broad_results/${treatment}_${lib}_rep1_q30_peaks_macs.bed \
		-b ${dirLoc}/MACS2_broad_results/${treatment}_${lib}_rep2_q30_peaks_macs.bed \
		${dirLoc}/MACS2_broad_results/${treatment}_${lib}_rep3_q30_peaks_macs.bed \
		-names rep2 rep3 -loj -sorted > ${dirLoc}/loj_overlaps/${treatment}_${lib}_rep1_loj_overlaps.bed

		bedtools intersect -a ${dirLoc}/MACS2_broad_results/${treatment}_${lib}_rep2_q30_peaks_macs.bed \
		-b ${dirLoc}/MACS2_broad_results/${treatment}_${lib}_rep1_q30_peaks_macs.bed \
		${dirLoc}/MACS2_broad_results/${treatment}_${lib}_rep3_q30_peaks_macs.bed \
		-names rep1 rep3 -loj -sorted > ${dirLoc}/loj_overlaps/${treatment}_${lib}_rep2_loj_overlaps.bed

		bedtools intersect -a ${dirLoc}/MACS2_broad_results/${treatment}_${lib}_rep3_q30_peaks_macs.bed \
		-b ${dirLoc}/MACS2_broad_results/${treatment}_${lib}_rep1_q30_peaks_macs.bed \
		${dirLoc}/MACS2_broad_results/${treatment}_${lib}_rep2_q30_peaks_macs.bed \
		-names rep1 rep2 -loj -sorted > ${dirLoc}/loj_overlaps/${treatment}_${lib}_rep3_loj_overlaps.bed
	done
done
