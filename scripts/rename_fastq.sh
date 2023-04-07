#!/bin/bash

dirLoc="/data/RUNX_KO_ChIPseq"

for i in rep1 rep2 rep3
do
	cp ${dirLoc}/fastq/fastq/${i}/CONT_H3K4me1.fastq.gz ${dirLoc}/fastq/CONT_H3K4me1_${i}.fastq.gz
	cp ${dirLoc}/fastq/fastq/${i}/CONT_H3K4me2.fastq.gz ${dirLoc}/fastq/CONT_H3K4me2_${i}.fastq.gz
	cp ${dirLoc}/fastq/fastq/${i}/CONT_H3K4me3.fastq.gz ${dirLoc}/fastq/CONT_H3K4me3_${i}.fastq.gz
	cp ${dirLoc}/fastq/fastq/${i}/CONT_H3K27ac.fastq.gz ${dirLoc}/fastq/CONT_H3K27ac_${i}.fastq.gz
	cp ${dirLoc}/fastq/fastq/${i}/CONT_H3K27me3.fastq.gz ${dirLoc}/fastq/CONT_H3K27me3_${i}.fastq.gz
	cp ${dirLoc}/fastq/fastq/${i}/CONT_input.fastq.gz ${dirLoc}/fastq/CONT_input_${i}.fastq.gz
	cp ${dirLoc}/fastq/fastq/${i}/POS_H3K4me1.fastq.gz ${dirLoc}/fastq/POS_H3K4me1_${i}.fastq.gz
	cp ${dirLoc}/fastq/fastq/${i}/POS_H3K4me2.fastq.gz ${dirLoc}/fastq/POS_H3K4me2_${i}.fastq.gz
	cp ${dirLoc}/fastq/fastq/${i}/POS_H3K4me3.fastq.gz ${dirLoc}/fastq/POS_H3K4me3_${i}.fastq.gz
	cp ${dirLoc}/fastq/fastq/${i}/POS_H3K27ac.fastq.gz ${dirLoc}/fastq/POS_H3K27ac_${i}.fastq.gz
	cp ${dirLoc}/fastq/fastq/${i}/POS_H3K27me3.fastq.gz ${dirLoc}/fastq/POS_H3K27me3_${i}.fastq.gz
	cp ${dirLoc}/fastq/fastq/${i}/POS_input.fastq.gz ${dirLoc}/fastq/POS_input_${i}.fastq.gz
done