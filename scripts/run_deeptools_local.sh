#!/bin/bash

dirLoc="/RUNX1T1/RUNX_KO_ChIPseq"

fq=()

for f in ${dirLoc}/aligned_bam/*.q30.sorted.bam
do
        fq+=( $(echo $f | cut -d '/' -f 10 | cut -d '.' -f 1) )
done

mkdir ${dirLoc}/coverage_bw

for f in ${fq[@]}
do
	lib=$f
	bamCoverage --bam ${dirLoc}/aligned_bam/${lib}.q30.sorted.bam \
		-o ${dirLoc}/coverage_bw/${lib}.bw \
		--binSize 10
done



