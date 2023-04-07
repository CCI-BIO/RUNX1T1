#!/bin/bash/
## Script for submitting chipseq_MACS.sh jobs to the hpc.
## Based on job submission template.


## Date: 26/09/2018
## Last updated: 08/09/2021

# Specify necessary parameters for scripts.
dirLoc="/data/RUNX_KO_ChIPseq"

fq=()

for f in ${dirLoc}/fastq/*.fastq.gz
do
        fq+=( $(echo $f | cut -d '/' -f 6 | cut -d '.' -f 1) )
done

#fq=(${fq[@]//*trimmed*})
fq=(${fq[@]//*input*})

mkdir ${dirLoc}/MACS2_broad_results
mkdir ${dirLoc}/MACS2_narrow_results

for f in ${fq[@]}
do
	if [[ $f =~ "H3K4me1" ]]
	then
		lib=$f
		script2Run=$dirLoc"/scripts/chipseq_MACS_broad_peaks.sh"
		echo "Submitting $lib for broad peak detection"
		qsub $script2Run -N $lib"_MACS2" -v lib="$lib",dirLoc="$dirLoc"
	fi

	if [[ $f =~ "H3K4me2" ]]
	then
		lib=$f
		script2Run=$dirLoc"/scripts/chipseq_MACS_narrow_peaks.sh"
		echo "Submitting $lib for narrow peak detection"
		#qsub $script2Run -N $lib"_MACS2" -v lib="$lib",dirLoc="$dirLoc"
	fi	

	if [[ $f =~ "H3K4me3" ]]
	then
		lib=$f
		script2Run=$dirLoc"/scripts/chipseq_MACS_narrow_peaks.sh"
		echo "Submitting $lib for narrow peak detection"
		#qsub $script2Run -N $lib"_MACS2" -v lib="$lib",dirLoc="$dirLoc"
	fi

	if [[ $f =~ "H3K27ac" ]]
	then
		lib=$f
		script2Run=$dirLoc"/scripts/chipseq_MACS_broad_peaks.sh"
		echo "Submitting $lib for broad peak detection"
		qsub $script2Run -N $lib"_MACS2" -v lib="$lib",dirLoc="$dirLoc"
	fi

	if [[ $f =~ "H3K27me3" ]]
	then
		lib=$f
		script2Run=$dirLoc"/scripts/chipseq_MACS_broad_peaks.sh"
		echo "Submitting $lib for broad peak detection"
		qsub $script2Run -N $lib"_MACS2" -v lib="$lib",dirLoc="$dirLoc"
	fi

done
