#!/bin/bash

export RAP_ID=def-stbil30
export JOB_MAIL="christophe.tav@gmail.com"
module load sra-toolkit/2.8.2-1

export OUTDIR="/home/chris11/scratch/sb_cofactor/MCF7/raw/chip-seq"

for sra in SRX490426S SRX490427 SRX490428 SRX490429 SRX490430 SRX490431
do
	echo "Downloading $sra"
	fastq-dump $sra -outdir $OUTDIR
done



