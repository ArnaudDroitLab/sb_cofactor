#!/bin/bash

module load sra-toolkit/2.8.2-1

export OUTDIR="/home/chris11/scratch/sb_cofactor/MCF7/raw/chip-seq"

#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --account=def-stbil30

for sra in SRX490426S SRX490427 SRX490428 SRX490429 SRX490430 SRX490431
do
	echo "Downloading $sra"
	fastq-dump $sra -outdir $OUTDIR
done



