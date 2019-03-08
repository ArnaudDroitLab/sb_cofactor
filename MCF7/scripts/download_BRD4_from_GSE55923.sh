#!/bin/bash

module load sra-toolkit/2.8.2-1

export OUTDIR="/home/chris11/scratch/sb_cofactor/MCF7/raw/chip-seq"

#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --account=def-stbil30

for sra in SRR1193526 SRR1193527 SRR1193528 SRR1193529 SRR1193530 SRR1193531 SRR1193562 SRR1193563
do
	echo "Downloading $sra"
	fastq-dump $sra -outdir $OUTDIR
done
