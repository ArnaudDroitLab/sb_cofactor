#!/bin/bash

export JOB_MAIL="christophe.tav@gmail.com"

output_mugqicdir="/home/chris11/scratch/sb_cofactor/K562/output/chip-pipeline-GRCh38"
mkdir -p $output_mugqicdir/tracks

bash mugqic_utils/generate_bigwig_chip.sh \
	--slurm \
	--mugqicdir $output_mugqicdir \
	--rapid def-stbil30
	
