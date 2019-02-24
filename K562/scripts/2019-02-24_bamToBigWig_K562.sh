#!/bin/bash

export JOB_MAIL="christophe.tav@gmail.com"

bash mugqic_utils/generate_bigwig_chip.sh \
	--slurm \
	--mugqicdir /home/chris11/scratch/sb_cofactor/K562/output/chip-pipeline-GRCh38 \
	--rapid def-stbil30
	
