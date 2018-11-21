#!/bin/bash

export JOB_MAIL="christophe.tav@gmail.com"

bash mugqic_utils/generate_bigwig_chip.sh \
	--slurm \
	--mugqicdir sb_cofactor_hr/GM12878/output/chip-pipeline-GRCh38 \
	--rapid def-stbil30
	
