#!/bin/bash

#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

#SBATCH --account=def-stbil30

module load mugqic/R_Bioconductor
Rscript scripts/metagene/metagene_cofactor/script_metagene.obj_GRbound30m_annotated.R
