#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

#-------------------------------------------------------------------------------
# RnaSeq SlurmScheduler Job Submission Bash script
# Version: 3.0.1-beta
# Created on: 2018-04-24T14:24:42
# Steps:
#   differential_expression: 1 job
#   TOTAL: 1 job
#-------------------------------------------------------------------------------

OUTPUT_DIR=/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38
JOB_OUTPUT_DIR=$OUTPUT_DIR/job_output
TIMESTAMP=`date +%FT%H.%M.%S`
JOB_LIST=$JOB_OUTPUT_DIR/RnaSeq_job_list_$TIMESTAMP
export CONFIG_FILES="/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini"
mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR


#-------------------------------------------------------------------------------
# STEP: differential_expression
#-------------------------------------------------------------------------------
STEP=differential_expression
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: differential_expression_1_JOB_ID: differential_expression
#-------------------------------------------------------------------------------
JOB_NAME=differential_expression
JOB_DEPENDENCIES=
JOB_DONE=job_output/differential_expression/differential_expression.2390a427c4bbd8b545ea3c42a141c9a8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'differential_expression.2390a427c4bbd8b545ea3c42a141c9a8.mugqic.done'
module load mugqic/mugqic_tools/2.1.9 mugqic/R_Bioconductor/3.4.2_3.6 && \
mkdir -p DGE && \
Rscript $R_TOOLS/edger.R \
  -d ../../raw/rna-seq/design.txt \
  -c DGE/rawCountMatrix.csv \
  -o DGE && \
Rscript $R_TOOLS/deseq.R \
  -d ../../raw/rna-seq/design.txt \
  -c DGE/rawCountMatrix.csv \
  -o DGE \
  
differential_expression.2390a427c4bbd8b545ea3c42a141c9a8.mugqic.done
)
differential_expression_1_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"differential_expression\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_1.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_1.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_1.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_3.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/ETOH_nonMamm_1.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/ETOH_nonMamm_2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/ETOH_shMED1_1.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/ETOH_shNIPBL_2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/ETOH_shNIPBL_1.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/ETOH_shMED1_2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/ETOH_shSMC1A_3.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/ETOH_shSMC1A_2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"differential_expression\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_1.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_1.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_1.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_3.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/ETOH_nonMamm_1.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/ETOH_nonMamm_2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/ETOH_shMED1_1.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/ETOH_shNIPBL_2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/ETOH_shNIPBL_1.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/ETOH_shMED1_2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/ETOH_shSMC1A_3.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/ETOH_shSMC1A_2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=10:00:0 --mem=4G -N 1 -n 1 | grep "[0-9]" | cut -d\  -f4)
echo "$differential_expression_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# Call home with pipeline statistics
#-------------------------------------------------------------------------------
wget "http://mugqic.hpc.mcgill.ca/cgi-bin/pipeline.cgi?hostname=cedar5.cedar.computecanada.ca&ip=206.12.124.6&pipeline=RnaSeq&steps=differential_expression&samples=16&AnonymizedList=dddd03cb6e65c1183657e5e5bf559f48,444c16ed9974bf0f01c9f5de6751c17b,15b56eae2cdf5203323e847cafca2478,9bb6d10c02e361fb068595ab4db4aed6,b55b936bea7aa8fa3185729524325737,cff800af6cafbbb70c28346717710b25,ba0cee9f1a8fbbefae102e42e845ca2b,277b756c8c3e10902132376d11c4df12,288bbdd32a25fd2689b78e7a576829ec,ac8379cdfc0119f8bb100d8925e216fd,6aa81e81b0a21275efefd2698d4b38d7,cdb870911b7cbb4b8a8795618f3e80e9,74f036ab16e84e94b1e03a92ba25555d,8befc711a91bdcd48cc37e62a68b2feb,e8f8bd8dfd005d92c279b40cd9c7e216,da8d2067113079941b7b9b4ded748341" --quiet --output-document=/dev/null

