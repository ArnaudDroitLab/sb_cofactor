#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

#-------------------------------------------------------------------------------
# RnaSeq SlurmScheduler Job Submission Bash script
# Version: 3.0.1-beta
# Created on: 2018-04-11T20:22:21
# Steps:
#   picard_sam_to_fastq: 0 job... skipping
#   trimmomatic: 4 jobs
#   merge_trimmomatic_stats: 1 job
#   star: 10 jobs
#   picard_merge_sam_files: 0 job... skipping
#   picard_sort_sam: 4 jobs
#   picard_mark_duplicates: 4 jobs
#   TOTAL: 23 jobs
#-------------------------------------------------------------------------------

OUTPUT_DIR=/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38
JOB_OUTPUT_DIR=$OUTPUT_DIR/job_output
TIMESTAMP=`date +%FT%H.%M.%S`
JOB_LIST=$JOB_OUTPUT_DIR/RnaSeq_job_list_$TIMESTAMP
export CONFIG_FILES="/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini"
mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR


#-------------------------------------------------------------------------------
# STEP: trimmomatic
#-------------------------------------------------------------------------------
STEP=trimmomatic
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: trimmomatic_1_JOB_ID: trimmomatic.HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2.359e70586ee19d7949bb067aa4ff2a17.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2.359e70586ee19d7949bb067aa4ff2a17.mugqic.done'
module load java/1.8.0_121 mugqic/trimmomatic/0.36 && \
mkdir -p trim/DEX_shNIPBL_2 && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR PE \
  -threads 6 \
  -phred33 \
  /project/6001942/Working_Directory/Eric/CofactorHR/A549/raw/rna-seq/HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2_R1.fastq.gz \
  /project/6001942/Working_Directory/Eric/CofactorHR/A549/raw/rna-seq/HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2_R2.fastq.gz \
  trim/DEX_shNIPBL_2/HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2.trim.pair1.fastq.gz \
  trim/DEX_shNIPBL_2/HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2.trim.single1.fastq.gz \
  trim/DEX_shNIPBL_2/HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2.trim.pair2.fastq.gz \
  trim/DEX_shNIPBL_2/HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2.trim.single2.fastq.gz \
  ILLUMINACLIP:/cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.0/bfx/adapters-truseq.fa:2:30:15:8:true \
  TRAILING:30 \
  MINLEN:32 \
  2> trim/DEX_shNIPBL_2/HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2.trim.log
trimmomatic.HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2.359e70586ee19d7949bb067aa4ff2a17.mugqic.done
)
trimmomatic_1_JOB_ID=$(echo "#! /bin/bash 
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
  -s \"trimmomatic\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"trimmomatic\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=24G -N 1 -n 6 | grep "[0-9]" | cut -d\  -f4)
echo "$trimmomatic_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: trimmomatic_2_JOB_ID: trimmomatic.HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2.48f771a71f5fc25a5cca9a963564e594.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2.48f771a71f5fc25a5cca9a963564e594.mugqic.done'
module load java/1.8.0_121 mugqic/trimmomatic/0.36 && \
mkdir -p trim/DEX_shMED1_2 && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR PE \
  -threads 6 \
  -phred33 \
  /project/6001942/Working_Directory/Eric/CofactorHR/A549/raw/rna-seq/HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2_R1.fastq.gz \
  /project/6001942/Working_Directory/Eric/CofactorHR/A549/raw/rna-seq/HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2_R2.fastq.gz \
  trim/DEX_shMED1_2/HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2.trim.pair1.fastq.gz \
  trim/DEX_shMED1_2/HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2.trim.single1.fastq.gz \
  trim/DEX_shMED1_2/HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2.trim.pair2.fastq.gz \
  trim/DEX_shMED1_2/HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2.trim.single2.fastq.gz \
  ILLUMINACLIP:/cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.0/bfx/adapters-truseq.fa:2:30:15:8:true \
  TRAILING:30 \
  MINLEN:32 \
  2> trim/DEX_shMED1_2/HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2.trim.log
trimmomatic.HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2.48f771a71f5fc25a5cca9a963564e594.mugqic.done
)
trimmomatic_2_JOB_ID=$(echo "#! /bin/bash 
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
  -s \"trimmomatic\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"trimmomatic\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=24G -N 1 -n 6 | grep "[0-9]" | cut -d\  -f4)
echo "$trimmomatic_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: trimmomatic_3_JOB_ID: trimmomatic.HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2.87c3ed68e5869bd64655ed98ba4c2162.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2.87c3ed68e5869bd64655ed98ba4c2162.mugqic.done'
module load java/1.8.0_121 mugqic/trimmomatic/0.36 && \
mkdir -p trim/DEX_shSMC1A_4 && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR PE \
  -threads 6 \
  -phred33 \
  /project/6001942/Working_Directory/Eric/CofactorHR/A549/raw/rna-seq/HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2_R1.fastq.gz \
  /project/6001942/Working_Directory/Eric/CofactorHR/A549/raw/rna-seq/HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2_R2.fastq.gz \
  trim/DEX_shSMC1A_4/HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2.trim.pair1.fastq.gz \
  trim/DEX_shSMC1A_4/HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2.trim.single1.fastq.gz \
  trim/DEX_shSMC1A_4/HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2.trim.pair2.fastq.gz \
  trim/DEX_shSMC1A_4/HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2.trim.single2.fastq.gz \
  ILLUMINACLIP:/cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.0/bfx/adapters-truseq.fa:2:30:15:8:true \
  TRAILING:30 \
  MINLEN:32 \
  2> trim/DEX_shSMC1A_4/HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2.trim.log
trimmomatic.HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2.87c3ed68e5869bd64655ed98ba4c2162.mugqic.done
)
trimmomatic_3_JOB_ID=$(echo "#! /bin/bash 
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
  -s \"trimmomatic\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"trimmomatic\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=24G -N 1 -n 6 | grep "[0-9]" | cut -d\  -f4)
echo "$trimmomatic_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: trimmomatic_4_JOB_ID: trimmomatic.HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2.4079b0bac62097b02378879ce3e6d2c5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2.4079b0bac62097b02378879ce3e6d2c5.mugqic.done'
module load java/1.8.0_121 mugqic/trimmomatic/0.36 && \
mkdir -p trim/DEX_nonMamm_2 && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR PE \
  -threads 6 \
  -phred33 \
  /project/6001942/Working_Directory/Eric/CofactorHR/A549/raw/rna-seq/HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2_R1.fastq.gz \
  /project/6001942/Working_Directory/Eric/CofactorHR/A549/raw/rna-seq/HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2_R2.fastq.gz \
  trim/DEX_nonMamm_2/HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2.trim.pair1.fastq.gz \
  trim/DEX_nonMamm_2/HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2.trim.single1.fastq.gz \
  trim/DEX_nonMamm_2/HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2.trim.pair2.fastq.gz \
  trim/DEX_nonMamm_2/HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2.trim.single2.fastq.gz \
  ILLUMINACLIP:/cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.0/bfx/adapters-truseq.fa:2:30:15:8:true \
  TRAILING:30 \
  MINLEN:32 \
  2> trim/DEX_nonMamm_2/HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2.trim.log
trimmomatic.HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2.4079b0bac62097b02378879ce3e6d2c5.mugqic.done
)
trimmomatic_4_JOB_ID=$(echo "#! /bin/bash 
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
  -s \"trimmomatic\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"trimmomatic\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=24G -N 1 -n 6 | grep "[0-9]" | cut -d\  -f4)
echo "$trimmomatic_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# STEP: merge_trimmomatic_stats
#-------------------------------------------------------------------------------
STEP=merge_trimmomatic_stats
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: merge_trimmomatic_stats_1_JOB_ID: merge_trimmomatic_stats
#-------------------------------------------------------------------------------
JOB_NAME=merge_trimmomatic_stats
JOB_DEPENDENCIES=$trimmomatic_1_JOB_ID:$trimmomatic_2_JOB_ID:$trimmomatic_3_JOB_ID:$trimmomatic_4_JOB_ID
JOB_DONE=job_output/merge_trimmomatic_stats/merge_trimmomatic_stats.1801b2bb6902c0d7c8374d34eab2a8e2.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'merge_trimmomatic_stats.1801b2bb6902c0d7c8374d34eab2a8e2.mugqic.done'
module load mugqic/pandoc/1.15.2 && \
mkdir -p metrics && \
echo 'Sample	Readset	Raw Paired Reads #	Surviving Paired Reads #	Surviving Paired Reads %' > metrics/trimReadsetTable.tsv && \
grep ^Input trim/DEX_shNIPBL_2/HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2.trim.log | \
perl -pe 's/^Input Read Pairs: (\d+).*Both Surviving: (\d+).*Forward Only Surviving: (\d+).*$/DEX_shNIPBL_2	HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/DEX_shMED1_2/HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2.trim.log | \
perl -pe 's/^Input Read Pairs: (\d+).*Both Surviving: (\d+).*Forward Only Surviving: (\d+).*$/DEX_shMED1_2	HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/DEX_shSMC1A_4/HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2.trim.log | \
perl -pe 's/^Input Read Pairs: (\d+).*Both Surviving: (\d+).*Forward Only Surviving: (\d+).*$/DEX_shSMC1A_4	HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/DEX_nonMamm_2/HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2.trim.log | \
perl -pe 's/^Input Read Pairs: (\d+).*Both Surviving: (\d+).*Forward Only Surviving: (\d+).*$/DEX_nonMamm_2	HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
cut -f1,3- metrics/trimReadsetTable.tsv | awk -F"	" '{OFS="	"; if (NR==1) {if ($2=="Raw Paired Reads #") {paired=1};print "Sample", "Raw Reads #", "Surviving Reads #", "Surviving %"} else {if (paired) {$2=$2*2; $3=$3*2}; raw[$1]+=$2; surviving[$1]+=$3}}END{for (sample in raw){print sample, raw[sample], surviving[sample], surviving[sample] / raw[sample] * 100}}' \
  > metrics/trimSampleTable.tsv && \
mkdir -p report && \
cp metrics/trimReadsetTable.tsv metrics/trimSampleTable.tsv report/ && \
trim_readset_table_md=`LC_NUMERIC=en_CA awk -F "	" '{OFS="|"; if (NR == 1) {$1 = $1; print $0; print "-----|-----|-----:|-----:|-----:"} else {print $1, $2, sprintf("%\47d", $3), sprintf("%\47d", $4), sprintf("%.1f", $5)}}' metrics/trimReadsetTable.tsv` && \
pandoc \
  /home/efournie/genpipes/bfx/report/Illumina.merge_trimmomatic_stats.md \
  --template /home/efournie/genpipes/bfx/report/Illumina.merge_trimmomatic_stats.md \
  --variable trailing_min_quality=30 \
  --variable min_length=32 \
  --variable read_type=Paired \
  --variable trim_readset_table="$trim_readset_table_md" \
  --to markdown \
  > report/Illumina.merge_trimmomatic_stats.md
merge_trimmomatic_stats.1801b2bb6902c0d7c8374d34eab2a8e2.mugqic.done
)
merge_trimmomatic_stats_1_JOB_ID=$(echo "#! /bin/bash 
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
  -s \"merge_trimmomatic_stats\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"merge_trimmomatic_stats\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$merge_trimmomatic_stats_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# STEP: star
#-------------------------------------------------------------------------------
STEP=star
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: star_1_JOB_ID: star_align.1.HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2
#-------------------------------------------------------------------------------
JOB_NAME=star_align.1.HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2
JOB_DEPENDENCIES=$trimmomatic_1_JOB_ID
JOB_DONE=job_output/star/star_align.1.HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2.6c8b70b729c6e407289561b0c7ccb23f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.1.HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2.6c8b70b729c6e407289561b0c7ccb23f.mugqic.done'
module load mugqic/star/2.5.3a && \
mkdir -p alignment_1stPass/DEX_shNIPBL_2/HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2 && \
STAR --runMode alignReads \
  --genomeDir /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/star_index/Ensembl87.sjdbOverhang99 \
  --readFilesIn \
    trim/DEX_shNIPBL_2/HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2.trim.pair1.fastq.gz \
    trim/DEX_shNIPBL_2/HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix alignment_1stPass/DEX_shNIPBL_2/HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2/ \
  --outSAMattrRGline ID:"HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2" 	PL:"ILLUMINA" 			SM:"DEX_shNIPBL_2" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 100000000000 \
  --limitIObufferSize 4000000000
star_align.1.HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2.6c8b70b729c6e407289561b0c7ccb23f.mugqic.done
)
star_1_JOB_ID=$(echo "#! /bin/bash 
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
  -s \"star\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"star\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=128G -N 1 -n 16 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$star_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: star_2_JOB_ID: star_align.1.HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2
#-------------------------------------------------------------------------------
JOB_NAME=star_align.1.HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2
JOB_DEPENDENCIES=$trimmomatic_2_JOB_ID
JOB_DONE=job_output/star/star_align.1.HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2.8ed8cea0f42c4242fdc4e877f1cbf6c0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.1.HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2.8ed8cea0f42c4242fdc4e877f1cbf6c0.mugqic.done'
module load mugqic/star/2.5.3a && \
mkdir -p alignment_1stPass/DEX_shMED1_2/HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2 && \
STAR --runMode alignReads \
  --genomeDir /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/star_index/Ensembl87.sjdbOverhang99 \
  --readFilesIn \
    trim/DEX_shMED1_2/HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2.trim.pair1.fastq.gz \
    trim/DEX_shMED1_2/HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix alignment_1stPass/DEX_shMED1_2/HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2/ \
  --outSAMattrRGline ID:"HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2" 	PL:"ILLUMINA" 			SM:"DEX_shMED1_2" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 100000000000 \
  --limitIObufferSize 4000000000
star_align.1.HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2.8ed8cea0f42c4242fdc4e877f1cbf6c0.mugqic.done
)
star_2_JOB_ID=$(echo "#! /bin/bash 
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
  -s \"star\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"star\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=128G -N 1 -n 16 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$star_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: star_3_JOB_ID: star_align.1.HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2
#-------------------------------------------------------------------------------
JOB_NAME=star_align.1.HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2
JOB_DEPENDENCIES=$trimmomatic_3_JOB_ID
JOB_DONE=job_output/star/star_align.1.HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2.ad9d1370e23371ef291c0eebb78f1003.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.1.HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2.ad9d1370e23371ef291c0eebb78f1003.mugqic.done'
module load mugqic/star/2.5.3a && \
mkdir -p alignment_1stPass/DEX_shSMC1A_4/HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2 && \
STAR --runMode alignReads \
  --genomeDir /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/star_index/Ensembl87.sjdbOverhang99 \
  --readFilesIn \
    trim/DEX_shSMC1A_4/HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2.trim.pair1.fastq.gz \
    trim/DEX_shSMC1A_4/HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix alignment_1stPass/DEX_shSMC1A_4/HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2/ \
  --outSAMattrRGline ID:"HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2" 	PL:"ILLUMINA" 			SM:"DEX_shSMC1A_4" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 100000000000 \
  --limitIObufferSize 4000000000
star_align.1.HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2.ad9d1370e23371ef291c0eebb78f1003.mugqic.done
)
star_3_JOB_ID=$(echo "#! /bin/bash 
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
  -s \"star\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"star\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=128G -N 1 -n 16 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$star_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: star_4_JOB_ID: star_align.1.HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2
#-------------------------------------------------------------------------------
JOB_NAME=star_align.1.HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2
JOB_DEPENDENCIES=$trimmomatic_4_JOB_ID
JOB_DONE=job_output/star/star_align.1.HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2.1b65e8960276dd91b240760cb5102985.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.1.HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2.1b65e8960276dd91b240760cb5102985.mugqic.done'
module load mugqic/star/2.5.3a && \
mkdir -p alignment_1stPass/DEX_nonMamm_2/HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2 && \
STAR --runMode alignReads \
  --genomeDir /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/star_index/Ensembl87.sjdbOverhang99 \
  --readFilesIn \
    trim/DEX_nonMamm_2/HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2.trim.pair1.fastq.gz \
    trim/DEX_nonMamm_2/HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix alignment_1stPass/DEX_nonMamm_2/HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2/ \
  --outSAMattrRGline ID:"HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2" 	PL:"ILLUMINA" 			SM:"DEX_nonMamm_2" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 100000000000 \
  --limitIObufferSize 4000000000
star_align.1.HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2.1b65e8960276dd91b240760cb5102985.mugqic.done
)
star_4_JOB_ID=$(echo "#! /bin/bash 
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
  -s \"star\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"star\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=128G -N 1 -n 16 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$star_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: star_5_JOB_ID: star_index.AllSamples
#-------------------------------------------------------------------------------
JOB_NAME=star_index.AllSamples
JOB_DEPENDENCIES=$star_1_JOB_ID:$star_2_JOB_ID:$star_3_JOB_ID:$star_4_JOB_ID
JOB_DONE=job_output/star/star_index.AllSamples.426238a228653d9344ee0d21721f12c7.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_index.AllSamples.426238a228653d9344ee0d21721f12c7.mugqic.done'
module load mugqic/star/2.5.3a && \
cat \
  alignment_1stPass/DEX_shNIPBL_2/HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2/SJ.out.tab \
  alignment_1stPass/DEX_shMED1_2/HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2/SJ.out.tab \
  alignment_1stPass/DEX_shSMC1A_4/HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2/SJ.out.tab \
  alignment_1stPass/DEX_nonMamm_2/HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2/SJ.out.tab | \
awk 'BEGIN {OFS="	"; strChar[0]="."; strChar[1]="+"; strChar[2]="-"} {if($5>0){print $1,$2,$3,strChar[$4]}}' | sort -k1,1h -k2,2n > alignment_1stPass/AllSamples.SJ.out.tab && \
mkdir -p reference.Merged && \
STAR --runMode genomeGenerate \
  --genomeDir reference.Merged \
  --genomeFastaFiles /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --runThreadN 16 \
  --limitGenomeGenerateRAM 100000000000 \
  --sjdbFileChrStartEnd alignment_1stPass/AllSamples.SJ.out.tab \
  --limitIObufferSize 1000000000 \
  --sjdbOverhang 99
star_index.AllSamples.426238a228653d9344ee0d21721f12c7.mugqic.done
)
star_5_JOB_ID=$(echo "#! /bin/bash 
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
  -s \"star\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"star\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=15:00:0 --mem=128G -N 1 -n 16 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$star_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: star_6_JOB_ID: star_align.2.HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2
#-------------------------------------------------------------------------------
JOB_NAME=star_align.2.HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2
JOB_DEPENDENCIES=$trimmomatic_1_JOB_ID:$star_5_JOB_ID
JOB_DONE=job_output/star/star_align.2.HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2.8ecff57d5be6a2c0eef193cfb1f6a8be.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.2.HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2.8ecff57d5be6a2c0eef193cfb1f6a8be.mugqic.done'
module load mugqic/star/2.5.3a && \
mkdir -p alignment/DEX_shNIPBL_2/HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2 && \
STAR --runMode alignReads \
  --genomeDir reference.Merged \
  --readFilesIn \
    trim/DEX_shNIPBL_2/HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2.trim.pair1.fastq.gz \
    trim/DEX_shNIPBL_2/HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix alignment/DEX_shNIPBL_2/HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2/ \
  --outSAMattrRGline ID:"HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2" 	PL:"ILLUMINA" 			SM:"DEX_shNIPBL_2" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 100000000000 \
  --limitBAMsortRAM 100000000000 \
  --limitIObufferSize 4000000000 \
  --outWigType wiggle read1_5p --outWigStrand Stranded --outWigReferencesPrefix chr \
  --chimSegmentMin 21 && \
ln -s -f HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2/Aligned.sortedByCoord.out.bam alignment/DEX_shNIPBL_2/DEX_shNIPBL_2.sorted.bam
star_align.2.HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2.8ecff57d5be6a2c0eef193cfb1f6a8be.mugqic.done
)
star_6_JOB_ID=$(echo "#! /bin/bash 
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
  -s \"star\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"star\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=128G -N 1 -n 16 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$star_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: star_7_JOB_ID: star_align.2.HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2
#-------------------------------------------------------------------------------
JOB_NAME=star_align.2.HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2
JOB_DEPENDENCIES=$trimmomatic_2_JOB_ID:$star_5_JOB_ID
JOB_DONE=job_output/star/star_align.2.HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2.6594585bda9a33609823ce56616d7fb9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.2.HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2.6594585bda9a33609823ce56616d7fb9.mugqic.done'
module load mugqic/star/2.5.3a && \
mkdir -p alignment/DEX_shMED1_2/HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2 && \
STAR --runMode alignReads \
  --genomeDir reference.Merged \
  --readFilesIn \
    trim/DEX_shMED1_2/HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2.trim.pair1.fastq.gz \
    trim/DEX_shMED1_2/HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix alignment/DEX_shMED1_2/HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2/ \
  --outSAMattrRGline ID:"HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2" 	PL:"ILLUMINA" 			SM:"DEX_shMED1_2" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 100000000000 \
  --limitBAMsortRAM 100000000000 \
  --limitIObufferSize 4000000000 \
  --outWigType wiggle read1_5p --outWigStrand Stranded --outWigReferencesPrefix chr \
  --chimSegmentMin 21 && \
ln -s -f HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2/Aligned.sortedByCoord.out.bam alignment/DEX_shMED1_2/DEX_shMED1_2.sorted.bam
star_align.2.HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2.6594585bda9a33609823ce56616d7fb9.mugqic.done
)
star_7_JOB_ID=$(echo "#! /bin/bash 
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
  -s \"star\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"star\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=128G -N 1 -n 16 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$star_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: star_8_JOB_ID: star_align.2.HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2
#-------------------------------------------------------------------------------
JOB_NAME=star_align.2.HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2
JOB_DEPENDENCIES=$trimmomatic_3_JOB_ID:$star_5_JOB_ID
JOB_DONE=job_output/star/star_align.2.HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2.679f18c9ea26a0352c8a4c202b14c469.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.2.HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2.679f18c9ea26a0352c8a4c202b14c469.mugqic.done'
module load mugqic/star/2.5.3a && \
mkdir -p alignment/DEX_shSMC1A_4/HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2 && \
STAR --runMode alignReads \
  --genomeDir reference.Merged \
  --readFilesIn \
    trim/DEX_shSMC1A_4/HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2.trim.pair1.fastq.gz \
    trim/DEX_shSMC1A_4/HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix alignment/DEX_shSMC1A_4/HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2/ \
  --outSAMattrRGline ID:"HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2" 	PL:"ILLUMINA" 			SM:"DEX_shSMC1A_4" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 100000000000 \
  --limitBAMsortRAM 100000000000 \
  --limitIObufferSize 4000000000 \
  --outWigType wiggle read1_5p --outWigStrand Stranded --outWigReferencesPrefix chr \
  --chimSegmentMin 21 && \
ln -s -f HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2/Aligned.sortedByCoord.out.bam alignment/DEX_shSMC1A_4/DEX_shSMC1A_4.sorted.bam
star_align.2.HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2.679f18c9ea26a0352c8a4c202b14c469.mugqic.done
)
star_8_JOB_ID=$(echo "#! /bin/bash 
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
  -s \"star\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"star\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=128G -N 1 -n 16 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$star_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: star_9_JOB_ID: star_align.2.HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2
#-------------------------------------------------------------------------------
JOB_NAME=star_align.2.HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2
JOB_DEPENDENCIES=$trimmomatic_4_JOB_ID:$star_5_JOB_ID
JOB_DONE=job_output/star/star_align.2.HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2.ec7dfcc78ec02487d7169d1202753a3f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.2.HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2.ec7dfcc78ec02487d7169d1202753a3f.mugqic.done'
module load mugqic/star/2.5.3a && \
mkdir -p alignment/DEX_nonMamm_2/HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2 && \
STAR --runMode alignReads \
  --genomeDir reference.Merged \
  --readFilesIn \
    trim/DEX_nonMamm_2/HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2.trim.pair1.fastq.gz \
    trim/DEX_nonMamm_2/HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix alignment/DEX_nonMamm_2/HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2/ \
  --outSAMattrRGline ID:"HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2" 	PL:"ILLUMINA" 			SM:"DEX_nonMamm_2" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 100000000000 \
  --limitBAMsortRAM 100000000000 \
  --limitIObufferSize 4000000000 \
  --outWigType wiggle read1_5p --outWigStrand Stranded --outWigReferencesPrefix chr \
  --chimSegmentMin 21 && \
ln -s -f HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2/Aligned.sortedByCoord.out.bam alignment/DEX_nonMamm_2/DEX_nonMamm_2.sorted.bam
star_align.2.HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2.ec7dfcc78ec02487d7169d1202753a3f.mugqic.done
)
star_9_JOB_ID=$(echo "#! /bin/bash 
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
  -s \"star\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"star\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=128G -N 1 -n 16 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$star_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: star_10_JOB_ID: star_report
#-------------------------------------------------------------------------------
JOB_NAME=star_report
JOB_DEPENDENCIES=$star_6_JOB_ID:$star_7_JOB_ID:$star_8_JOB_ID:$star_9_JOB_ID
JOB_DONE=job_output/star/star_report.20b949a56558f0d07a6baecd0297fd46.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_report.20b949a56558f0d07a6baecd0297fd46.mugqic.done'
module load mugqic/pandoc/1.15.2 && \
mkdir -p report && \
pandoc --to=markdown \
  --template /home/efournie/genpipes/bfx/report/RnaSeq.star.md \
  --variable scientific_name="Homo_sapiens" \
  --variable assembly="GRCh38" \
  /home/efournie/genpipes/bfx/report/RnaSeq.star.md \
  > report/RnaSeq.star.md
star_report.20b949a56558f0d07a6baecd0297fd46.mugqic.done
)
star_10_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE &&  $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$star_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# STEP: picard_sort_sam
#-------------------------------------------------------------------------------
STEP=picard_sort_sam
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: picard_sort_sam_1_JOB_ID: picard_sort_sam.DEX_shNIPBL_2
#-------------------------------------------------------------------------------
JOB_NAME=picard_sort_sam.DEX_shNIPBL_2
JOB_DEPENDENCIES=$star_6_JOB_ID
JOB_DONE=job_output/picard_sort_sam/picard_sort_sam.DEX_shNIPBL_2.b926075a49ce1786cc5b9b67137facce.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sort_sam.DEX_shNIPBL_2.b926075a49ce1786cc5b9b67137facce.mugqic.done'
module load java/1.8.0_121 mugqic/picard/2.9.0 && \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=4 -Dsamjdk.buffer_size=1048576 -Xmx40G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/DEX_shNIPBL_2/DEX_shNIPBL_2.sorted.bam \
 OUTPUT=alignment/DEX_shNIPBL_2/DEX_shNIPBL_2.QueryNameSorted.bam \
 SORT_ORDER=queryname \
 MAX_RECORDS_IN_RAM=5750000
picard_sort_sam.DEX_shNIPBL_2.b926075a49ce1786cc5b9b67137facce.mugqic.done
)
picard_sort_sam_1_JOB_ID=$(echo "#! /bin/bash 
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
  -s \"picard_sort_sam\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_sort_sam\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=48G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$picard_sort_sam_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: picard_sort_sam_2_JOB_ID: picard_sort_sam.DEX_shMED1_2
#-------------------------------------------------------------------------------
JOB_NAME=picard_sort_sam.DEX_shMED1_2
JOB_DEPENDENCIES=$star_7_JOB_ID
JOB_DONE=job_output/picard_sort_sam/picard_sort_sam.DEX_shMED1_2.a3aff1793dc437fb22660ac33ede1bd5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sort_sam.DEX_shMED1_2.a3aff1793dc437fb22660ac33ede1bd5.mugqic.done'
module load java/1.8.0_121 mugqic/picard/2.9.0 && \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=4 -Dsamjdk.buffer_size=1048576 -Xmx40G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/DEX_shMED1_2/DEX_shMED1_2.sorted.bam \
 OUTPUT=alignment/DEX_shMED1_2/DEX_shMED1_2.QueryNameSorted.bam \
 SORT_ORDER=queryname \
 MAX_RECORDS_IN_RAM=5750000
picard_sort_sam.DEX_shMED1_2.a3aff1793dc437fb22660ac33ede1bd5.mugqic.done
)
picard_sort_sam_2_JOB_ID=$(echo "#! /bin/bash 
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
  -s \"picard_sort_sam\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_sort_sam\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=48G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$picard_sort_sam_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: picard_sort_sam_3_JOB_ID: picard_sort_sam.DEX_shSMC1A_4
#-------------------------------------------------------------------------------
JOB_NAME=picard_sort_sam.DEX_shSMC1A_4
JOB_DEPENDENCIES=$star_8_JOB_ID
JOB_DONE=job_output/picard_sort_sam/picard_sort_sam.DEX_shSMC1A_4.021d4c59c9533a361aa9322c7769da04.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sort_sam.DEX_shSMC1A_4.021d4c59c9533a361aa9322c7769da04.mugqic.done'
module load java/1.8.0_121 mugqic/picard/2.9.0 && \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=4 -Dsamjdk.buffer_size=1048576 -Xmx40G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/DEX_shSMC1A_4/DEX_shSMC1A_4.sorted.bam \
 OUTPUT=alignment/DEX_shSMC1A_4/DEX_shSMC1A_4.QueryNameSorted.bam \
 SORT_ORDER=queryname \
 MAX_RECORDS_IN_RAM=5750000
picard_sort_sam.DEX_shSMC1A_4.021d4c59c9533a361aa9322c7769da04.mugqic.done
)
picard_sort_sam_3_JOB_ID=$(echo "#! /bin/bash 
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
  -s \"picard_sort_sam\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_sort_sam\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=48G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$picard_sort_sam_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: picard_sort_sam_4_JOB_ID: picard_sort_sam.DEX_nonMamm_2
#-------------------------------------------------------------------------------
JOB_NAME=picard_sort_sam.DEX_nonMamm_2
JOB_DEPENDENCIES=$star_9_JOB_ID
JOB_DONE=job_output/picard_sort_sam/picard_sort_sam.DEX_nonMamm_2.92fbf653003f66b8f917ca37fbdf4c68.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sort_sam.DEX_nonMamm_2.92fbf653003f66b8f917ca37fbdf4c68.mugqic.done'
module load java/1.8.0_121 mugqic/picard/2.9.0 && \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=4 -Dsamjdk.buffer_size=1048576 -Xmx40G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/DEX_nonMamm_2/DEX_nonMamm_2.sorted.bam \
 OUTPUT=alignment/DEX_nonMamm_2/DEX_nonMamm_2.QueryNameSorted.bam \
 SORT_ORDER=queryname \
 MAX_RECORDS_IN_RAM=5750000
picard_sort_sam.DEX_nonMamm_2.92fbf653003f66b8f917ca37fbdf4c68.mugqic.done
)
picard_sort_sam_4_JOB_ID=$(echo "#! /bin/bash 
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
  -s \"picard_sort_sam\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_sort_sam\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=48G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$picard_sort_sam_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# STEP: picard_mark_duplicates
#-------------------------------------------------------------------------------
STEP=picard_mark_duplicates
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_1_JOB_ID: picard_mark_duplicates.DEX_shNIPBL_2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.DEX_shNIPBL_2
JOB_DEPENDENCIES=$star_6_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.DEX_shNIPBL_2.f34787f1d3cdeda4a77108b09a2247c9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.DEX_shNIPBL_2.f34787f1d3cdeda4a77108b09a2247c9.mugqic.done'
module load java/1.8.0_121 mugqic/picard/2.9.0 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=4 -Dsamjdk.buffer_size=1048576 -Xmx14G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/DEX_shNIPBL_2/DEX_shNIPBL_2.sorted.bam \
 OUTPUT=alignment/DEX_shNIPBL_2/DEX_shNIPBL_2.sorted.mdup.bam \
 METRICS_FILE=alignment/DEX_shNIPBL_2/DEX_shNIPBL_2.sorted.mdup.metrics \
 MAX_RECORDS_IN_RAM=3500000
picard_mark_duplicates.DEX_shNIPBL_2.f34787f1d3cdeda4a77108b09a2247c9.mugqic.done
)
picard_mark_duplicates_1_JOB_ID=$(echo "#! /bin/bash 
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
  -s \"picard_mark_duplicates\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_mark_duplicates\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=48:00:0 --mem=20G -N 1 -n 5 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$picard_mark_duplicates_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_2_JOB_ID: picard_mark_duplicates.DEX_shMED1_2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.DEX_shMED1_2
JOB_DEPENDENCIES=$star_7_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.DEX_shMED1_2.33ac3b331204ec14c3095208b6a8222f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.DEX_shMED1_2.33ac3b331204ec14c3095208b6a8222f.mugqic.done'
module load java/1.8.0_121 mugqic/picard/2.9.0 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=4 -Dsamjdk.buffer_size=1048576 -Xmx14G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/DEX_shMED1_2/DEX_shMED1_2.sorted.bam \
 OUTPUT=alignment/DEX_shMED1_2/DEX_shMED1_2.sorted.mdup.bam \
 METRICS_FILE=alignment/DEX_shMED1_2/DEX_shMED1_2.sorted.mdup.metrics \
 MAX_RECORDS_IN_RAM=3500000
picard_mark_duplicates.DEX_shMED1_2.33ac3b331204ec14c3095208b6a8222f.mugqic.done
)
picard_mark_duplicates_2_JOB_ID=$(echo "#! /bin/bash 
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
  -s \"picard_mark_duplicates\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_mark_duplicates\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=48:00:0 --mem=20G -N 1 -n 5 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$picard_mark_duplicates_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_3_JOB_ID: picard_mark_duplicates.DEX_shSMC1A_4
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.DEX_shSMC1A_4
JOB_DEPENDENCIES=$star_8_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.DEX_shSMC1A_4.b17ffd1ddda3856df37db5d10fbbca4d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.DEX_shSMC1A_4.b17ffd1ddda3856df37db5d10fbbca4d.mugqic.done'
module load java/1.8.0_121 mugqic/picard/2.9.0 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=4 -Dsamjdk.buffer_size=1048576 -Xmx14G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/DEX_shSMC1A_4/DEX_shSMC1A_4.sorted.bam \
 OUTPUT=alignment/DEX_shSMC1A_4/DEX_shSMC1A_4.sorted.mdup.bam \
 METRICS_FILE=alignment/DEX_shSMC1A_4/DEX_shSMC1A_4.sorted.mdup.metrics \
 MAX_RECORDS_IN_RAM=3500000
picard_mark_duplicates.DEX_shSMC1A_4.b17ffd1ddda3856df37db5d10fbbca4d.mugqic.done
)
picard_mark_duplicates_3_JOB_ID=$(echo "#! /bin/bash 
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
  -s \"picard_mark_duplicates\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_mark_duplicates\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=48:00:0 --mem=20G -N 1 -n 5 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$picard_mark_duplicates_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_4_JOB_ID: picard_mark_duplicates.DEX_nonMamm_2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.DEX_nonMamm_2
JOB_DEPENDENCIES=$star_9_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.DEX_nonMamm_2.7293865370216ee1634bd5fe0b85ffc6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.DEX_nonMamm_2.7293865370216ee1634bd5fe0b85ffc6.mugqic.done'
module load java/1.8.0_121 mugqic/picard/2.9.0 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=4 -Dsamjdk.buffer_size=1048576 -Xmx14G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/DEX_nonMamm_2/DEX_nonMamm_2.sorted.bam \
 OUTPUT=alignment/DEX_nonMamm_2/DEX_nonMamm_2.sorted.mdup.bam \
 METRICS_FILE=alignment/DEX_nonMamm_2/DEX_nonMamm_2.sorted.mdup.metrics \
 MAX_RECORDS_IN_RAM=3500000
picard_mark_duplicates.DEX_nonMamm_2.7293865370216ee1634bd5fe0b85ffc6.mugqic.done
)
picard_mark_duplicates_4_JOB_ID=$(echo "#! /bin/bash 
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
  -s \"picard_mark_duplicates\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_mark_duplicates\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=48:00:0 --mem=20G -N 1 -n 5 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$picard_mark_duplicates_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# Call home with pipeline statistics
#-------------------------------------------------------------------------------
wget "http://mugqic.hpc.mcgill.ca/cgi-bin/pipeline.cgi?hostname=cedar5.cedar.computecanada.ca&ip=206.12.124.6&pipeline=RnaSeq&steps=picard_sam_to_fastq,trimmomatic,merge_trimmomatic_stats,star,picard_merge_sam_files,picard_sort_sam,picard_mark_duplicates&samples=4&AnonymizedList=288bbdd32a25fd2689b78e7a576829ec,8befc711a91bdcd48cc37e62a68b2feb,cdb870911b7cbb4b8a8795618f3e80e9,74f036ab16e84e94b1e03a92ba25555d" --quiet --output-document=/dev/null

