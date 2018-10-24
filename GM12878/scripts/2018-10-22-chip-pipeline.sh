#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

#-------------------------------------------------------------------------------
# ChipSeq SlurmScheduler Job Submission Bash script
# Version: 3.0.1-beta
# Created on: 2018-10-22T17:20:23
# Steps:
#   picard_sam_to_fastq: 10 jobs
#   trimmomatic: 10 jobs
#   merge_trimmomatic_stats: 1 job
#   bwa_mem_picard_sort_sam: 11 jobs
#   samtools_view_filter: 11 jobs
#   picard_merge_sam_files: 10 jobs
#   picard_mark_duplicates: 11 jobs
#   metrics: 2 jobs
#   homer_make_tag_directory: 10 jobs
#   qc_metrics: 1 job
#   homer_make_ucsc_file: 21 jobs
#   macs2_callpeak: 17 jobs
#   TOTAL: 115 jobs
#-------------------------------------------------------------------------------

OUTPUT_DIR=/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38
JOB_OUTPUT_DIR=$OUTPUT_DIR/job_output
TIMESTAMP=`date +%FT%H.%M.%S`
JOB_LIST=$JOB_OUTPUT_DIR/ChipSeq_job_list_$TIMESTAMP
export CONFIG_FILES="/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini"
mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR


#-------------------------------------------------------------------------------
# STEP: picard_sam_to_fastq
#-------------------------------------------------------------------------------
STEP=picard_sam_to_fastq
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_1_JOB_ID: picard_sam_to_fastq.HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.4aebe103211f7af7890fe06b8d687c27.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.4aebe103211f7af7890fe06b8d687c27.mugqic.done'
module load java/1.8.0_121 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/raw/HI.1385.006.Index_1.GM12878_WCE_IB_rep1.bam \
 FASTQ=/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/raw/HI.1385.006.Index_1.GM12878_WCE_IB_rep1.single.fastq.gz
picard_sam_to_fastq.HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.4aebe103211f7af7890fe06b8d687c27.mugqic.done
)
picard_sam_to_fastq_1_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_sam_to_fastq\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_WCE_rep1.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_sam_to_fastq\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_WCE_rep1.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=12G -N 1 -n 3 | grep "[0-9]" | cut -d\  -f4)
echo "$picard_sam_to_fastq_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_2_JOB_ID: picard_sam_to_fastq.HI.4552.003.Index_13.A549_ChIP_BRD4_NA_EtOH_0-01_1h_Rep2_R1
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.4552.003.Index_13.A549_ChIP_BRD4_NA_EtOH_0-01_1h_Rep2_R1
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.4552.003.Index_13.A549_ChIP_BRD4_NA_EtOH_0-01_1h_Rep2_R1.919395ac5b652fb381a215ba34b4ac23.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.4552.003.Index_13.A549_ChIP_BRD4_NA_EtOH_0-01_1h_Rep2_R1.919395ac5b652fb381a215ba34b4ac23.mugqic.done'
module load java/1.8.0_121 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/raw/HI.1385.006.Index_2.GM12878_NIPBL_IB_ChIP3_1_rep1.bam \
 FASTQ=/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/raw/HI.1385.006.Index_2.GM12878_NIPBL_IB_ChIP3_1_rep1.single.fastq.gz
picard_sam_to_fastq.HI.4552.003.Index_13.A549_ChIP_BRD4_NA_EtOH_0-01_1h_Rep2_R1.919395ac5b652fb381a215ba34b4ac23.mugqic.done
)
picard_sam_to_fastq_2_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_sam_to_fastq\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_NIPBL_rep1.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_sam_to_fastq\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_NIPBL_rep1.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=12G -N 1 -n 3 | grep "[0-9]" | cut -d\  -f4)
echo "$picard_sam_to_fastq_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_3_JOB_ID: picard_sam_to_fastq.HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.5e4b4f961fef956c5ba5f4f837c7eeed.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.5e4b4f961fef956c5ba5f4f837c7eeed.mugqic.done'
module load java/1.8.0_121 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/raw/HI.1385.006.Index_3.GM12878_MED1_IB_ChIP3_2_rep1.bam \
 FASTQ=/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/raw/HI.1385.006.Index_3.GM12878_MED1_IB_ChIP3_2_rep1.single.fastq.gz
picard_sam_to_fastq.HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.5e4b4f961fef956c5ba5f4f837c7eeed.mugqic.done
)
picard_sam_to_fastq_3_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_sam_to_fastq\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_MED1_rep1.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_sam_to_fastq\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_MED1_rep1.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=12G -N 1 -n 3 | grep "[0-9]" | cut -d\  -f4)
echo "$picard_sam_to_fastq_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_4_JOB_ID: picard_sam_to_fastq.HI.4552.003.Index_6.A549_ChIP_CDK9_NA_EtOH_0-01_1h_Rep2_R1
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.4552.003.Index_6.A549_ChIP_CDK9_NA_EtOH_0-01_1h_Rep2_R1
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.4552.003.Index_6.A549_ChIP_CDK9_NA_EtOH_0-01_1h_Rep2_R1.5309b0cf1bf8bbffa74ed5f4141edcbe.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.4552.003.Index_6.A549_ChIP_CDK9_NA_EtOH_0-01_1h_Rep2_R1.5309b0cf1bf8bbffa74ed5f4141edcbe.mugqic.done'
module load java/1.8.0_121 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/raw/HI.1385.006.Index_4.GM12878_SMC1_IB_ChIP3_3_rep1.bam \
 FASTQ=/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/raw/HI.1385.006.Index_4.GM12878_SMC1_IB_ChIP3_3_rep1.single.fastq.gz
picard_sam_to_fastq.HI.4552.003.Index_6.A549_ChIP_CDK9_NA_EtOH_0-01_1h_Rep2_R1.5309b0cf1bf8bbffa74ed5f4141edcbe.mugqic.done
)
picard_sam_to_fastq_4_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_sam_to_fastq\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_SMC1_rep1.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_sam_to_fastq\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_SMC1_rep1.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=12G -N 1 -n 3 | grep "[0-9]" | cut -d\  -f4)
echo "$picard_sam_to_fastq_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_5_JOB_ID: picard_sam_to_fastq.HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.6cbe80101d4abb086e15a49b8792fa68.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.6cbe80101d4abb086e15a49b8792fa68.mugqic.done'
module load java/1.8.0_121 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/raw/HI.1385.006.Index_5.GM12878_CDK9_IB_ChIP3_4_rep1.bam \
 FASTQ=/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/raw/HI.1385.006.Index_5.GM12878_CDK9_IB_ChIP3_4_rep1.single.fastq.gz
picard_sam_to_fastq.HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.6cbe80101d4abb086e15a49b8792fa68.mugqic.done
)
picard_sam_to_fastq_5_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_sam_to_fastq\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_CDK9_rep1.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_sam_to_fastq\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_CDK9_rep1.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=12G -N 1 -n 3 | grep "[0-9]" | cut -d\  -f4)
echo "$picard_sam_to_fastq_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_6_JOB_ID: picard_sam_to_fastq.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.4bb611380e1df217c5a08547a93db71c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.4bb611380e1df217c5a08547a93db71c.mugqic.done'
module load java/1.8.0_121 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/raw/HI.1997.006.Index_14.GM12878_WCE_IB_CHIP7_rep2.bam \
 FASTQ=/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/raw/HI.1997.006.Index_14.GM12878_WCE_IB_CHIP7_rep2.single.fastq.gz
picard_sam_to_fastq.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.4bb611380e1df217c5a08547a93db71c.mugqic.done
)
picard_sam_to_fastq_6_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_sam_to_fastq\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_WCE_rep2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_sam_to_fastq\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_WCE_rep2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=12G -N 1 -n 3 | grep "[0-9]" | cut -d\  -f4)
echo "$picard_sam_to_fastq_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_7_JOB_ID: picard_sam_to_fastq.HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam.7e1fb183a7ef10f88719fb796e4a25db.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam.7e1fb183a7ef10f88719fb796e4a25db.mugqic.done'
module load java/1.8.0_121 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/raw/HI.1997.006.Index_15.GM12878_SMC1_IB_CHIP7_rep2.bam \
 FASTQ=/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/raw/HI.1997.006.Index_15.GM12878_SMC1_IB_CHIP7_rep2.single.fastq.gz
picard_sam_to_fastq.HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam.7e1fb183a7ef10f88719fb796e4a25db.mugqic.done
)
picard_sam_to_fastq_7_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_sam_to_fastq\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_SMC1_rep2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_sam_to_fastq\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_SMC1_rep2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=12G -N 1 -n 3 | grep "[0-9]" | cut -d\  -f4)
echo "$picard_sam_to_fastq_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_8_JOB_ID: picard_sam_to_fastq.HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.55d0a49bcc7d95e698b3b880d5d0b5c1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.55d0a49bcc7d95e698b3b880d5d0b5c1.mugqic.done'
module load java/1.8.0_121 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/raw/HI.1997.006.Index_18.GM12878_CDK9_IB_CHIP7_rep2.bam \
 FASTQ=/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/raw/HI.1997.006.Index_18.GM12878_CDK9_IB_CHIP7_rep2.single.fastq.gz
picard_sam_to_fastq.HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.55d0a49bcc7d95e698b3b880d5d0b5c1.mugqic.done
)
picard_sam_to_fastq_8_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_sam_to_fastq\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_CDK9_rep2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_sam_to_fastq\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_CDK9_rep2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=12G -N 1 -n 3 | grep "[0-9]" | cut -d\  -f4)
echo "$picard_sam_to_fastq_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_9_JOB_ID: picard_sam_to_fastq.HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.7b916bfb8a621d8aea1adf839e0abc12.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.7b916bfb8a621d8aea1adf839e0abc12.mugqic.done'
module load java/1.8.0_121 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/raw/HI.1997.006.Index_6.GM12878_NIPBL_IB_CHIP7_rep2.bam \
 FASTQ=/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/raw/HI.1997.006.Index_6.GM12878_NIPBL_IB_CHIP7_rep2.single.fastq.gz
picard_sam_to_fastq.HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.7b916bfb8a621d8aea1adf839e0abc12.mugqic.done
)
picard_sam_to_fastq_9_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_sam_to_fastq\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_NIPBL_rep2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_sam_to_fastq\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_NIPBL_rep2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=12G -N 1 -n 3 | grep "[0-9]" | cut -d\  -f4)
echo "$picard_sam_to_fastq_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_10_JOB_ID: picard_sam_to_fastq.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.844da20a6f5061fb788b0fffe097804a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.844da20a6f5061fb788b0fffe097804a.mugqic.done'
module load java/1.8.0_121 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/raw/HI.1997.006.Index_7.GM12878_MED1_IB_CHIP7_rep2.bam \
 FASTQ=/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/raw/HI.1997.006.Index_7.GM12878_MED1_IB_CHIP7_rep2.single.fastq.gz
picard_sam_to_fastq.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.844da20a6f5061fb788b0fffe097804a.mugqic.done
)
picard_sam_to_fastq_10_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_sam_to_fastq\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_MED1_rep2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_sam_to_fastq\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_MED1_rep2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=12G -N 1 -n 3 | grep "[0-9]" | cut -d\  -f4)
echo "$picard_sam_to_fastq_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# STEP: trimmomatic
#-------------------------------------------------------------------------------
STEP=trimmomatic
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: trimmomatic_1_JOB_ID: trimmomatic.HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_1_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.90ef8c34b69ea9af40e26a04bc45dd3c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.90ef8c34b69ea9af40e26a04bc45dd3c.mugqic.done'
module load java/1.8.0_121 mugqic/trimmomatic/0.36 && \
mkdir -p trim/GM12878_WCE_rep1 && \
`cat > trim/GM12878_WCE_rep1/HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /project/6001942/Working_Directory/Eric/CofactorHR/GM12878/raw/HI.1385.006.Index_1.GM12878_WCE_IB_rep1.single.fastq.gz \
  trim/GM12878_WCE_rep1/HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/GM12878_WCE_rep1/HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/GM12878_WCE_rep1/HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.trim.log
trimmomatic.HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.90ef8c34b69ea9af40e26a04bc45dd3c.mugqic.done
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
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"trimmomatic\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_WCE_rep1.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"trimmomatic\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_WCE_rep1.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=24G -N 1 -n 6 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$trimmomatic_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: trimmomatic_2_JOB_ID: trimmomatic.HI.4552.003.Index_13.A549_ChIP_BRD4_NA_EtOH_0-01_1h_Rep2_R1
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.4552.003.Index_13.A549_ChIP_BRD4_NA_EtOH_0-01_1h_Rep2_R1
JOB_DEPENDENCIES=$picard_sam_to_fastq_2_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.4552.003.Index_13.A549_ChIP_BRD4_NA_EtOH_0-01_1h_Rep2_R1.f19c6d7fb91ba4c30f8aa79ec2f003ae.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.4552.003.Index_13.A549_ChIP_BRD4_NA_EtOH_0-01_1h_Rep2_R1.f19c6d7fb91ba4c30f8aa79ec2f003ae.mugqic.done'
module load java/1.8.0_121 mugqic/trimmomatic/0.36 && \
mkdir -p trim/GM12878_NIPBL_rep1 && \
`cat > trim/GM12878_NIPBL_rep1/HI.4552.003.Index_13.A549_ChIP_BRD4_NA_EtOH_0-01_1h_Rep2_R1.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /project/6001942/Working_Directory/Eric/CofactorHR/GM12878/raw/HI.1385.006.Index_2.GM12878_NIPBL_IB_ChIP3_1_rep1.single.fastq.gz \
  trim/GM12878_NIPBL_rep1/HI.4552.003.Index_13.A549_ChIP_BRD4_NA_EtOH_0-01_1h_Rep2_R1.trim.single.fastq.gz \
  ILLUMINACLIP:trim/GM12878_NIPBL_rep1/HI.4552.003.Index_13.A549_ChIP_BRD4_NA_EtOH_0-01_1h_Rep2_R1.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/GM12878_NIPBL_rep1/HI.4552.003.Index_13.A549_ChIP_BRD4_NA_EtOH_0-01_1h_Rep2_R1.trim.log
trimmomatic.HI.4552.003.Index_13.A549_ChIP_BRD4_NA_EtOH_0-01_1h_Rep2_R1.f19c6d7fb91ba4c30f8aa79ec2f003ae.mugqic.done
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
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"trimmomatic\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_NIPBL_rep1.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"trimmomatic\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_NIPBL_rep1.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=24G -N 1 -n 6 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$trimmomatic_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: trimmomatic_3_JOB_ID: trimmomatic.HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_3_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.54a8c5d24351cedb959b11baba9662f8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.54a8c5d24351cedb959b11baba9662f8.mugqic.done'
module load java/1.8.0_121 mugqic/trimmomatic/0.36 && \
mkdir -p trim/GM12878_MED1_rep1 && \
`cat > trim/GM12878_MED1_rep1/HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /project/6001942/Working_Directory/Eric/CofactorHR/GM12878/raw/HI.1385.006.Index_3.GM12878_MED1_IB_ChIP3_2_rep1.single.fastq.gz \
  trim/GM12878_MED1_rep1/HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/GM12878_MED1_rep1/HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/GM12878_MED1_rep1/HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.trim.log
trimmomatic.HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.54a8c5d24351cedb959b11baba9662f8.mugqic.done
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
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"trimmomatic\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_MED1_rep1.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"trimmomatic\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_MED1_rep1.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=24G -N 1 -n 6 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$trimmomatic_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: trimmomatic_4_JOB_ID: trimmomatic.HI.4552.003.Index_6.A549_ChIP_CDK9_NA_EtOH_0-01_1h_Rep2_R1
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.4552.003.Index_6.A549_ChIP_CDK9_NA_EtOH_0-01_1h_Rep2_R1
JOB_DEPENDENCIES=$picard_sam_to_fastq_4_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.4552.003.Index_6.A549_ChIP_CDK9_NA_EtOH_0-01_1h_Rep2_R1.260612c9c300cfdcfeb2b3e7edabe3b6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.4552.003.Index_6.A549_ChIP_CDK9_NA_EtOH_0-01_1h_Rep2_R1.260612c9c300cfdcfeb2b3e7edabe3b6.mugqic.done'
module load java/1.8.0_121 mugqic/trimmomatic/0.36 && \
mkdir -p trim/GM12878_SMC1_rep1 && \
`cat > trim/GM12878_SMC1_rep1/HI.4552.003.Index_6.A549_ChIP_CDK9_NA_EtOH_0-01_1h_Rep2_R1.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /project/6001942/Working_Directory/Eric/CofactorHR/GM12878/raw/HI.1385.006.Index_4.GM12878_SMC1_IB_ChIP3_3_rep1.single.fastq.gz \
  trim/GM12878_SMC1_rep1/HI.4552.003.Index_6.A549_ChIP_CDK9_NA_EtOH_0-01_1h_Rep2_R1.trim.single.fastq.gz \
  ILLUMINACLIP:trim/GM12878_SMC1_rep1/HI.4552.003.Index_6.A549_ChIP_CDK9_NA_EtOH_0-01_1h_Rep2_R1.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/GM12878_SMC1_rep1/HI.4552.003.Index_6.A549_ChIP_CDK9_NA_EtOH_0-01_1h_Rep2_R1.trim.log
trimmomatic.HI.4552.003.Index_6.A549_ChIP_CDK9_NA_EtOH_0-01_1h_Rep2_R1.260612c9c300cfdcfeb2b3e7edabe3b6.mugqic.done
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
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"trimmomatic\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_SMC1_rep1.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"trimmomatic\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_SMC1_rep1.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=24G -N 1 -n 6 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$trimmomatic_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: trimmomatic_5_JOB_ID: trimmomatic.HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_5_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.aafee298717dcad196b4a8ee29a8f4d9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.aafee298717dcad196b4a8ee29a8f4d9.mugqic.done'
module load java/1.8.0_121 mugqic/trimmomatic/0.36 && \
mkdir -p trim/GM12878_CDK9_rep1 && \
`cat > trim/GM12878_CDK9_rep1/HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /project/6001942/Working_Directory/Eric/CofactorHR/GM12878/raw/HI.1385.006.Index_5.GM12878_CDK9_IB_ChIP3_4_rep1.single.fastq.gz \
  trim/GM12878_CDK9_rep1/HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/GM12878_CDK9_rep1/HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/GM12878_CDK9_rep1/HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.trim.log
trimmomatic.HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.aafee298717dcad196b4a8ee29a8f4d9.mugqic.done
)
trimmomatic_5_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"trimmomatic\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_CDK9_rep1.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"trimmomatic\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_CDK9_rep1.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=24G -N 1 -n 6 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$trimmomatic_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: trimmomatic_6_JOB_ID: trimmomatic.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_6_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.b76d533c0040942e6a76d38323d2aeda.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.b76d533c0040942e6a76d38323d2aeda.mugqic.done'
module load java/1.8.0_121 mugqic/trimmomatic/0.36 && \
mkdir -p trim/GM12878_WCE_rep2 && \
`cat > trim/GM12878_WCE_rep2/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /project/6001942/Working_Directory/Eric/CofactorHR/GM12878/raw/HI.1997.006.Index_14.GM12878_WCE_IB_CHIP7_rep2.single.fastq.gz \
  trim/GM12878_WCE_rep2/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/GM12878_WCE_rep2/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/GM12878_WCE_rep2/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.trim.log
trimmomatic.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.b76d533c0040942e6a76d38323d2aeda.mugqic.done
)
trimmomatic_6_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"trimmomatic\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_WCE_rep2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"trimmomatic\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_WCE_rep2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=24G -N 1 -n 6 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$trimmomatic_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: trimmomatic_7_JOB_ID: trimmomatic.HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_7_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam.9b9d56c84b5183184ad8cc0def562d6c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam.9b9d56c84b5183184ad8cc0def562d6c.mugqic.done'
module load java/1.8.0_121 mugqic/trimmomatic/0.36 && \
mkdir -p trim/GM12878_SMC1_rep2 && \
`cat > trim/GM12878_SMC1_rep2/HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /project/6001942/Working_Directory/Eric/CofactorHR/GM12878/raw/HI.1997.006.Index_15.GM12878_SMC1_IB_CHIP7_rep2.single.fastq.gz \
  trim/GM12878_SMC1_rep2/HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/GM12878_SMC1_rep2/HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/GM12878_SMC1_rep2/HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam.trim.log
trimmomatic.HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam.9b9d56c84b5183184ad8cc0def562d6c.mugqic.done
)
trimmomatic_7_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"trimmomatic\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_SMC1_rep2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"trimmomatic\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_SMC1_rep2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=24G -N 1 -n 6 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$trimmomatic_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: trimmomatic_8_JOB_ID: trimmomatic.HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_8_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.766f7fad9578c1ab1ecc0bbb5fe48166.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.766f7fad9578c1ab1ecc0bbb5fe48166.mugqic.done'
module load java/1.8.0_121 mugqic/trimmomatic/0.36 && \
mkdir -p trim/GM12878_CDK9_rep2 && \
`cat > trim/GM12878_CDK9_rep2/HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /project/6001942/Working_Directory/Eric/CofactorHR/GM12878/raw/HI.1997.006.Index_18.GM12878_CDK9_IB_CHIP7_rep2.single.fastq.gz \
  trim/GM12878_CDK9_rep2/HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/GM12878_CDK9_rep2/HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/GM12878_CDK9_rep2/HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.trim.log
trimmomatic.HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.766f7fad9578c1ab1ecc0bbb5fe48166.mugqic.done
)
trimmomatic_8_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"trimmomatic\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_CDK9_rep2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"trimmomatic\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_CDK9_rep2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=24G -N 1 -n 6 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$trimmomatic_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: trimmomatic_9_JOB_ID: trimmomatic.HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_9_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.1ef2f68036e25b080bdaebb811ab725e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.1ef2f68036e25b080bdaebb811ab725e.mugqic.done'
module load java/1.8.0_121 mugqic/trimmomatic/0.36 && \
mkdir -p trim/GM12878_NIPBL_rep2 && \
`cat > trim/GM12878_NIPBL_rep2/HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /project/6001942/Working_Directory/Eric/CofactorHR/GM12878/raw/HI.1997.006.Index_6.GM12878_NIPBL_IB_CHIP7_rep2.single.fastq.gz \
  trim/GM12878_NIPBL_rep2/HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/GM12878_NIPBL_rep2/HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/GM12878_NIPBL_rep2/HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.trim.log
trimmomatic.HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.1ef2f68036e25b080bdaebb811ab725e.mugqic.done
)
trimmomatic_9_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"trimmomatic\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_NIPBL_rep2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"trimmomatic\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_NIPBL_rep2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=24G -N 1 -n 6 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$trimmomatic_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: trimmomatic_10_JOB_ID: trimmomatic.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_10_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.0703bbd2217959e52ba5f2c6c45c9d2c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.0703bbd2217959e52ba5f2c6c45c9d2c.mugqic.done'
module load java/1.8.0_121 mugqic/trimmomatic/0.36 && \
mkdir -p trim/GM12878_MED1_rep2 && \
`cat > trim/GM12878_MED1_rep2/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /project/6001942/Working_Directory/Eric/CofactorHR/GM12878/raw/HI.1997.006.Index_7.GM12878_MED1_IB_CHIP7_rep2.single.fastq.gz \
  trim/GM12878_MED1_rep2/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/GM12878_MED1_rep2/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/GM12878_MED1_rep2/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.trim.log
trimmomatic.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.0703bbd2217959e52ba5f2c6c45c9d2c.mugqic.done
)
trimmomatic_10_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"trimmomatic\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_MED1_rep2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"trimmomatic\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_MED1_rep2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=24G -N 1 -n 6 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$trimmomatic_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

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
JOB_DEPENDENCIES=$trimmomatic_1_JOB_ID:$trimmomatic_2_JOB_ID:$trimmomatic_3_JOB_ID:$trimmomatic_4_JOB_ID:$trimmomatic_5_JOB_ID:$trimmomatic_6_JOB_ID:$trimmomatic_7_JOB_ID:$trimmomatic_8_JOB_ID:$trimmomatic_9_JOB_ID:$trimmomatic_10_JOB_ID
JOB_DONE=job_output/merge_trimmomatic_stats/merge_trimmomatic_stats.1dc73ae5d192a6e676b52311cd848031.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'merge_trimmomatic_stats.1dc73ae5d192a6e676b52311cd848031.mugqic.done'
module load mugqic/pandoc/1.15.2 && \
mkdir -p metrics && \
echo 'Sample	Readset	Raw Single Reads #	Surviving Single Reads #	Surviving Single Reads %' > metrics/trimReadsetTable.tsv && \
grep ^Input trim/GM12878_WCE_rep1/HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/GM12878_WCE_rep1	HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/GM12878_NIPBL_rep1/HI.4552.003.Index_13.A549_ChIP_BRD4_NA_EtOH_0-01_1h_Rep2_R1.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/GM12878_NIPBL_rep1	HI.4552.003.Index_13.A549_ChIP_BRD4_NA_EtOH_0-01_1h_Rep2_R1	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/GM12878_MED1_rep1/HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/GM12878_MED1_rep1	HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/GM12878_SMC1_rep1/HI.4552.003.Index_6.A549_ChIP_CDK9_NA_EtOH_0-01_1h_Rep2_R1.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/GM12878_SMC1_rep1	HI.4552.003.Index_6.A549_ChIP_CDK9_NA_EtOH_0-01_1h_Rep2_R1	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/GM12878_CDK9_rep1/HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/GM12878_CDK9_rep1	HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/GM12878_WCE_rep2/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/GM12878_WCE_rep2	HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/GM12878_SMC1_rep2/HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/GM12878_SMC1_rep2	HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/GM12878_CDK9_rep2/HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/GM12878_CDK9_rep2	HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/GM12878_NIPBL_rep2/HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/GM12878_NIPBL_rep2	HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/GM12878_MED1_rep2/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/GM12878_MED1_rep2	HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam	\1	\2/' | \
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
  --variable min_length=50 \
  --variable read_type=Single \
  --variable trim_readset_table="$trim_readset_table_md" \
  --to markdown \
  > report/Illumina.merge_trimmomatic_stats.md
merge_trimmomatic_stats.1dc73ae5d192a6e676b52311cd848031.mugqic.done
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
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"merge_trimmomatic_stats\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_WCE_rep1.json,/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_NIPBL_rep1.json,/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_MED1_rep1.json,/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_SMC1_rep1.json,/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_CDK9_rep1.json,/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_WCE_rep2.json,/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_SMC1_rep2.json,/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_CDK9_rep2.json,/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_NIPBL_rep2.json,/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_MED1_rep2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"merge_trimmomatic_stats\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_WCE_rep1.json,/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_NIPBL_rep1.json,/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_MED1_rep1.json,/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_SMC1_rep1.json,/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_CDK9_rep1.json,/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_WCE_rep2.json,/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_SMC1_rep2.json,/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_CDK9_rep2.json,/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_NIPBL_rep2.json,/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_MED1_rep2.json\" \
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
# STEP: bwa_mem_picard_sort_sam
#-------------------------------------------------------------------------------
STEP=bwa_mem_picard_sort_sam
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_1_JOB_ID: bwa_mem_picard_sort_sam.HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_1_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.527724cfc6d3617aa335d921e347d3a4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.527724cfc6d3617aa335d921e347d3a4.mugqic.done'
module load mugqic/bwa/0.7.12 java/1.8.0_121 mugqic/picard/2.0.1 && \
mkdir -p alignment/GM12878_WCE_rep1/HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam	SM:GM12878_WCE_rep1	LB:GM12878_WCE_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/GM12878_WCE_rep1/HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/GM12878_WCE_rep1/HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam/HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.527724cfc6d3617aa335d921e347d3a4.mugqic.done
)
bwa_mem_picard_sort_sam_1_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"bwa_mem_picard_sort_sam\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_WCE_rep1.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"bwa_mem_picard_sort_sam\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_WCE_rep1.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=48G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$bwa_mem_picard_sort_sam_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_2_JOB_ID: bwa_mem_picard_sort_sam.HI.4552.003.Index_13.A549_ChIP_BRD4_NA_EtOH_0-01_1h_Rep2_R1
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.4552.003.Index_13.A549_ChIP_BRD4_NA_EtOH_0-01_1h_Rep2_R1
JOB_DEPENDENCIES=$trimmomatic_2_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.4552.003.Index_13.A549_ChIP_BRD4_NA_EtOH_0-01_1h_Rep2_R1.1ce3cf1cf880ad4d6f9e8abb66a30732.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.4552.003.Index_13.A549_ChIP_BRD4_NA_EtOH_0-01_1h_Rep2_R1.1ce3cf1cf880ad4d6f9e8abb66a30732.mugqic.done'
module load mugqic/bwa/0.7.12 java/1.8.0_121 mugqic/picard/2.0.1 && \
mkdir -p alignment/GM12878_NIPBL_rep1/HI.4552.003.Index_13.A549_ChIP_BRD4_NA_EtOH_0-01_1h_Rep2_R1 && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:HI.4552.003.Index_13.A549_ChIP_BRD4_NA_EtOH_0-01_1h_Rep2_R1	SM:GM12878_NIPBL_rep1	LB:GM12878_NIPBL_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/GM12878_NIPBL_rep1/HI.4552.003.Index_13.A549_ChIP_BRD4_NA_EtOH_0-01_1h_Rep2_R1.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/GM12878_NIPBL_rep1/HI.4552.003.Index_13.A549_ChIP_BRD4_NA_EtOH_0-01_1h_Rep2_R1/HI.4552.003.Index_13.A549_ChIP_BRD4_NA_EtOH_0-01_1h_Rep2_R1.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.HI.4552.003.Index_13.A549_ChIP_BRD4_NA_EtOH_0-01_1h_Rep2_R1.1ce3cf1cf880ad4d6f9e8abb66a30732.mugqic.done
)
bwa_mem_picard_sort_sam_2_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"bwa_mem_picard_sort_sam\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_NIPBL_rep1.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"bwa_mem_picard_sort_sam\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_NIPBL_rep1.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=48G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$bwa_mem_picard_sort_sam_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_3_JOB_ID: bwa_mem_picard_sort_sam.HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_3_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.18798aab2ee203fc54e247ca3877690e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.18798aab2ee203fc54e247ca3877690e.mugqic.done'
module load mugqic/bwa/0.7.12 java/1.8.0_121 mugqic/picard/2.0.1 && \
mkdir -p alignment/GM12878_MED1_rep1/HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam	SM:GM12878_MED1_rep1	LB:GM12878_MED1_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/GM12878_MED1_rep1/HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/GM12878_MED1_rep1/HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam/HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.18798aab2ee203fc54e247ca3877690e.mugqic.done
)
bwa_mem_picard_sort_sam_3_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"bwa_mem_picard_sort_sam\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_MED1_rep1.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"bwa_mem_picard_sort_sam\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_MED1_rep1.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=48G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$bwa_mem_picard_sort_sam_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_4_JOB_ID: bwa_mem_picard_sort_sam.HI.4552.003.Index_6.A549_ChIP_CDK9_NA_EtOH_0-01_1h_Rep2_R1
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.4552.003.Index_6.A549_ChIP_CDK9_NA_EtOH_0-01_1h_Rep2_R1
JOB_DEPENDENCIES=$trimmomatic_4_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.4552.003.Index_6.A549_ChIP_CDK9_NA_EtOH_0-01_1h_Rep2_R1.8d96d947a01a26666e7042f041eb6be6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.4552.003.Index_6.A549_ChIP_CDK9_NA_EtOH_0-01_1h_Rep2_R1.8d96d947a01a26666e7042f041eb6be6.mugqic.done'
module load mugqic/bwa/0.7.12 java/1.8.0_121 mugqic/picard/2.0.1 && \
mkdir -p alignment/GM12878_SMC1_rep1/HI.4552.003.Index_6.A549_ChIP_CDK9_NA_EtOH_0-01_1h_Rep2_R1 && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:HI.4552.003.Index_6.A549_ChIP_CDK9_NA_EtOH_0-01_1h_Rep2_R1	SM:GM12878_SMC1_rep1	LB:GM12878_SMC1_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/GM12878_SMC1_rep1/HI.4552.003.Index_6.A549_ChIP_CDK9_NA_EtOH_0-01_1h_Rep2_R1.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/GM12878_SMC1_rep1/HI.4552.003.Index_6.A549_ChIP_CDK9_NA_EtOH_0-01_1h_Rep2_R1/HI.4552.003.Index_6.A549_ChIP_CDK9_NA_EtOH_0-01_1h_Rep2_R1.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.HI.4552.003.Index_6.A549_ChIP_CDK9_NA_EtOH_0-01_1h_Rep2_R1.8d96d947a01a26666e7042f041eb6be6.mugqic.done
)
bwa_mem_picard_sort_sam_4_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"bwa_mem_picard_sort_sam\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_SMC1_rep1.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"bwa_mem_picard_sort_sam\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_SMC1_rep1.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=48G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$bwa_mem_picard_sort_sam_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_5_JOB_ID: bwa_mem_picard_sort_sam.HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_5_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.acf083f1a900e6e7e10d23594084597a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.acf083f1a900e6e7e10d23594084597a.mugqic.done'
module load mugqic/bwa/0.7.12 java/1.8.0_121 mugqic/picard/2.0.1 && \
mkdir -p alignment/GM12878_CDK9_rep1/HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam	SM:GM12878_CDK9_rep1	LB:GM12878_CDK9_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/GM12878_CDK9_rep1/HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/GM12878_CDK9_rep1/HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam/HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.acf083f1a900e6e7e10d23594084597a.mugqic.done
)
bwa_mem_picard_sort_sam_5_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"bwa_mem_picard_sort_sam\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_CDK9_rep1.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"bwa_mem_picard_sort_sam\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_CDK9_rep1.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=48G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$bwa_mem_picard_sort_sam_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_6_JOB_ID: bwa_mem_picard_sort_sam.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_6_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.c3d08ad4e32cc113144f8c179350e7ba.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.c3d08ad4e32cc113144f8c179350e7ba.mugqic.done'
module load mugqic/bwa/0.7.12 java/1.8.0_121 mugqic/picard/2.0.1 && \
mkdir -p alignment/GM12878_WCE_rep2/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam	SM:GM12878_WCE_rep2	LB:GM12878_WCE_rep2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/GM12878_WCE_rep2/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/GM12878_WCE_rep2/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.c3d08ad4e32cc113144f8c179350e7ba.mugqic.done
)
bwa_mem_picard_sort_sam_6_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"bwa_mem_picard_sort_sam\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_WCE_rep2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"bwa_mem_picard_sort_sam\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_WCE_rep2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=48G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$bwa_mem_picard_sort_sam_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_7_JOB_ID: bwa_mem_picard_sort_sam.HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_7_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam.973dd012f1e6af4887314d2662cb37ce.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam.973dd012f1e6af4887314d2662cb37ce.mugqic.done'
module load mugqic/bwa/0.7.12 java/1.8.0_121 mugqic/picard/2.0.1 && \
mkdir -p alignment/GM12878_SMC1_rep2/HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam	SM:GM12878_SMC1_rep2	LB:GM12878_SMC1_rep2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/GM12878_SMC1_rep2/HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/GM12878_SMC1_rep2/HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam/HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam.973dd012f1e6af4887314d2662cb37ce.mugqic.done
)
bwa_mem_picard_sort_sam_7_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"bwa_mem_picard_sort_sam\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_SMC1_rep2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"bwa_mem_picard_sort_sam\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_SMC1_rep2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=48G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$bwa_mem_picard_sort_sam_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_8_JOB_ID: bwa_mem_picard_sort_sam.HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_8_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.deb39ad775df6631c8f39d0c58cbed95.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.deb39ad775df6631c8f39d0c58cbed95.mugqic.done'
module load mugqic/bwa/0.7.12 java/1.8.0_121 mugqic/picard/2.0.1 && \
mkdir -p alignment/GM12878_CDK9_rep2/HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam	SM:GM12878_CDK9_rep2	LB:GM12878_CDK9_rep2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/GM12878_CDK9_rep2/HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/GM12878_CDK9_rep2/HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam/HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.deb39ad775df6631c8f39d0c58cbed95.mugqic.done
)
bwa_mem_picard_sort_sam_8_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"bwa_mem_picard_sort_sam\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_CDK9_rep2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"bwa_mem_picard_sort_sam\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_CDK9_rep2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=48G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$bwa_mem_picard_sort_sam_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_9_JOB_ID: bwa_mem_picard_sort_sam.HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_9_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.94d76c4bb595cd5c533c78617f053f37.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.94d76c4bb595cd5c533c78617f053f37.mugqic.done'
module load mugqic/bwa/0.7.12 java/1.8.0_121 mugqic/picard/2.0.1 && \
mkdir -p alignment/GM12878_NIPBL_rep2/HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam	SM:GM12878_NIPBL_rep2	LB:GM12878_NIPBL_rep2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/GM12878_NIPBL_rep2/HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/GM12878_NIPBL_rep2/HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam/HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.94d76c4bb595cd5c533c78617f053f37.mugqic.done
)
bwa_mem_picard_sort_sam_9_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"bwa_mem_picard_sort_sam\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_NIPBL_rep2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"bwa_mem_picard_sort_sam\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_NIPBL_rep2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=48G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$bwa_mem_picard_sort_sam_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_10_JOB_ID: bwa_mem_picard_sort_sam.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_10_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.dba898802218e953942c9f7745d5a7ed.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.dba898802218e953942c9f7745d5a7ed.mugqic.done'
module load mugqic/bwa/0.7.12 java/1.8.0_121 mugqic/picard/2.0.1 && \
mkdir -p alignment/GM12878_MED1_rep2/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam	SM:GM12878_MED1_rep2	LB:GM12878_MED1_rep2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/GM12878_MED1_rep2/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/GM12878_MED1_rep2/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.dba898802218e953942c9f7745d5a7ed.mugqic.done
)
bwa_mem_picard_sort_sam_10_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"bwa_mem_picard_sort_sam\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_MED1_rep2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"bwa_mem_picard_sort_sam\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_MED1_rep2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=48G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$bwa_mem_picard_sort_sam_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_11_JOB_ID: bwa_mem_picard_sort_sam_report
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam_report
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_1_JOB_ID:$bwa_mem_picard_sort_sam_2_JOB_ID:$bwa_mem_picard_sort_sam_3_JOB_ID:$bwa_mem_picard_sort_sam_4_JOB_ID:$bwa_mem_picard_sort_sam_5_JOB_ID:$bwa_mem_picard_sort_sam_6_JOB_ID:$bwa_mem_picard_sort_sam_7_JOB_ID:$bwa_mem_picard_sort_sam_8_JOB_ID:$bwa_mem_picard_sort_sam_9_JOB_ID:$bwa_mem_picard_sort_sam_10_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam_report.22d7ff6198a068bd5763589f2515f7f6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam_report.22d7ff6198a068bd5763589f2515f7f6.mugqic.done'
module load mugqic/pandoc/1.15.2 && \
mkdir -p report && \
pandoc --to=markdown \
  --template /home/efournie/genpipes/bfx/report/DnaSeq.bwa_mem_picard_sort_sam.md \
  --variable scientific_name="Homo_sapiens" \
  --variable assembly="GRCh38" \
  /home/efournie/genpipes/bfx/report/DnaSeq.bwa_mem_picard_sort_sam.md \
  > report/DnaSeq.bwa_mem_picard_sort_sam.md
bwa_mem_picard_sort_sam_report.22d7ff6198a068bd5763589f2515f7f6.mugqic.done
)
bwa_mem_picard_sort_sam_11_JOB_ID=$(echo "#! /bin/bash 
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
echo "$bwa_mem_picard_sort_sam_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# STEP: samtools_view_filter
#-------------------------------------------------------------------------------
STEP=samtools_view_filter
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_1_JOB_ID: samtools_view_filter.HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_1_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.d5483523510da5dddddb5e6591d73136.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.d5483523510da5dddddb5e6591d73136.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/GM12878_WCE_rep1/HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam/HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.sorted.bam \
  > alignment/GM12878_WCE_rep1/HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam/HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.d5483523510da5dddddb5e6591d73136.mugqic.done
)
samtools_view_filter_1_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"samtools_view_filter\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_WCE_rep1.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"samtools_view_filter\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_WCE_rep1.json\" \
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
echo "$samtools_view_filter_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_2_JOB_ID: samtools_view_filter.HI.4552.003.Index_13.A549_ChIP_BRD4_NA_EtOH_0-01_1h_Rep2_R1
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.4552.003.Index_13.A549_ChIP_BRD4_NA_EtOH_0-01_1h_Rep2_R1
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_2_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.4552.003.Index_13.A549_ChIP_BRD4_NA_EtOH_0-01_1h_Rep2_R1.7923dd95f4c49bed807872b6071e2913.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.4552.003.Index_13.A549_ChIP_BRD4_NA_EtOH_0-01_1h_Rep2_R1.7923dd95f4c49bed807872b6071e2913.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/GM12878_NIPBL_rep1/HI.4552.003.Index_13.A549_ChIP_BRD4_NA_EtOH_0-01_1h_Rep2_R1/HI.4552.003.Index_13.A549_ChIP_BRD4_NA_EtOH_0-01_1h_Rep2_R1.sorted.bam \
  > alignment/GM12878_NIPBL_rep1/HI.4552.003.Index_13.A549_ChIP_BRD4_NA_EtOH_0-01_1h_Rep2_R1/HI.4552.003.Index_13.A549_ChIP_BRD4_NA_EtOH_0-01_1h_Rep2_R1.sorted.filtered.bam
samtools_view_filter.HI.4552.003.Index_13.A549_ChIP_BRD4_NA_EtOH_0-01_1h_Rep2_R1.7923dd95f4c49bed807872b6071e2913.mugqic.done
)
samtools_view_filter_2_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"samtools_view_filter\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_NIPBL_rep1.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"samtools_view_filter\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_NIPBL_rep1.json\" \
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
echo "$samtools_view_filter_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_3_JOB_ID: samtools_view_filter.HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_3_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.7a6f7e6a197dad8dbe2adf9685a8736b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.7a6f7e6a197dad8dbe2adf9685a8736b.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/GM12878_MED1_rep1/HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam/HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.sorted.bam \
  > alignment/GM12878_MED1_rep1/HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam/HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.7a6f7e6a197dad8dbe2adf9685a8736b.mugqic.done
)
samtools_view_filter_3_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"samtools_view_filter\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_MED1_rep1.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"samtools_view_filter\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_MED1_rep1.json\" \
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
echo "$samtools_view_filter_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_4_JOB_ID: samtools_view_filter.HI.4552.003.Index_6.A549_ChIP_CDK9_NA_EtOH_0-01_1h_Rep2_R1
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.4552.003.Index_6.A549_ChIP_CDK9_NA_EtOH_0-01_1h_Rep2_R1
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_4_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.4552.003.Index_6.A549_ChIP_CDK9_NA_EtOH_0-01_1h_Rep2_R1.1d27525b932b824e3856a1efad0f653c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.4552.003.Index_6.A549_ChIP_CDK9_NA_EtOH_0-01_1h_Rep2_R1.1d27525b932b824e3856a1efad0f653c.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/GM12878_SMC1_rep1/HI.4552.003.Index_6.A549_ChIP_CDK9_NA_EtOH_0-01_1h_Rep2_R1/HI.4552.003.Index_6.A549_ChIP_CDK9_NA_EtOH_0-01_1h_Rep2_R1.sorted.bam \
  > alignment/GM12878_SMC1_rep1/HI.4552.003.Index_6.A549_ChIP_CDK9_NA_EtOH_0-01_1h_Rep2_R1/HI.4552.003.Index_6.A549_ChIP_CDK9_NA_EtOH_0-01_1h_Rep2_R1.sorted.filtered.bam
samtools_view_filter.HI.4552.003.Index_6.A549_ChIP_CDK9_NA_EtOH_0-01_1h_Rep2_R1.1d27525b932b824e3856a1efad0f653c.mugqic.done
)
samtools_view_filter_4_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"samtools_view_filter\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_SMC1_rep1.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"samtools_view_filter\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_SMC1_rep1.json\" \
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
echo "$samtools_view_filter_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_5_JOB_ID: samtools_view_filter.HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_5_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.1bdc3a2f9d2ecb27368f2e23392db08a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.1bdc3a2f9d2ecb27368f2e23392db08a.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/GM12878_CDK9_rep1/HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam/HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.sorted.bam \
  > alignment/GM12878_CDK9_rep1/HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam/HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.1bdc3a2f9d2ecb27368f2e23392db08a.mugqic.done
)
samtools_view_filter_5_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"samtools_view_filter\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_CDK9_rep1.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"samtools_view_filter\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_CDK9_rep1.json\" \
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
echo "$samtools_view_filter_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_6_JOB_ID: samtools_view_filter.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_6_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.9c8daae1c775a6807ede6e39b385d6c2.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.9c8daae1c775a6807ede6e39b385d6c2.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/GM12878_WCE_rep2/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.sorted.bam \
  > alignment/GM12878_WCE_rep2/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.9c8daae1c775a6807ede6e39b385d6c2.mugqic.done
)
samtools_view_filter_6_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"samtools_view_filter\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_WCE_rep2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"samtools_view_filter\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_WCE_rep2.json\" \
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
echo "$samtools_view_filter_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_7_JOB_ID: samtools_view_filter.HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_7_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam.1a1a86fd4a11997a55357a449608ce1b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam.1a1a86fd4a11997a55357a449608ce1b.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/GM12878_SMC1_rep2/HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam/HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam.sorted.bam \
  > alignment/GM12878_SMC1_rep2/HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam/HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam.1a1a86fd4a11997a55357a449608ce1b.mugqic.done
)
samtools_view_filter_7_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"samtools_view_filter\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_SMC1_rep2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"samtools_view_filter\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_SMC1_rep2.json\" \
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
echo "$samtools_view_filter_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_8_JOB_ID: samtools_view_filter.HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_8_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.bda21031428e3ce34cb2e782b26794c5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.bda21031428e3ce34cb2e782b26794c5.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/GM12878_CDK9_rep2/HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam/HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.sorted.bam \
  > alignment/GM12878_CDK9_rep2/HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam/HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.bda21031428e3ce34cb2e782b26794c5.mugqic.done
)
samtools_view_filter_8_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"samtools_view_filter\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_CDK9_rep2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"samtools_view_filter\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_CDK9_rep2.json\" \
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
echo "$samtools_view_filter_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_9_JOB_ID: samtools_view_filter.HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_9_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.fa9158dcda6197e38e71731dba8fbcec.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.fa9158dcda6197e38e71731dba8fbcec.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/GM12878_NIPBL_rep2/HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam/HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.sorted.bam \
  > alignment/GM12878_NIPBL_rep2/HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam/HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.fa9158dcda6197e38e71731dba8fbcec.mugqic.done
)
samtools_view_filter_9_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"samtools_view_filter\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_NIPBL_rep2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"samtools_view_filter\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_NIPBL_rep2.json\" \
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
echo "$samtools_view_filter_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_10_JOB_ID: samtools_view_filter.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_10_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.abdde763c869aec4c164ae64aa07658a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.abdde763c869aec4c164ae64aa07658a.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/GM12878_MED1_rep2/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.sorted.bam \
  > alignment/GM12878_MED1_rep2/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.abdde763c869aec4c164ae64aa07658a.mugqic.done
)
samtools_view_filter_10_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"samtools_view_filter\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_MED1_rep2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"samtools_view_filter\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_MED1_rep2.json\" \
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
echo "$samtools_view_filter_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_11_JOB_ID: samtools_view_filter_report
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter_report
JOB_DEPENDENCIES=$samtools_view_filter_1_JOB_ID:$samtools_view_filter_2_JOB_ID:$samtools_view_filter_3_JOB_ID:$samtools_view_filter_4_JOB_ID:$samtools_view_filter_5_JOB_ID:$samtools_view_filter_6_JOB_ID:$samtools_view_filter_7_JOB_ID:$samtools_view_filter_8_JOB_ID:$samtools_view_filter_9_JOB_ID:$samtools_view_filter_10_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter_report.adda516c06f797351611d7331668a33d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter_report.adda516c06f797351611d7331668a33d.mugqic.done'
module load mugqic/pandoc/1.15.2 && \
mkdir -p report && \
pandoc --to=markdown \
  --template /home/efournie/genpipes/bfx/report/ChipSeq.samtools_view_filter.md \
  --variable min_mapq="20" \
  /home/efournie/genpipes/bfx/report/ChipSeq.samtools_view_filter.md \
  > report/ChipSeq.samtools_view_filter.md
samtools_view_filter_report.adda516c06f797351611d7331668a33d.mugqic.done
)
samtools_view_filter_11_JOB_ID=$(echo "#! /bin/bash 
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
echo "$samtools_view_filter_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# STEP: picard_merge_sam_files
#-------------------------------------------------------------------------------
STEP=picard_merge_sam_files
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_1_JOB_ID: symlink_readset_sample_bam.GM12878_WCE_rep1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.GM12878_WCE_rep1
JOB_DEPENDENCIES=$samtools_view_filter_1_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.GM12878_WCE_rep1.fa0a6511844227171388d4913be8a0c5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.GM12878_WCE_rep1.fa0a6511844227171388d4913be8a0c5.mugqic.done'
mkdir -p alignment/GM12878_WCE_rep1 && \
ln -s -f HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam/HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.sorted.filtered.bam alignment/GM12878_WCE_rep1/GM12878_WCE_rep1.merged.bam
symlink_readset_sample_bam.GM12878_WCE_rep1.fa0a6511844227171388d4913be8a0c5.mugqic.done
)
picard_merge_sam_files_1_JOB_ID=$(echo "#! /bin/bash 
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
echo "$picard_merge_sam_files_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_2_JOB_ID: symlink_readset_sample_bam.GM12878_NIPBL_rep1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.GM12878_NIPBL_rep1
JOB_DEPENDENCIES=$samtools_view_filter_2_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.GM12878_NIPBL_rep1.5bd76157733893dd4c51c72bd7666944.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.GM12878_NIPBL_rep1.5bd76157733893dd4c51c72bd7666944.mugqic.done'
mkdir -p alignment/GM12878_NIPBL_rep1 && \
ln -s -f HI.4552.003.Index_13.A549_ChIP_BRD4_NA_EtOH_0-01_1h_Rep2_R1/HI.4552.003.Index_13.A549_ChIP_BRD4_NA_EtOH_0-01_1h_Rep2_R1.sorted.filtered.bam alignment/GM12878_NIPBL_rep1/GM12878_NIPBL_rep1.merged.bam
symlink_readset_sample_bam.GM12878_NIPBL_rep1.5bd76157733893dd4c51c72bd7666944.mugqic.done
)
picard_merge_sam_files_2_JOB_ID=$(echo "#! /bin/bash 
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
echo "$picard_merge_sam_files_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_3_JOB_ID: symlink_readset_sample_bam.GM12878_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.GM12878_MED1_rep1
JOB_DEPENDENCIES=$samtools_view_filter_3_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.GM12878_MED1_rep1.2abbbf4bf291a6999b96186984b761ab.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.GM12878_MED1_rep1.2abbbf4bf291a6999b96186984b761ab.mugqic.done'
mkdir -p alignment/GM12878_MED1_rep1 && \
ln -s -f HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam/HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.sorted.filtered.bam alignment/GM12878_MED1_rep1/GM12878_MED1_rep1.merged.bam
symlink_readset_sample_bam.GM12878_MED1_rep1.2abbbf4bf291a6999b96186984b761ab.mugqic.done
)
picard_merge_sam_files_3_JOB_ID=$(echo "#! /bin/bash 
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
echo "$picard_merge_sam_files_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_4_JOB_ID: symlink_readset_sample_bam.GM12878_SMC1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.GM12878_SMC1_rep1
JOB_DEPENDENCIES=$samtools_view_filter_4_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.GM12878_SMC1_rep1.50daba3dd1c33d62bc0e17312f780706.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.GM12878_SMC1_rep1.50daba3dd1c33d62bc0e17312f780706.mugqic.done'
mkdir -p alignment/GM12878_SMC1_rep1 && \
ln -s -f HI.4552.003.Index_6.A549_ChIP_CDK9_NA_EtOH_0-01_1h_Rep2_R1/HI.4552.003.Index_6.A549_ChIP_CDK9_NA_EtOH_0-01_1h_Rep2_R1.sorted.filtered.bam alignment/GM12878_SMC1_rep1/GM12878_SMC1_rep1.merged.bam
symlink_readset_sample_bam.GM12878_SMC1_rep1.50daba3dd1c33d62bc0e17312f780706.mugqic.done
)
picard_merge_sam_files_4_JOB_ID=$(echo "#! /bin/bash 
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
echo "$picard_merge_sam_files_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_5_JOB_ID: symlink_readset_sample_bam.GM12878_CDK9_rep1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.GM12878_CDK9_rep1
JOB_DEPENDENCIES=$samtools_view_filter_5_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.GM12878_CDK9_rep1.5e907b36fe9c96d84519e031417e8330.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.GM12878_CDK9_rep1.5e907b36fe9c96d84519e031417e8330.mugqic.done'
mkdir -p alignment/GM12878_CDK9_rep1 && \
ln -s -f HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam/HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.sorted.filtered.bam alignment/GM12878_CDK9_rep1/GM12878_CDK9_rep1.merged.bam
symlink_readset_sample_bam.GM12878_CDK9_rep1.5e907b36fe9c96d84519e031417e8330.mugqic.done
)
picard_merge_sam_files_5_JOB_ID=$(echo "#! /bin/bash 
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
echo "$picard_merge_sam_files_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_6_JOB_ID: symlink_readset_sample_bam.GM12878_WCE_rep2
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.GM12878_WCE_rep2
JOB_DEPENDENCIES=$samtools_view_filter_6_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.GM12878_WCE_rep2.769d6283a8a9eefb4c036e65dd5b2e28.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.GM12878_WCE_rep2.769d6283a8a9eefb4c036e65dd5b2e28.mugqic.done'
mkdir -p alignment/GM12878_WCE_rep2 && \
ln -s -f HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.sorted.filtered.bam alignment/GM12878_WCE_rep2/GM12878_WCE_rep2.merged.bam
symlink_readset_sample_bam.GM12878_WCE_rep2.769d6283a8a9eefb4c036e65dd5b2e28.mugqic.done
)
picard_merge_sam_files_6_JOB_ID=$(echo "#! /bin/bash 
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
echo "$picard_merge_sam_files_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_7_JOB_ID: symlink_readset_sample_bam.GM12878_SMC1_rep2
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.GM12878_SMC1_rep2
JOB_DEPENDENCIES=$samtools_view_filter_7_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.GM12878_SMC1_rep2.cde532e7e174ad3fa5e77249bc15f47d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.GM12878_SMC1_rep2.cde532e7e174ad3fa5e77249bc15f47d.mugqic.done'
mkdir -p alignment/GM12878_SMC1_rep2 && \
ln -s -f HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam/HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam.sorted.filtered.bam alignment/GM12878_SMC1_rep2/GM12878_SMC1_rep2.merged.bam
symlink_readset_sample_bam.GM12878_SMC1_rep2.cde532e7e174ad3fa5e77249bc15f47d.mugqic.done
)
picard_merge_sam_files_7_JOB_ID=$(echo "#! /bin/bash 
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
echo "$picard_merge_sam_files_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_8_JOB_ID: symlink_readset_sample_bam.GM12878_CDK9_rep2
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.GM12878_CDK9_rep2
JOB_DEPENDENCIES=$samtools_view_filter_8_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.GM12878_CDK9_rep2.cde2ff66a52611e5e92cee5b61357c13.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.GM12878_CDK9_rep2.cde2ff66a52611e5e92cee5b61357c13.mugqic.done'
mkdir -p alignment/GM12878_CDK9_rep2 && \
ln -s -f HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam/HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.sorted.filtered.bam alignment/GM12878_CDK9_rep2/GM12878_CDK9_rep2.merged.bam
symlink_readset_sample_bam.GM12878_CDK9_rep2.cde2ff66a52611e5e92cee5b61357c13.mugqic.done
)
picard_merge_sam_files_8_JOB_ID=$(echo "#! /bin/bash 
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
echo "$picard_merge_sam_files_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_9_JOB_ID: symlink_readset_sample_bam.GM12878_NIPBL_rep2
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.GM12878_NIPBL_rep2
JOB_DEPENDENCIES=$samtools_view_filter_9_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.GM12878_NIPBL_rep2.e55abba566fb2b8215e864fdac8d45f0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.GM12878_NIPBL_rep2.e55abba566fb2b8215e864fdac8d45f0.mugqic.done'
mkdir -p alignment/GM12878_NIPBL_rep2 && \
ln -s -f HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam/HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.sorted.filtered.bam alignment/GM12878_NIPBL_rep2/GM12878_NIPBL_rep2.merged.bam
symlink_readset_sample_bam.GM12878_NIPBL_rep2.e55abba566fb2b8215e864fdac8d45f0.mugqic.done
)
picard_merge_sam_files_9_JOB_ID=$(echo "#! /bin/bash 
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
echo "$picard_merge_sam_files_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_10_JOB_ID: symlink_readset_sample_bam.GM12878_MED1_rep2
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.GM12878_MED1_rep2
JOB_DEPENDENCIES=$samtools_view_filter_10_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.GM12878_MED1_rep2.13a1366bae5324beda8183e44d10145f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.GM12878_MED1_rep2.13a1366bae5324beda8183e44d10145f.mugqic.done'
mkdir -p alignment/GM12878_MED1_rep2 && \
ln -s -f HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.sorted.filtered.bam alignment/GM12878_MED1_rep2/GM12878_MED1_rep2.merged.bam
symlink_readset_sample_bam.GM12878_MED1_rep2.13a1366bae5324beda8183e44d10145f.mugqic.done
)
picard_merge_sam_files_10_JOB_ID=$(echo "#! /bin/bash 
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
echo "$picard_merge_sam_files_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# STEP: picard_mark_duplicates
#-------------------------------------------------------------------------------
STEP=picard_mark_duplicates
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_1_JOB_ID: picard_mark_duplicates.GM12878_WCE_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.GM12878_WCE_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_1_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.GM12878_WCE_rep1.f1d68ceb1a1467ab9019b564ef25ca32.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.GM12878_WCE_rep1.f1d68ceb1a1467ab9019b564ef25ca32.mugqic.done'
module load java/1.8.0_121 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/GM12878_WCE_rep1/GM12878_WCE_rep1.merged.bam \
 OUTPUT=alignment/GM12878_WCE_rep1/GM12878_WCE_rep1.sorted.dup.bam \
 METRICS_FILE=alignment/GM12878_WCE_rep1/GM12878_WCE_rep1.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.GM12878_WCE_rep1.f1d68ceb1a1467ab9019b564ef25ca32.mugqic.done
)
picard_mark_duplicates_1_JOB_ID=$(echo "#! /bin/bash 
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=48:00:0 --mem=8G -N 1 -n 2 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$picard_mark_duplicates_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_2_JOB_ID: picard_mark_duplicates.GM12878_NIPBL_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.GM12878_NIPBL_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_2_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.GM12878_NIPBL_rep1.8fdf49b017beba67d6fc08429a408f74.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.GM12878_NIPBL_rep1.8fdf49b017beba67d6fc08429a408f74.mugqic.done'
module load java/1.8.0_121 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/GM12878_NIPBL_rep1/GM12878_NIPBL_rep1.merged.bam \
 OUTPUT=alignment/GM12878_NIPBL_rep1/GM12878_NIPBL_rep1.sorted.dup.bam \
 METRICS_FILE=alignment/GM12878_NIPBL_rep1/GM12878_NIPBL_rep1.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.GM12878_NIPBL_rep1.8fdf49b017beba67d6fc08429a408f74.mugqic.done
)
picard_mark_duplicates_2_JOB_ID=$(echo "#! /bin/bash 
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=48:00:0 --mem=8G -N 1 -n 2 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$picard_mark_duplicates_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_3_JOB_ID: picard_mark_duplicates.GM12878_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.GM12878_MED1_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_3_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.GM12878_MED1_rep1.4507573b9c068d795a658ead7ee37751.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.GM12878_MED1_rep1.4507573b9c068d795a658ead7ee37751.mugqic.done'
module load java/1.8.0_121 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/GM12878_MED1_rep1/GM12878_MED1_rep1.merged.bam \
 OUTPUT=alignment/GM12878_MED1_rep1/GM12878_MED1_rep1.sorted.dup.bam \
 METRICS_FILE=alignment/GM12878_MED1_rep1/GM12878_MED1_rep1.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.GM12878_MED1_rep1.4507573b9c068d795a658ead7ee37751.mugqic.done
)
picard_mark_duplicates_3_JOB_ID=$(echo "#! /bin/bash 
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=48:00:0 --mem=8G -N 1 -n 2 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$picard_mark_duplicates_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_4_JOB_ID: picard_mark_duplicates.GM12878_SMC1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.GM12878_SMC1_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_4_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.GM12878_SMC1_rep1.b68c8ba7fc2373577bfb6ccdcd4d1c59.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.GM12878_SMC1_rep1.b68c8ba7fc2373577bfb6ccdcd4d1c59.mugqic.done'
module load java/1.8.0_121 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/GM12878_SMC1_rep1/GM12878_SMC1_rep1.merged.bam \
 OUTPUT=alignment/GM12878_SMC1_rep1/GM12878_SMC1_rep1.sorted.dup.bam \
 METRICS_FILE=alignment/GM12878_SMC1_rep1/GM12878_SMC1_rep1.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.GM12878_SMC1_rep1.b68c8ba7fc2373577bfb6ccdcd4d1c59.mugqic.done
)
picard_mark_duplicates_4_JOB_ID=$(echo "#! /bin/bash 
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=48:00:0 --mem=8G -N 1 -n 2 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$picard_mark_duplicates_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_5_JOB_ID: picard_mark_duplicates.GM12878_CDK9_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.GM12878_CDK9_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_5_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.GM12878_CDK9_rep1.fe1ed7208167dc7009d12e59041d0a81.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.GM12878_CDK9_rep1.fe1ed7208167dc7009d12e59041d0a81.mugqic.done'
module load java/1.8.0_121 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/GM12878_CDK9_rep1/GM12878_CDK9_rep1.merged.bam \
 OUTPUT=alignment/GM12878_CDK9_rep1/GM12878_CDK9_rep1.sorted.dup.bam \
 METRICS_FILE=alignment/GM12878_CDK9_rep1/GM12878_CDK9_rep1.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.GM12878_CDK9_rep1.fe1ed7208167dc7009d12e59041d0a81.mugqic.done
)
picard_mark_duplicates_5_JOB_ID=$(echo "#! /bin/bash 
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=48:00:0 --mem=8G -N 1 -n 2 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$picard_mark_duplicates_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_6_JOB_ID: picard_mark_duplicates.GM12878_WCE_rep2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.GM12878_WCE_rep2
JOB_DEPENDENCIES=$picard_merge_sam_files_6_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.GM12878_WCE_rep2.4db93b4a13a63b5d36578255f8c2e3ea.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.GM12878_WCE_rep2.4db93b4a13a63b5d36578255f8c2e3ea.mugqic.done'
module load java/1.8.0_121 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/GM12878_WCE_rep2/GM12878_WCE_rep2.merged.bam \
 OUTPUT=alignment/GM12878_WCE_rep2/GM12878_WCE_rep2.sorted.dup.bam \
 METRICS_FILE=alignment/GM12878_WCE_rep2/GM12878_WCE_rep2.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.GM12878_WCE_rep2.4db93b4a13a63b5d36578255f8c2e3ea.mugqic.done
)
picard_mark_duplicates_6_JOB_ID=$(echo "#! /bin/bash 
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=48:00:0 --mem=8G -N 1 -n 2 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$picard_mark_duplicates_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_7_JOB_ID: picard_mark_duplicates.GM12878_SMC1_rep2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.GM12878_SMC1_rep2
JOB_DEPENDENCIES=$picard_merge_sam_files_7_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.GM12878_SMC1_rep2.68900ae53962deda18bc690c5bcde262.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.GM12878_SMC1_rep2.68900ae53962deda18bc690c5bcde262.mugqic.done'
module load java/1.8.0_121 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/GM12878_SMC1_rep2/GM12878_SMC1_rep2.merged.bam \
 OUTPUT=alignment/GM12878_SMC1_rep2/GM12878_SMC1_rep2.sorted.dup.bam \
 METRICS_FILE=alignment/GM12878_SMC1_rep2/GM12878_SMC1_rep2.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.GM12878_SMC1_rep2.68900ae53962deda18bc690c5bcde262.mugqic.done
)
picard_mark_duplicates_7_JOB_ID=$(echo "#! /bin/bash 
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=48:00:0 --mem=8G -N 1 -n 2 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$picard_mark_duplicates_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_8_JOB_ID: picard_mark_duplicates.GM12878_CDK9_rep2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.GM12878_CDK9_rep2
JOB_DEPENDENCIES=$picard_merge_sam_files_8_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.GM12878_CDK9_rep2.f7192efbc65db97a3d4b40fcd9f36eed.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.GM12878_CDK9_rep2.f7192efbc65db97a3d4b40fcd9f36eed.mugqic.done'
module load java/1.8.0_121 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/GM12878_CDK9_rep2/GM12878_CDK9_rep2.merged.bam \
 OUTPUT=alignment/GM12878_CDK9_rep2/GM12878_CDK9_rep2.sorted.dup.bam \
 METRICS_FILE=alignment/GM12878_CDK9_rep2/GM12878_CDK9_rep2.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.GM12878_CDK9_rep2.f7192efbc65db97a3d4b40fcd9f36eed.mugqic.done
)
picard_mark_duplicates_8_JOB_ID=$(echo "#! /bin/bash 
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=48:00:0 --mem=8G -N 1 -n 2 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$picard_mark_duplicates_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_9_JOB_ID: picard_mark_duplicates.GM12878_NIPBL_rep2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.GM12878_NIPBL_rep2
JOB_DEPENDENCIES=$picard_merge_sam_files_9_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.GM12878_NIPBL_rep2.0f0b326a06cca86829f65bcb79c5f9b4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.GM12878_NIPBL_rep2.0f0b326a06cca86829f65bcb79c5f9b4.mugqic.done'
module load java/1.8.0_121 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/GM12878_NIPBL_rep2/GM12878_NIPBL_rep2.merged.bam \
 OUTPUT=alignment/GM12878_NIPBL_rep2/GM12878_NIPBL_rep2.sorted.dup.bam \
 METRICS_FILE=alignment/GM12878_NIPBL_rep2/GM12878_NIPBL_rep2.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.GM12878_NIPBL_rep2.0f0b326a06cca86829f65bcb79c5f9b4.mugqic.done
)
picard_mark_duplicates_9_JOB_ID=$(echo "#! /bin/bash 
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=48:00:0 --mem=8G -N 1 -n 2 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$picard_mark_duplicates_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_10_JOB_ID: picard_mark_duplicates.GM12878_MED1_rep2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.GM12878_MED1_rep2
JOB_DEPENDENCIES=$picard_merge_sam_files_10_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.GM12878_MED1_rep2.d59eebf0df894f591c3420404687832a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.GM12878_MED1_rep2.d59eebf0df894f591c3420404687832a.mugqic.done'
module load java/1.8.0_121 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/GM12878_MED1_rep2/GM12878_MED1_rep2.merged.bam \
 OUTPUT=alignment/GM12878_MED1_rep2/GM12878_MED1_rep2.sorted.dup.bam \
 METRICS_FILE=alignment/GM12878_MED1_rep2/GM12878_MED1_rep2.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.GM12878_MED1_rep2.d59eebf0df894f591c3420404687832a.mugqic.done
)
picard_mark_duplicates_10_JOB_ID=$(echo "#! /bin/bash 
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=48:00:0 --mem=8G -N 1 -n 2 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$picard_mark_duplicates_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_11_JOB_ID: picard_mark_duplicates_report
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates_report
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_2_JOB_ID:$picard_mark_duplicates_3_JOB_ID:$picard_mark_duplicates_4_JOB_ID:$picard_mark_duplicates_5_JOB_ID:$picard_mark_duplicates_6_JOB_ID:$picard_mark_duplicates_7_JOB_ID:$picard_mark_duplicates_8_JOB_ID:$picard_mark_duplicates_9_JOB_ID:$picard_mark_duplicates_10_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates_report.8bfd4b03ed7ee119ffcd84a77e108045.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates_report.8bfd4b03ed7ee119ffcd84a77e108045.mugqic.done'
mkdir -p report && \
cp \
  /home/efournie/genpipes/bfx/report/ChipSeq.picard_mark_duplicates.md \
  report/ChipSeq.picard_mark_duplicates.md
picard_mark_duplicates_report.8bfd4b03ed7ee119ffcd84a77e108045.mugqic.done
)
picard_mark_duplicates_11_JOB_ID=$(echo "#! /bin/bash 
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
echo "$picard_mark_duplicates_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# STEP: metrics
#-------------------------------------------------------------------------------
STEP=metrics
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: metrics_1_JOB_ID: metrics.flagstat
#-------------------------------------------------------------------------------
JOB_NAME=metrics.flagstat
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_2_JOB_ID:$picard_mark_duplicates_3_JOB_ID:$picard_mark_duplicates_4_JOB_ID:$picard_mark_duplicates_5_JOB_ID:$picard_mark_duplicates_6_JOB_ID:$picard_mark_duplicates_7_JOB_ID:$picard_mark_duplicates_8_JOB_ID:$picard_mark_duplicates_9_JOB_ID:$picard_mark_duplicates_10_JOB_ID
JOB_DONE=job_output/metrics/metrics.flagstat.4581caacd43a128044ab67b29d9d838e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'metrics.flagstat.4581caacd43a128044ab67b29d9d838e.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools flagstat \
  alignment/GM12878_WCE_rep1/GM12878_WCE_rep1.sorted.dup.bam \
  > alignment/GM12878_WCE_rep1/GM12878_WCE_rep1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/GM12878_NIPBL_rep1/GM12878_NIPBL_rep1.sorted.dup.bam \
  > alignment/GM12878_NIPBL_rep1/GM12878_NIPBL_rep1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/GM12878_MED1_rep1/GM12878_MED1_rep1.sorted.dup.bam \
  > alignment/GM12878_MED1_rep1/GM12878_MED1_rep1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/GM12878_SMC1_rep1/GM12878_SMC1_rep1.sorted.dup.bam \
  > alignment/GM12878_SMC1_rep1/GM12878_SMC1_rep1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/GM12878_CDK9_rep1/GM12878_CDK9_rep1.sorted.dup.bam \
  > alignment/GM12878_CDK9_rep1/GM12878_CDK9_rep1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/GM12878_WCE_rep2/GM12878_WCE_rep2.sorted.dup.bam \
  > alignment/GM12878_WCE_rep2/GM12878_WCE_rep2.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/GM12878_SMC1_rep2/GM12878_SMC1_rep2.sorted.dup.bam \
  > alignment/GM12878_SMC1_rep2/GM12878_SMC1_rep2.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/GM12878_CDK9_rep2/GM12878_CDK9_rep2.sorted.dup.bam \
  > alignment/GM12878_CDK9_rep2/GM12878_CDK9_rep2.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/GM12878_NIPBL_rep2/GM12878_NIPBL_rep2.sorted.dup.bam \
  > alignment/GM12878_NIPBL_rep2/GM12878_NIPBL_rep2.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/GM12878_MED1_rep2/GM12878_MED1_rep2.sorted.dup.bam \
  > alignment/GM12878_MED1_rep2/GM12878_MED1_rep2.sorted.dup.bam.flagstat
metrics.flagstat.4581caacd43a128044ab67b29d9d838e.mugqic.done
)
metrics_1_JOB_ID=$(echo "#! /bin/bash 
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
echo "$metrics_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: metrics_2_JOB_ID: metrics_report
#-------------------------------------------------------------------------------
JOB_NAME=metrics_report
JOB_DEPENDENCIES=$metrics_1_JOB_ID
JOB_DONE=job_output/metrics/metrics_report.d3e95205f8292b478abd0f1711e86a36.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'metrics_report.d3e95205f8292b478abd0f1711e86a36.mugqic.done'
module load mugqic/pandoc/1.15.2 && \
for sample in GM12878_WCE_rep1 GM12878_NIPBL_rep1 GM12878_MED1_rep1 GM12878_SMC1_rep1 GM12878_CDK9_rep1 GM12878_WCE_rep2 GM12878_SMC1_rep2 GM12878_CDK9_rep2 GM12878_NIPBL_rep2 GM12878_MED1_rep2
do
  flagstat_file=alignment/$sample/$sample.sorted.dup.bam.flagstat
  echo -e "$sample	`grep -P '^\d+ \+ \d+ mapped' $flagstat_file | grep -Po '^\d+'`	`grep -P '^\d+ \+ \d+ duplicate' $flagstat_file | grep -Po '^\d+'`"
done | \
awk -F"	" '{OFS="	"; print $0, $3 / $2 * 100}' | sed '1iSample	Aligned Filtered Reads	Duplicate Reads	Duplicate %' \
  > metrics/SampleMetrics.stats && \
mkdir -p report && \
if [[ -f metrics/trimSampleTable.tsv ]]
then
  awk -F "	" 'FNR==NR{trim_line[$1]=$0; surviving[$1]=$3; next}{OFS="	"; if ($1=="Sample") {print trim_line[$1], $2, "Aligned Filtered %", $3, $4} else {print trim_line[$1], $2, $2 / surviving[$1] * 100, $3, $4}}' metrics/trimSampleTable.tsv metrics/SampleMetrics.stats \
  > report/trimMemSampleTable.tsv
else
  cp metrics/SampleMetrics.stats report/trimMemSampleTable.tsv
fi && \
trim_mem_sample_table=`if [[ -f metrics/trimSampleTable.tsv ]] ; then LC_NUMERIC=en_CA awk -F "	" '{OFS="|"; if (NR == 1) {$1 = $1; print $0; print "-----|-----:|-----:|-----:|-----:|-----:|-----:|-----:"} else {print $1, sprintf("%\47d", $2), sprintf("%\47d", $3), sprintf("%.1f", $4), sprintf("%\47d", $5), sprintf("%.1f", $6), sprintf("%\47d", $7), sprintf("%.1f", $8)}}' report/trimMemSampleTable.tsv ; else LC_NUMERIC=en_CA awk -F "	" '{OFS="|"; if (NR == 1) {$1 = $1; print $0; print "-----|-----:|-----:|-----:"} else {print $1, sprintf("%\47d", $2), sprintf("%\47d", $3), sprintf("%.1f", $4)}}' report/trimMemSampleTable.tsv ; fi` && \
pandoc --to=markdown \
  --template /home/efournie/genpipes/bfx/report/ChipSeq.metrics.md \
  --variable trim_mem_sample_table="$trim_mem_sample_table" \
  /home/efournie/genpipes/bfx/report/ChipSeq.metrics.md \
  > report/ChipSeq.metrics.md

metrics_report.d3e95205f8292b478abd0f1711e86a36.mugqic.done
)
metrics_2_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"metrics\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_WCE_rep1.json,/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_NIPBL_rep1.json,/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_MED1_rep1.json,/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_SMC1_rep1.json,/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_CDK9_rep1.json,/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_WCE_rep2.json,/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_SMC1_rep2.json,/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_CDK9_rep2.json,/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_NIPBL_rep2.json,/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_MED1_rep2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"metrics\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_WCE_rep1.json,/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_NIPBL_rep1.json,/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_MED1_rep1.json,/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_SMC1_rep1.json,/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_CDK9_rep1.json,/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_WCE_rep2.json,/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_SMC1_rep2.json,/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_CDK9_rep2.json,/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_NIPBL_rep2.json,/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_MED1_rep2.json\" \
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
echo "$metrics_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# STEP: homer_make_tag_directory
#-------------------------------------------------------------------------------
STEP=homer_make_tag_directory
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_1_JOB_ID: homer_make_tag_directory.GM12878_WCE_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.GM12878_WCE_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.GM12878_WCE_rep1.7f792f57b40ce75f938e2908e25307bc.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.GM12878_WCE_rep1.7f792f57b40ce75f938e2908e25307bc.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/GM12878_WCE_rep1 \
            alignment/GM12878_WCE_rep1/GM12878_WCE_rep1.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.GM12878_WCE_rep1.7f792f57b40ce75f938e2908e25307bc.mugqic.done
)
homer_make_tag_directory_1_JOB_ID=$(echo "#! /bin/bash 
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
echo "$homer_make_tag_directory_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_2_JOB_ID: homer_make_tag_directory.GM12878_NIPBL_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.GM12878_NIPBL_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_2_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.GM12878_NIPBL_rep1.a6c9827b1901274b0495553165d56f0f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.GM12878_NIPBL_rep1.a6c9827b1901274b0495553165d56f0f.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/GM12878_NIPBL_rep1 \
            alignment/GM12878_NIPBL_rep1/GM12878_NIPBL_rep1.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.GM12878_NIPBL_rep1.a6c9827b1901274b0495553165d56f0f.mugqic.done
)
homer_make_tag_directory_2_JOB_ID=$(echo "#! /bin/bash 
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
echo "$homer_make_tag_directory_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_3_JOB_ID: homer_make_tag_directory.GM12878_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.GM12878_MED1_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_3_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.GM12878_MED1_rep1.b7e63408e802cf8cf249f291b3e08797.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.GM12878_MED1_rep1.b7e63408e802cf8cf249f291b3e08797.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/GM12878_MED1_rep1 \
            alignment/GM12878_MED1_rep1/GM12878_MED1_rep1.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.GM12878_MED1_rep1.b7e63408e802cf8cf249f291b3e08797.mugqic.done
)
homer_make_tag_directory_3_JOB_ID=$(echo "#! /bin/bash 
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
echo "$homer_make_tag_directory_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_4_JOB_ID: homer_make_tag_directory.GM12878_SMC1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.GM12878_SMC1_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_4_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.GM12878_SMC1_rep1.4a066a6095fb326ee7629b79a98dd337.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.GM12878_SMC1_rep1.4a066a6095fb326ee7629b79a98dd337.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/GM12878_SMC1_rep1 \
            alignment/GM12878_SMC1_rep1/GM12878_SMC1_rep1.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.GM12878_SMC1_rep1.4a066a6095fb326ee7629b79a98dd337.mugqic.done
)
homer_make_tag_directory_4_JOB_ID=$(echo "#! /bin/bash 
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
echo "$homer_make_tag_directory_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_5_JOB_ID: homer_make_tag_directory.GM12878_CDK9_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.GM12878_CDK9_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_5_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.GM12878_CDK9_rep1.2e4133b1c4fe6bbc4322661134d05b68.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.GM12878_CDK9_rep1.2e4133b1c4fe6bbc4322661134d05b68.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/GM12878_CDK9_rep1 \
            alignment/GM12878_CDK9_rep1/GM12878_CDK9_rep1.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.GM12878_CDK9_rep1.2e4133b1c4fe6bbc4322661134d05b68.mugqic.done
)
homer_make_tag_directory_5_JOB_ID=$(echo "#! /bin/bash 
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
echo "$homer_make_tag_directory_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_6_JOB_ID: homer_make_tag_directory.GM12878_WCE_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.GM12878_WCE_rep2
JOB_DEPENDENCIES=$picard_mark_duplicates_6_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.GM12878_WCE_rep2.480c48860f51acc56025a7e10bc1c68e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.GM12878_WCE_rep2.480c48860f51acc56025a7e10bc1c68e.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/GM12878_WCE_rep2 \
            alignment/GM12878_WCE_rep2/GM12878_WCE_rep2.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.GM12878_WCE_rep2.480c48860f51acc56025a7e10bc1c68e.mugqic.done
)
homer_make_tag_directory_6_JOB_ID=$(echo "#! /bin/bash 
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
echo "$homer_make_tag_directory_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_7_JOB_ID: homer_make_tag_directory.GM12878_SMC1_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.GM12878_SMC1_rep2
JOB_DEPENDENCIES=$picard_mark_duplicates_7_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.GM12878_SMC1_rep2.03f4fb297c6570672170c5dea35d27b0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.GM12878_SMC1_rep2.03f4fb297c6570672170c5dea35d27b0.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/GM12878_SMC1_rep2 \
            alignment/GM12878_SMC1_rep2/GM12878_SMC1_rep2.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.GM12878_SMC1_rep2.03f4fb297c6570672170c5dea35d27b0.mugqic.done
)
homer_make_tag_directory_7_JOB_ID=$(echo "#! /bin/bash 
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
echo "$homer_make_tag_directory_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_8_JOB_ID: homer_make_tag_directory.GM12878_CDK9_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.GM12878_CDK9_rep2
JOB_DEPENDENCIES=$picard_mark_duplicates_8_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.GM12878_CDK9_rep2.79a938a7c3960f3c23135e862aaea1e3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.GM12878_CDK9_rep2.79a938a7c3960f3c23135e862aaea1e3.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/GM12878_CDK9_rep2 \
            alignment/GM12878_CDK9_rep2/GM12878_CDK9_rep2.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.GM12878_CDK9_rep2.79a938a7c3960f3c23135e862aaea1e3.mugqic.done
)
homer_make_tag_directory_8_JOB_ID=$(echo "#! /bin/bash 
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
echo "$homer_make_tag_directory_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_9_JOB_ID: homer_make_tag_directory.GM12878_NIPBL_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.GM12878_NIPBL_rep2
JOB_DEPENDENCIES=$picard_mark_duplicates_9_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.GM12878_NIPBL_rep2.9274e6443313b30459341504b412061d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.GM12878_NIPBL_rep2.9274e6443313b30459341504b412061d.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/GM12878_NIPBL_rep2 \
            alignment/GM12878_NIPBL_rep2/GM12878_NIPBL_rep2.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.GM12878_NIPBL_rep2.9274e6443313b30459341504b412061d.mugqic.done
)
homer_make_tag_directory_9_JOB_ID=$(echo "#! /bin/bash 
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
echo "$homer_make_tag_directory_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_10_JOB_ID: homer_make_tag_directory.GM12878_MED1_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.GM12878_MED1_rep2
JOB_DEPENDENCIES=$picard_mark_duplicates_10_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.GM12878_MED1_rep2.ff941c14985fbe93a6476d01de6da4f9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.GM12878_MED1_rep2.ff941c14985fbe93a6476d01de6da4f9.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/GM12878_MED1_rep2 \
            alignment/GM12878_MED1_rep2/GM12878_MED1_rep2.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.GM12878_MED1_rep2.ff941c14985fbe93a6476d01de6da4f9.mugqic.done
)
homer_make_tag_directory_10_JOB_ID=$(echo "#! /bin/bash 
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
echo "$homer_make_tag_directory_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# STEP: qc_metrics
#-------------------------------------------------------------------------------
STEP=qc_metrics
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: qc_metrics_1_JOB_ID: qc_plots_R
#-------------------------------------------------------------------------------
JOB_NAME=qc_plots_R
JOB_DEPENDENCIES=$homer_make_tag_directory_1_JOB_ID:$homer_make_tag_directory_2_JOB_ID:$homer_make_tag_directory_3_JOB_ID:$homer_make_tag_directory_4_JOB_ID:$homer_make_tag_directory_5_JOB_ID:$homer_make_tag_directory_6_JOB_ID:$homer_make_tag_directory_7_JOB_ID:$homer_make_tag_directory_8_JOB_ID:$homer_make_tag_directory_9_JOB_ID:$homer_make_tag_directory_10_JOB_ID
JOB_DONE=job_output/qc_metrics/qc_plots_R.a6d02d1d9c228c22f84905e9bc6e9fcc.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'qc_plots_R.a6d02d1d9c228c22f84905e9bc6e9fcc.mugqic.done'
module load mugqic/mugqic_tools/2.1.9 mugqic/R_Bioconductor/3.2.3_3.2 && \
mkdir -p graphs && \
Rscript $R_TOOLS/chipSeqGenerateQCMetrics.R \
  ../../raw/design.txt \
  /project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38 && \
cp /home/efournie/genpipes/bfx/report/ChipSeq.qc_metrics.md report/ChipSeq.qc_metrics.md && \
for sample in GM12878_WCE_rep1 GM12878_NIPBL_rep1 GM12878_MED1_rep1 GM12878_SMC1_rep1 GM12878_CDK9_rep1 GM12878_WCE_rep2 GM12878_SMC1_rep2 GM12878_CDK9_rep2 GM12878_NIPBL_rep2 GM12878_MED1_rep2
do
  cp --parents graphs/${sample}_QC_Metrics.ps report/
  convert -rotate 90 graphs/${sample}_QC_Metrics.ps report/graphs/${sample}_QC_Metrics.png
  echo -e "----

![QC Metrics for Sample $sample ([download high-res image](graphs/${sample}_QC_Metrics.ps))](graphs/${sample}_QC_Metrics.png)
" \
  >> report/ChipSeq.qc_metrics.md
done
qc_plots_R.a6d02d1d9c228c22f84905e9bc6e9fcc.mugqic.done
)
qc_metrics_1_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"qc_metrics\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_WCE_rep1.json,/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_NIPBL_rep1.json,/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_MED1_rep1.json,/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_SMC1_rep1.json,/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_CDK9_rep1.json,/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_WCE_rep2.json,/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_SMC1_rep2.json,/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_CDK9_rep2.json,/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_NIPBL_rep2.json,/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_MED1_rep2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"qc_metrics\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_WCE_rep1.json,/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_NIPBL_rep1.json,/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_MED1_rep1.json,/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_SMC1_rep1.json,/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_CDK9_rep1.json,/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_WCE_rep2.json,/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_SMC1_rep2.json,/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_CDK9_rep2.json,/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_NIPBL_rep2.json,/project/6001942/Working_Directory/Eric/CofactorHR/GM12878/output/chip-pipeline-GRCh38/json/GM12878_MED1_rep2.json\" \
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
echo "$qc_metrics_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# STEP: homer_make_ucsc_file
#-------------------------------------------------------------------------------
STEP=homer_make_ucsc_file
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_1_JOB_ID: homer_make_ucsc_file.GM12878_WCE_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.GM12878_WCE_rep1
JOB_DEPENDENCIES=$homer_make_tag_directory_1_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.GM12878_WCE_rep1.67b2dcfdef8a280b02b232d2774287de.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.GM12878_WCE_rep1.67b2dcfdef8a280b02b232d2774287de.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 && \
mkdir -p tracks/GM12878_WCE_rep1 && \
makeUCSCfile \
        tags/GM12878_WCE_rep1 > tracks/GM12878_WCE_rep1/GM12878_WCE_rep1.ucsc.bedGraph && \
        gzip -c tracks/GM12878_WCE_rep1/GM12878_WCE_rep1.ucsc.bedGraph > tracks/GM12878_WCE_rep1/GM12878_WCE_rep1.ucsc.bedGraph.gz
homer_make_ucsc_file.GM12878_WCE_rep1.67b2dcfdef8a280b02b232d2774287de.mugqic.done
)
homer_make_ucsc_file_1_JOB_ID=$(echo "#! /bin/bash 
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
echo "$homer_make_ucsc_file_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_2_JOB_ID: homer_make_ucsc_file_bigWig.GM12878_WCE_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file_bigWig.GM12878_WCE_rep1
JOB_DEPENDENCIES=$homer_make_ucsc_file_1_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file_bigWig.GM12878_WCE_rep1.746c99a24909a98e6cf78a2b682a3b35.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file_bigWig.GM12878_WCE_rep1.746c99a24909a98e6cf78a2b682a3b35.mugqic.done'
module load mugqic/ucsc/v346 && \
mkdir -p tracks/GM12878_WCE_rep1/bigWig && \
export TMPDIR=${SLURM_TMPDIR} && \
cat tracks/GM12878_WCE_rep1/GM12878_WCE_rep1.ucsc.bedGraph | head -n 1 > tracks/GM12878_WCE_rep1/GM12878_WCE_rep1.ucsc.bedGraph.head.tmp && \
cat tracks/GM12878_WCE_rep1/GM12878_WCE_rep1.ucsc.bedGraph | awk ' NR > 1 ' | sort  --temporary-directory=${SLURM_TMPDIR} -k1,1 -k2,2n | \
awk '{if($0 !~ /^[A-W]/) print ""$0; else print $0}' | grep -vP "GL|lambda|pUC19" | sed 's/MT/chrM/' | \
awk '{printf "%s\t%d\t%d\t%4.4g\n", $1,$2,$3,$4}' > tracks/GM12878_WCE_rep1/GM12878_WCE_rep1.ucsc.bedGraph.body.tmp && \
cat tracks/GM12878_WCE_rep1/GM12878_WCE_rep1.ucsc.bedGraph.head.tmp tracks/GM12878_WCE_rep1/GM12878_WCE_rep1.ucsc.bedGraph.body.tmp > tracks/GM12878_WCE_rep1/GM12878_WCE_rep1.ucsc.bedGraph.sorted && \
rm tracks/GM12878_WCE_rep1/GM12878_WCE_rep1.ucsc.bedGraph.head.tmp tracks/GM12878_WCE_rep1/GM12878_WCE_rep1.ucsc.bedGraph.body.tmp && \
bedGraphToBigWig \
  tracks/GM12878_WCE_rep1/GM12878_WCE_rep1.ucsc.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  tracks/GM12878_WCE_rep1/bigWig/GM12878_WCE_rep1.bw
homer_make_ucsc_file_bigWig.GM12878_WCE_rep1.746c99a24909a98e6cf78a2b682a3b35.mugqic.done
)
homer_make_ucsc_file_2_JOB_ID=$(echo "#! /bin/bash 
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
echo "$homer_make_ucsc_file_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_3_JOB_ID: homer_make_ucsc_file.GM12878_NIPBL_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.GM12878_NIPBL_rep1
JOB_DEPENDENCIES=$homer_make_tag_directory_2_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.GM12878_NIPBL_rep1.748aecbeff8991bbb8b25cb16c4b3831.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.GM12878_NIPBL_rep1.748aecbeff8991bbb8b25cb16c4b3831.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 && \
mkdir -p tracks/GM12878_NIPBL_rep1 && \
makeUCSCfile \
        tags/GM12878_NIPBL_rep1 > tracks/GM12878_NIPBL_rep1/GM12878_NIPBL_rep1.ucsc.bedGraph && \
        gzip -c tracks/GM12878_NIPBL_rep1/GM12878_NIPBL_rep1.ucsc.bedGraph > tracks/GM12878_NIPBL_rep1/GM12878_NIPBL_rep1.ucsc.bedGraph.gz
homer_make_ucsc_file.GM12878_NIPBL_rep1.748aecbeff8991bbb8b25cb16c4b3831.mugqic.done
)
homer_make_ucsc_file_3_JOB_ID=$(echo "#! /bin/bash 
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
echo "$homer_make_ucsc_file_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_4_JOB_ID: homer_make_ucsc_file_bigWig.GM12878_NIPBL_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file_bigWig.GM12878_NIPBL_rep1
JOB_DEPENDENCIES=$homer_make_ucsc_file_3_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file_bigWig.GM12878_NIPBL_rep1.58b26bcfbb52f948d7260c9efeef8ec0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file_bigWig.GM12878_NIPBL_rep1.58b26bcfbb52f948d7260c9efeef8ec0.mugqic.done'
module load mugqic/ucsc/v346 && \
mkdir -p tracks/GM12878_NIPBL_rep1/bigWig && \
export TMPDIR=${SLURM_TMPDIR} && \
cat tracks/GM12878_NIPBL_rep1/GM12878_NIPBL_rep1.ucsc.bedGraph | head -n 1 > tracks/GM12878_NIPBL_rep1/GM12878_NIPBL_rep1.ucsc.bedGraph.head.tmp && \
cat tracks/GM12878_NIPBL_rep1/GM12878_NIPBL_rep1.ucsc.bedGraph | awk ' NR > 1 ' | sort  --temporary-directory=${SLURM_TMPDIR} -k1,1 -k2,2n | \
awk '{if($0 !~ /^[A-W]/) print ""$0; else print $0}' | grep -vP "GL|lambda|pUC19" | sed 's/MT/chrM/' | \
awk '{printf "%s\t%d\t%d\t%4.4g\n", $1,$2,$3,$4}' > tracks/GM12878_NIPBL_rep1/GM12878_NIPBL_rep1.ucsc.bedGraph.body.tmp && \
cat tracks/GM12878_NIPBL_rep1/GM12878_NIPBL_rep1.ucsc.bedGraph.head.tmp tracks/GM12878_NIPBL_rep1/GM12878_NIPBL_rep1.ucsc.bedGraph.body.tmp > tracks/GM12878_NIPBL_rep1/GM12878_NIPBL_rep1.ucsc.bedGraph.sorted && \
rm tracks/GM12878_NIPBL_rep1/GM12878_NIPBL_rep1.ucsc.bedGraph.head.tmp tracks/GM12878_NIPBL_rep1/GM12878_NIPBL_rep1.ucsc.bedGraph.body.tmp && \
bedGraphToBigWig \
  tracks/GM12878_NIPBL_rep1/GM12878_NIPBL_rep1.ucsc.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  tracks/GM12878_NIPBL_rep1/bigWig/GM12878_NIPBL_rep1.bw
homer_make_ucsc_file_bigWig.GM12878_NIPBL_rep1.58b26bcfbb52f948d7260c9efeef8ec0.mugqic.done
)
homer_make_ucsc_file_4_JOB_ID=$(echo "#! /bin/bash 
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
echo "$homer_make_ucsc_file_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_5_JOB_ID: homer_make_ucsc_file.GM12878_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.GM12878_MED1_rep1
JOB_DEPENDENCIES=$homer_make_tag_directory_3_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.GM12878_MED1_rep1.86c57cb268365d07e9d7849a198732c9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.GM12878_MED1_rep1.86c57cb268365d07e9d7849a198732c9.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 && \
mkdir -p tracks/GM12878_MED1_rep1 && \
makeUCSCfile \
        tags/GM12878_MED1_rep1 > tracks/GM12878_MED1_rep1/GM12878_MED1_rep1.ucsc.bedGraph && \
        gzip -c tracks/GM12878_MED1_rep1/GM12878_MED1_rep1.ucsc.bedGraph > tracks/GM12878_MED1_rep1/GM12878_MED1_rep1.ucsc.bedGraph.gz
homer_make_ucsc_file.GM12878_MED1_rep1.86c57cb268365d07e9d7849a198732c9.mugqic.done
)
homer_make_ucsc_file_5_JOB_ID=$(echo "#! /bin/bash 
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
echo "$homer_make_ucsc_file_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_6_JOB_ID: homer_make_ucsc_file_bigWig.GM12878_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file_bigWig.GM12878_MED1_rep1
JOB_DEPENDENCIES=$homer_make_ucsc_file_5_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file_bigWig.GM12878_MED1_rep1.cd2da6ccff55bfa96bf5a404581d2f54.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file_bigWig.GM12878_MED1_rep1.cd2da6ccff55bfa96bf5a404581d2f54.mugqic.done'
module load mugqic/ucsc/v346 && \
mkdir -p tracks/GM12878_MED1_rep1/bigWig && \
export TMPDIR=${SLURM_TMPDIR} && \
cat tracks/GM12878_MED1_rep1/GM12878_MED1_rep1.ucsc.bedGraph | head -n 1 > tracks/GM12878_MED1_rep1/GM12878_MED1_rep1.ucsc.bedGraph.head.tmp && \
cat tracks/GM12878_MED1_rep1/GM12878_MED1_rep1.ucsc.bedGraph | awk ' NR > 1 ' | sort  --temporary-directory=${SLURM_TMPDIR} -k1,1 -k2,2n | \
awk '{if($0 !~ /^[A-W]/) print ""$0; else print $0}' | grep -vP "GL|lambda|pUC19" | sed 's/MT/chrM/' | \
awk '{printf "%s\t%d\t%d\t%4.4g\n", $1,$2,$3,$4}' > tracks/GM12878_MED1_rep1/GM12878_MED1_rep1.ucsc.bedGraph.body.tmp && \
cat tracks/GM12878_MED1_rep1/GM12878_MED1_rep1.ucsc.bedGraph.head.tmp tracks/GM12878_MED1_rep1/GM12878_MED1_rep1.ucsc.bedGraph.body.tmp > tracks/GM12878_MED1_rep1/GM12878_MED1_rep1.ucsc.bedGraph.sorted && \
rm tracks/GM12878_MED1_rep1/GM12878_MED1_rep1.ucsc.bedGraph.head.tmp tracks/GM12878_MED1_rep1/GM12878_MED1_rep1.ucsc.bedGraph.body.tmp && \
bedGraphToBigWig \
  tracks/GM12878_MED1_rep1/GM12878_MED1_rep1.ucsc.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  tracks/GM12878_MED1_rep1/bigWig/GM12878_MED1_rep1.bw
homer_make_ucsc_file_bigWig.GM12878_MED1_rep1.cd2da6ccff55bfa96bf5a404581d2f54.mugqic.done
)
homer_make_ucsc_file_6_JOB_ID=$(echo "#! /bin/bash 
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
echo "$homer_make_ucsc_file_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_7_JOB_ID: homer_make_ucsc_file.GM12878_SMC1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.GM12878_SMC1_rep1
JOB_DEPENDENCIES=$homer_make_tag_directory_4_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.GM12878_SMC1_rep1.d86748c8f36c2bc5e1a46c993dae687e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.GM12878_SMC1_rep1.d86748c8f36c2bc5e1a46c993dae687e.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 && \
mkdir -p tracks/GM12878_SMC1_rep1 && \
makeUCSCfile \
        tags/GM12878_SMC1_rep1 > tracks/GM12878_SMC1_rep1/GM12878_SMC1_rep1.ucsc.bedGraph && \
        gzip -c tracks/GM12878_SMC1_rep1/GM12878_SMC1_rep1.ucsc.bedGraph > tracks/GM12878_SMC1_rep1/GM12878_SMC1_rep1.ucsc.bedGraph.gz
homer_make_ucsc_file.GM12878_SMC1_rep1.d86748c8f36c2bc5e1a46c993dae687e.mugqic.done
)
homer_make_ucsc_file_7_JOB_ID=$(echo "#! /bin/bash 
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
echo "$homer_make_ucsc_file_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_8_JOB_ID: homer_make_ucsc_file_bigWig.GM12878_SMC1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file_bigWig.GM12878_SMC1_rep1
JOB_DEPENDENCIES=$homer_make_ucsc_file_7_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file_bigWig.GM12878_SMC1_rep1.d803c0edccb0048dd07b5dfc97e56e9d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file_bigWig.GM12878_SMC1_rep1.d803c0edccb0048dd07b5dfc97e56e9d.mugqic.done'
module load mugqic/ucsc/v346 && \
mkdir -p tracks/GM12878_SMC1_rep1/bigWig && \
export TMPDIR=${SLURM_TMPDIR} && \
cat tracks/GM12878_SMC1_rep1/GM12878_SMC1_rep1.ucsc.bedGraph | head -n 1 > tracks/GM12878_SMC1_rep1/GM12878_SMC1_rep1.ucsc.bedGraph.head.tmp && \
cat tracks/GM12878_SMC1_rep1/GM12878_SMC1_rep1.ucsc.bedGraph | awk ' NR > 1 ' | sort  --temporary-directory=${SLURM_TMPDIR} -k1,1 -k2,2n | \
awk '{if($0 !~ /^[A-W]/) print ""$0; else print $0}' | grep -vP "GL|lambda|pUC19" | sed 's/MT/chrM/' | \
awk '{printf "%s\t%d\t%d\t%4.4g\n", $1,$2,$3,$4}' > tracks/GM12878_SMC1_rep1/GM12878_SMC1_rep1.ucsc.bedGraph.body.tmp && \
cat tracks/GM12878_SMC1_rep1/GM12878_SMC1_rep1.ucsc.bedGraph.head.tmp tracks/GM12878_SMC1_rep1/GM12878_SMC1_rep1.ucsc.bedGraph.body.tmp > tracks/GM12878_SMC1_rep1/GM12878_SMC1_rep1.ucsc.bedGraph.sorted && \
rm tracks/GM12878_SMC1_rep1/GM12878_SMC1_rep1.ucsc.bedGraph.head.tmp tracks/GM12878_SMC1_rep1/GM12878_SMC1_rep1.ucsc.bedGraph.body.tmp && \
bedGraphToBigWig \
  tracks/GM12878_SMC1_rep1/GM12878_SMC1_rep1.ucsc.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  tracks/GM12878_SMC1_rep1/bigWig/GM12878_SMC1_rep1.bw
homer_make_ucsc_file_bigWig.GM12878_SMC1_rep1.d803c0edccb0048dd07b5dfc97e56e9d.mugqic.done
)
homer_make_ucsc_file_8_JOB_ID=$(echo "#! /bin/bash 
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
echo "$homer_make_ucsc_file_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_9_JOB_ID: homer_make_ucsc_file.GM12878_CDK9_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.GM12878_CDK9_rep1
JOB_DEPENDENCIES=$homer_make_tag_directory_5_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.GM12878_CDK9_rep1.84020cce9f6d1345dcd1ae7d31aeeddd.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.GM12878_CDK9_rep1.84020cce9f6d1345dcd1ae7d31aeeddd.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 && \
mkdir -p tracks/GM12878_CDK9_rep1 && \
makeUCSCfile \
        tags/GM12878_CDK9_rep1 > tracks/GM12878_CDK9_rep1/GM12878_CDK9_rep1.ucsc.bedGraph && \
        gzip -c tracks/GM12878_CDK9_rep1/GM12878_CDK9_rep1.ucsc.bedGraph > tracks/GM12878_CDK9_rep1/GM12878_CDK9_rep1.ucsc.bedGraph.gz
homer_make_ucsc_file.GM12878_CDK9_rep1.84020cce9f6d1345dcd1ae7d31aeeddd.mugqic.done
)
homer_make_ucsc_file_9_JOB_ID=$(echo "#! /bin/bash 
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
echo "$homer_make_ucsc_file_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_10_JOB_ID: homer_make_ucsc_file_bigWig.GM12878_CDK9_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file_bigWig.GM12878_CDK9_rep1
JOB_DEPENDENCIES=$homer_make_ucsc_file_9_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file_bigWig.GM12878_CDK9_rep1.ffed4fa7d2b87a9098f28c67f5788d3b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file_bigWig.GM12878_CDK9_rep1.ffed4fa7d2b87a9098f28c67f5788d3b.mugqic.done'
module load mugqic/ucsc/v346 && \
mkdir -p tracks/GM12878_CDK9_rep1/bigWig && \
export TMPDIR=${SLURM_TMPDIR} && \
cat tracks/GM12878_CDK9_rep1/GM12878_CDK9_rep1.ucsc.bedGraph | head -n 1 > tracks/GM12878_CDK9_rep1/GM12878_CDK9_rep1.ucsc.bedGraph.head.tmp && \
cat tracks/GM12878_CDK9_rep1/GM12878_CDK9_rep1.ucsc.bedGraph | awk ' NR > 1 ' | sort  --temporary-directory=${SLURM_TMPDIR} -k1,1 -k2,2n | \
awk '{if($0 !~ /^[A-W]/) print ""$0; else print $0}' | grep -vP "GL|lambda|pUC19" | sed 's/MT/chrM/' | \
awk '{printf "%s\t%d\t%d\t%4.4g\n", $1,$2,$3,$4}' > tracks/GM12878_CDK9_rep1/GM12878_CDK9_rep1.ucsc.bedGraph.body.tmp && \
cat tracks/GM12878_CDK9_rep1/GM12878_CDK9_rep1.ucsc.bedGraph.head.tmp tracks/GM12878_CDK9_rep1/GM12878_CDK9_rep1.ucsc.bedGraph.body.tmp > tracks/GM12878_CDK9_rep1/GM12878_CDK9_rep1.ucsc.bedGraph.sorted && \
rm tracks/GM12878_CDK9_rep1/GM12878_CDK9_rep1.ucsc.bedGraph.head.tmp tracks/GM12878_CDK9_rep1/GM12878_CDK9_rep1.ucsc.bedGraph.body.tmp && \
bedGraphToBigWig \
  tracks/GM12878_CDK9_rep1/GM12878_CDK9_rep1.ucsc.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  tracks/GM12878_CDK9_rep1/bigWig/GM12878_CDK9_rep1.bw
homer_make_ucsc_file_bigWig.GM12878_CDK9_rep1.ffed4fa7d2b87a9098f28c67f5788d3b.mugqic.done
)
homer_make_ucsc_file_10_JOB_ID=$(echo "#! /bin/bash 
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
echo "$homer_make_ucsc_file_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_11_JOB_ID: homer_make_ucsc_file.GM12878_WCE_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.GM12878_WCE_rep2
JOB_DEPENDENCIES=$homer_make_tag_directory_6_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.GM12878_WCE_rep2.a433412e869f1beae3657129435f48c4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.GM12878_WCE_rep2.a433412e869f1beae3657129435f48c4.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 && \
mkdir -p tracks/GM12878_WCE_rep2 && \
makeUCSCfile \
        tags/GM12878_WCE_rep2 > tracks/GM12878_WCE_rep2/GM12878_WCE_rep2.ucsc.bedGraph && \
        gzip -c tracks/GM12878_WCE_rep2/GM12878_WCE_rep2.ucsc.bedGraph > tracks/GM12878_WCE_rep2/GM12878_WCE_rep2.ucsc.bedGraph.gz
homer_make_ucsc_file.GM12878_WCE_rep2.a433412e869f1beae3657129435f48c4.mugqic.done
)
homer_make_ucsc_file_11_JOB_ID=$(echo "#! /bin/bash 
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
echo "$homer_make_ucsc_file_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_12_JOB_ID: homer_make_ucsc_file_bigWig.GM12878_WCE_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file_bigWig.GM12878_WCE_rep2
JOB_DEPENDENCIES=$homer_make_ucsc_file_11_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file_bigWig.GM12878_WCE_rep2.ef1849478ae1c66bd9ce42c1c5df9d78.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file_bigWig.GM12878_WCE_rep2.ef1849478ae1c66bd9ce42c1c5df9d78.mugqic.done'
module load mugqic/ucsc/v346 && \
mkdir -p tracks/GM12878_WCE_rep2/bigWig && \
export TMPDIR=${SLURM_TMPDIR} && \
cat tracks/GM12878_WCE_rep2/GM12878_WCE_rep2.ucsc.bedGraph | head -n 1 > tracks/GM12878_WCE_rep2/GM12878_WCE_rep2.ucsc.bedGraph.head.tmp && \
cat tracks/GM12878_WCE_rep2/GM12878_WCE_rep2.ucsc.bedGraph | awk ' NR > 1 ' | sort  --temporary-directory=${SLURM_TMPDIR} -k1,1 -k2,2n | \
awk '{if($0 !~ /^[A-W]/) print ""$0; else print $0}' | grep -vP "GL|lambda|pUC19" | sed 's/MT/chrM/' | \
awk '{printf "%s\t%d\t%d\t%4.4g\n", $1,$2,$3,$4}' > tracks/GM12878_WCE_rep2/GM12878_WCE_rep2.ucsc.bedGraph.body.tmp && \
cat tracks/GM12878_WCE_rep2/GM12878_WCE_rep2.ucsc.bedGraph.head.tmp tracks/GM12878_WCE_rep2/GM12878_WCE_rep2.ucsc.bedGraph.body.tmp > tracks/GM12878_WCE_rep2/GM12878_WCE_rep2.ucsc.bedGraph.sorted && \
rm tracks/GM12878_WCE_rep2/GM12878_WCE_rep2.ucsc.bedGraph.head.tmp tracks/GM12878_WCE_rep2/GM12878_WCE_rep2.ucsc.bedGraph.body.tmp && \
bedGraphToBigWig \
  tracks/GM12878_WCE_rep2/GM12878_WCE_rep2.ucsc.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  tracks/GM12878_WCE_rep2/bigWig/GM12878_WCE_rep2.bw
homer_make_ucsc_file_bigWig.GM12878_WCE_rep2.ef1849478ae1c66bd9ce42c1c5df9d78.mugqic.done
)
homer_make_ucsc_file_12_JOB_ID=$(echo "#! /bin/bash 
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
echo "$homer_make_ucsc_file_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_13_JOB_ID: homer_make_ucsc_file.GM12878_SMC1_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.GM12878_SMC1_rep2
JOB_DEPENDENCIES=$homer_make_tag_directory_7_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.GM12878_SMC1_rep2.85e3cf9560819a24abb22bc2917c20bc.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.GM12878_SMC1_rep2.85e3cf9560819a24abb22bc2917c20bc.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 && \
mkdir -p tracks/GM12878_SMC1_rep2 && \
makeUCSCfile \
        tags/GM12878_SMC1_rep2 > tracks/GM12878_SMC1_rep2/GM12878_SMC1_rep2.ucsc.bedGraph && \
        gzip -c tracks/GM12878_SMC1_rep2/GM12878_SMC1_rep2.ucsc.bedGraph > tracks/GM12878_SMC1_rep2/GM12878_SMC1_rep2.ucsc.bedGraph.gz
homer_make_ucsc_file.GM12878_SMC1_rep2.85e3cf9560819a24abb22bc2917c20bc.mugqic.done
)
homer_make_ucsc_file_13_JOB_ID=$(echo "#! /bin/bash 
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
echo "$homer_make_ucsc_file_13_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_14_JOB_ID: homer_make_ucsc_file_bigWig.GM12878_SMC1_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file_bigWig.GM12878_SMC1_rep2
JOB_DEPENDENCIES=$homer_make_ucsc_file_13_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file_bigWig.GM12878_SMC1_rep2.6c63bc09f8f8b392ae354182a96a5a68.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file_bigWig.GM12878_SMC1_rep2.6c63bc09f8f8b392ae354182a96a5a68.mugqic.done'
module load mugqic/ucsc/v346 && \
mkdir -p tracks/GM12878_SMC1_rep2/bigWig && \
export TMPDIR=${SLURM_TMPDIR} && \
cat tracks/GM12878_SMC1_rep2/GM12878_SMC1_rep2.ucsc.bedGraph | head -n 1 > tracks/GM12878_SMC1_rep2/GM12878_SMC1_rep2.ucsc.bedGraph.head.tmp && \
cat tracks/GM12878_SMC1_rep2/GM12878_SMC1_rep2.ucsc.bedGraph | awk ' NR > 1 ' | sort  --temporary-directory=${SLURM_TMPDIR} -k1,1 -k2,2n | \
awk '{if($0 !~ /^[A-W]/) print ""$0; else print $0}' | grep -vP "GL|lambda|pUC19" | sed 's/MT/chrM/' | \
awk '{printf "%s\t%d\t%d\t%4.4g\n", $1,$2,$3,$4}' > tracks/GM12878_SMC1_rep2/GM12878_SMC1_rep2.ucsc.bedGraph.body.tmp && \
cat tracks/GM12878_SMC1_rep2/GM12878_SMC1_rep2.ucsc.bedGraph.head.tmp tracks/GM12878_SMC1_rep2/GM12878_SMC1_rep2.ucsc.bedGraph.body.tmp > tracks/GM12878_SMC1_rep2/GM12878_SMC1_rep2.ucsc.bedGraph.sorted && \
rm tracks/GM12878_SMC1_rep2/GM12878_SMC1_rep2.ucsc.bedGraph.head.tmp tracks/GM12878_SMC1_rep2/GM12878_SMC1_rep2.ucsc.bedGraph.body.tmp && \
bedGraphToBigWig \
  tracks/GM12878_SMC1_rep2/GM12878_SMC1_rep2.ucsc.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  tracks/GM12878_SMC1_rep2/bigWig/GM12878_SMC1_rep2.bw
homer_make_ucsc_file_bigWig.GM12878_SMC1_rep2.6c63bc09f8f8b392ae354182a96a5a68.mugqic.done
)
homer_make_ucsc_file_14_JOB_ID=$(echo "#! /bin/bash 
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
echo "$homer_make_ucsc_file_14_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_15_JOB_ID: homer_make_ucsc_file.GM12878_CDK9_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.GM12878_CDK9_rep2
JOB_DEPENDENCIES=$homer_make_tag_directory_8_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.GM12878_CDK9_rep2.23104d109369ce2209e32a538354412c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.GM12878_CDK9_rep2.23104d109369ce2209e32a538354412c.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 && \
mkdir -p tracks/GM12878_CDK9_rep2 && \
makeUCSCfile \
        tags/GM12878_CDK9_rep2 > tracks/GM12878_CDK9_rep2/GM12878_CDK9_rep2.ucsc.bedGraph && \
        gzip -c tracks/GM12878_CDK9_rep2/GM12878_CDK9_rep2.ucsc.bedGraph > tracks/GM12878_CDK9_rep2/GM12878_CDK9_rep2.ucsc.bedGraph.gz
homer_make_ucsc_file.GM12878_CDK9_rep2.23104d109369ce2209e32a538354412c.mugqic.done
)
homer_make_ucsc_file_15_JOB_ID=$(echo "#! /bin/bash 
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
echo "$homer_make_ucsc_file_15_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_16_JOB_ID: homer_make_ucsc_file_bigWig.GM12878_CDK9_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file_bigWig.GM12878_CDK9_rep2
JOB_DEPENDENCIES=$homer_make_ucsc_file_15_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file_bigWig.GM12878_CDK9_rep2.fbf08e2a104bc146aab82cc165657de4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file_bigWig.GM12878_CDK9_rep2.fbf08e2a104bc146aab82cc165657de4.mugqic.done'
module load mugqic/ucsc/v346 && \
mkdir -p tracks/GM12878_CDK9_rep2/bigWig && \
export TMPDIR=${SLURM_TMPDIR} && \
cat tracks/GM12878_CDK9_rep2/GM12878_CDK9_rep2.ucsc.bedGraph | head -n 1 > tracks/GM12878_CDK9_rep2/GM12878_CDK9_rep2.ucsc.bedGraph.head.tmp && \
cat tracks/GM12878_CDK9_rep2/GM12878_CDK9_rep2.ucsc.bedGraph | awk ' NR > 1 ' | sort  --temporary-directory=${SLURM_TMPDIR} -k1,1 -k2,2n | \
awk '{if($0 !~ /^[A-W]/) print ""$0; else print $0}' | grep -vP "GL|lambda|pUC19" | sed 's/MT/chrM/' | \
awk '{printf "%s\t%d\t%d\t%4.4g\n", $1,$2,$3,$4}' > tracks/GM12878_CDK9_rep2/GM12878_CDK9_rep2.ucsc.bedGraph.body.tmp && \
cat tracks/GM12878_CDK9_rep2/GM12878_CDK9_rep2.ucsc.bedGraph.head.tmp tracks/GM12878_CDK9_rep2/GM12878_CDK9_rep2.ucsc.bedGraph.body.tmp > tracks/GM12878_CDK9_rep2/GM12878_CDK9_rep2.ucsc.bedGraph.sorted && \
rm tracks/GM12878_CDK9_rep2/GM12878_CDK9_rep2.ucsc.bedGraph.head.tmp tracks/GM12878_CDK9_rep2/GM12878_CDK9_rep2.ucsc.bedGraph.body.tmp && \
bedGraphToBigWig \
  tracks/GM12878_CDK9_rep2/GM12878_CDK9_rep2.ucsc.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  tracks/GM12878_CDK9_rep2/bigWig/GM12878_CDK9_rep2.bw
homer_make_ucsc_file_bigWig.GM12878_CDK9_rep2.fbf08e2a104bc146aab82cc165657de4.mugqic.done
)
homer_make_ucsc_file_16_JOB_ID=$(echo "#! /bin/bash 
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
echo "$homer_make_ucsc_file_16_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_17_JOB_ID: homer_make_ucsc_file.GM12878_NIPBL_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.GM12878_NIPBL_rep2
JOB_DEPENDENCIES=$homer_make_tag_directory_9_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.GM12878_NIPBL_rep2.7f05db9526a01d446a9681a19b1b482b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.GM12878_NIPBL_rep2.7f05db9526a01d446a9681a19b1b482b.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 && \
mkdir -p tracks/GM12878_NIPBL_rep2 && \
makeUCSCfile \
        tags/GM12878_NIPBL_rep2 > tracks/GM12878_NIPBL_rep2/GM12878_NIPBL_rep2.ucsc.bedGraph && \
        gzip -c tracks/GM12878_NIPBL_rep2/GM12878_NIPBL_rep2.ucsc.bedGraph > tracks/GM12878_NIPBL_rep2/GM12878_NIPBL_rep2.ucsc.bedGraph.gz
homer_make_ucsc_file.GM12878_NIPBL_rep2.7f05db9526a01d446a9681a19b1b482b.mugqic.done
)
homer_make_ucsc_file_17_JOB_ID=$(echo "#! /bin/bash 
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
echo "$homer_make_ucsc_file_17_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_18_JOB_ID: homer_make_ucsc_file_bigWig.GM12878_NIPBL_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file_bigWig.GM12878_NIPBL_rep2
JOB_DEPENDENCIES=$homer_make_ucsc_file_17_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file_bigWig.GM12878_NIPBL_rep2.7a9eb3285bf230b826632202fa85cd5f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file_bigWig.GM12878_NIPBL_rep2.7a9eb3285bf230b826632202fa85cd5f.mugqic.done'
module load mugqic/ucsc/v346 && \
mkdir -p tracks/GM12878_NIPBL_rep2/bigWig && \
export TMPDIR=${SLURM_TMPDIR} && \
cat tracks/GM12878_NIPBL_rep2/GM12878_NIPBL_rep2.ucsc.bedGraph | head -n 1 > tracks/GM12878_NIPBL_rep2/GM12878_NIPBL_rep2.ucsc.bedGraph.head.tmp && \
cat tracks/GM12878_NIPBL_rep2/GM12878_NIPBL_rep2.ucsc.bedGraph | awk ' NR > 1 ' | sort  --temporary-directory=${SLURM_TMPDIR} -k1,1 -k2,2n | \
awk '{if($0 !~ /^[A-W]/) print ""$0; else print $0}' | grep -vP "GL|lambda|pUC19" | sed 's/MT/chrM/' | \
awk '{printf "%s\t%d\t%d\t%4.4g\n", $1,$2,$3,$4}' > tracks/GM12878_NIPBL_rep2/GM12878_NIPBL_rep2.ucsc.bedGraph.body.tmp && \
cat tracks/GM12878_NIPBL_rep2/GM12878_NIPBL_rep2.ucsc.bedGraph.head.tmp tracks/GM12878_NIPBL_rep2/GM12878_NIPBL_rep2.ucsc.bedGraph.body.tmp > tracks/GM12878_NIPBL_rep2/GM12878_NIPBL_rep2.ucsc.bedGraph.sorted && \
rm tracks/GM12878_NIPBL_rep2/GM12878_NIPBL_rep2.ucsc.bedGraph.head.tmp tracks/GM12878_NIPBL_rep2/GM12878_NIPBL_rep2.ucsc.bedGraph.body.tmp && \
bedGraphToBigWig \
  tracks/GM12878_NIPBL_rep2/GM12878_NIPBL_rep2.ucsc.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  tracks/GM12878_NIPBL_rep2/bigWig/GM12878_NIPBL_rep2.bw
homer_make_ucsc_file_bigWig.GM12878_NIPBL_rep2.7a9eb3285bf230b826632202fa85cd5f.mugqic.done
)
homer_make_ucsc_file_18_JOB_ID=$(echo "#! /bin/bash 
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
echo "$homer_make_ucsc_file_18_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_19_JOB_ID: homer_make_ucsc_file.GM12878_MED1_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.GM12878_MED1_rep2
JOB_DEPENDENCIES=$homer_make_tag_directory_10_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.GM12878_MED1_rep2.b9207e9b85679970a42178cca04d1668.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.GM12878_MED1_rep2.b9207e9b85679970a42178cca04d1668.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 && \
mkdir -p tracks/GM12878_MED1_rep2 && \
makeUCSCfile \
        tags/GM12878_MED1_rep2 > tracks/GM12878_MED1_rep2/GM12878_MED1_rep2.ucsc.bedGraph && \
        gzip -c tracks/GM12878_MED1_rep2/GM12878_MED1_rep2.ucsc.bedGraph > tracks/GM12878_MED1_rep2/GM12878_MED1_rep2.ucsc.bedGraph.gz
homer_make_ucsc_file.GM12878_MED1_rep2.b9207e9b85679970a42178cca04d1668.mugqic.done
)
homer_make_ucsc_file_19_JOB_ID=$(echo "#! /bin/bash 
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
echo "$homer_make_ucsc_file_19_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_20_JOB_ID: homer_make_ucsc_file_bigWig.GM12878_MED1_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file_bigWig.GM12878_MED1_rep2
JOB_DEPENDENCIES=$homer_make_ucsc_file_19_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file_bigWig.GM12878_MED1_rep2.f13f76e3ebdb56ad43097ddba96714cb.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file_bigWig.GM12878_MED1_rep2.f13f76e3ebdb56ad43097ddba96714cb.mugqic.done'
module load mugqic/ucsc/v346 && \
mkdir -p tracks/GM12878_MED1_rep2/bigWig && \
export TMPDIR=${SLURM_TMPDIR} && \
cat tracks/GM12878_MED1_rep2/GM12878_MED1_rep2.ucsc.bedGraph | head -n 1 > tracks/GM12878_MED1_rep2/GM12878_MED1_rep2.ucsc.bedGraph.head.tmp && \
cat tracks/GM12878_MED1_rep2/GM12878_MED1_rep2.ucsc.bedGraph | awk ' NR > 1 ' | sort  --temporary-directory=${SLURM_TMPDIR} -k1,1 -k2,2n | \
awk '{if($0 !~ /^[A-W]/) print ""$0; else print $0}' | grep -vP "GL|lambda|pUC19" | sed 's/MT/chrM/' | \
awk '{printf "%s\t%d\t%d\t%4.4g\n", $1,$2,$3,$4}' > tracks/GM12878_MED1_rep2/GM12878_MED1_rep2.ucsc.bedGraph.body.tmp && \
cat tracks/GM12878_MED1_rep2/GM12878_MED1_rep2.ucsc.bedGraph.head.tmp tracks/GM12878_MED1_rep2/GM12878_MED1_rep2.ucsc.bedGraph.body.tmp > tracks/GM12878_MED1_rep2/GM12878_MED1_rep2.ucsc.bedGraph.sorted && \
rm tracks/GM12878_MED1_rep2/GM12878_MED1_rep2.ucsc.bedGraph.head.tmp tracks/GM12878_MED1_rep2/GM12878_MED1_rep2.ucsc.bedGraph.body.tmp && \
bedGraphToBigWig \
  tracks/GM12878_MED1_rep2/GM12878_MED1_rep2.ucsc.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  tracks/GM12878_MED1_rep2/bigWig/GM12878_MED1_rep2.bw
homer_make_ucsc_file_bigWig.GM12878_MED1_rep2.f13f76e3ebdb56ad43097ddba96714cb.mugqic.done
)
homer_make_ucsc_file_20_JOB_ID=$(echo "#! /bin/bash 
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
echo "$homer_make_ucsc_file_20_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_21_JOB_ID: homer_make_ucsc_file_report
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file_report
JOB_DEPENDENCIES=$homer_make_ucsc_file_19_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file_report.90b1cb89146cadaf805ea66244074321.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file_report.90b1cb89146cadaf805ea66244074321.mugqic.done'
mkdir -p report && \
zip -r report/tracks.zip tracks/*/*.ucsc.bedGraph.gz && \
cp /home/efournie/genpipes/bfx/report/ChipSeq.homer_make_ucsc_file.md report/
homer_make_ucsc_file_report.90b1cb89146cadaf805ea66244074321.mugqic.done
)
homer_make_ucsc_file_21_JOB_ID=$(echo "#! /bin/bash 
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
echo "$homer_make_ucsc_file_21_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# STEP: macs2_callpeak
#-------------------------------------------------------------------------------
STEP=macs2_callpeak
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_1_JOB_ID: macs2_callpeak.NIPBL_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.NIPBL_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_2_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.NIPBL_rep1.4f8e760999a5a36f30de98b46a196bc4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.NIPBL_rep1.4f8e760999a5a36f30de98b46a196bc4.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/NIPBL_rep1 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/GM12878_NIPBL_rep1/GM12878_NIPBL_rep1.sorted.dup.bam \
  --control \
  alignment/GM12878_WCE_rep1/GM12878_WCE_rep1.sorted.dup.bam \
  --name peak_call/NIPBL_rep1/NIPBL_rep1 \
  >& peak_call/NIPBL_rep1/NIPBL_rep1.diag.macs.out
macs2_callpeak.NIPBL_rep1.4f8e760999a5a36f30de98b46a196bc4.mugqic.done
)
macs2_callpeak_1_JOB_ID=$(echo "#! /bin/bash 
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
echo "$macs2_callpeak_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_2_JOB_ID: macs2_callpeak_bigBed.NIPBL_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.NIPBL_rep1
JOB_DEPENDENCIES=$macs2_callpeak_1_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.NIPBL_rep1.9ce8705b7ae8da5514d24c2fd7dd101d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.NIPBL_rep1.9ce8705b7ae8da5514d24c2fd7dd101d.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/NIPBL_rep1/NIPBL_rep1_peaks.narrowPeak > peak_call/NIPBL_rep1/NIPBL_rep1_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/NIPBL_rep1/NIPBL_rep1_peaks.narrowPeak.bed \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/NIPBL_rep1/NIPBL_rep1_peaks.narrowPeak.bb
macs2_callpeak_bigBed.NIPBL_rep1.9ce8705b7ae8da5514d24c2fd7dd101d.mugqic.done
)
macs2_callpeak_2_JOB_ID=$(echo "#! /bin/bash 
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
echo "$macs2_callpeak_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_3_JOB_ID: macs2_callpeak.MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.MED1_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_3_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.MED1_rep1.82484d1210386e01b88763d8f9d0fd1f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.MED1_rep1.82484d1210386e01b88763d8f9d0fd1f.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/MED1_rep1 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/GM12878_MED1_rep1/GM12878_MED1_rep1.sorted.dup.bam \
  --control \
  alignment/GM12878_WCE_rep1/GM12878_WCE_rep1.sorted.dup.bam \
  --name peak_call/MED1_rep1/MED1_rep1 \
  >& peak_call/MED1_rep1/MED1_rep1.diag.macs.out
macs2_callpeak.MED1_rep1.82484d1210386e01b88763d8f9d0fd1f.mugqic.done
)
macs2_callpeak_3_JOB_ID=$(echo "#! /bin/bash 
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
echo "$macs2_callpeak_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_4_JOB_ID: macs2_callpeak_bigBed.MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.MED1_rep1
JOB_DEPENDENCIES=$macs2_callpeak_3_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.MED1_rep1.ac9c7da7d289164b576ff81fc48df095.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.MED1_rep1.ac9c7da7d289164b576ff81fc48df095.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/MED1_rep1/MED1_rep1_peaks.narrowPeak > peak_call/MED1_rep1/MED1_rep1_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/MED1_rep1/MED1_rep1_peaks.narrowPeak.bed \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/MED1_rep1/MED1_rep1_peaks.narrowPeak.bb
macs2_callpeak_bigBed.MED1_rep1.ac9c7da7d289164b576ff81fc48df095.mugqic.done
)
macs2_callpeak_4_JOB_ID=$(echo "#! /bin/bash 
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
echo "$macs2_callpeak_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_5_JOB_ID: macs2_callpeak.SMC1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.SMC1_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_4_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.SMC1_rep1.cd1c6cf306081e789e40adc24c12bb6e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.SMC1_rep1.cd1c6cf306081e789e40adc24c12bb6e.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/SMC1_rep1 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/GM12878_SMC1_rep1/GM12878_SMC1_rep1.sorted.dup.bam \
  --control \
  alignment/GM12878_WCE_rep1/GM12878_WCE_rep1.sorted.dup.bam \
  --name peak_call/SMC1_rep1/SMC1_rep1 \
  >& peak_call/SMC1_rep1/SMC1_rep1.diag.macs.out
macs2_callpeak.SMC1_rep1.cd1c6cf306081e789e40adc24c12bb6e.mugqic.done
)
macs2_callpeak_5_JOB_ID=$(echo "#! /bin/bash 
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
echo "$macs2_callpeak_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_6_JOB_ID: macs2_callpeak_bigBed.SMC1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.SMC1_rep1
JOB_DEPENDENCIES=$macs2_callpeak_5_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.SMC1_rep1.547310d830603267d4c2c50ddc01f216.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.SMC1_rep1.547310d830603267d4c2c50ddc01f216.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/SMC1_rep1/SMC1_rep1_peaks.narrowPeak > peak_call/SMC1_rep1/SMC1_rep1_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/SMC1_rep1/SMC1_rep1_peaks.narrowPeak.bed \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/SMC1_rep1/SMC1_rep1_peaks.narrowPeak.bb
macs2_callpeak_bigBed.SMC1_rep1.547310d830603267d4c2c50ddc01f216.mugqic.done
)
macs2_callpeak_6_JOB_ID=$(echo "#! /bin/bash 
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
echo "$macs2_callpeak_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_7_JOB_ID: macs2_callpeak.CDK9_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.CDK9_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_5_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.CDK9_rep1.d614e22e1656b9000769eafd7cfbac3c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.CDK9_rep1.d614e22e1656b9000769eafd7cfbac3c.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/CDK9_rep1 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/GM12878_CDK9_rep1/GM12878_CDK9_rep1.sorted.dup.bam \
  --control \
  alignment/GM12878_WCE_rep1/GM12878_WCE_rep1.sorted.dup.bam \
  --name peak_call/CDK9_rep1/CDK9_rep1 \
  >& peak_call/CDK9_rep1/CDK9_rep1.diag.macs.out
macs2_callpeak.CDK9_rep1.d614e22e1656b9000769eafd7cfbac3c.mugqic.done
)
macs2_callpeak_7_JOB_ID=$(echo "#! /bin/bash 
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
echo "$macs2_callpeak_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_8_JOB_ID: macs2_callpeak_bigBed.CDK9_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.CDK9_rep1
JOB_DEPENDENCIES=$macs2_callpeak_7_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.CDK9_rep1.c27659afecbde5d3031d9b14e421c361.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.CDK9_rep1.c27659afecbde5d3031d9b14e421c361.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/CDK9_rep1/CDK9_rep1_peaks.narrowPeak > peak_call/CDK9_rep1/CDK9_rep1_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/CDK9_rep1/CDK9_rep1_peaks.narrowPeak.bed \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/CDK9_rep1/CDK9_rep1_peaks.narrowPeak.bb
macs2_callpeak_bigBed.CDK9_rep1.c27659afecbde5d3031d9b14e421c361.mugqic.done
)
macs2_callpeak_8_JOB_ID=$(echo "#! /bin/bash 
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
echo "$macs2_callpeak_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_9_JOB_ID: macs2_callpeak.NIPBL_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.NIPBL_rep2
JOB_DEPENDENCIES=$picard_mark_duplicates_6_JOB_ID:$picard_mark_duplicates_9_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.NIPBL_rep2.42424b7dea0d0d248340da08a6a0be2d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.NIPBL_rep2.42424b7dea0d0d248340da08a6a0be2d.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/NIPBL_rep2 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/GM12878_NIPBL_rep2/GM12878_NIPBL_rep2.sorted.dup.bam \
  --control \
  alignment/GM12878_WCE_rep2/GM12878_WCE_rep2.sorted.dup.bam \
  --name peak_call/NIPBL_rep2/NIPBL_rep2 \
  >& peak_call/NIPBL_rep2/NIPBL_rep2.diag.macs.out
macs2_callpeak.NIPBL_rep2.42424b7dea0d0d248340da08a6a0be2d.mugqic.done
)
macs2_callpeak_9_JOB_ID=$(echo "#! /bin/bash 
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
echo "$macs2_callpeak_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_10_JOB_ID: macs2_callpeak_bigBed.NIPBL_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.NIPBL_rep2
JOB_DEPENDENCIES=$macs2_callpeak_9_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.NIPBL_rep2.fad8c73b685257a7ff06c1fa0f86a4bc.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.NIPBL_rep2.fad8c73b685257a7ff06c1fa0f86a4bc.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/NIPBL_rep2/NIPBL_rep2_peaks.narrowPeak > peak_call/NIPBL_rep2/NIPBL_rep2_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/NIPBL_rep2/NIPBL_rep2_peaks.narrowPeak.bed \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/NIPBL_rep2/NIPBL_rep2_peaks.narrowPeak.bb
macs2_callpeak_bigBed.NIPBL_rep2.fad8c73b685257a7ff06c1fa0f86a4bc.mugqic.done
)
macs2_callpeak_10_JOB_ID=$(echo "#! /bin/bash 
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
echo "$macs2_callpeak_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_11_JOB_ID: macs2_callpeak.MED1_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.MED1_rep2
JOB_DEPENDENCIES=$picard_mark_duplicates_6_JOB_ID:$picard_mark_duplicates_10_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.MED1_rep2.2ccca2c7d401527ef743291bad147110.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.MED1_rep2.2ccca2c7d401527ef743291bad147110.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/MED1_rep2 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/GM12878_MED1_rep2/GM12878_MED1_rep2.sorted.dup.bam \
  --control \
  alignment/GM12878_WCE_rep2/GM12878_WCE_rep2.sorted.dup.bam \
  --name peak_call/MED1_rep2/MED1_rep2 \
  >& peak_call/MED1_rep2/MED1_rep2.diag.macs.out
macs2_callpeak.MED1_rep2.2ccca2c7d401527ef743291bad147110.mugqic.done
)
macs2_callpeak_11_JOB_ID=$(echo "#! /bin/bash 
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
echo "$macs2_callpeak_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_12_JOB_ID: macs2_callpeak_bigBed.MED1_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.MED1_rep2
JOB_DEPENDENCIES=$macs2_callpeak_11_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.MED1_rep2.f4c8a0c6c2a4e13bedd9f92a49b4f153.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.MED1_rep2.f4c8a0c6c2a4e13bedd9f92a49b4f153.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/MED1_rep2/MED1_rep2_peaks.narrowPeak > peak_call/MED1_rep2/MED1_rep2_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/MED1_rep2/MED1_rep2_peaks.narrowPeak.bed \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/MED1_rep2/MED1_rep2_peaks.narrowPeak.bb
macs2_callpeak_bigBed.MED1_rep2.f4c8a0c6c2a4e13bedd9f92a49b4f153.mugqic.done
)
macs2_callpeak_12_JOB_ID=$(echo "#! /bin/bash 
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
echo "$macs2_callpeak_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_13_JOB_ID: macs2_callpeak.SMC1_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.SMC1_rep2
JOB_DEPENDENCIES=$picard_mark_duplicates_6_JOB_ID:$picard_mark_duplicates_7_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.SMC1_rep2.0418656737b898e663e5764a3b82c9ee.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.SMC1_rep2.0418656737b898e663e5764a3b82c9ee.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/SMC1_rep2 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/GM12878_SMC1_rep2/GM12878_SMC1_rep2.sorted.dup.bam \
  --control \
  alignment/GM12878_WCE_rep2/GM12878_WCE_rep2.sorted.dup.bam \
  --name peak_call/SMC1_rep2/SMC1_rep2 \
  >& peak_call/SMC1_rep2/SMC1_rep2.diag.macs.out
macs2_callpeak.SMC1_rep2.0418656737b898e663e5764a3b82c9ee.mugqic.done
)
macs2_callpeak_13_JOB_ID=$(echo "#! /bin/bash 
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
echo "$macs2_callpeak_13_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_14_JOB_ID: macs2_callpeak_bigBed.SMC1_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.SMC1_rep2
JOB_DEPENDENCIES=$macs2_callpeak_13_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.SMC1_rep2.ca8a4309b0986f4495f4f5984308822f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.SMC1_rep2.ca8a4309b0986f4495f4f5984308822f.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/SMC1_rep2/SMC1_rep2_peaks.narrowPeak > peak_call/SMC1_rep2/SMC1_rep2_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/SMC1_rep2/SMC1_rep2_peaks.narrowPeak.bed \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/SMC1_rep2/SMC1_rep2_peaks.narrowPeak.bb
macs2_callpeak_bigBed.SMC1_rep2.ca8a4309b0986f4495f4f5984308822f.mugqic.done
)
macs2_callpeak_14_JOB_ID=$(echo "#! /bin/bash 
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
echo "$macs2_callpeak_14_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_15_JOB_ID: macs2_callpeak.CDK9_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.CDK9_rep2
JOB_DEPENDENCIES=$picard_mark_duplicates_6_JOB_ID:$picard_mark_duplicates_8_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.CDK9_rep2.a89e0a606b55a9a11872b747634098e7.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.CDK9_rep2.a89e0a606b55a9a11872b747634098e7.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/CDK9_rep2 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/GM12878_CDK9_rep2/GM12878_CDK9_rep2.sorted.dup.bam \
  --control \
  alignment/GM12878_WCE_rep2/GM12878_WCE_rep2.sorted.dup.bam \
  --name peak_call/CDK9_rep2/CDK9_rep2 \
  >& peak_call/CDK9_rep2/CDK9_rep2.diag.macs.out
macs2_callpeak.CDK9_rep2.a89e0a606b55a9a11872b747634098e7.mugqic.done
)
macs2_callpeak_15_JOB_ID=$(echo "#! /bin/bash 
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
echo "$macs2_callpeak_15_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_16_JOB_ID: macs2_callpeak_bigBed.CDK9_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.CDK9_rep2
JOB_DEPENDENCIES=$macs2_callpeak_15_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.CDK9_rep2.7264128ca8a16ae478c6d04df7081ccc.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.CDK9_rep2.7264128ca8a16ae478c6d04df7081ccc.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/CDK9_rep2/CDK9_rep2_peaks.narrowPeak > peak_call/CDK9_rep2/CDK9_rep2_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/CDK9_rep2/CDK9_rep2_peaks.narrowPeak.bed \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/CDK9_rep2/CDK9_rep2_peaks.narrowPeak.bb
macs2_callpeak_bigBed.CDK9_rep2.7264128ca8a16ae478c6d04df7081ccc.mugqic.done
)
macs2_callpeak_16_JOB_ID=$(echo "#! /bin/bash 
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
echo "$macs2_callpeak_16_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_17_JOB_ID: macs2_callpeak_report
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_report
JOB_DEPENDENCIES=$macs2_callpeak_1_JOB_ID:$macs2_callpeak_3_JOB_ID:$macs2_callpeak_5_JOB_ID:$macs2_callpeak_7_JOB_ID:$macs2_callpeak_9_JOB_ID:$macs2_callpeak_11_JOB_ID:$macs2_callpeak_13_JOB_ID:$macs2_callpeak_15_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_report.d954ea8cc0fd27fb5113bb5ebff2aac1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_report.d954ea8cc0fd27fb5113bb5ebff2aac1.mugqic.done'
mkdir -p report && \
cp /home/efournie/genpipes/bfx/report/ChipSeq.macs2_callpeak.md report/ && \
for contrast in NIPBL_rep1 MED1_rep1 SMC1_rep1 CDK9_rep1 NIPBL_rep2 MED1_rep2 SMC1_rep2 CDK9_rep2
do
  cp -a --parents peak_call/$contrast/ report/ && \
  echo -e "* [Peak Calls File for Design $contrast](peak_call/$contrast/${contrast}_peaks.xls)" \
  >> report/ChipSeq.macs2_callpeak.md
done
macs2_callpeak_report.d954ea8cc0fd27fb5113bb5ebff2aac1.mugqic.done
)
macs2_callpeak_17_JOB_ID=$(echo "#! /bin/bash 
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
echo "$macs2_callpeak_17_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# Call home with pipeline statistics
#-------------------------------------------------------------------------------
wget "http://mugqic.hpc.mcgill.ca/cgi-bin/pipeline.cgi?hostname=cedar1.cedar.computecanada.ca&ip=206.12.124.2&pipeline=ChipSeq&steps=picard_sam_to_fastq,trimmomatic,merge_trimmomatic_stats,bwa_mem_picard_sort_sam,samtools_view_filter,picard_merge_sam_files,picard_mark_duplicates,metrics,homer_make_tag_directory,qc_metrics,homer_make_ucsc_file,macs2_callpeak&samples=10&AnonymizedList=8e3a32c6841b3667a36f592251a4d299,d21af0b56344f0d1f55c38834c06b9a3,749b4ef76c1e98f10fdd9d3d43eaf70b,d77cb5d5e2b76f57a420535616889fc6,2a8ac303d25ae0d1d4993388cd470647,4320e1eed1a76d7b846a76989532a7cc,be36b6d5a09b2d5ce19dd9766b74acd0,3d9042ca2e59082607fd856ed0c535ab,ed16bbf6d36f654ec401abbc3a8b897b,569c3651629e2af9ac4fb059cfc7cf0d" --quiet --output-document=/dev/null

