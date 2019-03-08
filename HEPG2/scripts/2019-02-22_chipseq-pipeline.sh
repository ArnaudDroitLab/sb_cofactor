#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

#-------------------------------------------------------------------------------
# ChipSeq SlurmScheduler Job Submission Bash script
# Version: 3.1.2-beta
# Created on: 2019-02-22T10:09:10
# Steps:
#   picard_sam_to_fastq: 9 jobs
#   trimmomatic: 13 jobs
#   merge_trimmomatic_stats: 1 job
#   bwa_mem_picard_sort_sam: 14 jobs
#   samtools_view_filter: 14 jobs
#   picard_merge_sam_files: 13 jobs
#   picard_mark_duplicates: 14 jobs
#   metrics: 2 jobs
#   homer_make_tag_directory: 13 jobs
#   qc_metrics: 1 job
#   macs2_callpeak: 19 jobs
#   TOTAL: 113 jobs
#-------------------------------------------------------------------------------

OUTPUT_DIR=/scratch/chris11/sb_cofactor/HEPG2/output/chip-pipeline-GRCh38
JOB_OUTPUT_DIR=$OUTPUT_DIR/job_output
TIMESTAMP=`date +%FT%H.%M.%S`
JOB_LIST=$JOB_OUTPUT_DIR/ChipSeq_job_list_$TIMESTAMP
export CONFIG_FILES="/cvmfs/soft.mugqic/CentOS6/software/genpipes/genpipes-3.1.2/pipelines/chipseq/chipseq.base.ini,/cvmfs/soft.mugqic/CentOS6/software/genpipes/genpipes-3.1.2/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini,input/cedar_assembly_bugfix.txt"
mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR


#-------------------------------------------------------------------------------
# STEP: picard_sam_to_fastq
#-------------------------------------------------------------------------------
STEP=picard_sam_to_fastq
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_1_JOB_ID: picard_sam_to_fastq.HI.1308.004.Index_18.HEPG2_WCE_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.1308.004.Index_18.HEPG2_WCE_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.1308.004.Index_18.HEPG2_WCE_rep1.bam.22abddfa29e97ea09590065dcb75b494.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.1308.004.Index_18.HEPG2_WCE_rep1.bam.22abddfa29e97ea09590065dcb75b494.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/scratch/chris11/sb_cofactor/HEPG2/raw/chip-seq/HI.1308.004.Index_18.HEPG2_WCE_rep1.bam \
 FASTQ=/scratch/chris11/sb_cofactor/HEPG2/raw/chip-seq/HI.1308.004.Index_18.HEPG2_WCE_rep1.single.fastq.gz
picard_sam_to_fastq.HI.1308.004.Index_18.HEPG2_WCE_rep1.bam.22abddfa29e97ea09590065dcb75b494.mugqic.done
)
picard_sam_to_fastq_1_JOB_ID=$(echo "#! /bin/bash
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=12G -N 1 -n 3 | grep "[0-9]" | cut -d\  -f4)
echo "$picard_sam_to_fastq_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_2_JOB_ID: picard_sam_to_fastq.HI.1308.004.Index_19.HEPG2_NIPBL_ChIP67_1G_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.1308.004.Index_19.HEPG2_NIPBL_ChIP67_1G_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.1308.004.Index_19.HEPG2_NIPBL_ChIP67_1G_rep1.bam.88181392c82bb2be3ae39d16992246ba.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.1308.004.Index_19.HEPG2_NIPBL_ChIP67_1G_rep1.bam.88181392c82bb2be3ae39d16992246ba.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/scratch/chris11/sb_cofactor/HEPG2/raw/chip-seq/HI.1308.004.Index_19.HEPG2_NIPBL_ChIP67_1G_rep1.bam \
 FASTQ=/scratch/chris11/sb_cofactor/HEPG2/raw/chip-seq/HI.1308.004.Index_19.HEPG2_NIPBL_ChIP67_1G_rep1.single.fastq.gz
picard_sam_to_fastq.HI.1308.004.Index_19.HEPG2_NIPBL_ChIP67_1G_rep1.bam.88181392c82bb2be3ae39d16992246ba.mugqic.done
)
picard_sam_to_fastq_2_JOB_ID=$(echo "#! /bin/bash
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=12G -N 1 -n 3 | grep "[0-9]" | cut -d\  -f4)
echo "$picard_sam_to_fastq_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_3_JOB_ID: picard_sam_to_fastq.HI.1308.004.Index_22.HEPG2_Cdk9_ChIP67_5G_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.1308.004.Index_22.HEPG2_Cdk9_ChIP67_5G_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.1308.004.Index_22.HEPG2_Cdk9_ChIP67_5G_rep1.bam.ae0a29a1b6886594c59d31ca4c89259d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.1308.004.Index_22.HEPG2_Cdk9_ChIP67_5G_rep1.bam.ae0a29a1b6886594c59d31ca4c89259d.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/scratch/chris11/sb_cofactor/HEPG2/raw/chip-seq/HI.1308.004.Index_22.HEPG2_Cdk9_ChIP67_5G_rep1.bam \
 FASTQ=/scratch/chris11/sb_cofactor/HEPG2/raw/chip-seq/HI.1308.004.Index_22.HEPG2_Cdk9_ChIP67_5G_rep1.single.fastq.gz
picard_sam_to_fastq.HI.1308.004.Index_22.HEPG2_Cdk9_ChIP67_5G_rep1.bam.ae0a29a1b6886594c59d31ca4c89259d.mugqic.done
)
picard_sam_to_fastq_3_JOB_ID=$(echo "#! /bin/bash
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=12G -N 1 -n 3 | grep "[0-9]" | cut -d\  -f4)
echo "$picard_sam_to_fastq_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_4_JOB_ID: picard_sam_to_fastq.HI.1308.004.Index_20.HEPG2_Med1_ChIP67_2G_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.1308.004.Index_20.HEPG2_Med1_ChIP67_2G_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.1308.004.Index_20.HEPG2_Med1_ChIP67_2G_rep1.bam.e1c7a6124371e2529e1ee1f81fb6f130.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.1308.004.Index_20.HEPG2_Med1_ChIP67_2G_rep1.bam.e1c7a6124371e2529e1ee1f81fb6f130.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/scratch/chris11/sb_cofactor/HEPG2/raw/chip-seq/HI.1308.004.Index_20.HEPG2_Med1_ChIP67_2G_rep1.bam \
 FASTQ=/scratch/chris11/sb_cofactor/HEPG2/raw/chip-seq/HI.1308.004.Index_20.HEPG2_Med1_ChIP67_2G_rep1.single.fastq.gz
picard_sam_to_fastq.HI.1308.004.Index_20.HEPG2_Med1_ChIP67_2G_rep1.bam.e1c7a6124371e2529e1ee1f81fb6f130.mugqic.done
)
picard_sam_to_fastq_4_JOB_ID=$(echo "#! /bin/bash
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=12G -N 1 -n 3 | grep "[0-9]" | cut -d\  -f4)
echo "$picard_sam_to_fastq_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_5_JOB_ID: picard_sam_to_fastq.HI.1308.004.Index_21.HEPG2_Smc1_ChIP67_3G_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.1308.004.Index_21.HEPG2_Smc1_ChIP67_3G_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.1308.004.Index_21.HEPG2_Smc1_ChIP67_3G_rep1.bam.5ce0840d6bd1a787fdced6e5590f21eb.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.1308.004.Index_21.HEPG2_Smc1_ChIP67_3G_rep1.bam.5ce0840d6bd1a787fdced6e5590f21eb.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/scratch/chris11/sb_cofactor/HEPG2/raw/chip-seq/HI.1308.004.Index_21.HEPG2_Smc1_ChIP67_3G_rep1.bam \
 FASTQ=/scratch/chris11/sb_cofactor/HEPG2/raw/chip-seq/HI.1308.004.Index_21.HEPG2_Smc1_ChIP67_3G_rep1.single.fastq.gz
picard_sam_to_fastq.HI.1308.004.Index_21.HEPG2_Smc1_ChIP67_3G_rep1.bam.5ce0840d6bd1a787fdced6e5590f21eb.mugqic.done
)
picard_sam_to_fastq_5_JOB_ID=$(echo "#! /bin/bash
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=12G -N 1 -n 3 | grep "[0-9]" | cut -d\  -f4)
echo "$picard_sam_to_fastq_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_6_JOB_ID: picard_sam_to_fastq.HI.2449.006.Index_8.HEPG2_ctrl_WCE_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.2449.006.Index_8.HEPG2_ctrl_WCE_MF6_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.2449.006.Index_8.HEPG2_ctrl_WCE_MF6_rep1.bam.0796a582e8c01faa74d76273c9d3298f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.2449.006.Index_8.HEPG2_ctrl_WCE_MF6_rep1.bam.0796a582e8c01faa74d76273c9d3298f.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/scratch/chris11/sb_cofactor/HEPG2/raw/chip-seq/HI.2449.006.Index_8.HEPG2_ctrl_WCE_MF6_rep1.bam \
 FASTQ=/scratch/chris11/sb_cofactor/HEPG2/raw/chip-seq/HI.2449.006.Index_8.HEPG2_ctrl_WCE_MF6_rep1.single.fastq.gz
picard_sam_to_fastq.HI.2449.006.Index_8.HEPG2_ctrl_WCE_MF6_rep1.bam.0796a582e8c01faa74d76273c9d3298f.mugqic.done
)
picard_sam_to_fastq_6_JOB_ID=$(echo "#! /bin/bash
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=12G -N 1 -n 3 | grep "[0-9]" | cut -d\  -f4)
echo "$picard_sam_to_fastq_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_7_JOB_ID: picard_sam_to_fastq.HI.2449.006.Index_9.HEPG2_ctrl_SMC1A_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.2449.006.Index_9.HEPG2_ctrl_SMC1A_MF6_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.2449.006.Index_9.HEPG2_ctrl_SMC1A_MF6_rep1.bam.f04ccce7fafa0664a5260eaf14e11d81.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.2449.006.Index_9.HEPG2_ctrl_SMC1A_MF6_rep1.bam.f04ccce7fafa0664a5260eaf14e11d81.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/scratch/chris11/sb_cofactor/HEPG2/raw/chip-seq/HI.2449.006.Index_9.HEPG2_ctrl_SMC1A_MF6_rep1.bam \
 FASTQ=/scratch/chris11/sb_cofactor/HEPG2/raw/chip-seq/HI.2449.006.Index_9.HEPG2_ctrl_SMC1A_MF6_rep1.single.fastq.gz
picard_sam_to_fastq.HI.2449.006.Index_9.HEPG2_ctrl_SMC1A_MF6_rep1.bam.f04ccce7fafa0664a5260eaf14e11d81.mugqic.done
)
picard_sam_to_fastq_7_JOB_ID=$(echo "#! /bin/bash
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=12G -N 1 -n 3 | grep "[0-9]" | cut -d\  -f4)
echo "$picard_sam_to_fastq_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_8_JOB_ID: picard_sam_to_fastq.HI.2449.006.Index_1.HEPG2_ctrl_NIPBL_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.2449.006.Index_1.HEPG2_ctrl_NIPBL_MF6_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.2449.006.Index_1.HEPG2_ctrl_NIPBL_MF6_rep1.bam.8168cfb170c9f236e6da2bfc45bfc29c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.2449.006.Index_1.HEPG2_ctrl_NIPBL_MF6_rep1.bam.8168cfb170c9f236e6da2bfc45bfc29c.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/scratch/chris11/sb_cofactor/HEPG2/raw/chip-seq/HI.2449.006.Index_1.HEPG2_ctrl_NIPBL_MF6_rep1.bam \
 FASTQ=/scratch/chris11/sb_cofactor/HEPG2/raw/chip-seq/HI.2449.006.Index_1.HEPG2_ctrl_NIPBL_MF6_rep1.single.fastq.gz
picard_sam_to_fastq.HI.2449.006.Index_1.HEPG2_ctrl_NIPBL_MF6_rep1.bam.8168cfb170c9f236e6da2bfc45bfc29c.mugqic.done
)
picard_sam_to_fastq_8_JOB_ID=$(echo "#! /bin/bash
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=12G -N 1 -n 3 | grep "[0-9]" | cut -d\  -f4)
echo "$picard_sam_to_fastq_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_9_JOB_ID: picard_sam_to_fastq.HI.2449.006.Index_3.HEPG2_ctrl_MED1_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.2449.006.Index_3.HEPG2_ctrl_MED1_MF6_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.2449.006.Index_3.HEPG2_ctrl_MED1_MF6_rep1.bam.cb65d140ef02ccf572c268e6bfa4c958.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.2449.006.Index_3.HEPG2_ctrl_MED1_MF6_rep1.bam.cb65d140ef02ccf572c268e6bfa4c958.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/scratch/chris11/sb_cofactor/HEPG2/raw/chip-seq/HI.2449.006.Index_3.HEPG2_ctrl_MED1_MF6_rep1.bam \
 FASTQ=/scratch/chris11/sb_cofactor/HEPG2/raw/chip-seq/HI.2449.006.Index_3.HEPG2_ctrl_MED1_MF6_rep1.single.fastq.gz
picard_sam_to_fastq.HI.2449.006.Index_3.HEPG2_ctrl_MED1_MF6_rep1.bam.cb65d140ef02ccf572c268e6bfa4c958.mugqic.done
)
picard_sam_to_fastq_9_JOB_ID=$(echo "#! /bin/bash
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=12G -N 1 -n 3 | grep "[0-9]" | cut -d\  -f4)
echo "$picard_sam_to_fastq_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# STEP: trimmomatic
#-------------------------------------------------------------------------------
STEP=trimmomatic
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: trimmomatic_1_JOB_ID: trimmomatic.HI.1308.004.Index_18.HEPG2_WCE_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.1308.004.Index_18.HEPG2_WCE_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_1_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.1308.004.Index_18.HEPG2_WCE_rep1.bam.7f3c280af93db39019581d9de64279f0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.1308.004.Index_18.HEPG2_WCE_rep1.bam.7f3c280af93db39019581d9de64279f0.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/HEPG2_WCE_rep1 && \
`cat > trim/HEPG2_WCE_rep1/HI.1308.004.Index_18.HEPG2_WCE_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /scratch/chris11/sb_cofactor/HEPG2/raw/chip-seq/HI.1308.004.Index_18.HEPG2_WCE_rep1.single.fastq.gz \
  trim/HEPG2_WCE_rep1/HI.1308.004.Index_18.HEPG2_WCE_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/HEPG2_WCE_rep1/HI.1308.004.Index_18.HEPG2_WCE_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/HEPG2_WCE_rep1/HI.1308.004.Index_18.HEPG2_WCE_rep1.bam.trim.log
trimmomatic.HI.1308.004.Index_18.HEPG2_WCE_rep1.bam.7f3c280af93db39019581d9de64279f0.mugqic.done
)
trimmomatic_1_JOB_ID=$(echo "#! /bin/bash
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=24G -N 1 -n 6 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$trimmomatic_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: trimmomatic_2_JOB_ID: trimmomatic.HI.1308.004.Index_19.HEPG2_NIPBL_ChIP67_1G_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.1308.004.Index_19.HEPG2_NIPBL_ChIP67_1G_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_2_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.1308.004.Index_19.HEPG2_NIPBL_ChIP67_1G_rep1.bam.3d3d1800d038675c929b82bd2618930a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.1308.004.Index_19.HEPG2_NIPBL_ChIP67_1G_rep1.bam.3d3d1800d038675c929b82bd2618930a.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/HEPG2_NIPBL_rep1 && \
`cat > trim/HEPG2_NIPBL_rep1/HI.1308.004.Index_19.HEPG2_NIPBL_ChIP67_1G_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /scratch/chris11/sb_cofactor/HEPG2/raw/chip-seq/HI.1308.004.Index_19.HEPG2_NIPBL_ChIP67_1G_rep1.single.fastq.gz \
  trim/HEPG2_NIPBL_rep1/HI.1308.004.Index_19.HEPG2_NIPBL_ChIP67_1G_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/HEPG2_NIPBL_rep1/HI.1308.004.Index_19.HEPG2_NIPBL_ChIP67_1G_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/HEPG2_NIPBL_rep1/HI.1308.004.Index_19.HEPG2_NIPBL_ChIP67_1G_rep1.bam.trim.log
trimmomatic.HI.1308.004.Index_19.HEPG2_NIPBL_ChIP67_1G_rep1.bam.3d3d1800d038675c929b82bd2618930a.mugqic.done
)
trimmomatic_2_JOB_ID=$(echo "#! /bin/bash
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=24G -N 1 -n 6 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$trimmomatic_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: trimmomatic_3_JOB_ID: trimmomatic.HI.1308.004.Index_22.HEPG2_Cdk9_ChIP67_5G_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.1308.004.Index_22.HEPG2_Cdk9_ChIP67_5G_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_3_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.1308.004.Index_22.HEPG2_Cdk9_ChIP67_5G_rep1.bam.3cd71626b943b3f03d40c8c39cdc501a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.1308.004.Index_22.HEPG2_Cdk9_ChIP67_5G_rep1.bam.3cd71626b943b3f03d40c8c39cdc501a.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/HEPG2_CDK9_rep1 && \
`cat > trim/HEPG2_CDK9_rep1/HI.1308.004.Index_22.HEPG2_Cdk9_ChIP67_5G_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /scratch/chris11/sb_cofactor/HEPG2/raw/chip-seq/HI.1308.004.Index_22.HEPG2_Cdk9_ChIP67_5G_rep1.single.fastq.gz \
  trim/HEPG2_CDK9_rep1/HI.1308.004.Index_22.HEPG2_Cdk9_ChIP67_5G_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/HEPG2_CDK9_rep1/HI.1308.004.Index_22.HEPG2_Cdk9_ChIP67_5G_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/HEPG2_CDK9_rep1/HI.1308.004.Index_22.HEPG2_Cdk9_ChIP67_5G_rep1.bam.trim.log
trimmomatic.HI.1308.004.Index_22.HEPG2_Cdk9_ChIP67_5G_rep1.bam.3cd71626b943b3f03d40c8c39cdc501a.mugqic.done
)
trimmomatic_3_JOB_ID=$(echo "#! /bin/bash
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=24G -N 1 -n 6 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$trimmomatic_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: trimmomatic_4_JOB_ID: trimmomatic.HI.1308.004.Index_20.HEPG2_Med1_ChIP67_2G_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.1308.004.Index_20.HEPG2_Med1_ChIP67_2G_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_4_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.1308.004.Index_20.HEPG2_Med1_ChIP67_2G_rep1.bam.f7b9e1f9ebba763b58e83392ae2f9308.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.1308.004.Index_20.HEPG2_Med1_ChIP67_2G_rep1.bam.f7b9e1f9ebba763b58e83392ae2f9308.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/HEPG2_MED1_rep1 && \
`cat > trim/HEPG2_MED1_rep1/HI.1308.004.Index_20.HEPG2_Med1_ChIP67_2G_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /scratch/chris11/sb_cofactor/HEPG2/raw/chip-seq/HI.1308.004.Index_20.HEPG2_Med1_ChIP67_2G_rep1.single.fastq.gz \
  trim/HEPG2_MED1_rep1/HI.1308.004.Index_20.HEPG2_Med1_ChIP67_2G_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/HEPG2_MED1_rep1/HI.1308.004.Index_20.HEPG2_Med1_ChIP67_2G_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/HEPG2_MED1_rep1/HI.1308.004.Index_20.HEPG2_Med1_ChIP67_2G_rep1.bam.trim.log
trimmomatic.HI.1308.004.Index_20.HEPG2_Med1_ChIP67_2G_rep1.bam.f7b9e1f9ebba763b58e83392ae2f9308.mugqic.done
)
trimmomatic_4_JOB_ID=$(echo "#! /bin/bash
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=24G -N 1 -n 6 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$trimmomatic_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: trimmomatic_5_JOB_ID: trimmomatic.HI.1308.004.Index_21.HEPG2_Smc1_ChIP67_3G_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.1308.004.Index_21.HEPG2_Smc1_ChIP67_3G_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_5_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.1308.004.Index_21.HEPG2_Smc1_ChIP67_3G_rep1.bam.59590a6edaf8479d3cdb2b34b77b04df.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.1308.004.Index_21.HEPG2_Smc1_ChIP67_3G_rep1.bam.59590a6edaf8479d3cdb2b34b77b04df.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/HEPG2_SMC1A_rep1 && \
`cat > trim/HEPG2_SMC1A_rep1/HI.1308.004.Index_21.HEPG2_Smc1_ChIP67_3G_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /scratch/chris11/sb_cofactor/HEPG2/raw/chip-seq/HI.1308.004.Index_21.HEPG2_Smc1_ChIP67_3G_rep1.single.fastq.gz \
  trim/HEPG2_SMC1A_rep1/HI.1308.004.Index_21.HEPG2_Smc1_ChIP67_3G_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/HEPG2_SMC1A_rep1/HI.1308.004.Index_21.HEPG2_Smc1_ChIP67_3G_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/HEPG2_SMC1A_rep1/HI.1308.004.Index_21.HEPG2_Smc1_ChIP67_3G_rep1.bam.trim.log
trimmomatic.HI.1308.004.Index_21.HEPG2_Smc1_ChIP67_3G_rep1.bam.59590a6edaf8479d3cdb2b34b77b04df.mugqic.done
)
trimmomatic_5_JOB_ID=$(echo "#! /bin/bash
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=24G -N 1 -n 6 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$trimmomatic_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: trimmomatic_6_JOB_ID: trimmomatic.HI.2449.006.Index_8.HEPG2_ctrl_WCE_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.2449.006.Index_8.HEPG2_ctrl_WCE_MF6_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_6_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.2449.006.Index_8.HEPG2_ctrl_WCE_MF6_rep1.bam.466bf729c089235efa6ee08b7577d267.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.2449.006.Index_8.HEPG2_ctrl_WCE_MF6_rep1.bam.466bf729c089235efa6ee08b7577d267.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/HEPG2_WCE_rep2 && \
`cat > trim/HEPG2_WCE_rep2/HI.2449.006.Index_8.HEPG2_ctrl_WCE_MF6_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /scratch/chris11/sb_cofactor/HEPG2/raw/chip-seq/HI.2449.006.Index_8.HEPG2_ctrl_WCE_MF6_rep1.single.fastq.gz \
  trim/HEPG2_WCE_rep2/HI.2449.006.Index_8.HEPG2_ctrl_WCE_MF6_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/HEPG2_WCE_rep2/HI.2449.006.Index_8.HEPG2_ctrl_WCE_MF6_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/HEPG2_WCE_rep2/HI.2449.006.Index_8.HEPG2_ctrl_WCE_MF6_rep1.bam.trim.log
trimmomatic.HI.2449.006.Index_8.HEPG2_ctrl_WCE_MF6_rep1.bam.466bf729c089235efa6ee08b7577d267.mugqic.done
)
trimmomatic_6_JOB_ID=$(echo "#! /bin/bash
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=24G -N 1 -n 6 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$trimmomatic_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: trimmomatic_7_JOB_ID: trimmomatic.HI.2449.006.Index_9.HEPG2_ctrl_SMC1A_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.2449.006.Index_9.HEPG2_ctrl_SMC1A_MF6_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_7_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.2449.006.Index_9.HEPG2_ctrl_SMC1A_MF6_rep1.bam.98682d96a0c3f02b2a137a622a9fb70c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.2449.006.Index_9.HEPG2_ctrl_SMC1A_MF6_rep1.bam.98682d96a0c3f02b2a137a622a9fb70c.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/HEPG2_SMC1A_rep2 && \
`cat > trim/HEPG2_SMC1A_rep2/HI.2449.006.Index_9.HEPG2_ctrl_SMC1A_MF6_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /scratch/chris11/sb_cofactor/HEPG2/raw/chip-seq/HI.2449.006.Index_9.HEPG2_ctrl_SMC1A_MF6_rep1.single.fastq.gz \
  trim/HEPG2_SMC1A_rep2/HI.2449.006.Index_9.HEPG2_ctrl_SMC1A_MF6_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/HEPG2_SMC1A_rep2/HI.2449.006.Index_9.HEPG2_ctrl_SMC1A_MF6_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/HEPG2_SMC1A_rep2/HI.2449.006.Index_9.HEPG2_ctrl_SMC1A_MF6_rep1.bam.trim.log
trimmomatic.HI.2449.006.Index_9.HEPG2_ctrl_SMC1A_MF6_rep1.bam.98682d96a0c3f02b2a137a622a9fb70c.mugqic.done
)
trimmomatic_7_JOB_ID=$(echo "#! /bin/bash
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=24G -N 1 -n 6 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$trimmomatic_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: trimmomatic_8_JOB_ID: trimmomatic.HI.2449.006.Index_1.HEPG2_ctrl_NIPBL_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.2449.006.Index_1.HEPG2_ctrl_NIPBL_MF6_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_8_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.2449.006.Index_1.HEPG2_ctrl_NIPBL_MF6_rep1.bam.8e5353537a1712ee82a8a1ed5df62f3c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.2449.006.Index_1.HEPG2_ctrl_NIPBL_MF6_rep1.bam.8e5353537a1712ee82a8a1ed5df62f3c.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/HEPG2_NIPBL_rep2 && \
`cat > trim/HEPG2_NIPBL_rep2/HI.2449.006.Index_1.HEPG2_ctrl_NIPBL_MF6_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /scratch/chris11/sb_cofactor/HEPG2/raw/chip-seq/HI.2449.006.Index_1.HEPG2_ctrl_NIPBL_MF6_rep1.single.fastq.gz \
  trim/HEPG2_NIPBL_rep2/HI.2449.006.Index_1.HEPG2_ctrl_NIPBL_MF6_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/HEPG2_NIPBL_rep2/HI.2449.006.Index_1.HEPG2_ctrl_NIPBL_MF6_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/HEPG2_NIPBL_rep2/HI.2449.006.Index_1.HEPG2_ctrl_NIPBL_MF6_rep1.bam.trim.log
trimmomatic.HI.2449.006.Index_1.HEPG2_ctrl_NIPBL_MF6_rep1.bam.8e5353537a1712ee82a8a1ed5df62f3c.mugqic.done
)
trimmomatic_8_JOB_ID=$(echo "#! /bin/bash
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=24G -N 1 -n 6 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$trimmomatic_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: trimmomatic_9_JOB_ID: trimmomatic.HI.2449.006.Index_3.HEPG2_ctrl_MED1_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.2449.006.Index_3.HEPG2_ctrl_MED1_MF6_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_9_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.2449.006.Index_3.HEPG2_ctrl_MED1_MF6_rep1.bam.1ebd26f8a8fc2ee88500aa4c5981d6d8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.2449.006.Index_3.HEPG2_ctrl_MED1_MF6_rep1.bam.1ebd26f8a8fc2ee88500aa4c5981d6d8.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/HEPG2_MED1_rep2 && \
`cat > trim/HEPG2_MED1_rep2/HI.2449.006.Index_3.HEPG2_ctrl_MED1_MF6_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /scratch/chris11/sb_cofactor/HEPG2/raw/chip-seq/HI.2449.006.Index_3.HEPG2_ctrl_MED1_MF6_rep1.single.fastq.gz \
  trim/HEPG2_MED1_rep2/HI.2449.006.Index_3.HEPG2_ctrl_MED1_MF6_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/HEPG2_MED1_rep2/HI.2449.006.Index_3.HEPG2_ctrl_MED1_MF6_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/HEPG2_MED1_rep2/HI.2449.006.Index_3.HEPG2_ctrl_MED1_MF6_rep1.bam.trim.log
trimmomatic.HI.2449.006.Index_3.HEPG2_ctrl_MED1_MF6_rep1.bam.1ebd26f8a8fc2ee88500aa4c5981d6d8.mugqic.done
)
trimmomatic_9_JOB_ID=$(echo "#! /bin/bash
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=24G -N 1 -n 6 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$trimmomatic_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: trimmomatic_10_JOB_ID: trimmomatic.ENCFF953UBF.fastq
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.ENCFF953UBF.fastq
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.ENCFF953UBF.fastq.d0ffccb1a4a51b69bf4813201f821cdf.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.ENCFF953UBF.fastq.d0ffccb1a4a51b69bf4813201f821cdf.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/HEPG2_BRD4_rep1 && \
`cat > trim/HEPG2_BRD4_rep1/ENCFF953UBF.fastq.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /scratch/chris11/sb_cofactor/HEPG2/raw/chip-seq/ENCFF953UBF.fastq \
  trim/HEPG2_BRD4_rep1/ENCFF953UBF.fastq.trim.single.fastq.gz \
  ILLUMINACLIP:trim/HEPG2_BRD4_rep1/ENCFF953UBF.fastq.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/HEPG2_BRD4_rep1/ENCFF953UBF.fastq.trim.log
trimmomatic.ENCFF953UBF.fastq.d0ffccb1a4a51b69bf4813201f821cdf.mugqic.done
)
trimmomatic_10_JOB_ID=$(echo "#! /bin/bash
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=24G -N 1 -n 6 | grep "[0-9]" | cut -d\  -f4)
echo "$trimmomatic_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: trimmomatic_11_JOB_ID: trimmomatic.ENCFF946FEX.fastq
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.ENCFF946FEX.fastq
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.ENCFF946FEX.fastq.935f192d8ba2fc599c4ec83d1c18c1f1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.ENCFF946FEX.fastq.935f192d8ba2fc599c4ec83d1c18c1f1.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/HEPG2_BRD4_rep2 && \
`cat > trim/HEPG2_BRD4_rep2/ENCFF946FEX.fastq.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /scratch/chris11/sb_cofactor/HEPG2/raw/chip-seq/ENCFF946FEX.fastq \
  trim/HEPG2_BRD4_rep2/ENCFF946FEX.fastq.trim.single.fastq.gz \
  ILLUMINACLIP:trim/HEPG2_BRD4_rep2/ENCFF946FEX.fastq.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/HEPG2_BRD4_rep2/ENCFF946FEX.fastq.trim.log
trimmomatic.ENCFF946FEX.fastq.935f192d8ba2fc599c4ec83d1c18c1f1.mugqic.done
)
trimmomatic_11_JOB_ID=$(echo "#! /bin/bash
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=24G -N 1 -n 6 | grep "[0-9]" | cut -d\  -f4)
echo "$trimmomatic_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: trimmomatic_12_JOB_ID: trimmomatic.ENCFF474ZPM.fastq
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.ENCFF474ZPM.fastq
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.ENCFF474ZPM.fastq.b27c0f071d93fde834de1beeebc9bdcb.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.ENCFF474ZPM.fastq.b27c0f071d93fde834de1beeebc9bdcb.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/HEPG2_BRD4_WCE_rep1 && \
`cat > trim/HEPG2_BRD4_WCE_rep1/ENCFF474ZPM.fastq.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /scratch/chris11/sb_cofactor/HEPG2/raw/chip-seq/ENCFF474ZPM.fastq \
  trim/HEPG2_BRD4_WCE_rep1/ENCFF474ZPM.fastq.trim.single.fastq.gz \
  ILLUMINACLIP:trim/HEPG2_BRD4_WCE_rep1/ENCFF474ZPM.fastq.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/HEPG2_BRD4_WCE_rep1/ENCFF474ZPM.fastq.trim.log
trimmomatic.ENCFF474ZPM.fastq.b27c0f071d93fde834de1beeebc9bdcb.mugqic.done
)
trimmomatic_12_JOB_ID=$(echo "#! /bin/bash
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=24G -N 1 -n 6 | grep "[0-9]" | cut -d\  -f4)
echo "$trimmomatic_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: trimmomatic_13_JOB_ID: trimmomatic.ENCFF653JWL.fastq
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.ENCFF653JWL.fastq
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.ENCFF653JWL.fastq.a2dfa35c3cffebea4755f1c7ee002f94.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.ENCFF653JWL.fastq.a2dfa35c3cffebea4755f1c7ee002f94.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/HEPG2_BRD4_WCE_rep2 && \
`cat > trim/HEPG2_BRD4_WCE_rep2/ENCFF653JWL.fastq.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /scratch/chris11/sb_cofactor/HEPG2/raw/chip-seq/ENCFF653JWL.fastq \
  trim/HEPG2_BRD4_WCE_rep2/ENCFF653JWL.fastq.trim.single.fastq.gz \
  ILLUMINACLIP:trim/HEPG2_BRD4_WCE_rep2/ENCFF653JWL.fastq.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/HEPG2_BRD4_WCE_rep2/ENCFF653JWL.fastq.trim.log
trimmomatic.ENCFF653JWL.fastq.a2dfa35c3cffebea4755f1c7ee002f94.mugqic.done
)
trimmomatic_13_JOB_ID=$(echo "#! /bin/bash
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=24G -N 1 -n 6 | grep "[0-9]" | cut -d\  -f4)
echo "$trimmomatic_13_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# STEP: merge_trimmomatic_stats
#-------------------------------------------------------------------------------
STEP=merge_trimmomatic_stats
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: merge_trimmomatic_stats_1_JOB_ID: merge_trimmomatic_stats
#-------------------------------------------------------------------------------
JOB_NAME=merge_trimmomatic_stats
JOB_DEPENDENCIES=$trimmomatic_1_JOB_ID:$trimmomatic_2_JOB_ID:$trimmomatic_3_JOB_ID:$trimmomatic_4_JOB_ID:$trimmomatic_5_JOB_ID:$trimmomatic_6_JOB_ID:$trimmomatic_7_JOB_ID:$trimmomatic_8_JOB_ID:$trimmomatic_9_JOB_ID:$trimmomatic_10_JOB_ID:$trimmomatic_11_JOB_ID:$trimmomatic_12_JOB_ID:$trimmomatic_13_JOB_ID
JOB_DONE=job_output/merge_trimmomatic_stats/merge_trimmomatic_stats.b08357e409aca59c54e3b0e2beba8ff3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'merge_trimmomatic_stats.b08357e409aca59c54e3b0e2beba8ff3.mugqic.done'
module load mugqic/pandoc/1.15.2 && \
mkdir -p metrics && \
echo 'Sample	Readset	Raw Single Reads #	Surviving Single Reads #	Surviving Single Reads %' > metrics/trimReadsetTable.tsv && \
grep ^Input trim/HEPG2_WCE_rep1/HI.1308.004.Index_18.HEPG2_WCE_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/HEPG2_WCE_rep1	HI.1308.004.Index_18.HEPG2_WCE_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/HEPG2_NIPBL_rep1/HI.1308.004.Index_19.HEPG2_NIPBL_ChIP67_1G_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/HEPG2_NIPBL_rep1	HI.1308.004.Index_19.HEPG2_NIPBL_ChIP67_1G_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/HEPG2_CDK9_rep1/HI.1308.004.Index_22.HEPG2_Cdk9_ChIP67_5G_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/HEPG2_CDK9_rep1	HI.1308.004.Index_22.HEPG2_Cdk9_ChIP67_5G_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/HEPG2_MED1_rep1/HI.1308.004.Index_20.HEPG2_Med1_ChIP67_2G_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/HEPG2_MED1_rep1	HI.1308.004.Index_20.HEPG2_Med1_ChIP67_2G_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/HEPG2_SMC1A_rep1/HI.1308.004.Index_21.HEPG2_Smc1_ChIP67_3G_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/HEPG2_SMC1A_rep1	HI.1308.004.Index_21.HEPG2_Smc1_ChIP67_3G_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/HEPG2_WCE_rep2/HI.2449.006.Index_8.HEPG2_ctrl_WCE_MF6_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/HEPG2_WCE_rep2	HI.2449.006.Index_8.HEPG2_ctrl_WCE_MF6_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/HEPG2_SMC1A_rep2/HI.2449.006.Index_9.HEPG2_ctrl_SMC1A_MF6_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/HEPG2_SMC1A_rep2	HI.2449.006.Index_9.HEPG2_ctrl_SMC1A_MF6_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/HEPG2_NIPBL_rep2/HI.2449.006.Index_1.HEPG2_ctrl_NIPBL_MF6_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/HEPG2_NIPBL_rep2	HI.2449.006.Index_1.HEPG2_ctrl_NIPBL_MF6_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/HEPG2_MED1_rep2/HI.2449.006.Index_3.HEPG2_ctrl_MED1_MF6_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/HEPG2_MED1_rep2	HI.2449.006.Index_3.HEPG2_ctrl_MED1_MF6_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/HEPG2_BRD4_rep1/ENCFF953UBF.fastq.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/HEPG2_BRD4_rep1	ENCFF953UBF.fastq	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/HEPG2_BRD4_rep2/ENCFF946FEX.fastq.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/HEPG2_BRD4_rep2	ENCFF946FEX.fastq	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/HEPG2_BRD4_WCE_rep1/ENCFF474ZPM.fastq.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/HEPG2_BRD4_WCE_rep1	ENCFF474ZPM.fastq	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/HEPG2_BRD4_WCE_rep2/ENCFF653JWL.fastq.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/HEPG2_BRD4_WCE_rep2	ENCFF653JWL.fastq	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
cut -f1,3- metrics/trimReadsetTable.tsv | awk -F"	" '{OFS="	"; if (NR==1) {if ($2=="Raw Paired Reads #") {paired=1};print "Sample", "Raw Reads #", "Surviving Reads #", "Surviving %"} else {if (paired) {$2=$2*2; $3=$3*2}; raw[$1]+=$2; surviving[$1]+=$3}}END{for (sample in raw){print sample, raw[sample], surviving[sample], surviving[sample] / raw[sample] * 100}}' \
  > metrics/trimSampleTable.tsv && \
mkdir -p report && \
cp metrics/trimReadsetTable.tsv metrics/trimSampleTable.tsv report/ && \
trim_readset_table_md=`LC_NUMERIC=en_CA awk -F "	" '{OFS="|"; if (NR == 1) {$1 = $1; print $0; print "-----|-----|-----:|-----:|-----:"} else {print $1, $2, sprintf("%\47d", $3), sprintf("%\47d", $4), sprintf("%.1f", $5)}}' metrics/trimReadsetTable.tsv` && \
pandoc \
  /cvmfs/soft.mugqic/CentOS6/software/genpipes/genpipes-3.1.2/bfx/report/Illumina.merge_trimmomatic_stats.md \
  --template /cvmfs/soft.mugqic/CentOS6/software/genpipes/genpipes-3.1.2/bfx/report/Illumina.merge_trimmomatic_stats.md \
  --variable trailing_min_quality=30 \
  --variable min_length=50 \
  --variable read_type=Single \
  --variable trim_readset_table="$trim_readset_table_md" \
  --to markdown \
  > report/Illumina.merge_trimmomatic_stats.md
merge_trimmomatic_stats.b08357e409aca59c54e3b0e2beba8ff3.mugqic.done
)
merge_trimmomatic_stats_1_JOB_ID=$(echo "#! /bin/bash
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
echo "$merge_trimmomatic_stats_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# STEP: bwa_mem_picard_sort_sam
#-------------------------------------------------------------------------------
STEP=bwa_mem_picard_sort_sam
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_1_JOB_ID: bwa_mem_picard_sort_sam.HI.1308.004.Index_18.HEPG2_WCE_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.1308.004.Index_18.HEPG2_WCE_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_1_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.1308.004.Index_18.HEPG2_WCE_rep1.bam.f9f421419eb8233fbceef7b31eaac83a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.1308.004.Index_18.HEPG2_WCE_rep1.bam.f9f421419eb8233fbceef7b31eaac83a.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p alignment/HEPG2_WCE_rep1/HI.1308.004.Index_18.HEPG2_WCE_rep1.bam && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:HI.1308.004.Index_18.HEPG2_WCE_rep1.bam	SM:HEPG2_WCE_rep1	LB:HEPG2_WCE_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/HEPG2_WCE_rep1/HI.1308.004.Index_18.HEPG2_WCE_rep1.bam.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/HEPG2_WCE_rep1/HI.1308.004.Index_18.HEPG2_WCE_rep1.bam/HI.1308.004.Index_18.HEPG2_WCE_rep1.bam.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.HI.1308.004.Index_18.HEPG2_WCE_rep1.bam.f9f421419eb8233fbceef7b31eaac83a.mugqic.done
)
bwa_mem_picard_sort_sam_1_JOB_ID=$(echo "#! /bin/bash
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=64G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$bwa_mem_picard_sort_sam_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_2_JOB_ID: bwa_mem_picard_sort_sam.HI.1308.004.Index_19.HEPG2_NIPBL_ChIP67_1G_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.1308.004.Index_19.HEPG2_NIPBL_ChIP67_1G_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_2_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.1308.004.Index_19.HEPG2_NIPBL_ChIP67_1G_rep1.bam.4c4841d01336436188cfbfdb61119756.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.1308.004.Index_19.HEPG2_NIPBL_ChIP67_1G_rep1.bam.4c4841d01336436188cfbfdb61119756.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p alignment/HEPG2_NIPBL_rep1/HI.1308.004.Index_19.HEPG2_NIPBL_ChIP67_1G_rep1.bam && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:HI.1308.004.Index_19.HEPG2_NIPBL_ChIP67_1G_rep1.bam	SM:HEPG2_NIPBL_rep1	LB:HEPG2_NIPBL_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/HEPG2_NIPBL_rep1/HI.1308.004.Index_19.HEPG2_NIPBL_ChIP67_1G_rep1.bam.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/HEPG2_NIPBL_rep1/HI.1308.004.Index_19.HEPG2_NIPBL_ChIP67_1G_rep1.bam/HI.1308.004.Index_19.HEPG2_NIPBL_ChIP67_1G_rep1.bam.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.HI.1308.004.Index_19.HEPG2_NIPBL_ChIP67_1G_rep1.bam.4c4841d01336436188cfbfdb61119756.mugqic.done
)
bwa_mem_picard_sort_sam_2_JOB_ID=$(echo "#! /bin/bash
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=64G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$bwa_mem_picard_sort_sam_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_3_JOB_ID: bwa_mem_picard_sort_sam.HI.1308.004.Index_22.HEPG2_Cdk9_ChIP67_5G_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.1308.004.Index_22.HEPG2_Cdk9_ChIP67_5G_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_3_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.1308.004.Index_22.HEPG2_Cdk9_ChIP67_5G_rep1.bam.de7ea2e9da32baff426e30c34127d1f3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.1308.004.Index_22.HEPG2_Cdk9_ChIP67_5G_rep1.bam.de7ea2e9da32baff426e30c34127d1f3.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p alignment/HEPG2_CDK9_rep1/HI.1308.004.Index_22.HEPG2_Cdk9_ChIP67_5G_rep1.bam && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:HI.1308.004.Index_22.HEPG2_Cdk9_ChIP67_5G_rep1.bam	SM:HEPG2_CDK9_rep1	LB:HEPG2_CDK9_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/HEPG2_CDK9_rep1/HI.1308.004.Index_22.HEPG2_Cdk9_ChIP67_5G_rep1.bam.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/HEPG2_CDK9_rep1/HI.1308.004.Index_22.HEPG2_Cdk9_ChIP67_5G_rep1.bam/HI.1308.004.Index_22.HEPG2_Cdk9_ChIP67_5G_rep1.bam.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.HI.1308.004.Index_22.HEPG2_Cdk9_ChIP67_5G_rep1.bam.de7ea2e9da32baff426e30c34127d1f3.mugqic.done
)
bwa_mem_picard_sort_sam_3_JOB_ID=$(echo "#! /bin/bash
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=64G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$bwa_mem_picard_sort_sam_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_4_JOB_ID: bwa_mem_picard_sort_sam.HI.1308.004.Index_20.HEPG2_Med1_ChIP67_2G_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.1308.004.Index_20.HEPG2_Med1_ChIP67_2G_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_4_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.1308.004.Index_20.HEPG2_Med1_ChIP67_2G_rep1.bam.0d16ae641ff09d52cd5298929699b5ef.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.1308.004.Index_20.HEPG2_Med1_ChIP67_2G_rep1.bam.0d16ae641ff09d52cd5298929699b5ef.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p alignment/HEPG2_MED1_rep1/HI.1308.004.Index_20.HEPG2_Med1_ChIP67_2G_rep1.bam && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:HI.1308.004.Index_20.HEPG2_Med1_ChIP67_2G_rep1.bam	SM:HEPG2_MED1_rep1	LB:HEPG2_MED1_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/HEPG2_MED1_rep1/HI.1308.004.Index_20.HEPG2_Med1_ChIP67_2G_rep1.bam.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/HEPG2_MED1_rep1/HI.1308.004.Index_20.HEPG2_Med1_ChIP67_2G_rep1.bam/HI.1308.004.Index_20.HEPG2_Med1_ChIP67_2G_rep1.bam.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.HI.1308.004.Index_20.HEPG2_Med1_ChIP67_2G_rep1.bam.0d16ae641ff09d52cd5298929699b5ef.mugqic.done
)
bwa_mem_picard_sort_sam_4_JOB_ID=$(echo "#! /bin/bash
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=64G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$bwa_mem_picard_sort_sam_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_5_JOB_ID: bwa_mem_picard_sort_sam.HI.1308.004.Index_21.HEPG2_Smc1_ChIP67_3G_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.1308.004.Index_21.HEPG2_Smc1_ChIP67_3G_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_5_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.1308.004.Index_21.HEPG2_Smc1_ChIP67_3G_rep1.bam.0bb413e0e071bad47ebb903dfe4a0a08.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.1308.004.Index_21.HEPG2_Smc1_ChIP67_3G_rep1.bam.0bb413e0e071bad47ebb903dfe4a0a08.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p alignment/HEPG2_SMC1A_rep1/HI.1308.004.Index_21.HEPG2_Smc1_ChIP67_3G_rep1.bam && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:HI.1308.004.Index_21.HEPG2_Smc1_ChIP67_3G_rep1.bam	SM:HEPG2_SMC1A_rep1	LB:HEPG2_SMC1A_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/HEPG2_SMC1A_rep1/HI.1308.004.Index_21.HEPG2_Smc1_ChIP67_3G_rep1.bam.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/HEPG2_SMC1A_rep1/HI.1308.004.Index_21.HEPG2_Smc1_ChIP67_3G_rep1.bam/HI.1308.004.Index_21.HEPG2_Smc1_ChIP67_3G_rep1.bam.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.HI.1308.004.Index_21.HEPG2_Smc1_ChIP67_3G_rep1.bam.0bb413e0e071bad47ebb903dfe4a0a08.mugqic.done
)
bwa_mem_picard_sort_sam_5_JOB_ID=$(echo "#! /bin/bash
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=64G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$bwa_mem_picard_sort_sam_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_6_JOB_ID: bwa_mem_picard_sort_sam.HI.2449.006.Index_8.HEPG2_ctrl_WCE_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.2449.006.Index_8.HEPG2_ctrl_WCE_MF6_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_6_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.2449.006.Index_8.HEPG2_ctrl_WCE_MF6_rep1.bam.c93a23d32784fd5631f12f474a71886e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.2449.006.Index_8.HEPG2_ctrl_WCE_MF6_rep1.bam.c93a23d32784fd5631f12f474a71886e.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p alignment/HEPG2_WCE_rep2/HI.2449.006.Index_8.HEPG2_ctrl_WCE_MF6_rep1.bam && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:HI.2449.006.Index_8.HEPG2_ctrl_WCE_MF6_rep1.bam	SM:HEPG2_WCE_rep2	LB:HEPG2_WCE_rep2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/HEPG2_WCE_rep2/HI.2449.006.Index_8.HEPG2_ctrl_WCE_MF6_rep1.bam.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/HEPG2_WCE_rep2/HI.2449.006.Index_8.HEPG2_ctrl_WCE_MF6_rep1.bam/HI.2449.006.Index_8.HEPG2_ctrl_WCE_MF6_rep1.bam.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.HI.2449.006.Index_8.HEPG2_ctrl_WCE_MF6_rep1.bam.c93a23d32784fd5631f12f474a71886e.mugqic.done
)
bwa_mem_picard_sort_sam_6_JOB_ID=$(echo "#! /bin/bash
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=64G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$bwa_mem_picard_sort_sam_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_7_JOB_ID: bwa_mem_picard_sort_sam.HI.2449.006.Index_9.HEPG2_ctrl_SMC1A_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.2449.006.Index_9.HEPG2_ctrl_SMC1A_MF6_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_7_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.2449.006.Index_9.HEPG2_ctrl_SMC1A_MF6_rep1.bam.f39d29fec17c8446b11ff960d324d5fa.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.2449.006.Index_9.HEPG2_ctrl_SMC1A_MF6_rep1.bam.f39d29fec17c8446b11ff960d324d5fa.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p alignment/HEPG2_SMC1A_rep2/HI.2449.006.Index_9.HEPG2_ctrl_SMC1A_MF6_rep1.bam && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:HI.2449.006.Index_9.HEPG2_ctrl_SMC1A_MF6_rep1.bam	SM:HEPG2_SMC1A_rep2	LB:HEPG2_SMC1A_rep2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/HEPG2_SMC1A_rep2/HI.2449.006.Index_9.HEPG2_ctrl_SMC1A_MF6_rep1.bam.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/HEPG2_SMC1A_rep2/HI.2449.006.Index_9.HEPG2_ctrl_SMC1A_MF6_rep1.bam/HI.2449.006.Index_9.HEPG2_ctrl_SMC1A_MF6_rep1.bam.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.HI.2449.006.Index_9.HEPG2_ctrl_SMC1A_MF6_rep1.bam.f39d29fec17c8446b11ff960d324d5fa.mugqic.done
)
bwa_mem_picard_sort_sam_7_JOB_ID=$(echo "#! /bin/bash
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=64G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$bwa_mem_picard_sort_sam_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_8_JOB_ID: bwa_mem_picard_sort_sam.HI.2449.006.Index_1.HEPG2_ctrl_NIPBL_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.2449.006.Index_1.HEPG2_ctrl_NIPBL_MF6_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_8_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.2449.006.Index_1.HEPG2_ctrl_NIPBL_MF6_rep1.bam.47fb631df72b1ae1f91ae00b51c889fc.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.2449.006.Index_1.HEPG2_ctrl_NIPBL_MF6_rep1.bam.47fb631df72b1ae1f91ae00b51c889fc.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p alignment/HEPG2_NIPBL_rep2/HI.2449.006.Index_1.HEPG2_ctrl_NIPBL_MF6_rep1.bam && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:HI.2449.006.Index_1.HEPG2_ctrl_NIPBL_MF6_rep1.bam	SM:HEPG2_NIPBL_rep2	LB:HEPG2_NIPBL_rep2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/HEPG2_NIPBL_rep2/HI.2449.006.Index_1.HEPG2_ctrl_NIPBL_MF6_rep1.bam.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/HEPG2_NIPBL_rep2/HI.2449.006.Index_1.HEPG2_ctrl_NIPBL_MF6_rep1.bam/HI.2449.006.Index_1.HEPG2_ctrl_NIPBL_MF6_rep1.bam.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.HI.2449.006.Index_1.HEPG2_ctrl_NIPBL_MF6_rep1.bam.47fb631df72b1ae1f91ae00b51c889fc.mugqic.done
)
bwa_mem_picard_sort_sam_8_JOB_ID=$(echo "#! /bin/bash
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=64G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$bwa_mem_picard_sort_sam_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_9_JOB_ID: bwa_mem_picard_sort_sam.HI.2449.006.Index_3.HEPG2_ctrl_MED1_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.2449.006.Index_3.HEPG2_ctrl_MED1_MF6_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_9_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.2449.006.Index_3.HEPG2_ctrl_MED1_MF6_rep1.bam.e7dca4809fc9d16c6be8f2319073d21e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.2449.006.Index_3.HEPG2_ctrl_MED1_MF6_rep1.bam.e7dca4809fc9d16c6be8f2319073d21e.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p alignment/HEPG2_MED1_rep2/HI.2449.006.Index_3.HEPG2_ctrl_MED1_MF6_rep1.bam && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:HI.2449.006.Index_3.HEPG2_ctrl_MED1_MF6_rep1.bam	SM:HEPG2_MED1_rep2	LB:HEPG2_MED1_rep2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/HEPG2_MED1_rep2/HI.2449.006.Index_3.HEPG2_ctrl_MED1_MF6_rep1.bam.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/HEPG2_MED1_rep2/HI.2449.006.Index_3.HEPG2_ctrl_MED1_MF6_rep1.bam/HI.2449.006.Index_3.HEPG2_ctrl_MED1_MF6_rep1.bam.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.HI.2449.006.Index_3.HEPG2_ctrl_MED1_MF6_rep1.bam.e7dca4809fc9d16c6be8f2319073d21e.mugqic.done
)
bwa_mem_picard_sort_sam_9_JOB_ID=$(echo "#! /bin/bash
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=64G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$bwa_mem_picard_sort_sam_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_10_JOB_ID: bwa_mem_picard_sort_sam.ENCFF953UBF.fastq
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.ENCFF953UBF.fastq
JOB_DEPENDENCIES=$trimmomatic_10_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.ENCFF953UBF.fastq.7f45e512676e89a9d31762e4e0f44cd8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.ENCFF953UBF.fastq.7f45e512676e89a9d31762e4e0f44cd8.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p alignment/HEPG2_BRD4_rep1/ENCFF953UBF.fastq && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:ENCFF953UBF.fastq	SM:HEPG2_BRD4_rep1	LB:HEPG2_BRD4_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/HEPG2_BRD4_rep1/ENCFF953UBF.fastq.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/HEPG2_BRD4_rep1/ENCFF953UBF.fastq/ENCFF953UBF.fastq.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.ENCFF953UBF.fastq.7f45e512676e89a9d31762e4e0f44cd8.mugqic.done
)
bwa_mem_picard_sort_sam_10_JOB_ID=$(echo "#! /bin/bash
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=64G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$bwa_mem_picard_sort_sam_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_11_JOB_ID: bwa_mem_picard_sort_sam.ENCFF946FEX.fastq
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.ENCFF946FEX.fastq
JOB_DEPENDENCIES=$trimmomatic_11_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.ENCFF946FEX.fastq.901884b42dd5d16cb28f62010fd26f32.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.ENCFF946FEX.fastq.901884b42dd5d16cb28f62010fd26f32.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p alignment/HEPG2_BRD4_rep2/ENCFF946FEX.fastq && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:ENCFF946FEX.fastq	SM:HEPG2_BRD4_rep2	LB:HEPG2_BRD4_rep2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/HEPG2_BRD4_rep2/ENCFF946FEX.fastq.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/HEPG2_BRD4_rep2/ENCFF946FEX.fastq/ENCFF946FEX.fastq.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.ENCFF946FEX.fastq.901884b42dd5d16cb28f62010fd26f32.mugqic.done
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=64G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$bwa_mem_picard_sort_sam_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_12_JOB_ID: bwa_mem_picard_sort_sam.ENCFF474ZPM.fastq
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.ENCFF474ZPM.fastq
JOB_DEPENDENCIES=$trimmomatic_12_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.ENCFF474ZPM.fastq.eaacace0cc458c84841ad67ffe505b6b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.ENCFF474ZPM.fastq.eaacace0cc458c84841ad67ffe505b6b.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p alignment/HEPG2_BRD4_WCE_rep1/ENCFF474ZPM.fastq && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:ENCFF474ZPM.fastq	SM:HEPG2_BRD4_WCE_rep1	LB:HEPG2_BRD4_WCE_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/HEPG2_BRD4_WCE_rep1/ENCFF474ZPM.fastq.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/HEPG2_BRD4_WCE_rep1/ENCFF474ZPM.fastq/ENCFF474ZPM.fastq.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.ENCFF474ZPM.fastq.eaacace0cc458c84841ad67ffe505b6b.mugqic.done
)
bwa_mem_picard_sort_sam_12_JOB_ID=$(echo "#! /bin/bash
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=64G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$bwa_mem_picard_sort_sam_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_13_JOB_ID: bwa_mem_picard_sort_sam.ENCFF653JWL.fastq
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.ENCFF653JWL.fastq
JOB_DEPENDENCIES=$trimmomatic_13_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.ENCFF653JWL.fastq.cd4175cb2e09863758791baf299455cd.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.ENCFF653JWL.fastq.cd4175cb2e09863758791baf299455cd.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p alignment/HEPG2_BRD4_WCE_rep2/ENCFF653JWL.fastq && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:ENCFF653JWL.fastq	SM:HEPG2_BRD4_WCE_rep2	LB:HEPG2_BRD4_WCE_rep2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/HEPG2_BRD4_WCE_rep2/ENCFF653JWL.fastq.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/HEPG2_BRD4_WCE_rep2/ENCFF653JWL.fastq/ENCFF653JWL.fastq.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.ENCFF653JWL.fastq.cd4175cb2e09863758791baf299455cd.mugqic.done
)
bwa_mem_picard_sort_sam_13_JOB_ID=$(echo "#! /bin/bash
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=64G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$bwa_mem_picard_sort_sam_13_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_14_JOB_ID: bwa_mem_picard_sort_sam_report
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam_report
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_1_JOB_ID:$bwa_mem_picard_sort_sam_2_JOB_ID:$bwa_mem_picard_sort_sam_3_JOB_ID:$bwa_mem_picard_sort_sam_4_JOB_ID:$bwa_mem_picard_sort_sam_5_JOB_ID:$bwa_mem_picard_sort_sam_6_JOB_ID:$bwa_mem_picard_sort_sam_7_JOB_ID:$bwa_mem_picard_sort_sam_8_JOB_ID:$bwa_mem_picard_sort_sam_9_JOB_ID:$bwa_mem_picard_sort_sam_10_JOB_ID:$bwa_mem_picard_sort_sam_11_JOB_ID:$bwa_mem_picard_sort_sam_12_JOB_ID:$bwa_mem_picard_sort_sam_13_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam_report.6d966dfaed922f59ab242900bce11284.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam_report.6d966dfaed922f59ab242900bce11284.mugqic.done'
module load mugqic/pandoc/1.15.2 && \
mkdir -p report && \
pandoc --to=markdown \
  --template /cvmfs/soft.mugqic/CentOS6/software/genpipes/genpipes-3.1.2/bfx/report/DnaSeq.bwa_mem_picard_sort_sam.md \
  --variable scientific_name="Homo_sapiens" \
  --variable assembly="GRCh38" \
  /cvmfs/soft.mugqic/CentOS6/software/genpipes/genpipes-3.1.2/bfx/report/DnaSeq.bwa_mem_picard_sort_sam.md \
  > report/DnaSeq.bwa_mem_picard_sort_sam.md
bwa_mem_picard_sort_sam_report.6d966dfaed922f59ab242900bce11284.mugqic.done
)
bwa_mem_picard_sort_sam_14_JOB_ID=$(echo "#! /bin/bash
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
echo "$bwa_mem_picard_sort_sam_14_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# STEP: samtools_view_filter
#-------------------------------------------------------------------------------
STEP=samtools_view_filter
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_1_JOB_ID: samtools_view_filter.HI.1308.004.Index_18.HEPG2_WCE_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.1308.004.Index_18.HEPG2_WCE_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_1_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.1308.004.Index_18.HEPG2_WCE_rep1.bam.74b75b17601847d8a94642a767bb9bbc.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.1308.004.Index_18.HEPG2_WCE_rep1.bam.74b75b17601847d8a94642a767bb9bbc.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/HEPG2_WCE_rep1/HI.1308.004.Index_18.HEPG2_WCE_rep1.bam/HI.1308.004.Index_18.HEPG2_WCE_rep1.bam.sorted.bam \
  > alignment/HEPG2_WCE_rep1/HI.1308.004.Index_18.HEPG2_WCE_rep1.bam/HI.1308.004.Index_18.HEPG2_WCE_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.1308.004.Index_18.HEPG2_WCE_rep1.bam.74b75b17601847d8a94642a767bb9bbc.mugqic.done
)
samtools_view_filter_1_JOB_ID=$(echo "#! /bin/bash
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
echo "$samtools_view_filter_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_2_JOB_ID: samtools_view_filter.HI.1308.004.Index_19.HEPG2_NIPBL_ChIP67_1G_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.1308.004.Index_19.HEPG2_NIPBL_ChIP67_1G_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_2_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.1308.004.Index_19.HEPG2_NIPBL_ChIP67_1G_rep1.bam.cf723c52d617b878489b5e35cec84e93.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.1308.004.Index_19.HEPG2_NIPBL_ChIP67_1G_rep1.bam.cf723c52d617b878489b5e35cec84e93.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/HEPG2_NIPBL_rep1/HI.1308.004.Index_19.HEPG2_NIPBL_ChIP67_1G_rep1.bam/HI.1308.004.Index_19.HEPG2_NIPBL_ChIP67_1G_rep1.bam.sorted.bam \
  > alignment/HEPG2_NIPBL_rep1/HI.1308.004.Index_19.HEPG2_NIPBL_ChIP67_1G_rep1.bam/HI.1308.004.Index_19.HEPG2_NIPBL_ChIP67_1G_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.1308.004.Index_19.HEPG2_NIPBL_ChIP67_1G_rep1.bam.cf723c52d617b878489b5e35cec84e93.mugqic.done
)
samtools_view_filter_2_JOB_ID=$(echo "#! /bin/bash
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
echo "$samtools_view_filter_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_3_JOB_ID: samtools_view_filter.HI.1308.004.Index_22.HEPG2_Cdk9_ChIP67_5G_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.1308.004.Index_22.HEPG2_Cdk9_ChIP67_5G_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_3_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.1308.004.Index_22.HEPG2_Cdk9_ChIP67_5G_rep1.bam.8193a1ba3f26fc0162aea8cf1de3809e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.1308.004.Index_22.HEPG2_Cdk9_ChIP67_5G_rep1.bam.8193a1ba3f26fc0162aea8cf1de3809e.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/HEPG2_CDK9_rep1/HI.1308.004.Index_22.HEPG2_Cdk9_ChIP67_5G_rep1.bam/HI.1308.004.Index_22.HEPG2_Cdk9_ChIP67_5G_rep1.bam.sorted.bam \
  > alignment/HEPG2_CDK9_rep1/HI.1308.004.Index_22.HEPG2_Cdk9_ChIP67_5G_rep1.bam/HI.1308.004.Index_22.HEPG2_Cdk9_ChIP67_5G_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.1308.004.Index_22.HEPG2_Cdk9_ChIP67_5G_rep1.bam.8193a1ba3f26fc0162aea8cf1de3809e.mugqic.done
)
samtools_view_filter_3_JOB_ID=$(echo "#! /bin/bash
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
echo "$samtools_view_filter_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_4_JOB_ID: samtools_view_filter.HI.1308.004.Index_20.HEPG2_Med1_ChIP67_2G_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.1308.004.Index_20.HEPG2_Med1_ChIP67_2G_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_4_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.1308.004.Index_20.HEPG2_Med1_ChIP67_2G_rep1.bam.a72280e766b5ef74f412e81481be1768.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.1308.004.Index_20.HEPG2_Med1_ChIP67_2G_rep1.bam.a72280e766b5ef74f412e81481be1768.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/HEPG2_MED1_rep1/HI.1308.004.Index_20.HEPG2_Med1_ChIP67_2G_rep1.bam/HI.1308.004.Index_20.HEPG2_Med1_ChIP67_2G_rep1.bam.sorted.bam \
  > alignment/HEPG2_MED1_rep1/HI.1308.004.Index_20.HEPG2_Med1_ChIP67_2G_rep1.bam/HI.1308.004.Index_20.HEPG2_Med1_ChIP67_2G_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.1308.004.Index_20.HEPG2_Med1_ChIP67_2G_rep1.bam.a72280e766b5ef74f412e81481be1768.mugqic.done
)
samtools_view_filter_4_JOB_ID=$(echo "#! /bin/bash
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
echo "$samtools_view_filter_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_5_JOB_ID: samtools_view_filter.HI.1308.004.Index_21.HEPG2_Smc1_ChIP67_3G_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.1308.004.Index_21.HEPG2_Smc1_ChIP67_3G_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_5_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.1308.004.Index_21.HEPG2_Smc1_ChIP67_3G_rep1.bam.eb4c3e8de6f345a9ac829c033d6c3c85.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.1308.004.Index_21.HEPG2_Smc1_ChIP67_3G_rep1.bam.eb4c3e8de6f345a9ac829c033d6c3c85.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/HEPG2_SMC1A_rep1/HI.1308.004.Index_21.HEPG2_Smc1_ChIP67_3G_rep1.bam/HI.1308.004.Index_21.HEPG2_Smc1_ChIP67_3G_rep1.bam.sorted.bam \
  > alignment/HEPG2_SMC1A_rep1/HI.1308.004.Index_21.HEPG2_Smc1_ChIP67_3G_rep1.bam/HI.1308.004.Index_21.HEPG2_Smc1_ChIP67_3G_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.1308.004.Index_21.HEPG2_Smc1_ChIP67_3G_rep1.bam.eb4c3e8de6f345a9ac829c033d6c3c85.mugqic.done
)
samtools_view_filter_5_JOB_ID=$(echo "#! /bin/bash
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
echo "$samtools_view_filter_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_6_JOB_ID: samtools_view_filter.HI.2449.006.Index_8.HEPG2_ctrl_WCE_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.2449.006.Index_8.HEPG2_ctrl_WCE_MF6_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_6_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.2449.006.Index_8.HEPG2_ctrl_WCE_MF6_rep1.bam.d0cc91ffaf4c5f916c725f853289181f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.2449.006.Index_8.HEPG2_ctrl_WCE_MF6_rep1.bam.d0cc91ffaf4c5f916c725f853289181f.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/HEPG2_WCE_rep2/HI.2449.006.Index_8.HEPG2_ctrl_WCE_MF6_rep1.bam/HI.2449.006.Index_8.HEPG2_ctrl_WCE_MF6_rep1.bam.sorted.bam \
  > alignment/HEPG2_WCE_rep2/HI.2449.006.Index_8.HEPG2_ctrl_WCE_MF6_rep1.bam/HI.2449.006.Index_8.HEPG2_ctrl_WCE_MF6_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.2449.006.Index_8.HEPG2_ctrl_WCE_MF6_rep1.bam.d0cc91ffaf4c5f916c725f853289181f.mugqic.done
)
samtools_view_filter_6_JOB_ID=$(echo "#! /bin/bash
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
echo "$samtools_view_filter_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_7_JOB_ID: samtools_view_filter.HI.2449.006.Index_9.HEPG2_ctrl_SMC1A_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.2449.006.Index_9.HEPG2_ctrl_SMC1A_MF6_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_7_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.2449.006.Index_9.HEPG2_ctrl_SMC1A_MF6_rep1.bam.1b2d834da539d35fa49755b284b1524c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.2449.006.Index_9.HEPG2_ctrl_SMC1A_MF6_rep1.bam.1b2d834da539d35fa49755b284b1524c.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/HEPG2_SMC1A_rep2/HI.2449.006.Index_9.HEPG2_ctrl_SMC1A_MF6_rep1.bam/HI.2449.006.Index_9.HEPG2_ctrl_SMC1A_MF6_rep1.bam.sorted.bam \
  > alignment/HEPG2_SMC1A_rep2/HI.2449.006.Index_9.HEPG2_ctrl_SMC1A_MF6_rep1.bam/HI.2449.006.Index_9.HEPG2_ctrl_SMC1A_MF6_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.2449.006.Index_9.HEPG2_ctrl_SMC1A_MF6_rep1.bam.1b2d834da539d35fa49755b284b1524c.mugqic.done
)
samtools_view_filter_7_JOB_ID=$(echo "#! /bin/bash
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
echo "$samtools_view_filter_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_8_JOB_ID: samtools_view_filter.HI.2449.006.Index_1.HEPG2_ctrl_NIPBL_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.2449.006.Index_1.HEPG2_ctrl_NIPBL_MF6_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_8_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.2449.006.Index_1.HEPG2_ctrl_NIPBL_MF6_rep1.bam.cc0b6709da434aa52a8e2a3661ed14e6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.2449.006.Index_1.HEPG2_ctrl_NIPBL_MF6_rep1.bam.cc0b6709da434aa52a8e2a3661ed14e6.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/HEPG2_NIPBL_rep2/HI.2449.006.Index_1.HEPG2_ctrl_NIPBL_MF6_rep1.bam/HI.2449.006.Index_1.HEPG2_ctrl_NIPBL_MF6_rep1.bam.sorted.bam \
  > alignment/HEPG2_NIPBL_rep2/HI.2449.006.Index_1.HEPG2_ctrl_NIPBL_MF6_rep1.bam/HI.2449.006.Index_1.HEPG2_ctrl_NIPBL_MF6_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.2449.006.Index_1.HEPG2_ctrl_NIPBL_MF6_rep1.bam.cc0b6709da434aa52a8e2a3661ed14e6.mugqic.done
)
samtools_view_filter_8_JOB_ID=$(echo "#! /bin/bash
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
echo "$samtools_view_filter_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_9_JOB_ID: samtools_view_filter.HI.2449.006.Index_3.HEPG2_ctrl_MED1_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.2449.006.Index_3.HEPG2_ctrl_MED1_MF6_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_9_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.2449.006.Index_3.HEPG2_ctrl_MED1_MF6_rep1.bam.db0df42fc755017e0743b6982eb87d34.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.2449.006.Index_3.HEPG2_ctrl_MED1_MF6_rep1.bam.db0df42fc755017e0743b6982eb87d34.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/HEPG2_MED1_rep2/HI.2449.006.Index_3.HEPG2_ctrl_MED1_MF6_rep1.bam/HI.2449.006.Index_3.HEPG2_ctrl_MED1_MF6_rep1.bam.sorted.bam \
  > alignment/HEPG2_MED1_rep2/HI.2449.006.Index_3.HEPG2_ctrl_MED1_MF6_rep1.bam/HI.2449.006.Index_3.HEPG2_ctrl_MED1_MF6_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.2449.006.Index_3.HEPG2_ctrl_MED1_MF6_rep1.bam.db0df42fc755017e0743b6982eb87d34.mugqic.done
)
samtools_view_filter_9_JOB_ID=$(echo "#! /bin/bash
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
echo "$samtools_view_filter_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_10_JOB_ID: samtools_view_filter.ENCFF953UBF.fastq
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.ENCFF953UBF.fastq
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_10_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.ENCFF953UBF.fastq.319af3484449160342095215398fa5ef.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.ENCFF953UBF.fastq.319af3484449160342095215398fa5ef.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/HEPG2_BRD4_rep1/ENCFF953UBF.fastq/ENCFF953UBF.fastq.sorted.bam \
  > alignment/HEPG2_BRD4_rep1/ENCFF953UBF.fastq/ENCFF953UBF.fastq.sorted.filtered.bam
samtools_view_filter.ENCFF953UBF.fastq.319af3484449160342095215398fa5ef.mugqic.done
)
samtools_view_filter_10_JOB_ID=$(echo "#! /bin/bash
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
echo "$samtools_view_filter_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_11_JOB_ID: samtools_view_filter.ENCFF946FEX.fastq
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.ENCFF946FEX.fastq
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_11_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.ENCFF946FEX.fastq.41139531fc862707c4f6e9943bd40d8c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.ENCFF946FEX.fastq.41139531fc862707c4f6e9943bd40d8c.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/HEPG2_BRD4_rep2/ENCFF946FEX.fastq/ENCFF946FEX.fastq.sorted.bam \
  > alignment/HEPG2_BRD4_rep2/ENCFF946FEX.fastq/ENCFF946FEX.fastq.sorted.filtered.bam
samtools_view_filter.ENCFF946FEX.fastq.41139531fc862707c4f6e9943bd40d8c.mugqic.done
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

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_12_JOB_ID: samtools_view_filter.ENCFF474ZPM.fastq
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.ENCFF474ZPM.fastq
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_12_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.ENCFF474ZPM.fastq.2de9d51a2465dd406d994b7ab89f859f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.ENCFF474ZPM.fastq.2de9d51a2465dd406d994b7ab89f859f.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/HEPG2_BRD4_WCE_rep1/ENCFF474ZPM.fastq/ENCFF474ZPM.fastq.sorted.bam \
  > alignment/HEPG2_BRD4_WCE_rep1/ENCFF474ZPM.fastq/ENCFF474ZPM.fastq.sorted.filtered.bam
samtools_view_filter.ENCFF474ZPM.fastq.2de9d51a2465dd406d994b7ab89f859f.mugqic.done
)
samtools_view_filter_12_JOB_ID=$(echo "#! /bin/bash
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
echo "$samtools_view_filter_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_13_JOB_ID: samtools_view_filter.ENCFF653JWL.fastq
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.ENCFF653JWL.fastq
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_13_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.ENCFF653JWL.fastq.b2c6668a20b84d349403b7c2bf468511.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.ENCFF653JWL.fastq.b2c6668a20b84d349403b7c2bf468511.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/HEPG2_BRD4_WCE_rep2/ENCFF653JWL.fastq/ENCFF653JWL.fastq.sorted.bam \
  > alignment/HEPG2_BRD4_WCE_rep2/ENCFF653JWL.fastq/ENCFF653JWL.fastq.sorted.filtered.bam
samtools_view_filter.ENCFF653JWL.fastq.b2c6668a20b84d349403b7c2bf468511.mugqic.done
)
samtools_view_filter_13_JOB_ID=$(echo "#! /bin/bash
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
echo "$samtools_view_filter_13_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_14_JOB_ID: samtools_view_filter_report
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter_report
JOB_DEPENDENCIES=$samtools_view_filter_1_JOB_ID:$samtools_view_filter_2_JOB_ID:$samtools_view_filter_3_JOB_ID:$samtools_view_filter_4_JOB_ID:$samtools_view_filter_5_JOB_ID:$samtools_view_filter_6_JOB_ID:$samtools_view_filter_7_JOB_ID:$samtools_view_filter_8_JOB_ID:$samtools_view_filter_9_JOB_ID:$samtools_view_filter_10_JOB_ID:$samtools_view_filter_11_JOB_ID:$samtools_view_filter_12_JOB_ID:$samtools_view_filter_13_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter_report.8d06991c14c95c9c183f78f67effc8fb.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter_report.8d06991c14c95c9c183f78f67effc8fb.mugqic.done'
module load mugqic/pandoc/1.15.2 && \
mkdir -p report && \
pandoc --to=markdown \
  --template /cvmfs/soft.mugqic/CentOS6/software/genpipes/genpipes-3.1.2/bfx/report/ChipSeq.samtools_view_filter.md \
  --variable min_mapq="20" \
  /cvmfs/soft.mugqic/CentOS6/software/genpipes/genpipes-3.1.2/bfx/report/ChipSeq.samtools_view_filter.md \
  > report/ChipSeq.samtools_view_filter.md
samtools_view_filter_report.8d06991c14c95c9c183f78f67effc8fb.mugqic.done
)
samtools_view_filter_14_JOB_ID=$(echo "#! /bin/bash
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
echo "$samtools_view_filter_14_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# STEP: picard_merge_sam_files
#-------------------------------------------------------------------------------
STEP=picard_merge_sam_files
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_1_JOB_ID: symlink_readset_sample_bam.HEPG2_WCE_rep1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.HEPG2_WCE_rep1
JOB_DEPENDENCIES=$samtools_view_filter_1_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.HEPG2_WCE_rep1.8fc3df0f833d3519d972cc54852e80e4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.HEPG2_WCE_rep1.8fc3df0f833d3519d972cc54852e80e4.mugqic.done'
mkdir -p alignment/HEPG2_WCE_rep1 && \
ln -s -f HI.1308.004.Index_18.HEPG2_WCE_rep1.bam/HI.1308.004.Index_18.HEPG2_WCE_rep1.bam.sorted.filtered.bam alignment/HEPG2_WCE_rep1/HEPG2_WCE_rep1.merged.bam
symlink_readset_sample_bam.HEPG2_WCE_rep1.8fc3df0f833d3519d972cc54852e80e4.mugqic.done
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

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_2_JOB_ID: symlink_readset_sample_bam.HEPG2_NIPBL_rep1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.HEPG2_NIPBL_rep1
JOB_DEPENDENCIES=$samtools_view_filter_2_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.HEPG2_NIPBL_rep1.0dd72ecf25cdb89b66f441e67d701607.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.HEPG2_NIPBL_rep1.0dd72ecf25cdb89b66f441e67d701607.mugqic.done'
mkdir -p alignment/HEPG2_NIPBL_rep1 && \
ln -s -f HI.1308.004.Index_19.HEPG2_NIPBL_ChIP67_1G_rep1.bam/HI.1308.004.Index_19.HEPG2_NIPBL_ChIP67_1G_rep1.bam.sorted.filtered.bam alignment/HEPG2_NIPBL_rep1/HEPG2_NIPBL_rep1.merged.bam
symlink_readset_sample_bam.HEPG2_NIPBL_rep1.0dd72ecf25cdb89b66f441e67d701607.mugqic.done
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

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_3_JOB_ID: symlink_readset_sample_bam.HEPG2_CDK9_rep1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.HEPG2_CDK9_rep1
JOB_DEPENDENCIES=$samtools_view_filter_3_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.HEPG2_CDK9_rep1.400315f8a391b42f37e7c77e55ec53d5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.HEPG2_CDK9_rep1.400315f8a391b42f37e7c77e55ec53d5.mugqic.done'
mkdir -p alignment/HEPG2_CDK9_rep1 && \
ln -s -f HI.1308.004.Index_22.HEPG2_Cdk9_ChIP67_5G_rep1.bam/HI.1308.004.Index_22.HEPG2_Cdk9_ChIP67_5G_rep1.bam.sorted.filtered.bam alignment/HEPG2_CDK9_rep1/HEPG2_CDK9_rep1.merged.bam
symlink_readset_sample_bam.HEPG2_CDK9_rep1.400315f8a391b42f37e7c77e55ec53d5.mugqic.done
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

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_4_JOB_ID: symlink_readset_sample_bam.HEPG2_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.HEPG2_MED1_rep1
JOB_DEPENDENCIES=$samtools_view_filter_4_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.HEPG2_MED1_rep1.7f4fbd1dcbc0c6dc6ffe17b1c4a6e6e2.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.HEPG2_MED1_rep1.7f4fbd1dcbc0c6dc6ffe17b1c4a6e6e2.mugqic.done'
mkdir -p alignment/HEPG2_MED1_rep1 && \
ln -s -f HI.1308.004.Index_20.HEPG2_Med1_ChIP67_2G_rep1.bam/HI.1308.004.Index_20.HEPG2_Med1_ChIP67_2G_rep1.bam.sorted.filtered.bam alignment/HEPG2_MED1_rep1/HEPG2_MED1_rep1.merged.bam
symlink_readset_sample_bam.HEPG2_MED1_rep1.7f4fbd1dcbc0c6dc6ffe17b1c4a6e6e2.mugqic.done
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

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_5_JOB_ID: symlink_readset_sample_bam.HEPG2_SMC1A_rep1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.HEPG2_SMC1A_rep1
JOB_DEPENDENCIES=$samtools_view_filter_5_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.HEPG2_SMC1A_rep1.4e43692495f54b0798dc8273949c394e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.HEPG2_SMC1A_rep1.4e43692495f54b0798dc8273949c394e.mugqic.done'
mkdir -p alignment/HEPG2_SMC1A_rep1 && \
ln -s -f HI.1308.004.Index_21.HEPG2_Smc1_ChIP67_3G_rep1.bam/HI.1308.004.Index_21.HEPG2_Smc1_ChIP67_3G_rep1.bam.sorted.filtered.bam alignment/HEPG2_SMC1A_rep1/HEPG2_SMC1A_rep1.merged.bam
symlink_readset_sample_bam.HEPG2_SMC1A_rep1.4e43692495f54b0798dc8273949c394e.mugqic.done
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

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_6_JOB_ID: symlink_readset_sample_bam.HEPG2_WCE_rep2
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.HEPG2_WCE_rep2
JOB_DEPENDENCIES=$samtools_view_filter_6_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.HEPG2_WCE_rep2.187961afd58f69fb3a98ba2751a2f4cd.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.HEPG2_WCE_rep2.187961afd58f69fb3a98ba2751a2f4cd.mugqic.done'
mkdir -p alignment/HEPG2_WCE_rep2 && \
ln -s -f HI.2449.006.Index_8.HEPG2_ctrl_WCE_MF6_rep1.bam/HI.2449.006.Index_8.HEPG2_ctrl_WCE_MF6_rep1.bam.sorted.filtered.bam alignment/HEPG2_WCE_rep2/HEPG2_WCE_rep2.merged.bam
symlink_readset_sample_bam.HEPG2_WCE_rep2.187961afd58f69fb3a98ba2751a2f4cd.mugqic.done
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

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_7_JOB_ID: symlink_readset_sample_bam.HEPG2_SMC1A_rep2
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.HEPG2_SMC1A_rep2
JOB_DEPENDENCIES=$samtools_view_filter_7_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.HEPG2_SMC1A_rep2.3cde51969e27346fdb2b03213533a9ad.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.HEPG2_SMC1A_rep2.3cde51969e27346fdb2b03213533a9ad.mugqic.done'
mkdir -p alignment/HEPG2_SMC1A_rep2 && \
ln -s -f HI.2449.006.Index_9.HEPG2_ctrl_SMC1A_MF6_rep1.bam/HI.2449.006.Index_9.HEPG2_ctrl_SMC1A_MF6_rep1.bam.sorted.filtered.bam alignment/HEPG2_SMC1A_rep2/HEPG2_SMC1A_rep2.merged.bam
symlink_readset_sample_bam.HEPG2_SMC1A_rep2.3cde51969e27346fdb2b03213533a9ad.mugqic.done
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

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_8_JOB_ID: symlink_readset_sample_bam.HEPG2_NIPBL_rep2
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.HEPG2_NIPBL_rep2
JOB_DEPENDENCIES=$samtools_view_filter_8_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.HEPG2_NIPBL_rep2.8844a6ea0c58f45e1040216beeeb4c60.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.HEPG2_NIPBL_rep2.8844a6ea0c58f45e1040216beeeb4c60.mugqic.done'
mkdir -p alignment/HEPG2_NIPBL_rep2 && \
ln -s -f HI.2449.006.Index_1.HEPG2_ctrl_NIPBL_MF6_rep1.bam/HI.2449.006.Index_1.HEPG2_ctrl_NIPBL_MF6_rep1.bam.sorted.filtered.bam alignment/HEPG2_NIPBL_rep2/HEPG2_NIPBL_rep2.merged.bam
symlink_readset_sample_bam.HEPG2_NIPBL_rep2.8844a6ea0c58f45e1040216beeeb4c60.mugqic.done
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

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_9_JOB_ID: symlink_readset_sample_bam.HEPG2_MED1_rep2
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.HEPG2_MED1_rep2
JOB_DEPENDENCIES=$samtools_view_filter_9_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.HEPG2_MED1_rep2.bfa3b92e3a544bb006302eea06596ab3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.HEPG2_MED1_rep2.bfa3b92e3a544bb006302eea06596ab3.mugqic.done'
mkdir -p alignment/HEPG2_MED1_rep2 && \
ln -s -f HI.2449.006.Index_3.HEPG2_ctrl_MED1_MF6_rep1.bam/HI.2449.006.Index_3.HEPG2_ctrl_MED1_MF6_rep1.bam.sorted.filtered.bam alignment/HEPG2_MED1_rep2/HEPG2_MED1_rep2.merged.bam
symlink_readset_sample_bam.HEPG2_MED1_rep2.bfa3b92e3a544bb006302eea06596ab3.mugqic.done
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

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_10_JOB_ID: symlink_readset_sample_bam.HEPG2_BRD4_rep1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.HEPG2_BRD4_rep1
JOB_DEPENDENCIES=$samtools_view_filter_10_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.HEPG2_BRD4_rep1.0e0c141673c627a9e79759fcbb72ada6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.HEPG2_BRD4_rep1.0e0c141673c627a9e79759fcbb72ada6.mugqic.done'
mkdir -p alignment/HEPG2_BRD4_rep1 && \
ln -s -f ENCFF953UBF.fastq/ENCFF953UBF.fastq.sorted.filtered.bam alignment/HEPG2_BRD4_rep1/HEPG2_BRD4_rep1.merged.bam
symlink_readset_sample_bam.HEPG2_BRD4_rep1.0e0c141673c627a9e79759fcbb72ada6.mugqic.done
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

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_11_JOB_ID: symlink_readset_sample_bam.HEPG2_BRD4_rep2
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.HEPG2_BRD4_rep2
JOB_DEPENDENCIES=$samtools_view_filter_11_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.HEPG2_BRD4_rep2.1e318cc9bf2d970b21a46c5e7045d2e3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.HEPG2_BRD4_rep2.1e318cc9bf2d970b21a46c5e7045d2e3.mugqic.done'
mkdir -p alignment/HEPG2_BRD4_rep2 && \
ln -s -f ENCFF946FEX.fastq/ENCFF946FEX.fastq.sorted.filtered.bam alignment/HEPG2_BRD4_rep2/HEPG2_BRD4_rep2.merged.bam
symlink_readset_sample_bam.HEPG2_BRD4_rep2.1e318cc9bf2d970b21a46c5e7045d2e3.mugqic.done
)
picard_merge_sam_files_11_JOB_ID=$(echo "#! /bin/bash
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
echo "$picard_merge_sam_files_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_12_JOB_ID: symlink_readset_sample_bam.HEPG2_BRD4_WCE_rep1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.HEPG2_BRD4_WCE_rep1
JOB_DEPENDENCIES=$samtools_view_filter_12_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.HEPG2_BRD4_WCE_rep1.dab835fe312f5c480433de90873678b5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.HEPG2_BRD4_WCE_rep1.dab835fe312f5c480433de90873678b5.mugqic.done'
mkdir -p alignment/HEPG2_BRD4_WCE_rep1 && \
ln -s -f ENCFF474ZPM.fastq/ENCFF474ZPM.fastq.sorted.filtered.bam alignment/HEPG2_BRD4_WCE_rep1/HEPG2_BRD4_WCE_rep1.merged.bam
symlink_readset_sample_bam.HEPG2_BRD4_WCE_rep1.dab835fe312f5c480433de90873678b5.mugqic.done
)
picard_merge_sam_files_12_JOB_ID=$(echo "#! /bin/bash
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
echo "$picard_merge_sam_files_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_13_JOB_ID: symlink_readset_sample_bam.HEPG2_BRD4_WCE_rep2
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.HEPG2_BRD4_WCE_rep2
JOB_DEPENDENCIES=$samtools_view_filter_13_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.HEPG2_BRD4_WCE_rep2.af857b1562a2e89cb8b32924180948e8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.HEPG2_BRD4_WCE_rep2.af857b1562a2e89cb8b32924180948e8.mugqic.done'
mkdir -p alignment/HEPG2_BRD4_WCE_rep2 && \
ln -s -f ENCFF653JWL.fastq/ENCFF653JWL.fastq.sorted.filtered.bam alignment/HEPG2_BRD4_WCE_rep2/HEPG2_BRD4_WCE_rep2.merged.bam
symlink_readset_sample_bam.HEPG2_BRD4_WCE_rep2.af857b1562a2e89cb8b32924180948e8.mugqic.done
)
picard_merge_sam_files_13_JOB_ID=$(echo "#! /bin/bash
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
echo "$picard_merge_sam_files_13_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# STEP: picard_mark_duplicates
#-------------------------------------------------------------------------------
STEP=picard_mark_duplicates
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_1_JOB_ID: picard_mark_duplicates.HEPG2_WCE_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.HEPG2_WCE_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_1_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.HEPG2_WCE_rep1.d0880c1a815de80bd71657d1d49eae17.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.HEPG2_WCE_rep1.d0880c1a815de80bd71657d1d49eae17.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/HEPG2_WCE_rep1/HEPG2_WCE_rep1.merged.bam \
 OUTPUT=alignment/HEPG2_WCE_rep1/HEPG2_WCE_rep1.sorted.dup.bam \
 METRICS_FILE=alignment/HEPG2_WCE_rep1/HEPG2_WCE_rep1.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.HEPG2_WCE_rep1.d0880c1a815de80bd71657d1d49eae17.mugqic.done
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

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_2_JOB_ID: picard_mark_duplicates.HEPG2_NIPBL_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.HEPG2_NIPBL_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_2_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.HEPG2_NIPBL_rep1.4fd78db1a0484415d46d38cf7accebd7.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.HEPG2_NIPBL_rep1.4fd78db1a0484415d46d38cf7accebd7.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/HEPG2_NIPBL_rep1/HEPG2_NIPBL_rep1.merged.bam \
 OUTPUT=alignment/HEPG2_NIPBL_rep1/HEPG2_NIPBL_rep1.sorted.dup.bam \
 METRICS_FILE=alignment/HEPG2_NIPBL_rep1/HEPG2_NIPBL_rep1.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.HEPG2_NIPBL_rep1.4fd78db1a0484415d46d38cf7accebd7.mugqic.done
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

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_3_JOB_ID: picard_mark_duplicates.HEPG2_CDK9_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.HEPG2_CDK9_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_3_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.HEPG2_CDK9_rep1.69f0dc6e2979f8754f7e09fc1963cffa.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.HEPG2_CDK9_rep1.69f0dc6e2979f8754f7e09fc1963cffa.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/HEPG2_CDK9_rep1/HEPG2_CDK9_rep1.merged.bam \
 OUTPUT=alignment/HEPG2_CDK9_rep1/HEPG2_CDK9_rep1.sorted.dup.bam \
 METRICS_FILE=alignment/HEPG2_CDK9_rep1/HEPG2_CDK9_rep1.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.HEPG2_CDK9_rep1.69f0dc6e2979f8754f7e09fc1963cffa.mugqic.done
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

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_4_JOB_ID: picard_mark_duplicates.HEPG2_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.HEPG2_MED1_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_4_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.HEPG2_MED1_rep1.9ca83f6c9fd04f44a3355e0d1314bf77.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.HEPG2_MED1_rep1.9ca83f6c9fd04f44a3355e0d1314bf77.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/HEPG2_MED1_rep1/HEPG2_MED1_rep1.merged.bam \
 OUTPUT=alignment/HEPG2_MED1_rep1/HEPG2_MED1_rep1.sorted.dup.bam \
 METRICS_FILE=alignment/HEPG2_MED1_rep1/HEPG2_MED1_rep1.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.HEPG2_MED1_rep1.9ca83f6c9fd04f44a3355e0d1314bf77.mugqic.done
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

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_5_JOB_ID: picard_mark_duplicates.HEPG2_SMC1A_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.HEPG2_SMC1A_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_5_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.HEPG2_SMC1A_rep1.a9e8e2428918c9785d861cab9136ef70.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.HEPG2_SMC1A_rep1.a9e8e2428918c9785d861cab9136ef70.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/HEPG2_SMC1A_rep1/HEPG2_SMC1A_rep1.merged.bam \
 OUTPUT=alignment/HEPG2_SMC1A_rep1/HEPG2_SMC1A_rep1.sorted.dup.bam \
 METRICS_FILE=alignment/HEPG2_SMC1A_rep1/HEPG2_SMC1A_rep1.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.HEPG2_SMC1A_rep1.a9e8e2428918c9785d861cab9136ef70.mugqic.done
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

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_6_JOB_ID: picard_mark_duplicates.HEPG2_WCE_rep2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.HEPG2_WCE_rep2
JOB_DEPENDENCIES=$picard_merge_sam_files_6_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.HEPG2_WCE_rep2.4208d08642a70f97e8302eff0d1d3a91.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.HEPG2_WCE_rep2.4208d08642a70f97e8302eff0d1d3a91.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/HEPG2_WCE_rep2/HEPG2_WCE_rep2.merged.bam \
 OUTPUT=alignment/HEPG2_WCE_rep2/HEPG2_WCE_rep2.sorted.dup.bam \
 METRICS_FILE=alignment/HEPG2_WCE_rep2/HEPG2_WCE_rep2.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.HEPG2_WCE_rep2.4208d08642a70f97e8302eff0d1d3a91.mugqic.done
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

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_7_JOB_ID: picard_mark_duplicates.HEPG2_SMC1A_rep2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.HEPG2_SMC1A_rep2
JOB_DEPENDENCIES=$picard_merge_sam_files_7_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.HEPG2_SMC1A_rep2.e993d673ae03b1af60beacf44463575a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.HEPG2_SMC1A_rep2.e993d673ae03b1af60beacf44463575a.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/HEPG2_SMC1A_rep2/HEPG2_SMC1A_rep2.merged.bam \
 OUTPUT=alignment/HEPG2_SMC1A_rep2/HEPG2_SMC1A_rep2.sorted.dup.bam \
 METRICS_FILE=alignment/HEPG2_SMC1A_rep2/HEPG2_SMC1A_rep2.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.HEPG2_SMC1A_rep2.e993d673ae03b1af60beacf44463575a.mugqic.done
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

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_8_JOB_ID: picard_mark_duplicates.HEPG2_NIPBL_rep2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.HEPG2_NIPBL_rep2
JOB_DEPENDENCIES=$picard_merge_sam_files_8_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.HEPG2_NIPBL_rep2.1ac8e1205ab16c691e79fbb700fb02cb.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.HEPG2_NIPBL_rep2.1ac8e1205ab16c691e79fbb700fb02cb.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/HEPG2_NIPBL_rep2/HEPG2_NIPBL_rep2.merged.bam \
 OUTPUT=alignment/HEPG2_NIPBL_rep2/HEPG2_NIPBL_rep2.sorted.dup.bam \
 METRICS_FILE=alignment/HEPG2_NIPBL_rep2/HEPG2_NIPBL_rep2.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.HEPG2_NIPBL_rep2.1ac8e1205ab16c691e79fbb700fb02cb.mugqic.done
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

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_9_JOB_ID: picard_mark_duplicates.HEPG2_MED1_rep2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.HEPG2_MED1_rep2
JOB_DEPENDENCIES=$picard_merge_sam_files_9_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.HEPG2_MED1_rep2.2b832349893d896ea71a36ba00fd2671.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.HEPG2_MED1_rep2.2b832349893d896ea71a36ba00fd2671.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/HEPG2_MED1_rep2/HEPG2_MED1_rep2.merged.bam \
 OUTPUT=alignment/HEPG2_MED1_rep2/HEPG2_MED1_rep2.sorted.dup.bam \
 METRICS_FILE=alignment/HEPG2_MED1_rep2/HEPG2_MED1_rep2.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.HEPG2_MED1_rep2.2b832349893d896ea71a36ba00fd2671.mugqic.done
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

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_10_JOB_ID: picard_mark_duplicates.HEPG2_BRD4_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.HEPG2_BRD4_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_10_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.HEPG2_BRD4_rep1.50a7734eab4601260580d384e78ff60c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.HEPG2_BRD4_rep1.50a7734eab4601260580d384e78ff60c.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/HEPG2_BRD4_rep1/HEPG2_BRD4_rep1.merged.bam \
 OUTPUT=alignment/HEPG2_BRD4_rep1/HEPG2_BRD4_rep1.sorted.dup.bam \
 METRICS_FILE=alignment/HEPG2_BRD4_rep1/HEPG2_BRD4_rep1.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.HEPG2_BRD4_rep1.50a7734eab4601260580d384e78ff60c.mugqic.done
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

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_11_JOB_ID: picard_mark_duplicates.HEPG2_BRD4_rep2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.HEPG2_BRD4_rep2
JOB_DEPENDENCIES=$picard_merge_sam_files_11_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.HEPG2_BRD4_rep2.1b5a3f91b8eef7256d8b1b0baab41b2b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.HEPG2_BRD4_rep2.1b5a3f91b8eef7256d8b1b0baab41b2b.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/HEPG2_BRD4_rep2/HEPG2_BRD4_rep2.merged.bam \
 OUTPUT=alignment/HEPG2_BRD4_rep2/HEPG2_BRD4_rep2.sorted.dup.bam \
 METRICS_FILE=alignment/HEPG2_BRD4_rep2/HEPG2_BRD4_rep2.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.HEPG2_BRD4_rep2.1b5a3f91b8eef7256d8b1b0baab41b2b.mugqic.done
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=48:00:0 --mem=8G -N 1 -n 2 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$picard_mark_duplicates_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_12_JOB_ID: picard_mark_duplicates.HEPG2_BRD4_WCE_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.HEPG2_BRD4_WCE_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_12_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.HEPG2_BRD4_WCE_rep1.357dcfe28b5d472e0bb74940f09650b9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.HEPG2_BRD4_WCE_rep1.357dcfe28b5d472e0bb74940f09650b9.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/HEPG2_BRD4_WCE_rep1/HEPG2_BRD4_WCE_rep1.merged.bam \
 OUTPUT=alignment/HEPG2_BRD4_WCE_rep1/HEPG2_BRD4_WCE_rep1.sorted.dup.bam \
 METRICS_FILE=alignment/HEPG2_BRD4_WCE_rep1/HEPG2_BRD4_WCE_rep1.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.HEPG2_BRD4_WCE_rep1.357dcfe28b5d472e0bb74940f09650b9.mugqic.done
)
picard_mark_duplicates_12_JOB_ID=$(echo "#! /bin/bash
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
echo "$picard_mark_duplicates_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_13_JOB_ID: picard_mark_duplicates.HEPG2_BRD4_WCE_rep2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.HEPG2_BRD4_WCE_rep2
JOB_DEPENDENCIES=$picard_merge_sam_files_13_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.HEPG2_BRD4_WCE_rep2.e950c60bd00861a8d1b018d2853a5d83.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.HEPG2_BRD4_WCE_rep2.e950c60bd00861a8d1b018d2853a5d83.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/HEPG2_BRD4_WCE_rep2/HEPG2_BRD4_WCE_rep2.merged.bam \
 OUTPUT=alignment/HEPG2_BRD4_WCE_rep2/HEPG2_BRD4_WCE_rep2.sorted.dup.bam \
 METRICS_FILE=alignment/HEPG2_BRD4_WCE_rep2/HEPG2_BRD4_WCE_rep2.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.HEPG2_BRD4_WCE_rep2.e950c60bd00861a8d1b018d2853a5d83.mugqic.done
)
picard_mark_duplicates_13_JOB_ID=$(echo "#! /bin/bash
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
echo "$picard_mark_duplicates_13_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_14_JOB_ID: picard_mark_duplicates_report
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates_report
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_2_JOB_ID:$picard_mark_duplicates_3_JOB_ID:$picard_mark_duplicates_4_JOB_ID:$picard_mark_duplicates_5_JOB_ID:$picard_mark_duplicates_6_JOB_ID:$picard_mark_duplicates_7_JOB_ID:$picard_mark_duplicates_8_JOB_ID:$picard_mark_duplicates_9_JOB_ID:$picard_mark_duplicates_10_JOB_ID:$picard_mark_duplicates_11_JOB_ID:$picard_mark_duplicates_12_JOB_ID:$picard_mark_duplicates_13_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates_report.0d6756c3b0268d30a7f6ff885850a9fd.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates_report.0d6756c3b0268d30a7f6ff885850a9fd.mugqic.done'
mkdir -p report && \
cp \
  /cvmfs/soft.mugqic/CentOS6/software/genpipes/genpipes-3.1.2/bfx/report/ChipSeq.picard_mark_duplicates.md \
  report/ChipSeq.picard_mark_duplicates.md
picard_mark_duplicates_report.0d6756c3b0268d30a7f6ff885850a9fd.mugqic.done
)
picard_mark_duplicates_14_JOB_ID=$(echo "#! /bin/bash
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
echo "$picard_mark_duplicates_14_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# STEP: metrics
#-------------------------------------------------------------------------------
STEP=metrics
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: metrics_1_JOB_ID: metrics.flagstat
#-------------------------------------------------------------------------------
JOB_NAME=metrics.flagstat
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_2_JOB_ID:$picard_mark_duplicates_3_JOB_ID:$picard_mark_duplicates_4_JOB_ID:$picard_mark_duplicates_5_JOB_ID:$picard_mark_duplicates_6_JOB_ID:$picard_mark_duplicates_7_JOB_ID:$picard_mark_duplicates_8_JOB_ID:$picard_mark_duplicates_9_JOB_ID:$picard_mark_duplicates_10_JOB_ID:$picard_mark_duplicates_11_JOB_ID:$picard_mark_duplicates_12_JOB_ID:$picard_mark_duplicates_13_JOB_ID
JOB_DONE=job_output/metrics/metrics.flagstat.619d0ac05ecd2113a4f9a312fe6b157a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'metrics.flagstat.619d0ac05ecd2113a4f9a312fe6b157a.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools flagstat \
  alignment/HEPG2_WCE_rep1/HEPG2_WCE_rep1.sorted.dup.bam \
  > alignment/HEPG2_WCE_rep1/HEPG2_WCE_rep1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/HEPG2_NIPBL_rep1/HEPG2_NIPBL_rep1.sorted.dup.bam \
  > alignment/HEPG2_NIPBL_rep1/HEPG2_NIPBL_rep1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/HEPG2_CDK9_rep1/HEPG2_CDK9_rep1.sorted.dup.bam \
  > alignment/HEPG2_CDK9_rep1/HEPG2_CDK9_rep1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/HEPG2_MED1_rep1/HEPG2_MED1_rep1.sorted.dup.bam \
  > alignment/HEPG2_MED1_rep1/HEPG2_MED1_rep1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/HEPG2_SMC1A_rep1/HEPG2_SMC1A_rep1.sorted.dup.bam \
  > alignment/HEPG2_SMC1A_rep1/HEPG2_SMC1A_rep1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/HEPG2_WCE_rep2/HEPG2_WCE_rep2.sorted.dup.bam \
  > alignment/HEPG2_WCE_rep2/HEPG2_WCE_rep2.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/HEPG2_SMC1A_rep2/HEPG2_SMC1A_rep2.sorted.dup.bam \
  > alignment/HEPG2_SMC1A_rep2/HEPG2_SMC1A_rep2.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/HEPG2_NIPBL_rep2/HEPG2_NIPBL_rep2.sorted.dup.bam \
  > alignment/HEPG2_NIPBL_rep2/HEPG2_NIPBL_rep2.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/HEPG2_MED1_rep2/HEPG2_MED1_rep2.sorted.dup.bam \
  > alignment/HEPG2_MED1_rep2/HEPG2_MED1_rep2.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/HEPG2_BRD4_rep1/HEPG2_BRD4_rep1.sorted.dup.bam \
  > alignment/HEPG2_BRD4_rep1/HEPG2_BRD4_rep1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/HEPG2_BRD4_rep2/HEPG2_BRD4_rep2.sorted.dup.bam \
  > alignment/HEPG2_BRD4_rep2/HEPG2_BRD4_rep2.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/HEPG2_BRD4_WCE_rep1/HEPG2_BRD4_WCE_rep1.sorted.dup.bam \
  > alignment/HEPG2_BRD4_WCE_rep1/HEPG2_BRD4_WCE_rep1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/HEPG2_BRD4_WCE_rep2/HEPG2_BRD4_WCE_rep2.sorted.dup.bam \
  > alignment/HEPG2_BRD4_WCE_rep2/HEPG2_BRD4_WCE_rep2.sorted.dup.bam.flagstat
metrics.flagstat.619d0ac05ecd2113a4f9a312fe6b157a.mugqic.done
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

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: metrics_2_JOB_ID: metrics_report
#-------------------------------------------------------------------------------
JOB_NAME=metrics_report
JOB_DEPENDENCIES=$metrics_1_JOB_ID
JOB_DONE=job_output/metrics/metrics_report.52393feda35c09442d7470ea02357fd6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'metrics_report.52393feda35c09442d7470ea02357fd6.mugqic.done'
module load mugqic/pandoc/1.15.2 && \
for sample in HEPG2_WCE_rep1 HEPG2_NIPBL_rep1 HEPG2_CDK9_rep1 HEPG2_MED1_rep1 HEPG2_SMC1A_rep1 HEPG2_WCE_rep2 HEPG2_SMC1A_rep2 HEPG2_NIPBL_rep2 HEPG2_MED1_rep2 HEPG2_BRD4_rep1 HEPG2_BRD4_rep2 HEPG2_BRD4_WCE_rep1 HEPG2_BRD4_WCE_rep2
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
  --template /cvmfs/soft.mugqic/CentOS6/software/genpipes/genpipes-3.1.2/bfx/report/ChipSeq.metrics.md \
  --variable trim_mem_sample_table="$trim_mem_sample_table" \
  /cvmfs/soft.mugqic/CentOS6/software/genpipes/genpipes-3.1.2/bfx/report/ChipSeq.metrics.md \
  > report/ChipSeq.metrics.md

metrics_report.52393feda35c09442d7470ea02357fd6.mugqic.done
)
metrics_2_JOB_ID=$(echo "#! /bin/bash
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
echo "$metrics_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# STEP: homer_make_tag_directory
#-------------------------------------------------------------------------------
STEP=homer_make_tag_directory
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_1_JOB_ID: homer_make_tag_directory.HEPG2_WCE_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.HEPG2_WCE_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.HEPG2_WCE_rep1.69f60c3356f7d8d5429db96a20964300.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.HEPG2_WCE_rep1.69f60c3356f7d8d5429db96a20964300.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/HEPG2_WCE_rep1 \
            alignment/HEPG2_WCE_rep1/HEPG2_WCE_rep1.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.HEPG2_WCE_rep1.69f60c3356f7d8d5429db96a20964300.mugqic.done
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

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_2_JOB_ID: homer_make_tag_directory.HEPG2_NIPBL_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.HEPG2_NIPBL_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_2_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.HEPG2_NIPBL_rep1.c0e719612b06df9a1c29359f23aa75a0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.HEPG2_NIPBL_rep1.c0e719612b06df9a1c29359f23aa75a0.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/HEPG2_NIPBL_rep1 \
            alignment/HEPG2_NIPBL_rep1/HEPG2_NIPBL_rep1.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.HEPG2_NIPBL_rep1.c0e719612b06df9a1c29359f23aa75a0.mugqic.done
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

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_3_JOB_ID: homer_make_tag_directory.HEPG2_CDK9_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.HEPG2_CDK9_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_3_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.HEPG2_CDK9_rep1.c53aa96b0b0902480ad8b6cac234dd29.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.HEPG2_CDK9_rep1.c53aa96b0b0902480ad8b6cac234dd29.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/HEPG2_CDK9_rep1 \
            alignment/HEPG2_CDK9_rep1/HEPG2_CDK9_rep1.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.HEPG2_CDK9_rep1.c53aa96b0b0902480ad8b6cac234dd29.mugqic.done
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

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_4_JOB_ID: homer_make_tag_directory.HEPG2_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.HEPG2_MED1_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_4_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.HEPG2_MED1_rep1.46006a8b420a8503b08aa003367f2cc5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.HEPG2_MED1_rep1.46006a8b420a8503b08aa003367f2cc5.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/HEPG2_MED1_rep1 \
            alignment/HEPG2_MED1_rep1/HEPG2_MED1_rep1.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.HEPG2_MED1_rep1.46006a8b420a8503b08aa003367f2cc5.mugqic.done
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

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_5_JOB_ID: homer_make_tag_directory.HEPG2_SMC1A_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.HEPG2_SMC1A_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_5_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.HEPG2_SMC1A_rep1.d540559f13e91dfafe2b61ef1fd3c5e7.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.HEPG2_SMC1A_rep1.d540559f13e91dfafe2b61ef1fd3c5e7.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/HEPG2_SMC1A_rep1 \
            alignment/HEPG2_SMC1A_rep1/HEPG2_SMC1A_rep1.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.HEPG2_SMC1A_rep1.d540559f13e91dfafe2b61ef1fd3c5e7.mugqic.done
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

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_6_JOB_ID: homer_make_tag_directory.HEPG2_WCE_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.HEPG2_WCE_rep2
JOB_DEPENDENCIES=$picard_mark_duplicates_6_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.HEPG2_WCE_rep2.bb6f5d68b15d5bd5950dd72648655ec2.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.HEPG2_WCE_rep2.bb6f5d68b15d5bd5950dd72648655ec2.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/HEPG2_WCE_rep2 \
            alignment/HEPG2_WCE_rep2/HEPG2_WCE_rep2.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.HEPG2_WCE_rep2.bb6f5d68b15d5bd5950dd72648655ec2.mugqic.done
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

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_7_JOB_ID: homer_make_tag_directory.HEPG2_SMC1A_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.HEPG2_SMC1A_rep2
JOB_DEPENDENCIES=$picard_mark_duplicates_7_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.HEPG2_SMC1A_rep2.fa62e88bf17dd5cc38e4e362de99fba6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.HEPG2_SMC1A_rep2.fa62e88bf17dd5cc38e4e362de99fba6.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/HEPG2_SMC1A_rep2 \
            alignment/HEPG2_SMC1A_rep2/HEPG2_SMC1A_rep2.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.HEPG2_SMC1A_rep2.fa62e88bf17dd5cc38e4e362de99fba6.mugqic.done
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

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_8_JOB_ID: homer_make_tag_directory.HEPG2_NIPBL_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.HEPG2_NIPBL_rep2
JOB_DEPENDENCIES=$picard_mark_duplicates_8_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.HEPG2_NIPBL_rep2.aa54b5137f482c1b0352d0256204bb62.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.HEPG2_NIPBL_rep2.aa54b5137f482c1b0352d0256204bb62.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/HEPG2_NIPBL_rep2 \
            alignment/HEPG2_NIPBL_rep2/HEPG2_NIPBL_rep2.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.HEPG2_NIPBL_rep2.aa54b5137f482c1b0352d0256204bb62.mugqic.done
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

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_9_JOB_ID: homer_make_tag_directory.HEPG2_MED1_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.HEPG2_MED1_rep2
JOB_DEPENDENCIES=$picard_mark_duplicates_9_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.HEPG2_MED1_rep2.7de322ad9023916a2c8e442890967f5a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.HEPG2_MED1_rep2.7de322ad9023916a2c8e442890967f5a.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/HEPG2_MED1_rep2 \
            alignment/HEPG2_MED1_rep2/HEPG2_MED1_rep2.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.HEPG2_MED1_rep2.7de322ad9023916a2c8e442890967f5a.mugqic.done
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

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_10_JOB_ID: homer_make_tag_directory.HEPG2_BRD4_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.HEPG2_BRD4_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_10_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.HEPG2_BRD4_rep1.eb2b54b85e846ae49c28ff504e1339c3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.HEPG2_BRD4_rep1.eb2b54b85e846ae49c28ff504e1339c3.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/HEPG2_BRD4_rep1 \
            alignment/HEPG2_BRD4_rep1/HEPG2_BRD4_rep1.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.HEPG2_BRD4_rep1.eb2b54b85e846ae49c28ff504e1339c3.mugqic.done
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

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_11_JOB_ID: homer_make_tag_directory.HEPG2_BRD4_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.HEPG2_BRD4_rep2
JOB_DEPENDENCIES=$picard_mark_duplicates_11_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.HEPG2_BRD4_rep2.709380be87caf0639b99486f4ddb8e15.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.HEPG2_BRD4_rep2.709380be87caf0639b99486f4ddb8e15.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/HEPG2_BRD4_rep2 \
            alignment/HEPG2_BRD4_rep2/HEPG2_BRD4_rep2.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.HEPG2_BRD4_rep2.709380be87caf0639b99486f4ddb8e15.mugqic.done
)
homer_make_tag_directory_11_JOB_ID=$(echo "#! /bin/bash
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
echo "$homer_make_tag_directory_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_12_JOB_ID: homer_make_tag_directory.HEPG2_BRD4_WCE_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.HEPG2_BRD4_WCE_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_12_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.HEPG2_BRD4_WCE_rep1.fb657a6f1359e3e3ed4b354c9d8fe4bb.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.HEPG2_BRD4_WCE_rep1.fb657a6f1359e3e3ed4b354c9d8fe4bb.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/HEPG2_BRD4_WCE_rep1 \
            alignment/HEPG2_BRD4_WCE_rep1/HEPG2_BRD4_WCE_rep1.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.HEPG2_BRD4_WCE_rep1.fb657a6f1359e3e3ed4b354c9d8fe4bb.mugqic.done
)
homer_make_tag_directory_12_JOB_ID=$(echo "#! /bin/bash
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
echo "$homer_make_tag_directory_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_13_JOB_ID: homer_make_tag_directory.HEPG2_BRD4_WCE_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.HEPG2_BRD4_WCE_rep2
JOB_DEPENDENCIES=$picard_mark_duplicates_13_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.HEPG2_BRD4_WCE_rep2.dc550df112ddca804ceeb18c29ff0071.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.HEPG2_BRD4_WCE_rep2.dc550df112ddca804ceeb18c29ff0071.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/HEPG2_BRD4_WCE_rep2 \
            alignment/HEPG2_BRD4_WCE_rep2/HEPG2_BRD4_WCE_rep2.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.HEPG2_BRD4_WCE_rep2.dc550df112ddca804ceeb18c29ff0071.mugqic.done
)
homer_make_tag_directory_13_JOB_ID=$(echo "#! /bin/bash
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
echo "$homer_make_tag_directory_13_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# STEP: qc_metrics
#-------------------------------------------------------------------------------
STEP=qc_metrics
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: qc_metrics_1_JOB_ID: qc_plots_R
#-------------------------------------------------------------------------------
JOB_NAME=qc_plots_R
JOB_DEPENDENCIES=$homer_make_tag_directory_1_JOB_ID:$homer_make_tag_directory_2_JOB_ID:$homer_make_tag_directory_3_JOB_ID:$homer_make_tag_directory_4_JOB_ID:$homer_make_tag_directory_5_JOB_ID:$homer_make_tag_directory_6_JOB_ID:$homer_make_tag_directory_7_JOB_ID:$homer_make_tag_directory_8_JOB_ID:$homer_make_tag_directory_9_JOB_ID:$homer_make_tag_directory_10_JOB_ID:$homer_make_tag_directory_11_JOB_ID:$homer_make_tag_directory_12_JOB_ID:$homer_make_tag_directory_13_JOB_ID
JOB_DONE=job_output/qc_metrics/qc_plots_R.a394b513977fba299a749c1b645ad824.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'qc_plots_R.a394b513977fba299a749c1b645ad824.mugqic.done'
module load mugqic/mugqic_tools/2.1.9 mugqic/R_Bioconductor/3.5.0_3.7 && \
mkdir -p graphs && \
Rscript $R_TOOLS/chipSeqGenerateQCMetrics.R \
  ../../raw/chip-seq/design.txt \
  /scratch/chris11/sb_cofactor/HEPG2/output/chip-pipeline-GRCh38 && \
cp /cvmfs/soft.mugqic/CentOS6/software/genpipes/genpipes-3.1.2/bfx/report/ChipSeq.qc_metrics.md report/ChipSeq.qc_metrics.md && \
for sample in HEPG2_WCE_rep1 HEPG2_NIPBL_rep1 HEPG2_CDK9_rep1 HEPG2_MED1_rep1 HEPG2_SMC1A_rep1 HEPG2_WCE_rep2 HEPG2_SMC1A_rep2 HEPG2_NIPBL_rep2 HEPG2_MED1_rep2 HEPG2_BRD4_rep1 HEPG2_BRD4_rep2 HEPG2_BRD4_WCE_rep1 HEPG2_BRD4_WCE_rep2
do
  cp --parents graphs/${sample}_QC_Metrics.ps report/
  convert -rotate 90 graphs/${sample}_QC_Metrics.ps report/graphs/${sample}_QC_Metrics.png
  echo -e "----

![QC Metrics for Sample $sample ([download high-res image](graphs/${sample}_QC_Metrics.ps))](graphs/${sample}_QC_Metrics.png)
" \
  >> report/ChipSeq.qc_metrics.md
done
qc_plots_R.a394b513977fba299a749c1b645ad824.mugqic.done
)
qc_metrics_1_JOB_ID=$(echo "#! /bin/bash
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
echo "$qc_metrics_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# STEP: macs2_callpeak
#-------------------------------------------------------------------------------
STEP=macs2_callpeak
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_1_JOB_ID: macs2_callpeak.HEPG2_NIPBL_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.HEPG2_NIPBL_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_2_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.HEPG2_NIPBL_rep1.996bbfc77fe7c8bf6b1cc11dcc524b18.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.HEPG2_NIPBL_rep1.996bbfc77fe7c8bf6b1cc11dcc524b18.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/HEPG2_NIPBL_rep1 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/HEPG2_NIPBL_rep1/HEPG2_NIPBL_rep1.sorted.dup.bam \
  --control \
  alignment/HEPG2_WCE_rep1/HEPG2_WCE_rep1.sorted.dup.bam \
  --name peak_call/HEPG2_NIPBL_rep1/HEPG2_NIPBL_rep1 \
  >& peak_call/HEPG2_NIPBL_rep1/HEPG2_NIPBL_rep1.diag.macs.out
macs2_callpeak.HEPG2_NIPBL_rep1.996bbfc77fe7c8bf6b1cc11dcc524b18.mugqic.done
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

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_2_JOB_ID: macs2_callpeak_bigBed.HEPG2_NIPBL_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.HEPG2_NIPBL_rep1
JOB_DEPENDENCIES=$macs2_callpeak_1_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.HEPG2_NIPBL_rep1.a8564ac9b74295c8ade443b117617f65.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.HEPG2_NIPBL_rep1.a8564ac9b74295c8ade443b117617f65.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/HEPG2_NIPBL_rep1/HEPG2_NIPBL_rep1_peaks.narrowPeak > peak_call/HEPG2_NIPBL_rep1/HEPG2_NIPBL_rep1_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/HEPG2_NIPBL_rep1/HEPG2_NIPBL_rep1_peaks.narrowPeak.bed \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/HEPG2_NIPBL_rep1/HEPG2_NIPBL_rep1_peaks.narrowPeak.bb
macs2_callpeak_bigBed.HEPG2_NIPBL_rep1.a8564ac9b74295c8ade443b117617f65.mugqic.done
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

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_3_JOB_ID: macs2_callpeak.HEPG2_CDK9_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.HEPG2_CDK9_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_3_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.HEPG2_CDK9_rep1.2ec3b0b7acce98829deefa09a16331ab.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.HEPG2_CDK9_rep1.2ec3b0b7acce98829deefa09a16331ab.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/HEPG2_CDK9_rep1 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/HEPG2_CDK9_rep1/HEPG2_CDK9_rep1.sorted.dup.bam \
  --control \
  alignment/HEPG2_WCE_rep1/HEPG2_WCE_rep1.sorted.dup.bam \
  --name peak_call/HEPG2_CDK9_rep1/HEPG2_CDK9_rep1 \
  >& peak_call/HEPG2_CDK9_rep1/HEPG2_CDK9_rep1.diag.macs.out
macs2_callpeak.HEPG2_CDK9_rep1.2ec3b0b7acce98829deefa09a16331ab.mugqic.done
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

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_4_JOB_ID: macs2_callpeak_bigBed.HEPG2_CDK9_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.HEPG2_CDK9_rep1
JOB_DEPENDENCIES=$macs2_callpeak_3_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.HEPG2_CDK9_rep1.d516ca0ad86fe26d0548a67bd6ac78a8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.HEPG2_CDK9_rep1.d516ca0ad86fe26d0548a67bd6ac78a8.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/HEPG2_CDK9_rep1/HEPG2_CDK9_rep1_peaks.narrowPeak > peak_call/HEPG2_CDK9_rep1/HEPG2_CDK9_rep1_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/HEPG2_CDK9_rep1/HEPG2_CDK9_rep1_peaks.narrowPeak.bed \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/HEPG2_CDK9_rep1/HEPG2_CDK9_rep1_peaks.narrowPeak.bb
macs2_callpeak_bigBed.HEPG2_CDK9_rep1.d516ca0ad86fe26d0548a67bd6ac78a8.mugqic.done
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

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_5_JOB_ID: macs2_callpeak.HEPG2_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.HEPG2_MED1_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_4_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.HEPG2_MED1_rep1.ed2cfc802c8c94391dc1a8c02894348c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.HEPG2_MED1_rep1.ed2cfc802c8c94391dc1a8c02894348c.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/HEPG2_MED1_rep1 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/HEPG2_MED1_rep1/HEPG2_MED1_rep1.sorted.dup.bam \
  --control \
  alignment/HEPG2_WCE_rep1/HEPG2_WCE_rep1.sorted.dup.bam \
  --name peak_call/HEPG2_MED1_rep1/HEPG2_MED1_rep1 \
  >& peak_call/HEPG2_MED1_rep1/HEPG2_MED1_rep1.diag.macs.out
macs2_callpeak.HEPG2_MED1_rep1.ed2cfc802c8c94391dc1a8c02894348c.mugqic.done
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

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_6_JOB_ID: macs2_callpeak_bigBed.HEPG2_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.HEPG2_MED1_rep1
JOB_DEPENDENCIES=$macs2_callpeak_5_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.HEPG2_MED1_rep1.ee5c81733c6edcd4cd7ea6e437b166fd.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.HEPG2_MED1_rep1.ee5c81733c6edcd4cd7ea6e437b166fd.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/HEPG2_MED1_rep1/HEPG2_MED1_rep1_peaks.narrowPeak > peak_call/HEPG2_MED1_rep1/HEPG2_MED1_rep1_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/HEPG2_MED1_rep1/HEPG2_MED1_rep1_peaks.narrowPeak.bed \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/HEPG2_MED1_rep1/HEPG2_MED1_rep1_peaks.narrowPeak.bb
macs2_callpeak_bigBed.HEPG2_MED1_rep1.ee5c81733c6edcd4cd7ea6e437b166fd.mugqic.done
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

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_7_JOB_ID: macs2_callpeak.HEPG2_SMC1A_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.HEPG2_SMC1A_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_5_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.HEPG2_SMC1A_rep1.73196d53928769f9bdc2826dcd84bdba.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.HEPG2_SMC1A_rep1.73196d53928769f9bdc2826dcd84bdba.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/HEPG2_SMC1A_rep1 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/HEPG2_SMC1A_rep1/HEPG2_SMC1A_rep1.sorted.dup.bam \
  --control \
  alignment/HEPG2_WCE_rep1/HEPG2_WCE_rep1.sorted.dup.bam \
  --name peak_call/HEPG2_SMC1A_rep1/HEPG2_SMC1A_rep1 \
  >& peak_call/HEPG2_SMC1A_rep1/HEPG2_SMC1A_rep1.diag.macs.out
macs2_callpeak.HEPG2_SMC1A_rep1.73196d53928769f9bdc2826dcd84bdba.mugqic.done
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

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_8_JOB_ID: macs2_callpeak_bigBed.HEPG2_SMC1A_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.HEPG2_SMC1A_rep1
JOB_DEPENDENCIES=$macs2_callpeak_7_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.HEPG2_SMC1A_rep1.359cfb036f70b80802b2ff4f77a8ba15.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.HEPG2_SMC1A_rep1.359cfb036f70b80802b2ff4f77a8ba15.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/HEPG2_SMC1A_rep1/HEPG2_SMC1A_rep1_peaks.narrowPeak > peak_call/HEPG2_SMC1A_rep1/HEPG2_SMC1A_rep1_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/HEPG2_SMC1A_rep1/HEPG2_SMC1A_rep1_peaks.narrowPeak.bed \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/HEPG2_SMC1A_rep1/HEPG2_SMC1A_rep1_peaks.narrowPeak.bb
macs2_callpeak_bigBed.HEPG2_SMC1A_rep1.359cfb036f70b80802b2ff4f77a8ba15.mugqic.done
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

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_9_JOB_ID: macs2_callpeak.HEPG2_SMC1A_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.HEPG2_SMC1A_rep2
JOB_DEPENDENCIES=$picard_mark_duplicates_6_JOB_ID:$picard_mark_duplicates_7_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.HEPG2_SMC1A_rep2.9436cb4640a199569c048f93b99b8057.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.HEPG2_SMC1A_rep2.9436cb4640a199569c048f93b99b8057.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/HEPG2_SMC1A_rep2 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/HEPG2_SMC1A_rep2/HEPG2_SMC1A_rep2.sorted.dup.bam \
  --control \
  alignment/HEPG2_WCE_rep2/HEPG2_WCE_rep2.sorted.dup.bam \
  --name peak_call/HEPG2_SMC1A_rep2/HEPG2_SMC1A_rep2 \
  >& peak_call/HEPG2_SMC1A_rep2/HEPG2_SMC1A_rep2.diag.macs.out
macs2_callpeak.HEPG2_SMC1A_rep2.9436cb4640a199569c048f93b99b8057.mugqic.done
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

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_10_JOB_ID: macs2_callpeak_bigBed.HEPG2_SMC1A_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.HEPG2_SMC1A_rep2
JOB_DEPENDENCIES=$macs2_callpeak_9_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.HEPG2_SMC1A_rep2.af4d33b89aef78aba67c4909a37b0e7b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.HEPG2_SMC1A_rep2.af4d33b89aef78aba67c4909a37b0e7b.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/HEPG2_SMC1A_rep2/HEPG2_SMC1A_rep2_peaks.narrowPeak > peak_call/HEPG2_SMC1A_rep2/HEPG2_SMC1A_rep2_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/HEPG2_SMC1A_rep2/HEPG2_SMC1A_rep2_peaks.narrowPeak.bed \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/HEPG2_SMC1A_rep2/HEPG2_SMC1A_rep2_peaks.narrowPeak.bb
macs2_callpeak_bigBed.HEPG2_SMC1A_rep2.af4d33b89aef78aba67c4909a37b0e7b.mugqic.done
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

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_11_JOB_ID: macs2_callpeak.HEPG2_NIPBL_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.HEPG2_NIPBL_rep2
JOB_DEPENDENCIES=$picard_mark_duplicates_6_JOB_ID:$picard_mark_duplicates_8_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.HEPG2_NIPBL_rep2.1041e8a7eebc608de5a9c43e19cdeb45.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.HEPG2_NIPBL_rep2.1041e8a7eebc608de5a9c43e19cdeb45.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/HEPG2_NIPBL_rep2 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/HEPG2_NIPBL_rep2/HEPG2_NIPBL_rep2.sorted.dup.bam \
  --control \
  alignment/HEPG2_WCE_rep2/HEPG2_WCE_rep2.sorted.dup.bam \
  --name peak_call/HEPG2_NIPBL_rep2/HEPG2_NIPBL_rep2 \
  >& peak_call/HEPG2_NIPBL_rep2/HEPG2_NIPBL_rep2.diag.macs.out
macs2_callpeak.HEPG2_NIPBL_rep2.1041e8a7eebc608de5a9c43e19cdeb45.mugqic.done
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

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_12_JOB_ID: macs2_callpeak_bigBed.HEPG2_NIPBL_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.HEPG2_NIPBL_rep2
JOB_DEPENDENCIES=$macs2_callpeak_11_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.HEPG2_NIPBL_rep2.d0f2db7c4446227650485165129c5e65.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.HEPG2_NIPBL_rep2.d0f2db7c4446227650485165129c5e65.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/HEPG2_NIPBL_rep2/HEPG2_NIPBL_rep2_peaks.narrowPeak > peak_call/HEPG2_NIPBL_rep2/HEPG2_NIPBL_rep2_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/HEPG2_NIPBL_rep2/HEPG2_NIPBL_rep2_peaks.narrowPeak.bed \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/HEPG2_NIPBL_rep2/HEPG2_NIPBL_rep2_peaks.narrowPeak.bb
macs2_callpeak_bigBed.HEPG2_NIPBL_rep2.d0f2db7c4446227650485165129c5e65.mugqic.done
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

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_13_JOB_ID: macs2_callpeak.HEPG2_MED1_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.HEPG2_MED1_rep2
JOB_DEPENDENCIES=$picard_mark_duplicates_6_JOB_ID:$picard_mark_duplicates_9_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.HEPG2_MED1_rep2.977d8477ee23f48603af15f0ad3cfcba.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.HEPG2_MED1_rep2.977d8477ee23f48603af15f0ad3cfcba.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/HEPG2_MED1_rep2 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/HEPG2_MED1_rep2/HEPG2_MED1_rep2.sorted.dup.bam \
  --control \
  alignment/HEPG2_WCE_rep2/HEPG2_WCE_rep2.sorted.dup.bam \
  --name peak_call/HEPG2_MED1_rep2/HEPG2_MED1_rep2 \
  >& peak_call/HEPG2_MED1_rep2/HEPG2_MED1_rep2.diag.macs.out
macs2_callpeak.HEPG2_MED1_rep2.977d8477ee23f48603af15f0ad3cfcba.mugqic.done
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

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_14_JOB_ID: macs2_callpeak_bigBed.HEPG2_MED1_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.HEPG2_MED1_rep2
JOB_DEPENDENCIES=$macs2_callpeak_13_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.HEPG2_MED1_rep2.6532b9b6e3247747a935943d2b8f0829.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.HEPG2_MED1_rep2.6532b9b6e3247747a935943d2b8f0829.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/HEPG2_MED1_rep2/HEPG2_MED1_rep2_peaks.narrowPeak > peak_call/HEPG2_MED1_rep2/HEPG2_MED1_rep2_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/HEPG2_MED1_rep2/HEPG2_MED1_rep2_peaks.narrowPeak.bed \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/HEPG2_MED1_rep2/HEPG2_MED1_rep2_peaks.narrowPeak.bb
macs2_callpeak_bigBed.HEPG2_MED1_rep2.6532b9b6e3247747a935943d2b8f0829.mugqic.done
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

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_15_JOB_ID: macs2_callpeak.HEPG2_BRD4_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.HEPG2_BRD4_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_10_JOB_ID:$picard_mark_duplicates_12_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.HEPG2_BRD4_rep1.d4d1d7d9162bd57911c7a313f81b3e90.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.HEPG2_BRD4_rep1.d4d1d7d9162bd57911c7a313f81b3e90.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/HEPG2_BRD4_rep1 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/HEPG2_BRD4_rep1/HEPG2_BRD4_rep1.sorted.dup.bam \
  --control \
  alignment/HEPG2_BRD4_WCE_rep1/HEPG2_BRD4_WCE_rep1.sorted.dup.bam \
  --name peak_call/HEPG2_BRD4_rep1/HEPG2_BRD4_rep1 \
  >& peak_call/HEPG2_BRD4_rep1/HEPG2_BRD4_rep1.diag.macs.out
macs2_callpeak.HEPG2_BRD4_rep1.d4d1d7d9162bd57911c7a313f81b3e90.mugqic.done
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

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_16_JOB_ID: macs2_callpeak_bigBed.HEPG2_BRD4_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.HEPG2_BRD4_rep1
JOB_DEPENDENCIES=$macs2_callpeak_15_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.HEPG2_BRD4_rep1.773cd7fe7beaa2b6c253947488ccb8cf.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.HEPG2_BRD4_rep1.773cd7fe7beaa2b6c253947488ccb8cf.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/HEPG2_BRD4_rep1/HEPG2_BRD4_rep1_peaks.narrowPeak > peak_call/HEPG2_BRD4_rep1/HEPG2_BRD4_rep1_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/HEPG2_BRD4_rep1/HEPG2_BRD4_rep1_peaks.narrowPeak.bed \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/HEPG2_BRD4_rep1/HEPG2_BRD4_rep1_peaks.narrowPeak.bb
macs2_callpeak_bigBed.HEPG2_BRD4_rep1.773cd7fe7beaa2b6c253947488ccb8cf.mugqic.done
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

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_17_JOB_ID: macs2_callpeak.HEPG2_BRD4_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.HEPG2_BRD4_rep2
JOB_DEPENDENCIES=$picard_mark_duplicates_11_JOB_ID:$picard_mark_duplicates_13_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.HEPG2_BRD4_rep2.0f8751f623f69fcd31fccd4a4fbc72e8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.HEPG2_BRD4_rep2.0f8751f623f69fcd31fccd4a4fbc72e8.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/HEPG2_BRD4_rep2 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/HEPG2_BRD4_rep2/HEPG2_BRD4_rep2.sorted.dup.bam \
  --control \
  alignment/HEPG2_BRD4_WCE_rep2/HEPG2_BRD4_WCE_rep2.sorted.dup.bam \
  --name peak_call/HEPG2_BRD4_rep2/HEPG2_BRD4_rep2 \
  >& peak_call/HEPG2_BRD4_rep2/HEPG2_BRD4_rep2.diag.macs.out
macs2_callpeak.HEPG2_BRD4_rep2.0f8751f623f69fcd31fccd4a4fbc72e8.mugqic.done
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

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_18_JOB_ID: macs2_callpeak_bigBed.HEPG2_BRD4_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.HEPG2_BRD4_rep2
JOB_DEPENDENCIES=$macs2_callpeak_17_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.HEPG2_BRD4_rep2.7058646fc70b5e0d9525ed0d3e7efe63.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.HEPG2_BRD4_rep2.7058646fc70b5e0d9525ed0d3e7efe63.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/HEPG2_BRD4_rep2/HEPG2_BRD4_rep2_peaks.narrowPeak > peak_call/HEPG2_BRD4_rep2/HEPG2_BRD4_rep2_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/HEPG2_BRD4_rep2/HEPG2_BRD4_rep2_peaks.narrowPeak.bed \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/HEPG2_BRD4_rep2/HEPG2_BRD4_rep2_peaks.narrowPeak.bb
macs2_callpeak_bigBed.HEPG2_BRD4_rep2.7058646fc70b5e0d9525ed0d3e7efe63.mugqic.done
)
macs2_callpeak_18_JOB_ID=$(echo "#! /bin/bash
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
echo "$macs2_callpeak_18_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_19_JOB_ID: macs2_callpeak_report
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_report
JOB_DEPENDENCIES=$macs2_callpeak_1_JOB_ID:$macs2_callpeak_3_JOB_ID:$macs2_callpeak_5_JOB_ID:$macs2_callpeak_7_JOB_ID:$macs2_callpeak_9_JOB_ID:$macs2_callpeak_11_JOB_ID:$macs2_callpeak_13_JOB_ID:$macs2_callpeak_15_JOB_ID:$macs2_callpeak_17_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_report.06607885865c410bbecd02edf0ad70e9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_report.06607885865c410bbecd02edf0ad70e9.mugqic.done'
mkdir -p report && \
cp /cvmfs/soft.mugqic/CentOS6/software/genpipes/genpipes-3.1.2/bfx/report/ChipSeq.macs2_callpeak.md report/ && \
for contrast in HEPG2_NIPBL_rep1 HEPG2_CDK9_rep1 HEPG2_MED1_rep1 HEPG2_SMC1A_rep1 HEPG2_SMC1A_rep2 HEPG2_NIPBL_rep2 HEPG2_MED1_rep2 HEPG2_BRD4_rep1 HEPG2_BRD4_rep2
do
  cp -a --parents peak_call/$contrast/ report/ && \
  echo -e "* [Peak Calls File for Design $contrast](peak_call/$contrast/${contrast}_peaks.xls)" \
  >> report/ChipSeq.macs2_callpeak.md
done
macs2_callpeak_report.06607885865c410bbecd02edf0ad70e9.mugqic.done
)
macs2_callpeak_19_JOB_ID=$(echo "#! /bin/bash
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
echo "$macs2_callpeak_19_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# Call home with pipeline statistics
#-------------------------------------------------------------------------------
LOG_MD5=$(echo $USER-'206.12.124.6-ChipSeq-HEPG2_BRD4_rep2.ENCFF946FEX.fastq,HEPG2_BRD4_rep1.ENCFF953UBF.fastq,HEPG2_CDK9_rep1.HI.1308.004.Index_22.HEPG2_Cdk9_ChIP67_5G_rep1.bam,HEPG2_NIPBL_rep2.HI.2449.006.Index_1.HEPG2_ctrl_NIPBL_MF6_rep1.bam,HEPG2_NIPBL_rep1.HI.1308.004.Index_19.HEPG2_NIPBL_ChIP67_1G_rep1.bam,HEPG2_WCE_rep1.HI.1308.004.Index_18.HEPG2_WCE_rep1.bam,HEPG2_WCE_rep2.HI.2449.006.Index_8.HEPG2_ctrl_WCE_MF6_rep1.bam,HEPG2_SMC1A_rep1.HI.1308.004.Index_21.HEPG2_Smc1_ChIP67_3G_rep1.bam,HEPG2_MED1_rep2.HI.2449.006.Index_3.HEPG2_ctrl_MED1_MF6_rep1.bam,HEPG2_MED1_rep1.HI.1308.004.Index_20.HEPG2_Med1_ChIP67_2G_rep1.bam,HEPG2_SMC1A_rep2.HI.2449.006.Index_9.HEPG2_ctrl_SMC1A_MF6_rep1.bam,HEPG2_BRD4_WCE_rep1.ENCFF474ZPM.fastq,HEPG2_BRD4_WCE_rep2.ENCFF653JWL.fastq' | md5sum | awk '{ print $1 }')
wget "http://mugqic.hpc.mcgill.ca/cgi-bin/pipeline.cgi?hostname=cedar5.cedar.computecanada.ca&ip=206.12.124.6&pipeline=ChipSeq&steps=picard_sam_to_fastq,trimmomatic,merge_trimmomatic_stats,bwa_mem_picard_sort_sam,samtools_view_filter,picard_merge_sam_files,picard_mark_duplicates,metrics,homer_make_tag_directory,qc_metrics,macs2_callpeak&samples=13&md5=$LOG_MD5" --quiet --output-document=/dev/null

