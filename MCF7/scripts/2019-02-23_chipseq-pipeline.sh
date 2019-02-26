#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

#-------------------------------------------------------------------------------
# ChipSeq SlurmScheduler Job Submission Bash script
# Version: 3.1.2-beta
# Created on: 2019-02-23T12:44:18
# Steps:
#   picard_sam_to_fastq: 22 jobs
#   trimmomatic: 30 jobs
#   merge_trimmomatic_stats: 1 job
#   bwa_mem_picard_sort_sam: 31 jobs
#   samtools_view_filter: 31 jobs
#   picard_merge_sam_files: 30 jobs
#   picard_mark_duplicates: 31 jobs
#   metrics: 2 jobs
#   homer_make_tag_directory: 30 jobs
#   qc_metrics: 1 job
#   macs2_callpeak: 41 jobs
#   TOTAL: 250 jobs
#-------------------------------------------------------------------------------

OUTPUT_DIR=/scratch/chris11/sb_cofactor/MCF7/output/chip-pipeline-GRCh38
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
# JOB: picard_sam_to_fastq_1_JOB_ID: picard_sam_to_fastq.HI.1943.008.Index_2.MCF7_CTRL_WCE_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.1943.008.Index_2.MCF7_CTRL_WCE_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.1943.008.Index_2.MCF7_CTRL_WCE_ChIP8_GB_rep1.bam.2a01d875d999f1abc7762443f25214a5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.1943.008.Index_2.MCF7_CTRL_WCE_ChIP8_GB_rep1.bam.2a01d875d999f1abc7762443f25214a5.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.1943.008.Index_2.MCF7_CTRL_WCE_ChIP8_GB_rep1.bam \
 FASTQ=/scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.1943.008.Index_2.MCF7_CTRL_WCE_ChIP8_GB_rep1.single.fastq.gz
picard_sam_to_fastq.HI.1943.008.Index_2.MCF7_CTRL_WCE_ChIP8_GB_rep1.bam.2a01d875d999f1abc7762443f25214a5.mugqic.done
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
# JOB: picard_sam_to_fastq_2_JOB_ID: picard_sam_to_fastq.HI.1997.004.Index_1.MCF7_CTRL_NIPBL_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.1997.004.Index_1.MCF7_CTRL_NIPBL_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.1997.004.Index_1.MCF7_CTRL_NIPBL_ChIP8_GB_rep1.bam.7253513a6b61fb7c94d72915923cfe25.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.1997.004.Index_1.MCF7_CTRL_NIPBL_ChIP8_GB_rep1.bam.7253513a6b61fb7c94d72915923cfe25.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.1997.004.Index_1.MCF7_CTRL_NIPBL_ChIP8_GB_rep1.bam \
 FASTQ=/scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.1997.004.Index_1.MCF7_CTRL_NIPBL_ChIP8_GB_rep1.single.fastq.gz
picard_sam_to_fastq.HI.1997.004.Index_1.MCF7_CTRL_NIPBL_ChIP8_GB_rep1.bam.7253513a6b61fb7c94d72915923cfe25.mugqic.done
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
# JOB: picard_sam_to_fastq_3_JOB_ID: picard_sam_to_fastq.HI.1997.004.Index_4.MCF7_CTRL_SMC1_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.1997.004.Index_4.MCF7_CTRL_SMC1_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.1997.004.Index_4.MCF7_CTRL_SMC1_ChIP8_GB_rep1.bam.d19eb16095dea86eb21001d6a45b474c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.1997.004.Index_4.MCF7_CTRL_SMC1_ChIP8_GB_rep1.bam.d19eb16095dea86eb21001d6a45b474c.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.1997.004.Index_4.MCF7_CTRL_SMC1_ChIP8_GB_rep1.bam \
 FASTQ=/scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.1997.004.Index_4.MCF7_CTRL_SMC1_ChIP8_GB_rep1.single.fastq.gz
picard_sam_to_fastq.HI.1997.004.Index_4.MCF7_CTRL_SMC1_ChIP8_GB_rep1.bam.d19eb16095dea86eb21001d6a45b474c.mugqic.done
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
# JOB: picard_sam_to_fastq_4_JOB_ID: picard_sam_to_fastq.HI.1943.008.Index_7.MCF7_CTRL_MED1_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.1943.008.Index_7.MCF7_CTRL_MED1_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.1943.008.Index_7.MCF7_CTRL_MED1_ChIP8_GB_rep1.bam.f3ac6e974513ce2cb369eef5c73b77a4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.1943.008.Index_7.MCF7_CTRL_MED1_ChIP8_GB_rep1.bam.f3ac6e974513ce2cb369eef5c73b77a4.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.1943.008.Index_7.MCF7_CTRL_MED1_ChIP8_GB_rep1.bam \
 FASTQ=/scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.1943.008.Index_7.MCF7_CTRL_MED1_ChIP8_GB_rep1.single.fastq.gz
picard_sam_to_fastq.HI.1943.008.Index_7.MCF7_CTRL_MED1_ChIP8_GB_rep1.bam.f3ac6e974513ce2cb369eef5c73b77a4.mugqic.done
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
# JOB: picard_sam_to_fastq_5_JOB_ID: picard_sam_to_fastq.HI.1997.005.Index_10.MCF7_CTRL_POL2_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.1997.005.Index_10.MCF7_CTRL_POL2_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.1997.005.Index_10.MCF7_CTRL_POL2_ChIP8_GB_rep1.bam.73c2cf90b47e264be7277f5c791442d6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.1997.005.Index_10.MCF7_CTRL_POL2_ChIP8_GB_rep1.bam.73c2cf90b47e264be7277f5c791442d6.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.1997.005.Index_10.MCF7_CTRL_POL2_ChIP8_GB_rep1.bam \
 FASTQ=/scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.1997.005.Index_10.MCF7_CTRL_POL2_ChIP8_GB_rep1.single.fastq.gz
picard_sam_to_fastq.HI.1997.005.Index_10.MCF7_CTRL_POL2_ChIP8_GB_rep1.bam.73c2cf90b47e264be7277f5c791442d6.mugqic.done
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
# JOB: picard_sam_to_fastq_6_JOB_ID: picard_sam_to_fastq.HI.1997.005.Index_25.MCF7_CTRL_ERA_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.1997.005.Index_25.MCF7_CTRL_ERA_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.1997.005.Index_25.MCF7_CTRL_ERA_ChIP8_GB_rep1.bam.3342f93eb0088dc47d175b2f033ee23f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.1997.005.Index_25.MCF7_CTRL_ERA_ChIP8_GB_rep1.bam.3342f93eb0088dc47d175b2f033ee23f.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.1997.005.Index_25.MCF7_CTRL_ERA_ChIP8_GB_rep1.bam \
 FASTQ=/scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.1997.005.Index_25.MCF7_CTRL_ERA_ChIP8_GB_rep1.single.fastq.gz
picard_sam_to_fastq.HI.1997.005.Index_25.MCF7_CTRL_ERA_ChIP8_GB_rep1.bam.3342f93eb0088dc47d175b2f033ee23f.mugqic.done
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
# JOB: picard_sam_to_fastq_7_JOB_ID: picard_sam_to_fastq.HI.1943.008.Index_13.MCF7_E2_WCE_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.1943.008.Index_13.MCF7_E2_WCE_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.1943.008.Index_13.MCF7_E2_WCE_ChIP8_GB_rep1.bam.66cd486a0652415a599b8174876b250b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.1943.008.Index_13.MCF7_E2_WCE_ChIP8_GB_rep1.bam.66cd486a0652415a599b8174876b250b.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.1943.008.Index_13.MCF7_E2_WCE_ChIP8_GB_rep1.bam \
 FASTQ=/scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.1943.008.Index_13.MCF7_E2_WCE_ChIP8_GB_rep1.single.fastq.gz
picard_sam_to_fastq.HI.1943.008.Index_13.MCF7_E2_WCE_ChIP8_GB_rep1.bam.66cd486a0652415a599b8174876b250b.mugqic.done
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
# JOB: picard_sam_to_fastq_8_JOB_ID: picard_sam_to_fastq.HI.1997.004.Index_3.MCF7_E2_NIPBL_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.1997.004.Index_3.MCF7_E2_NIPBL_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.1997.004.Index_3.MCF7_E2_NIPBL_ChIP8_GB_rep1.bam.0ab4b5dffd9e8ce7a572369956c52a5d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.1997.004.Index_3.MCF7_E2_NIPBL_ChIP8_GB_rep1.bam.0ab4b5dffd9e8ce7a572369956c52a5d.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.1997.004.Index_3.MCF7_E2_NIPBL_ChIP8_GB_rep1.bam \
 FASTQ=/scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.1997.004.Index_3.MCF7_E2_NIPBL_ChIP8_GB_rep1.single.fastq.gz
picard_sam_to_fastq.HI.1997.004.Index_3.MCF7_E2_NIPBL_ChIP8_GB_rep1.bam.0ab4b5dffd9e8ce7a572369956c52a5d.mugqic.done
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
# JOB: picard_sam_to_fastq_9_JOB_ID: picard_sam_to_fastq.HI.1997.004.Index_5.MCF7_E2_SMC1_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.1997.004.Index_5.MCF7_E2_SMC1_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.1997.004.Index_5.MCF7_E2_SMC1_ChIP8_GB_rep1.bam.571385ec11d495f3dec53e1b8b3a4bf2.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.1997.004.Index_5.MCF7_E2_SMC1_ChIP8_GB_rep1.bam.571385ec11d495f3dec53e1b8b3a4bf2.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.1997.004.Index_5.MCF7_E2_SMC1_ChIP8_GB_rep1.bam \
 FASTQ=/scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.1997.004.Index_5.MCF7_E2_SMC1_ChIP8_GB_rep1.single.fastq.gz
picard_sam_to_fastq.HI.1997.004.Index_5.MCF7_E2_SMC1_ChIP8_GB_rep1.bam.571385ec11d495f3dec53e1b8b3a4bf2.mugqic.done
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
# JOB: picard_sam_to_fastq_10_JOB_ID: picard_sam_to_fastq.HI.1943.008.Index_18.MCF7_E2_MED1_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.1943.008.Index_18.MCF7_E2_MED1_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.1943.008.Index_18.MCF7_E2_MED1_ChIP8_GB_rep1.bam.7bd755a9263ec89cd0bb2359ab11ed82.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.1943.008.Index_18.MCF7_E2_MED1_ChIP8_GB_rep1.bam.7bd755a9263ec89cd0bb2359ab11ed82.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.1943.008.Index_18.MCF7_E2_MED1_ChIP8_GB_rep1.bam \
 FASTQ=/scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.1943.008.Index_18.MCF7_E2_MED1_ChIP8_GB_rep1.single.fastq.gz
picard_sam_to_fastq.HI.1943.008.Index_18.MCF7_E2_MED1_ChIP8_GB_rep1.bam.7bd755a9263ec89cd0bb2359ab11ed82.mugqic.done
)
picard_sam_to_fastq_10_JOB_ID=$(echo "#! /bin/bash
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
echo "$picard_sam_to_fastq_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_11_JOB_ID: picard_sam_to_fastq.HI.1997.005.Index_11.MCF7_E2_POL2_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.1997.005.Index_11.MCF7_E2_POL2_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.1997.005.Index_11.MCF7_E2_POL2_ChIP8_GB_rep1.bam.66c7666ad54056173e4f5768a87c734d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.1997.005.Index_11.MCF7_E2_POL2_ChIP8_GB_rep1.bam.66c7666ad54056173e4f5768a87c734d.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.1997.005.Index_11.MCF7_E2_POL2_ChIP8_GB_rep1.bam \
 FASTQ=/scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.1997.005.Index_11.MCF7_E2_POL2_ChIP8_GB_rep1.single.fastq.gz
picard_sam_to_fastq.HI.1997.005.Index_11.MCF7_E2_POL2_ChIP8_GB_rep1.bam.66c7666ad54056173e4f5768a87c734d.mugqic.done
)
picard_sam_to_fastq_11_JOB_ID=$(echo "#! /bin/bash
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
echo "$picard_sam_to_fastq_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_12_JOB_ID: picard_sam_to_fastq.HI.1997.005.Index_21.MCF7_E2_ERA_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.1997.005.Index_21.MCF7_E2_ERA_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.1997.005.Index_21.MCF7_E2_ERA_ChIP8_GB_rep1.bam.f512385a91c8cb9c6221036acbeab785.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.1997.005.Index_21.MCF7_E2_ERA_ChIP8_GB_rep1.bam.f512385a91c8cb9c6221036acbeab785.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.1997.005.Index_21.MCF7_E2_ERA_ChIP8_GB_rep1.bam \
 FASTQ=/scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.1997.005.Index_21.MCF7_E2_ERA_ChIP8_GB_rep1.single.fastq.gz
picard_sam_to_fastq.HI.1997.005.Index_21.MCF7_E2_ERA_ChIP8_GB_rep1.bam.f512385a91c8cb9c6221036acbeab785.mugqic.done
)
picard_sam_to_fastq_12_JOB_ID=$(echo "#! /bin/bash
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
echo "$picard_sam_to_fastq_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_13_JOB_ID: picard_sam_to_fastq.HI.2791.003.Index_5.MCF7_WCE_C_rep2.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.2791.003.Index_5.MCF7_WCE_C_rep2.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.2791.003.Index_5.MCF7_WCE_C_rep2.bam.d1796096bb1f84975a2e11924fc25778.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.2791.003.Index_5.MCF7_WCE_C_rep2.bam.d1796096bb1f84975a2e11924fc25778.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.2791.003.Index_5.MCF7_WCE_C_rep2.bam \
 FASTQ=/scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.2791.003.Index_5.MCF7_WCE_C_rep2.single.fastq.gz
picard_sam_to_fastq.HI.2791.003.Index_5.MCF7_WCE_C_rep2.bam.d1796096bb1f84975a2e11924fc25778.mugqic.done
)
picard_sam_to_fastq_13_JOB_ID=$(echo "#! /bin/bash
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
echo "$picard_sam_to_fastq_13_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_14_JOB_ID: picard_sam_to_fastq.HI.2791.003.Index_9.MCF7_1_NIPBL_C_rep2.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.2791.003.Index_9.MCF7_1_NIPBL_C_rep2.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.2791.003.Index_9.MCF7_1_NIPBL_C_rep2.bam.1dad0898d3164fae593741f1a987e83b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.2791.003.Index_9.MCF7_1_NIPBL_C_rep2.bam.1dad0898d3164fae593741f1a987e83b.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.2791.003.Index_9.MCF7_1_NIPBL_C_rep2.bam \
 FASTQ=/scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.2791.003.Index_9.MCF7_1_NIPBL_C_rep2.single.fastq.gz
picard_sam_to_fastq.HI.2791.003.Index_9.MCF7_1_NIPBL_C_rep2.bam.1dad0898d3164fae593741f1a987e83b.mugqic.done
)
picard_sam_to_fastq_14_JOB_ID=$(echo "#! /bin/bash
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
echo "$picard_sam_to_fastq_14_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_15_JOB_ID: picard_sam_to_fastq.HI.2791.003.Index_11.MCF7_4_Smc1_C_rep2.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.2791.003.Index_11.MCF7_4_Smc1_C_rep2.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.2791.003.Index_11.MCF7_4_Smc1_C_rep2.bam.3f7b5dd4cd58d4e2bef74b32f33013b4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.2791.003.Index_11.MCF7_4_Smc1_C_rep2.bam.3f7b5dd4cd58d4e2bef74b32f33013b4.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.2791.003.Index_11.MCF7_4_Smc1_C_rep2.bam \
 FASTQ=/scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.2791.003.Index_11.MCF7_4_Smc1_C_rep2.single.fastq.gz
picard_sam_to_fastq.HI.2791.003.Index_11.MCF7_4_Smc1_C_rep2.bam.3f7b5dd4cd58d4e2bef74b32f33013b4.mugqic.done
)
picard_sam_to_fastq_15_JOB_ID=$(echo "#! /bin/bash
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
echo "$picard_sam_to_fastq_15_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_16_JOB_ID: picard_sam_to_fastq.HI.2791.003.Index_8.MCF7_2_MED1_C_rep2.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.2791.003.Index_8.MCF7_2_MED1_C_rep2.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.2791.003.Index_8.MCF7_2_MED1_C_rep2.bam.10e8e9844168fd3968aa1fa6485c5cf8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.2791.003.Index_8.MCF7_2_MED1_C_rep2.bam.10e8e9844168fd3968aa1fa6485c5cf8.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.2791.003.Index_8.MCF7_2_MED1_C_rep2.bam \
 FASTQ=/scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.2791.003.Index_8.MCF7_2_MED1_C_rep2.single.fastq.gz
picard_sam_to_fastq.HI.2791.003.Index_8.MCF7_2_MED1_C_rep2.bam.10e8e9844168fd3968aa1fa6485c5cf8.mugqic.done
)
picard_sam_to_fastq_16_JOB_ID=$(echo "#! /bin/bash
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
echo "$picard_sam_to_fastq_16_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_17_JOB_ID: picard_sam_to_fastq.HI.2791.003.Index_10.MCF7_3_Era_C_rep2.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.2791.003.Index_10.MCF7_3_Era_C_rep2.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.2791.003.Index_10.MCF7_3_Era_C_rep2.bam.796ff1c3872767d5d63f8ff36503231d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.2791.003.Index_10.MCF7_3_Era_C_rep2.bam.796ff1c3872767d5d63f8ff36503231d.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.2791.003.Index_10.MCF7_3_Era_C_rep2.bam \
 FASTQ=/scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.2791.003.Index_10.MCF7_3_Era_C_rep2.single.fastq.gz
picard_sam_to_fastq.HI.2791.003.Index_10.MCF7_3_Era_C_rep2.bam.796ff1c3872767d5d63f8ff36503231d.mugqic.done
)
picard_sam_to_fastq_17_JOB_ID=$(echo "#! /bin/bash
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
echo "$picard_sam_to_fastq_17_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_18_JOB_ID: picard_sam_to_fastq.HI.2791.003.Index_4.MCF7_WCE_E_rep2.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.2791.003.Index_4.MCF7_WCE_E_rep2.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.2791.003.Index_4.MCF7_WCE_E_rep2.bam.7bba4cdbe53414308b5b92a36fcae74a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.2791.003.Index_4.MCF7_WCE_E_rep2.bam.7bba4cdbe53414308b5b92a36fcae74a.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.2791.003.Index_4.MCF7_WCE_E_rep2.bam \
 FASTQ=/scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.2791.003.Index_4.MCF7_WCE_E_rep2.single.fastq.gz
picard_sam_to_fastq.HI.2791.003.Index_4.MCF7_WCE_E_rep2.bam.7bba4cdbe53414308b5b92a36fcae74a.mugqic.done
)
picard_sam_to_fastq_18_JOB_ID=$(echo "#! /bin/bash
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
echo "$picard_sam_to_fastq_18_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_19_JOB_ID: picard_sam_to_fastq.HI.2791.003.Index_12.MCF7_1_NIPBL_E_rep2.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.2791.003.Index_12.MCF7_1_NIPBL_E_rep2.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.2791.003.Index_12.MCF7_1_NIPBL_E_rep2.bam.c347ced6edd59a5dcb12531aeaa0f8c6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.2791.003.Index_12.MCF7_1_NIPBL_E_rep2.bam.c347ced6edd59a5dcb12531aeaa0f8c6.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.2791.003.Index_12.MCF7_1_NIPBL_E_rep2.bam \
 FASTQ=/scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.2791.003.Index_12.MCF7_1_NIPBL_E_rep2.single.fastq.gz
picard_sam_to_fastq.HI.2791.003.Index_12.MCF7_1_NIPBL_E_rep2.bam.c347ced6edd59a5dcb12531aeaa0f8c6.mugqic.done
)
picard_sam_to_fastq_19_JOB_ID=$(echo "#! /bin/bash
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
echo "$picard_sam_to_fastq_19_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_20_JOB_ID: picard_sam_to_fastq.HI.2791.003.Index_3.MCF7_4_Smc1_E_rep2.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.2791.003.Index_3.MCF7_4_Smc1_E_rep2.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.2791.003.Index_3.MCF7_4_Smc1_E_rep2.bam.b8e581cf8cde79422280a5e4e4f09fa1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.2791.003.Index_3.MCF7_4_Smc1_E_rep2.bam.b8e581cf8cde79422280a5e4e4f09fa1.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.2791.003.Index_3.MCF7_4_Smc1_E_rep2.bam \
 FASTQ=/scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.2791.003.Index_3.MCF7_4_Smc1_E_rep2.single.fastq.gz
picard_sam_to_fastq.HI.2791.003.Index_3.MCF7_4_Smc1_E_rep2.bam.b8e581cf8cde79422280a5e4e4f09fa1.mugqic.done
)
picard_sam_to_fastq_20_JOB_ID=$(echo "#! /bin/bash
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
echo "$picard_sam_to_fastq_20_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_21_JOB_ID: picard_sam_to_fastq.HI.2791.003.Index_19.MCF7_2_MED1_E_rep2.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.2791.003.Index_19.MCF7_2_MED1_E_rep2.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.2791.003.Index_19.MCF7_2_MED1_E_rep2.bam.998569ad4b13ad1674b03115b77cd9e8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.2791.003.Index_19.MCF7_2_MED1_E_rep2.bam.998569ad4b13ad1674b03115b77cd9e8.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.2791.003.Index_19.MCF7_2_MED1_E_rep2.bam \
 FASTQ=/scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.2791.003.Index_19.MCF7_2_MED1_E_rep2.single.fastq.gz
picard_sam_to_fastq.HI.2791.003.Index_19.MCF7_2_MED1_E_rep2.bam.998569ad4b13ad1674b03115b77cd9e8.mugqic.done
)
picard_sam_to_fastq_21_JOB_ID=$(echo "#! /bin/bash
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
echo "$picard_sam_to_fastq_21_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_22_JOB_ID: picard_sam_to_fastq.HI.2791.003.Index_1.MCF7_3_Era_E_rep2.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.2791.003.Index_1.MCF7_3_Era_E_rep2.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.2791.003.Index_1.MCF7_3_Era_E_rep2.bam.708d94f3a890fe335413cc96c99b64b9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.2791.003.Index_1.MCF7_3_Era_E_rep2.bam.708d94f3a890fe335413cc96c99b64b9.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.2791.003.Index_1.MCF7_3_Era_E_rep2.bam \
 FASTQ=/scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.2791.003.Index_1.MCF7_3_Era_E_rep2.single.fastq.gz
picard_sam_to_fastq.HI.2791.003.Index_1.MCF7_3_Era_E_rep2.bam.708d94f3a890fe335413cc96c99b64b9.mugqic.done
)
picard_sam_to_fastq_22_JOB_ID=$(echo "#! /bin/bash
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
echo "$picard_sam_to_fastq_22_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# STEP: trimmomatic
#-------------------------------------------------------------------------------
STEP=trimmomatic
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: trimmomatic_1_JOB_ID: trimmomatic.HI.1943.008.Index_2.MCF7_CTRL_WCE_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.1943.008.Index_2.MCF7_CTRL_WCE_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_1_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.1943.008.Index_2.MCF7_CTRL_WCE_ChIP8_GB_rep1.bam.f569060a08a06a1e99d905ba5e457919.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.1943.008.Index_2.MCF7_CTRL_WCE_ChIP8_GB_rep1.bam.f569060a08a06a1e99d905ba5e457919.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/MCF7_CTRL_WCE_rep1 && \
`cat > trim/MCF7_CTRL_WCE_rep1/HI.1943.008.Index_2.MCF7_CTRL_WCE_ChIP8_GB_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.1943.008.Index_2.MCF7_CTRL_WCE_ChIP8_GB_rep1.single.fastq.gz \
  trim/MCF7_CTRL_WCE_rep1/HI.1943.008.Index_2.MCF7_CTRL_WCE_ChIP8_GB_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/MCF7_CTRL_WCE_rep1/HI.1943.008.Index_2.MCF7_CTRL_WCE_ChIP8_GB_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/MCF7_CTRL_WCE_rep1/HI.1943.008.Index_2.MCF7_CTRL_WCE_ChIP8_GB_rep1.bam.trim.log
trimmomatic.HI.1943.008.Index_2.MCF7_CTRL_WCE_ChIP8_GB_rep1.bam.f569060a08a06a1e99d905ba5e457919.mugqic.done
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
# JOB: trimmomatic_2_JOB_ID: trimmomatic.HI.1997.004.Index_1.MCF7_CTRL_NIPBL_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.1997.004.Index_1.MCF7_CTRL_NIPBL_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_2_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.1997.004.Index_1.MCF7_CTRL_NIPBL_ChIP8_GB_rep1.bam.3394efca60da2deab687fca6349f3feb.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.1997.004.Index_1.MCF7_CTRL_NIPBL_ChIP8_GB_rep1.bam.3394efca60da2deab687fca6349f3feb.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/MCF7_CTRL_NIPBL_rep1 && \
`cat > trim/MCF7_CTRL_NIPBL_rep1/HI.1997.004.Index_1.MCF7_CTRL_NIPBL_ChIP8_GB_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.1997.004.Index_1.MCF7_CTRL_NIPBL_ChIP8_GB_rep1.single.fastq.gz \
  trim/MCF7_CTRL_NIPBL_rep1/HI.1997.004.Index_1.MCF7_CTRL_NIPBL_ChIP8_GB_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/MCF7_CTRL_NIPBL_rep1/HI.1997.004.Index_1.MCF7_CTRL_NIPBL_ChIP8_GB_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/MCF7_CTRL_NIPBL_rep1/HI.1997.004.Index_1.MCF7_CTRL_NIPBL_ChIP8_GB_rep1.bam.trim.log
trimmomatic.HI.1997.004.Index_1.MCF7_CTRL_NIPBL_ChIP8_GB_rep1.bam.3394efca60da2deab687fca6349f3feb.mugqic.done
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
# JOB: trimmomatic_3_JOB_ID: trimmomatic.HI.1997.004.Index_4.MCF7_CTRL_SMC1_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.1997.004.Index_4.MCF7_CTRL_SMC1_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_3_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.1997.004.Index_4.MCF7_CTRL_SMC1_ChIP8_GB_rep1.bam.43d92dc55e596e6e0422110b69df1e90.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.1997.004.Index_4.MCF7_CTRL_SMC1_ChIP8_GB_rep1.bam.43d92dc55e596e6e0422110b69df1e90.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/MCF7_CTRL_SMC1A_rep1 && \
`cat > trim/MCF7_CTRL_SMC1A_rep1/HI.1997.004.Index_4.MCF7_CTRL_SMC1_ChIP8_GB_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.1997.004.Index_4.MCF7_CTRL_SMC1_ChIP8_GB_rep1.single.fastq.gz \
  trim/MCF7_CTRL_SMC1A_rep1/HI.1997.004.Index_4.MCF7_CTRL_SMC1_ChIP8_GB_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/MCF7_CTRL_SMC1A_rep1/HI.1997.004.Index_4.MCF7_CTRL_SMC1_ChIP8_GB_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/MCF7_CTRL_SMC1A_rep1/HI.1997.004.Index_4.MCF7_CTRL_SMC1_ChIP8_GB_rep1.bam.trim.log
trimmomatic.HI.1997.004.Index_4.MCF7_CTRL_SMC1_ChIP8_GB_rep1.bam.43d92dc55e596e6e0422110b69df1e90.mugqic.done
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
# JOB: trimmomatic_4_JOB_ID: trimmomatic.HI.1943.008.Index_7.MCF7_CTRL_MED1_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.1943.008.Index_7.MCF7_CTRL_MED1_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_4_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.1943.008.Index_7.MCF7_CTRL_MED1_ChIP8_GB_rep1.bam.02f22ebe4eead9d490a53057ebee82ef.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.1943.008.Index_7.MCF7_CTRL_MED1_ChIP8_GB_rep1.bam.02f22ebe4eead9d490a53057ebee82ef.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/MCF7_CTRL_MED1_rep1 && \
`cat > trim/MCF7_CTRL_MED1_rep1/HI.1943.008.Index_7.MCF7_CTRL_MED1_ChIP8_GB_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.1943.008.Index_7.MCF7_CTRL_MED1_ChIP8_GB_rep1.single.fastq.gz \
  trim/MCF7_CTRL_MED1_rep1/HI.1943.008.Index_7.MCF7_CTRL_MED1_ChIP8_GB_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/MCF7_CTRL_MED1_rep1/HI.1943.008.Index_7.MCF7_CTRL_MED1_ChIP8_GB_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/MCF7_CTRL_MED1_rep1/HI.1943.008.Index_7.MCF7_CTRL_MED1_ChIP8_GB_rep1.bam.trim.log
trimmomatic.HI.1943.008.Index_7.MCF7_CTRL_MED1_ChIP8_GB_rep1.bam.02f22ebe4eead9d490a53057ebee82ef.mugqic.done
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
# JOB: trimmomatic_5_JOB_ID: trimmomatic.HI.1997.005.Index_10.MCF7_CTRL_POL2_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.1997.005.Index_10.MCF7_CTRL_POL2_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_5_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.1997.005.Index_10.MCF7_CTRL_POL2_ChIP8_GB_rep1.bam.1d416d9def7e97cb7433102d92e9bf77.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.1997.005.Index_10.MCF7_CTRL_POL2_ChIP8_GB_rep1.bam.1d416d9def7e97cb7433102d92e9bf77.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/MCF7_CTRL_POL2_rep1 && \
`cat > trim/MCF7_CTRL_POL2_rep1/HI.1997.005.Index_10.MCF7_CTRL_POL2_ChIP8_GB_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.1997.005.Index_10.MCF7_CTRL_POL2_ChIP8_GB_rep1.single.fastq.gz \
  trim/MCF7_CTRL_POL2_rep1/HI.1997.005.Index_10.MCF7_CTRL_POL2_ChIP8_GB_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/MCF7_CTRL_POL2_rep1/HI.1997.005.Index_10.MCF7_CTRL_POL2_ChIP8_GB_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/MCF7_CTRL_POL2_rep1/HI.1997.005.Index_10.MCF7_CTRL_POL2_ChIP8_GB_rep1.bam.trim.log
trimmomatic.HI.1997.005.Index_10.MCF7_CTRL_POL2_ChIP8_GB_rep1.bam.1d416d9def7e97cb7433102d92e9bf77.mugqic.done
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
# JOB: trimmomatic_6_JOB_ID: trimmomatic.HI.1997.005.Index_25.MCF7_CTRL_ERA_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.1997.005.Index_25.MCF7_CTRL_ERA_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_6_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.1997.005.Index_25.MCF7_CTRL_ERA_ChIP8_GB_rep1.bam.88addc681abb1890a2dd0887221fa105.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.1997.005.Index_25.MCF7_CTRL_ERA_ChIP8_GB_rep1.bam.88addc681abb1890a2dd0887221fa105.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/MCF7_CTRL_ERA_rep1 && \
`cat > trim/MCF7_CTRL_ERA_rep1/HI.1997.005.Index_25.MCF7_CTRL_ERA_ChIP8_GB_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.1997.005.Index_25.MCF7_CTRL_ERA_ChIP8_GB_rep1.single.fastq.gz \
  trim/MCF7_CTRL_ERA_rep1/HI.1997.005.Index_25.MCF7_CTRL_ERA_ChIP8_GB_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/MCF7_CTRL_ERA_rep1/HI.1997.005.Index_25.MCF7_CTRL_ERA_ChIP8_GB_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/MCF7_CTRL_ERA_rep1/HI.1997.005.Index_25.MCF7_CTRL_ERA_ChIP8_GB_rep1.bam.trim.log
trimmomatic.HI.1997.005.Index_25.MCF7_CTRL_ERA_ChIP8_GB_rep1.bam.88addc681abb1890a2dd0887221fa105.mugqic.done
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
# JOB: trimmomatic_7_JOB_ID: trimmomatic.HI.1943.008.Index_13.MCF7_E2_WCE_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.1943.008.Index_13.MCF7_E2_WCE_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_7_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.1943.008.Index_13.MCF7_E2_WCE_ChIP8_GB_rep1.bam.37510c3a8a1d2a0d3e53bef54a71783b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.1943.008.Index_13.MCF7_E2_WCE_ChIP8_GB_rep1.bam.37510c3a8a1d2a0d3e53bef54a71783b.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/MCF7_E2_WCE_rep1 && \
`cat > trim/MCF7_E2_WCE_rep1/HI.1943.008.Index_13.MCF7_E2_WCE_ChIP8_GB_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.1943.008.Index_13.MCF7_E2_WCE_ChIP8_GB_rep1.single.fastq.gz \
  trim/MCF7_E2_WCE_rep1/HI.1943.008.Index_13.MCF7_E2_WCE_ChIP8_GB_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/MCF7_E2_WCE_rep1/HI.1943.008.Index_13.MCF7_E2_WCE_ChIP8_GB_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/MCF7_E2_WCE_rep1/HI.1943.008.Index_13.MCF7_E2_WCE_ChIP8_GB_rep1.bam.trim.log
trimmomatic.HI.1943.008.Index_13.MCF7_E2_WCE_ChIP8_GB_rep1.bam.37510c3a8a1d2a0d3e53bef54a71783b.mugqic.done
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
# JOB: trimmomatic_8_JOB_ID: trimmomatic.HI.1997.004.Index_3.MCF7_E2_NIPBL_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.1997.004.Index_3.MCF7_E2_NIPBL_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_8_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.1997.004.Index_3.MCF7_E2_NIPBL_ChIP8_GB_rep1.bam.84f963840a7b7aac1eb8bade14e20df1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.1997.004.Index_3.MCF7_E2_NIPBL_ChIP8_GB_rep1.bam.84f963840a7b7aac1eb8bade14e20df1.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/MCF7_E2_NIPBL_rep1 && \
`cat > trim/MCF7_E2_NIPBL_rep1/HI.1997.004.Index_3.MCF7_E2_NIPBL_ChIP8_GB_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.1997.004.Index_3.MCF7_E2_NIPBL_ChIP8_GB_rep1.single.fastq.gz \
  trim/MCF7_E2_NIPBL_rep1/HI.1997.004.Index_3.MCF7_E2_NIPBL_ChIP8_GB_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/MCF7_E2_NIPBL_rep1/HI.1997.004.Index_3.MCF7_E2_NIPBL_ChIP8_GB_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/MCF7_E2_NIPBL_rep1/HI.1997.004.Index_3.MCF7_E2_NIPBL_ChIP8_GB_rep1.bam.trim.log
trimmomatic.HI.1997.004.Index_3.MCF7_E2_NIPBL_ChIP8_GB_rep1.bam.84f963840a7b7aac1eb8bade14e20df1.mugqic.done
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
# JOB: trimmomatic_9_JOB_ID: trimmomatic.HI.1997.004.Index_5.MCF7_E2_SMC1_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.1997.004.Index_5.MCF7_E2_SMC1_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_9_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.1997.004.Index_5.MCF7_E2_SMC1_ChIP8_GB_rep1.bam.bd2e9505eefc30c3348412b81c0efd79.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.1997.004.Index_5.MCF7_E2_SMC1_ChIP8_GB_rep1.bam.bd2e9505eefc30c3348412b81c0efd79.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/MCF7_E2_SMC1A_rep1 && \
`cat > trim/MCF7_E2_SMC1A_rep1/HI.1997.004.Index_5.MCF7_E2_SMC1_ChIP8_GB_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.1997.004.Index_5.MCF7_E2_SMC1_ChIP8_GB_rep1.single.fastq.gz \
  trim/MCF7_E2_SMC1A_rep1/HI.1997.004.Index_5.MCF7_E2_SMC1_ChIP8_GB_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/MCF7_E2_SMC1A_rep1/HI.1997.004.Index_5.MCF7_E2_SMC1_ChIP8_GB_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/MCF7_E2_SMC1A_rep1/HI.1997.004.Index_5.MCF7_E2_SMC1_ChIP8_GB_rep1.bam.trim.log
trimmomatic.HI.1997.004.Index_5.MCF7_E2_SMC1_ChIP8_GB_rep1.bam.bd2e9505eefc30c3348412b81c0efd79.mugqic.done
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
# JOB: trimmomatic_10_JOB_ID: trimmomatic.HI.1943.008.Index_18.MCF7_E2_MED1_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.1943.008.Index_18.MCF7_E2_MED1_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_10_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.1943.008.Index_18.MCF7_E2_MED1_ChIP8_GB_rep1.bam.ebcc4a07dfdae6974cbc14b723b85324.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.1943.008.Index_18.MCF7_E2_MED1_ChIP8_GB_rep1.bam.ebcc4a07dfdae6974cbc14b723b85324.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/MCF7_E2_MED1_rep1 && \
`cat > trim/MCF7_E2_MED1_rep1/HI.1943.008.Index_18.MCF7_E2_MED1_ChIP8_GB_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.1943.008.Index_18.MCF7_E2_MED1_ChIP8_GB_rep1.single.fastq.gz \
  trim/MCF7_E2_MED1_rep1/HI.1943.008.Index_18.MCF7_E2_MED1_ChIP8_GB_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/MCF7_E2_MED1_rep1/HI.1943.008.Index_18.MCF7_E2_MED1_ChIP8_GB_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/MCF7_E2_MED1_rep1/HI.1943.008.Index_18.MCF7_E2_MED1_ChIP8_GB_rep1.bam.trim.log
trimmomatic.HI.1943.008.Index_18.MCF7_E2_MED1_ChIP8_GB_rep1.bam.ebcc4a07dfdae6974cbc14b723b85324.mugqic.done
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=24G -N 1 -n 6 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$trimmomatic_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: trimmomatic_11_JOB_ID: trimmomatic.HI.1997.005.Index_11.MCF7_E2_POL2_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.1997.005.Index_11.MCF7_E2_POL2_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_11_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.1997.005.Index_11.MCF7_E2_POL2_ChIP8_GB_rep1.bam.18b90e8340399e0cb2e293c2f98e5df0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.1997.005.Index_11.MCF7_E2_POL2_ChIP8_GB_rep1.bam.18b90e8340399e0cb2e293c2f98e5df0.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/MCF7_E2_POL2_rep1 && \
`cat > trim/MCF7_E2_POL2_rep1/HI.1997.005.Index_11.MCF7_E2_POL2_ChIP8_GB_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.1997.005.Index_11.MCF7_E2_POL2_ChIP8_GB_rep1.single.fastq.gz \
  trim/MCF7_E2_POL2_rep1/HI.1997.005.Index_11.MCF7_E2_POL2_ChIP8_GB_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/MCF7_E2_POL2_rep1/HI.1997.005.Index_11.MCF7_E2_POL2_ChIP8_GB_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/MCF7_E2_POL2_rep1/HI.1997.005.Index_11.MCF7_E2_POL2_ChIP8_GB_rep1.bam.trim.log
trimmomatic.HI.1997.005.Index_11.MCF7_E2_POL2_ChIP8_GB_rep1.bam.18b90e8340399e0cb2e293c2f98e5df0.mugqic.done
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=24G -N 1 -n 6 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$trimmomatic_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: trimmomatic_12_JOB_ID: trimmomatic.HI.1997.005.Index_21.MCF7_E2_ERA_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.1997.005.Index_21.MCF7_E2_ERA_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_12_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.1997.005.Index_21.MCF7_E2_ERA_ChIP8_GB_rep1.bam.fd99b74a020a6ffd72ff37976e2561fa.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.1997.005.Index_21.MCF7_E2_ERA_ChIP8_GB_rep1.bam.fd99b74a020a6ffd72ff37976e2561fa.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/MCF7_E2_ERA_rep1 && \
`cat > trim/MCF7_E2_ERA_rep1/HI.1997.005.Index_21.MCF7_E2_ERA_ChIP8_GB_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.1997.005.Index_21.MCF7_E2_ERA_ChIP8_GB_rep1.single.fastq.gz \
  trim/MCF7_E2_ERA_rep1/HI.1997.005.Index_21.MCF7_E2_ERA_ChIP8_GB_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/MCF7_E2_ERA_rep1/HI.1997.005.Index_21.MCF7_E2_ERA_ChIP8_GB_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/MCF7_E2_ERA_rep1/HI.1997.005.Index_21.MCF7_E2_ERA_ChIP8_GB_rep1.bam.trim.log
trimmomatic.HI.1997.005.Index_21.MCF7_E2_ERA_ChIP8_GB_rep1.bam.fd99b74a020a6ffd72ff37976e2561fa.mugqic.done
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=24G -N 1 -n 6 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$trimmomatic_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: trimmomatic_13_JOB_ID: trimmomatic.HI.2791.003.Index_5.MCF7_WCE_C_rep2.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.2791.003.Index_5.MCF7_WCE_C_rep2.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_13_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.2791.003.Index_5.MCF7_WCE_C_rep2.bam.574ed734986b12b24aa6cd8d7511a7c7.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.2791.003.Index_5.MCF7_WCE_C_rep2.bam.574ed734986b12b24aa6cd8d7511a7c7.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/MCF7_CTRL_WCE_rep2 && \
`cat > trim/MCF7_CTRL_WCE_rep2/HI.2791.003.Index_5.MCF7_WCE_C_rep2.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.2791.003.Index_5.MCF7_WCE_C_rep2.single.fastq.gz \
  trim/MCF7_CTRL_WCE_rep2/HI.2791.003.Index_5.MCF7_WCE_C_rep2.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/MCF7_CTRL_WCE_rep2/HI.2791.003.Index_5.MCF7_WCE_C_rep2.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/MCF7_CTRL_WCE_rep2/HI.2791.003.Index_5.MCF7_WCE_C_rep2.bam.trim.log
trimmomatic.HI.2791.003.Index_5.MCF7_WCE_C_rep2.bam.574ed734986b12b24aa6cd8d7511a7c7.mugqic.done
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=24G -N 1 -n 6 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$trimmomatic_13_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: trimmomatic_14_JOB_ID: trimmomatic.HI.2791.003.Index_9.MCF7_1_NIPBL_C_rep2.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.2791.003.Index_9.MCF7_1_NIPBL_C_rep2.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_14_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.2791.003.Index_9.MCF7_1_NIPBL_C_rep2.bam.650e2fe99938e41a634d7c38eeaae209.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.2791.003.Index_9.MCF7_1_NIPBL_C_rep2.bam.650e2fe99938e41a634d7c38eeaae209.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/MCF7_CTRL_NIPBL_rep2 && \
`cat > trim/MCF7_CTRL_NIPBL_rep2/HI.2791.003.Index_9.MCF7_1_NIPBL_C_rep2.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.2791.003.Index_9.MCF7_1_NIPBL_C_rep2.single.fastq.gz \
  trim/MCF7_CTRL_NIPBL_rep2/HI.2791.003.Index_9.MCF7_1_NIPBL_C_rep2.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/MCF7_CTRL_NIPBL_rep2/HI.2791.003.Index_9.MCF7_1_NIPBL_C_rep2.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/MCF7_CTRL_NIPBL_rep2/HI.2791.003.Index_9.MCF7_1_NIPBL_C_rep2.bam.trim.log
trimmomatic.HI.2791.003.Index_9.MCF7_1_NIPBL_C_rep2.bam.650e2fe99938e41a634d7c38eeaae209.mugqic.done
)
trimmomatic_14_JOB_ID=$(echo "#! /bin/bash
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
echo "$trimmomatic_14_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: trimmomatic_15_JOB_ID: trimmomatic.HI.2791.003.Index_11.MCF7_4_Smc1_C_rep2.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.2791.003.Index_11.MCF7_4_Smc1_C_rep2.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_15_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.2791.003.Index_11.MCF7_4_Smc1_C_rep2.bam.f3feef01b40d397c1c959a8ac7de041d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.2791.003.Index_11.MCF7_4_Smc1_C_rep2.bam.f3feef01b40d397c1c959a8ac7de041d.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/MCF7_CTRL_SMC1A_rep2 && \
`cat > trim/MCF7_CTRL_SMC1A_rep2/HI.2791.003.Index_11.MCF7_4_Smc1_C_rep2.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.2791.003.Index_11.MCF7_4_Smc1_C_rep2.single.fastq.gz \
  trim/MCF7_CTRL_SMC1A_rep2/HI.2791.003.Index_11.MCF7_4_Smc1_C_rep2.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/MCF7_CTRL_SMC1A_rep2/HI.2791.003.Index_11.MCF7_4_Smc1_C_rep2.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/MCF7_CTRL_SMC1A_rep2/HI.2791.003.Index_11.MCF7_4_Smc1_C_rep2.bam.trim.log
trimmomatic.HI.2791.003.Index_11.MCF7_4_Smc1_C_rep2.bam.f3feef01b40d397c1c959a8ac7de041d.mugqic.done
)
trimmomatic_15_JOB_ID=$(echo "#! /bin/bash
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
echo "$trimmomatic_15_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: trimmomatic_16_JOB_ID: trimmomatic.HI.2791.003.Index_8.MCF7_2_MED1_C_rep2.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.2791.003.Index_8.MCF7_2_MED1_C_rep2.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_16_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.2791.003.Index_8.MCF7_2_MED1_C_rep2.bam.c4fef4fafd622254fc07d76a480b9de4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.2791.003.Index_8.MCF7_2_MED1_C_rep2.bam.c4fef4fafd622254fc07d76a480b9de4.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/MCF7_CTRL_MED1_rep2 && \
`cat > trim/MCF7_CTRL_MED1_rep2/HI.2791.003.Index_8.MCF7_2_MED1_C_rep2.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.2791.003.Index_8.MCF7_2_MED1_C_rep2.single.fastq.gz \
  trim/MCF7_CTRL_MED1_rep2/HI.2791.003.Index_8.MCF7_2_MED1_C_rep2.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/MCF7_CTRL_MED1_rep2/HI.2791.003.Index_8.MCF7_2_MED1_C_rep2.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/MCF7_CTRL_MED1_rep2/HI.2791.003.Index_8.MCF7_2_MED1_C_rep2.bam.trim.log
trimmomatic.HI.2791.003.Index_8.MCF7_2_MED1_C_rep2.bam.c4fef4fafd622254fc07d76a480b9de4.mugqic.done
)
trimmomatic_16_JOB_ID=$(echo "#! /bin/bash
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
echo "$trimmomatic_16_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: trimmomatic_17_JOB_ID: trimmomatic.HI.2791.003.Index_10.MCF7_3_Era_C_rep2.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.2791.003.Index_10.MCF7_3_Era_C_rep2.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_17_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.2791.003.Index_10.MCF7_3_Era_C_rep2.bam.c4e61091ef5071ce0291320989281ee8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.2791.003.Index_10.MCF7_3_Era_C_rep2.bam.c4e61091ef5071ce0291320989281ee8.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/MCF7_CTRL_ERA_rep2 && \
`cat > trim/MCF7_CTRL_ERA_rep2/HI.2791.003.Index_10.MCF7_3_Era_C_rep2.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.2791.003.Index_10.MCF7_3_Era_C_rep2.single.fastq.gz \
  trim/MCF7_CTRL_ERA_rep2/HI.2791.003.Index_10.MCF7_3_Era_C_rep2.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/MCF7_CTRL_ERA_rep2/HI.2791.003.Index_10.MCF7_3_Era_C_rep2.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/MCF7_CTRL_ERA_rep2/HI.2791.003.Index_10.MCF7_3_Era_C_rep2.bam.trim.log
trimmomatic.HI.2791.003.Index_10.MCF7_3_Era_C_rep2.bam.c4e61091ef5071ce0291320989281ee8.mugqic.done
)
trimmomatic_17_JOB_ID=$(echo "#! /bin/bash
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
echo "$trimmomatic_17_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: trimmomatic_18_JOB_ID: trimmomatic.HI.2791.003.Index_4.MCF7_WCE_E_rep2.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.2791.003.Index_4.MCF7_WCE_E_rep2.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_18_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.2791.003.Index_4.MCF7_WCE_E_rep2.bam.b82b121b689951c5c09531c68ed08536.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.2791.003.Index_4.MCF7_WCE_E_rep2.bam.b82b121b689951c5c09531c68ed08536.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/MCF7_E2_WCE_rep2 && \
`cat > trim/MCF7_E2_WCE_rep2/HI.2791.003.Index_4.MCF7_WCE_E_rep2.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.2791.003.Index_4.MCF7_WCE_E_rep2.single.fastq.gz \
  trim/MCF7_E2_WCE_rep2/HI.2791.003.Index_4.MCF7_WCE_E_rep2.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/MCF7_E2_WCE_rep2/HI.2791.003.Index_4.MCF7_WCE_E_rep2.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/MCF7_E2_WCE_rep2/HI.2791.003.Index_4.MCF7_WCE_E_rep2.bam.trim.log
trimmomatic.HI.2791.003.Index_4.MCF7_WCE_E_rep2.bam.b82b121b689951c5c09531c68ed08536.mugqic.done
)
trimmomatic_18_JOB_ID=$(echo "#! /bin/bash
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
echo "$trimmomatic_18_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: trimmomatic_19_JOB_ID: trimmomatic.HI.2791.003.Index_12.MCF7_1_NIPBL_E_rep2.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.2791.003.Index_12.MCF7_1_NIPBL_E_rep2.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_19_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.2791.003.Index_12.MCF7_1_NIPBL_E_rep2.bam.7cc3a2226419cff35c7574ef998ea6af.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.2791.003.Index_12.MCF7_1_NIPBL_E_rep2.bam.7cc3a2226419cff35c7574ef998ea6af.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/MCF7_E2_NIPBL_rep2 && \
`cat > trim/MCF7_E2_NIPBL_rep2/HI.2791.003.Index_12.MCF7_1_NIPBL_E_rep2.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.2791.003.Index_12.MCF7_1_NIPBL_E_rep2.single.fastq.gz \
  trim/MCF7_E2_NIPBL_rep2/HI.2791.003.Index_12.MCF7_1_NIPBL_E_rep2.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/MCF7_E2_NIPBL_rep2/HI.2791.003.Index_12.MCF7_1_NIPBL_E_rep2.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/MCF7_E2_NIPBL_rep2/HI.2791.003.Index_12.MCF7_1_NIPBL_E_rep2.bam.trim.log
trimmomatic.HI.2791.003.Index_12.MCF7_1_NIPBL_E_rep2.bam.7cc3a2226419cff35c7574ef998ea6af.mugqic.done
)
trimmomatic_19_JOB_ID=$(echo "#! /bin/bash
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
echo "$trimmomatic_19_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: trimmomatic_20_JOB_ID: trimmomatic.HI.2791.003.Index_3.MCF7_4_Smc1_E_rep2.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.2791.003.Index_3.MCF7_4_Smc1_E_rep2.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_20_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.2791.003.Index_3.MCF7_4_Smc1_E_rep2.bam.22f209b22a9e88eb7a853beed3d75781.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.2791.003.Index_3.MCF7_4_Smc1_E_rep2.bam.22f209b22a9e88eb7a853beed3d75781.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/MCF7_E2_SMC1A_rep2 && \
`cat > trim/MCF7_E2_SMC1A_rep2/HI.2791.003.Index_3.MCF7_4_Smc1_E_rep2.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.2791.003.Index_3.MCF7_4_Smc1_E_rep2.single.fastq.gz \
  trim/MCF7_E2_SMC1A_rep2/HI.2791.003.Index_3.MCF7_4_Smc1_E_rep2.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/MCF7_E2_SMC1A_rep2/HI.2791.003.Index_3.MCF7_4_Smc1_E_rep2.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/MCF7_E2_SMC1A_rep2/HI.2791.003.Index_3.MCF7_4_Smc1_E_rep2.bam.trim.log
trimmomatic.HI.2791.003.Index_3.MCF7_4_Smc1_E_rep2.bam.22f209b22a9e88eb7a853beed3d75781.mugqic.done
)
trimmomatic_20_JOB_ID=$(echo "#! /bin/bash
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
echo "$trimmomatic_20_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: trimmomatic_21_JOB_ID: trimmomatic.HI.2791.003.Index_19.MCF7_2_MED1_E_rep2.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.2791.003.Index_19.MCF7_2_MED1_E_rep2.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_21_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.2791.003.Index_19.MCF7_2_MED1_E_rep2.bam.fdac30f66cc37814f80f08c031f487f7.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.2791.003.Index_19.MCF7_2_MED1_E_rep2.bam.fdac30f66cc37814f80f08c031f487f7.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/MCF7_E2_MED1_rep2 && \
`cat > trim/MCF7_E2_MED1_rep2/HI.2791.003.Index_19.MCF7_2_MED1_E_rep2.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.2791.003.Index_19.MCF7_2_MED1_E_rep2.single.fastq.gz \
  trim/MCF7_E2_MED1_rep2/HI.2791.003.Index_19.MCF7_2_MED1_E_rep2.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/MCF7_E2_MED1_rep2/HI.2791.003.Index_19.MCF7_2_MED1_E_rep2.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/MCF7_E2_MED1_rep2/HI.2791.003.Index_19.MCF7_2_MED1_E_rep2.bam.trim.log
trimmomatic.HI.2791.003.Index_19.MCF7_2_MED1_E_rep2.bam.fdac30f66cc37814f80f08c031f487f7.mugqic.done
)
trimmomatic_21_JOB_ID=$(echo "#! /bin/bash
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
echo "$trimmomatic_21_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: trimmomatic_22_JOB_ID: trimmomatic.HI.2791.003.Index_1.MCF7_3_Era_E_rep2.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.2791.003.Index_1.MCF7_3_Era_E_rep2.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_22_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.2791.003.Index_1.MCF7_3_Era_E_rep2.bam.dbebed080d5706e2777011c22a38daa1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.2791.003.Index_1.MCF7_3_Era_E_rep2.bam.dbebed080d5706e2777011c22a38daa1.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/MCF7_E2_ERA_rep2 && \
`cat > trim/MCF7_E2_ERA_rep2/HI.2791.003.Index_1.MCF7_3_Era_E_rep2.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/HI.2791.003.Index_1.MCF7_3_Era_E_rep2.single.fastq.gz \
  trim/MCF7_E2_ERA_rep2/HI.2791.003.Index_1.MCF7_3_Era_E_rep2.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/MCF7_E2_ERA_rep2/HI.2791.003.Index_1.MCF7_3_Era_E_rep2.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/MCF7_E2_ERA_rep2/HI.2791.003.Index_1.MCF7_3_Era_E_rep2.bam.trim.log
trimmomatic.HI.2791.003.Index_1.MCF7_3_Era_E_rep2.bam.dbebed080d5706e2777011c22a38daa1.mugqic.done
)
trimmomatic_22_JOB_ID=$(echo "#! /bin/bash
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
echo "$trimmomatic_22_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: trimmomatic_23_JOB_ID: trimmomatic.SRR1193526.fastq
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.SRR1193526.fastq
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.SRR1193526.fastq.68610d0e95dead85720c41f1832c2984.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.SRR1193526.fastq.68610d0e95dead85720c41f1832c2984.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/MCF7_CTRL_BRD4_rep1 && \
`cat > trim/MCF7_CTRL_BRD4_rep1/SRR1193526.fastq.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/SRR1193526.fastq \
  trim/MCF7_CTRL_BRD4_rep1/SRR1193526.fastq.trim.single.fastq.gz \
  ILLUMINACLIP:trim/MCF7_CTRL_BRD4_rep1/SRR1193526.fastq.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/MCF7_CTRL_BRD4_rep1/SRR1193526.fastq.trim.log
trimmomatic.SRR1193526.fastq.68610d0e95dead85720c41f1832c2984.mugqic.done
)
trimmomatic_23_JOB_ID=$(echo "#! /bin/bash
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
echo "$trimmomatic_23_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: trimmomatic_24_JOB_ID: trimmomatic.SRR1193527.fastq
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.SRR1193527.fastq
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.SRR1193527.fastq.508d225fa16401f3494bf65c493f20bf.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.SRR1193527.fastq.508d225fa16401f3494bf65c493f20bf.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/MCF7_CTRL_BRD4_rep2 && \
`cat > trim/MCF7_CTRL_BRD4_rep2/SRR1193527.fastq.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/SRR1193527.fastq \
  trim/MCF7_CTRL_BRD4_rep2/SRR1193527.fastq.trim.single.fastq.gz \
  ILLUMINACLIP:trim/MCF7_CTRL_BRD4_rep2/SRR1193527.fastq.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/MCF7_CTRL_BRD4_rep2/SRR1193527.fastq.trim.log
trimmomatic.SRR1193527.fastq.508d225fa16401f3494bf65c493f20bf.mugqic.done
)
trimmomatic_24_JOB_ID=$(echo "#! /bin/bash
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
echo "$trimmomatic_24_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: trimmomatic_25_JOB_ID: trimmomatic.SRR1193528.fastq
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.SRR1193528.fastq
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.SRR1193528.fastq.4e3ca86811c97052c443f52e4723bf41.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.SRR1193528.fastq.4e3ca86811c97052c443f52e4723bf41.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/MCF7_CTRL_BRD4_rep3 && \
`cat > trim/MCF7_CTRL_BRD4_rep3/SRR1193528.fastq.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/SRR1193528.fastq \
  trim/MCF7_CTRL_BRD4_rep3/SRR1193528.fastq.trim.single.fastq.gz \
  ILLUMINACLIP:trim/MCF7_CTRL_BRD4_rep3/SRR1193528.fastq.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/MCF7_CTRL_BRD4_rep3/SRR1193528.fastq.trim.log
trimmomatic.SRR1193528.fastq.4e3ca86811c97052c443f52e4723bf41.mugqic.done
)
trimmomatic_25_JOB_ID=$(echo "#! /bin/bash
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
echo "$trimmomatic_25_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: trimmomatic_26_JOB_ID: trimmomatic.SRR1193529.fastq
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.SRR1193529.fastq
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.SRR1193529.fastq.41fbe085efe08da9bd94d1d0f800982f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.SRR1193529.fastq.41fbe085efe08da9bd94d1d0f800982f.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/MCF7_E2_BRD4_rep1 && \
`cat > trim/MCF7_E2_BRD4_rep1/SRR1193529.fastq.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/SRR1193529.fastq \
  trim/MCF7_E2_BRD4_rep1/SRR1193529.fastq.trim.single.fastq.gz \
  ILLUMINACLIP:trim/MCF7_E2_BRD4_rep1/SRR1193529.fastq.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/MCF7_E2_BRD4_rep1/SRR1193529.fastq.trim.log
trimmomatic.SRR1193529.fastq.41fbe085efe08da9bd94d1d0f800982f.mugqic.done
)
trimmomatic_26_JOB_ID=$(echo "#! /bin/bash
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
echo "$trimmomatic_26_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: trimmomatic_27_JOB_ID: trimmomatic.SRR1193530.fastq
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.SRR1193530.fastq
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.SRR1193530.fastq.fbf460cd61ed0a8a3741305b991d37d8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.SRR1193530.fastq.fbf460cd61ed0a8a3741305b991d37d8.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/MCF7_E2_BRD4_rep2 && \
`cat > trim/MCF7_E2_BRD4_rep2/SRR1193530.fastq.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/SRR1193530.fastq \
  trim/MCF7_E2_BRD4_rep2/SRR1193530.fastq.trim.single.fastq.gz \
  ILLUMINACLIP:trim/MCF7_E2_BRD4_rep2/SRR1193530.fastq.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/MCF7_E2_BRD4_rep2/SRR1193530.fastq.trim.log
trimmomatic.SRR1193530.fastq.fbf460cd61ed0a8a3741305b991d37d8.mugqic.done
)
trimmomatic_27_JOB_ID=$(echo "#! /bin/bash
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
echo "$trimmomatic_27_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: trimmomatic_28_JOB_ID: trimmomatic.SRR1193531.fastq
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.SRR1193531.fastq
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.SRR1193531.fastq.b015f123cc81752d7c9a4626c3b9934e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.SRR1193531.fastq.b015f123cc81752d7c9a4626c3b9934e.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/MCF7_E2_BRD4_rep3 && \
`cat > trim/MCF7_E2_BRD4_rep3/SRR1193531.fastq.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/SRR1193531.fastq \
  trim/MCF7_E2_BRD4_rep3/SRR1193531.fastq.trim.single.fastq.gz \
  ILLUMINACLIP:trim/MCF7_E2_BRD4_rep3/SRR1193531.fastq.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/MCF7_E2_BRD4_rep3/SRR1193531.fastq.trim.log
trimmomatic.SRR1193531.fastq.b015f123cc81752d7c9a4626c3b9934e.mugqic.done
)
trimmomatic_28_JOB_ID=$(echo "#! /bin/bash
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
echo "$trimmomatic_28_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: trimmomatic_29_JOB_ID: trimmomatic.SRR1193562.fastq
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.SRR1193562.fastq
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.SRR1193562.fastq.5591aa4bcb0ebb46f871d350e7850eaa.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.SRR1193562.fastq.5591aa4bcb0ebb46f871d350e7850eaa.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/MCF7_CTRL_WCE && \
`cat > trim/MCF7_CTRL_WCE/SRR1193562.fastq.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/SRR1193562.fastq \
  trim/MCF7_CTRL_WCE/SRR1193562.fastq.trim.single.fastq.gz \
  ILLUMINACLIP:trim/MCF7_CTRL_WCE/SRR1193562.fastq.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/MCF7_CTRL_WCE/SRR1193562.fastq.trim.log
trimmomatic.SRR1193562.fastq.5591aa4bcb0ebb46f871d350e7850eaa.mugqic.done
)
trimmomatic_29_JOB_ID=$(echo "#! /bin/bash
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
echo "$trimmomatic_29_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: trimmomatic_30_JOB_ID: trimmomatic.SRR1193563.fastq
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.SRR1193563.fastq
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.SRR1193563.fastq.4ab1835df71333bd21b3d25bb2f0aa6a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.SRR1193563.fastq.4ab1835df71333bd21b3d25bb2f0aa6a.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/MCF7_E2_WCE && \
`cat > trim/MCF7_E2_WCE/SRR1193563.fastq.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /scratch/chris11/sb_cofactor/MCF7/raw/chip-seq/SRR1193563.fastq \
  trim/MCF7_E2_WCE/SRR1193563.fastq.trim.single.fastq.gz \
  ILLUMINACLIP:trim/MCF7_E2_WCE/SRR1193563.fastq.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/MCF7_E2_WCE/SRR1193563.fastq.trim.log
trimmomatic.SRR1193563.fastq.4ab1835df71333bd21b3d25bb2f0aa6a.mugqic.done
)
trimmomatic_30_JOB_ID=$(echo "#! /bin/bash
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
echo "$trimmomatic_30_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

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
JOB_DEPENDENCIES=$trimmomatic_1_JOB_ID:$trimmomatic_2_JOB_ID:$trimmomatic_3_JOB_ID:$trimmomatic_4_JOB_ID:$trimmomatic_5_JOB_ID:$trimmomatic_6_JOB_ID:$trimmomatic_7_JOB_ID:$trimmomatic_8_JOB_ID:$trimmomatic_9_JOB_ID:$trimmomatic_10_JOB_ID:$trimmomatic_11_JOB_ID:$trimmomatic_12_JOB_ID:$trimmomatic_13_JOB_ID:$trimmomatic_14_JOB_ID:$trimmomatic_15_JOB_ID:$trimmomatic_16_JOB_ID:$trimmomatic_17_JOB_ID:$trimmomatic_18_JOB_ID:$trimmomatic_19_JOB_ID:$trimmomatic_20_JOB_ID:$trimmomatic_21_JOB_ID:$trimmomatic_22_JOB_ID:$trimmomatic_23_JOB_ID:$trimmomatic_24_JOB_ID:$trimmomatic_25_JOB_ID:$trimmomatic_26_JOB_ID:$trimmomatic_27_JOB_ID:$trimmomatic_28_JOB_ID:$trimmomatic_29_JOB_ID:$trimmomatic_30_JOB_ID
JOB_DONE=job_output/merge_trimmomatic_stats/merge_trimmomatic_stats.6554c2633306d7fa5f61133d3c160fef.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'merge_trimmomatic_stats.6554c2633306d7fa5f61133d3c160fef.mugqic.done'
module load mugqic/pandoc/1.15.2 && \
mkdir -p metrics && \
echo 'Sample	Readset	Raw Single Reads #	Surviving Single Reads #	Surviving Single Reads %' > metrics/trimReadsetTable.tsv && \
grep ^Input trim/MCF7_CTRL_WCE_rep1/HI.1943.008.Index_2.MCF7_CTRL_WCE_ChIP8_GB_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/MCF7_CTRL_WCE_rep1	HI.1943.008.Index_2.MCF7_CTRL_WCE_ChIP8_GB_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/MCF7_CTRL_NIPBL_rep1/HI.1997.004.Index_1.MCF7_CTRL_NIPBL_ChIP8_GB_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/MCF7_CTRL_NIPBL_rep1	HI.1997.004.Index_1.MCF7_CTRL_NIPBL_ChIP8_GB_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/MCF7_CTRL_SMC1A_rep1/HI.1997.004.Index_4.MCF7_CTRL_SMC1_ChIP8_GB_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/MCF7_CTRL_SMC1A_rep1	HI.1997.004.Index_4.MCF7_CTRL_SMC1_ChIP8_GB_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/MCF7_CTRL_MED1_rep1/HI.1943.008.Index_7.MCF7_CTRL_MED1_ChIP8_GB_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/MCF7_CTRL_MED1_rep1	HI.1943.008.Index_7.MCF7_CTRL_MED1_ChIP8_GB_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/MCF7_CTRL_POL2_rep1/HI.1997.005.Index_10.MCF7_CTRL_POL2_ChIP8_GB_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/MCF7_CTRL_POL2_rep1	HI.1997.005.Index_10.MCF7_CTRL_POL2_ChIP8_GB_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/MCF7_CTRL_ERA_rep1/HI.1997.005.Index_25.MCF7_CTRL_ERA_ChIP8_GB_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/MCF7_CTRL_ERA_rep1	HI.1997.005.Index_25.MCF7_CTRL_ERA_ChIP8_GB_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/MCF7_E2_WCE_rep1/HI.1943.008.Index_13.MCF7_E2_WCE_ChIP8_GB_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/MCF7_E2_WCE_rep1	HI.1943.008.Index_13.MCF7_E2_WCE_ChIP8_GB_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/MCF7_E2_NIPBL_rep1/HI.1997.004.Index_3.MCF7_E2_NIPBL_ChIP8_GB_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/MCF7_E2_NIPBL_rep1	HI.1997.004.Index_3.MCF7_E2_NIPBL_ChIP8_GB_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/MCF7_E2_SMC1A_rep1/HI.1997.004.Index_5.MCF7_E2_SMC1_ChIP8_GB_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/MCF7_E2_SMC1A_rep1	HI.1997.004.Index_5.MCF7_E2_SMC1_ChIP8_GB_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/MCF7_E2_MED1_rep1/HI.1943.008.Index_18.MCF7_E2_MED1_ChIP8_GB_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/MCF7_E2_MED1_rep1	HI.1943.008.Index_18.MCF7_E2_MED1_ChIP8_GB_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/MCF7_E2_POL2_rep1/HI.1997.005.Index_11.MCF7_E2_POL2_ChIP8_GB_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/MCF7_E2_POL2_rep1	HI.1997.005.Index_11.MCF7_E2_POL2_ChIP8_GB_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/MCF7_E2_ERA_rep1/HI.1997.005.Index_21.MCF7_E2_ERA_ChIP8_GB_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/MCF7_E2_ERA_rep1	HI.1997.005.Index_21.MCF7_E2_ERA_ChIP8_GB_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/MCF7_CTRL_WCE_rep2/HI.2791.003.Index_5.MCF7_WCE_C_rep2.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/MCF7_CTRL_WCE_rep2	HI.2791.003.Index_5.MCF7_WCE_C_rep2.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/MCF7_CTRL_NIPBL_rep2/HI.2791.003.Index_9.MCF7_1_NIPBL_C_rep2.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/MCF7_CTRL_NIPBL_rep2	HI.2791.003.Index_9.MCF7_1_NIPBL_C_rep2.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/MCF7_CTRL_SMC1A_rep2/HI.2791.003.Index_11.MCF7_4_Smc1_C_rep2.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/MCF7_CTRL_SMC1A_rep2	HI.2791.003.Index_11.MCF7_4_Smc1_C_rep2.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/MCF7_CTRL_MED1_rep2/HI.2791.003.Index_8.MCF7_2_MED1_C_rep2.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/MCF7_CTRL_MED1_rep2	HI.2791.003.Index_8.MCF7_2_MED1_C_rep2.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/MCF7_CTRL_ERA_rep2/HI.2791.003.Index_10.MCF7_3_Era_C_rep2.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/MCF7_CTRL_ERA_rep2	HI.2791.003.Index_10.MCF7_3_Era_C_rep2.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/MCF7_E2_WCE_rep2/HI.2791.003.Index_4.MCF7_WCE_E_rep2.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/MCF7_E2_WCE_rep2	HI.2791.003.Index_4.MCF7_WCE_E_rep2.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/MCF7_E2_NIPBL_rep2/HI.2791.003.Index_12.MCF7_1_NIPBL_E_rep2.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/MCF7_E2_NIPBL_rep2	HI.2791.003.Index_12.MCF7_1_NIPBL_E_rep2.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/MCF7_E2_SMC1A_rep2/HI.2791.003.Index_3.MCF7_4_Smc1_E_rep2.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/MCF7_E2_SMC1A_rep2	HI.2791.003.Index_3.MCF7_4_Smc1_E_rep2.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/MCF7_E2_MED1_rep2/HI.2791.003.Index_19.MCF7_2_MED1_E_rep2.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/MCF7_E2_MED1_rep2	HI.2791.003.Index_19.MCF7_2_MED1_E_rep2.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/MCF7_E2_ERA_rep2/HI.2791.003.Index_1.MCF7_3_Era_E_rep2.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/MCF7_E2_ERA_rep2	HI.2791.003.Index_1.MCF7_3_Era_E_rep2.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/MCF7_CTRL_BRD4_rep1/SRR1193526.fastq.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/MCF7_CTRL_BRD4_rep1	SRR1193526.fastq	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/MCF7_CTRL_BRD4_rep2/SRR1193527.fastq.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/MCF7_CTRL_BRD4_rep2	SRR1193527.fastq	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/MCF7_CTRL_BRD4_rep3/SRR1193528.fastq.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/MCF7_CTRL_BRD4_rep3	SRR1193528.fastq	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/MCF7_E2_BRD4_rep1/SRR1193529.fastq.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/MCF7_E2_BRD4_rep1	SRR1193529.fastq	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/MCF7_E2_BRD4_rep2/SRR1193530.fastq.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/MCF7_E2_BRD4_rep2	SRR1193530.fastq	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/MCF7_E2_BRD4_rep3/SRR1193531.fastq.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/MCF7_E2_BRD4_rep3	SRR1193531.fastq	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/MCF7_CTRL_WCE/SRR1193562.fastq.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/MCF7_CTRL_WCE	SRR1193562.fastq	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/MCF7_E2_WCE/SRR1193563.fastq.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/MCF7_E2_WCE	SRR1193563.fastq	\1	\2/' | \
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
merge_trimmomatic_stats.6554c2633306d7fa5f61133d3c160fef.mugqic.done
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
# JOB: bwa_mem_picard_sort_sam_1_JOB_ID: bwa_mem_picard_sort_sam.HI.1943.008.Index_2.MCF7_CTRL_WCE_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.1943.008.Index_2.MCF7_CTRL_WCE_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_1_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.1943.008.Index_2.MCF7_CTRL_WCE_ChIP8_GB_rep1.bam.2381ff135764e0769ce3d38afce031df.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.1943.008.Index_2.MCF7_CTRL_WCE_ChIP8_GB_rep1.bam.2381ff135764e0769ce3d38afce031df.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p alignment/MCF7_CTRL_WCE_rep1/HI.1943.008.Index_2.MCF7_CTRL_WCE_ChIP8_GB_rep1.bam && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:HI.1943.008.Index_2.MCF7_CTRL_WCE_ChIP8_GB_rep1.bam	SM:MCF7_CTRL_WCE_rep1	LB:MCF7_CTRL_WCE_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/MCF7_CTRL_WCE_rep1/HI.1943.008.Index_2.MCF7_CTRL_WCE_ChIP8_GB_rep1.bam.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/MCF7_CTRL_WCE_rep1/HI.1943.008.Index_2.MCF7_CTRL_WCE_ChIP8_GB_rep1.bam/HI.1943.008.Index_2.MCF7_CTRL_WCE_ChIP8_GB_rep1.bam.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.HI.1943.008.Index_2.MCF7_CTRL_WCE_ChIP8_GB_rep1.bam.2381ff135764e0769ce3d38afce031df.mugqic.done
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
# JOB: bwa_mem_picard_sort_sam_2_JOB_ID: bwa_mem_picard_sort_sam.HI.1997.004.Index_1.MCF7_CTRL_NIPBL_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.1997.004.Index_1.MCF7_CTRL_NIPBL_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_2_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.1997.004.Index_1.MCF7_CTRL_NIPBL_ChIP8_GB_rep1.bam.c2b6107a3d18c13b92aaf2461d750944.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.1997.004.Index_1.MCF7_CTRL_NIPBL_ChIP8_GB_rep1.bam.c2b6107a3d18c13b92aaf2461d750944.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p alignment/MCF7_CTRL_NIPBL_rep1/HI.1997.004.Index_1.MCF7_CTRL_NIPBL_ChIP8_GB_rep1.bam && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:HI.1997.004.Index_1.MCF7_CTRL_NIPBL_ChIP8_GB_rep1.bam	SM:MCF7_CTRL_NIPBL_rep1	LB:MCF7_CTRL_NIPBL_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/MCF7_CTRL_NIPBL_rep1/HI.1997.004.Index_1.MCF7_CTRL_NIPBL_ChIP8_GB_rep1.bam.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/MCF7_CTRL_NIPBL_rep1/HI.1997.004.Index_1.MCF7_CTRL_NIPBL_ChIP8_GB_rep1.bam/HI.1997.004.Index_1.MCF7_CTRL_NIPBL_ChIP8_GB_rep1.bam.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.HI.1997.004.Index_1.MCF7_CTRL_NIPBL_ChIP8_GB_rep1.bam.c2b6107a3d18c13b92aaf2461d750944.mugqic.done
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
# JOB: bwa_mem_picard_sort_sam_3_JOB_ID: bwa_mem_picard_sort_sam.HI.1997.004.Index_4.MCF7_CTRL_SMC1_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.1997.004.Index_4.MCF7_CTRL_SMC1_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_3_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.1997.004.Index_4.MCF7_CTRL_SMC1_ChIP8_GB_rep1.bam.1eede0abc4a4d9e960e0b4399e5ccb9c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.1997.004.Index_4.MCF7_CTRL_SMC1_ChIP8_GB_rep1.bam.1eede0abc4a4d9e960e0b4399e5ccb9c.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p alignment/MCF7_CTRL_SMC1A_rep1/HI.1997.004.Index_4.MCF7_CTRL_SMC1_ChIP8_GB_rep1.bam && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:HI.1997.004.Index_4.MCF7_CTRL_SMC1_ChIP8_GB_rep1.bam	SM:MCF7_CTRL_SMC1A_rep1	LB:MCF7_CTRL_SMC1A_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/MCF7_CTRL_SMC1A_rep1/HI.1997.004.Index_4.MCF7_CTRL_SMC1_ChIP8_GB_rep1.bam.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/MCF7_CTRL_SMC1A_rep1/HI.1997.004.Index_4.MCF7_CTRL_SMC1_ChIP8_GB_rep1.bam/HI.1997.004.Index_4.MCF7_CTRL_SMC1_ChIP8_GB_rep1.bam.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.HI.1997.004.Index_4.MCF7_CTRL_SMC1_ChIP8_GB_rep1.bam.1eede0abc4a4d9e960e0b4399e5ccb9c.mugqic.done
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
# JOB: bwa_mem_picard_sort_sam_4_JOB_ID: bwa_mem_picard_sort_sam.HI.1943.008.Index_7.MCF7_CTRL_MED1_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.1943.008.Index_7.MCF7_CTRL_MED1_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_4_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.1943.008.Index_7.MCF7_CTRL_MED1_ChIP8_GB_rep1.bam.60ac2852f23f0d2f4ae9100092266bee.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.1943.008.Index_7.MCF7_CTRL_MED1_ChIP8_GB_rep1.bam.60ac2852f23f0d2f4ae9100092266bee.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p alignment/MCF7_CTRL_MED1_rep1/HI.1943.008.Index_7.MCF7_CTRL_MED1_ChIP8_GB_rep1.bam && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:HI.1943.008.Index_7.MCF7_CTRL_MED1_ChIP8_GB_rep1.bam	SM:MCF7_CTRL_MED1_rep1	LB:MCF7_CTRL_MED1_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/MCF7_CTRL_MED1_rep1/HI.1943.008.Index_7.MCF7_CTRL_MED1_ChIP8_GB_rep1.bam.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/MCF7_CTRL_MED1_rep1/HI.1943.008.Index_7.MCF7_CTRL_MED1_ChIP8_GB_rep1.bam/HI.1943.008.Index_7.MCF7_CTRL_MED1_ChIP8_GB_rep1.bam.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.HI.1943.008.Index_7.MCF7_CTRL_MED1_ChIP8_GB_rep1.bam.60ac2852f23f0d2f4ae9100092266bee.mugqic.done
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
# JOB: bwa_mem_picard_sort_sam_5_JOB_ID: bwa_mem_picard_sort_sam.HI.1997.005.Index_10.MCF7_CTRL_POL2_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.1997.005.Index_10.MCF7_CTRL_POL2_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_5_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.1997.005.Index_10.MCF7_CTRL_POL2_ChIP8_GB_rep1.bam.26aecb25e22494e99d742c41f4c8a5f1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.1997.005.Index_10.MCF7_CTRL_POL2_ChIP8_GB_rep1.bam.26aecb25e22494e99d742c41f4c8a5f1.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p alignment/MCF7_CTRL_POL2_rep1/HI.1997.005.Index_10.MCF7_CTRL_POL2_ChIP8_GB_rep1.bam && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:HI.1997.005.Index_10.MCF7_CTRL_POL2_ChIP8_GB_rep1.bam	SM:MCF7_CTRL_POL2_rep1	LB:MCF7_CTRL_POL2_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/MCF7_CTRL_POL2_rep1/HI.1997.005.Index_10.MCF7_CTRL_POL2_ChIP8_GB_rep1.bam.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/MCF7_CTRL_POL2_rep1/HI.1997.005.Index_10.MCF7_CTRL_POL2_ChIP8_GB_rep1.bam/HI.1997.005.Index_10.MCF7_CTRL_POL2_ChIP8_GB_rep1.bam.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.HI.1997.005.Index_10.MCF7_CTRL_POL2_ChIP8_GB_rep1.bam.26aecb25e22494e99d742c41f4c8a5f1.mugqic.done
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
# JOB: bwa_mem_picard_sort_sam_6_JOB_ID: bwa_mem_picard_sort_sam.HI.1997.005.Index_25.MCF7_CTRL_ERA_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.1997.005.Index_25.MCF7_CTRL_ERA_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_6_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.1997.005.Index_25.MCF7_CTRL_ERA_ChIP8_GB_rep1.bam.e6533d52bb1bebf50ecb92502ac923dc.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.1997.005.Index_25.MCF7_CTRL_ERA_ChIP8_GB_rep1.bam.e6533d52bb1bebf50ecb92502ac923dc.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p alignment/MCF7_CTRL_ERA_rep1/HI.1997.005.Index_25.MCF7_CTRL_ERA_ChIP8_GB_rep1.bam && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:HI.1997.005.Index_25.MCF7_CTRL_ERA_ChIP8_GB_rep1.bam	SM:MCF7_CTRL_ERA_rep1	LB:MCF7_CTRL_ERA_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/MCF7_CTRL_ERA_rep1/HI.1997.005.Index_25.MCF7_CTRL_ERA_ChIP8_GB_rep1.bam.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/MCF7_CTRL_ERA_rep1/HI.1997.005.Index_25.MCF7_CTRL_ERA_ChIP8_GB_rep1.bam/HI.1997.005.Index_25.MCF7_CTRL_ERA_ChIP8_GB_rep1.bam.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.HI.1997.005.Index_25.MCF7_CTRL_ERA_ChIP8_GB_rep1.bam.e6533d52bb1bebf50ecb92502ac923dc.mugqic.done
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
# JOB: bwa_mem_picard_sort_sam_7_JOB_ID: bwa_mem_picard_sort_sam.HI.1943.008.Index_13.MCF7_E2_WCE_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.1943.008.Index_13.MCF7_E2_WCE_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_7_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.1943.008.Index_13.MCF7_E2_WCE_ChIP8_GB_rep1.bam.ed6c9e3e0ddd9935ff7d0f50d4db1322.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.1943.008.Index_13.MCF7_E2_WCE_ChIP8_GB_rep1.bam.ed6c9e3e0ddd9935ff7d0f50d4db1322.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p alignment/MCF7_E2_WCE_rep1/HI.1943.008.Index_13.MCF7_E2_WCE_ChIP8_GB_rep1.bam && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:HI.1943.008.Index_13.MCF7_E2_WCE_ChIP8_GB_rep1.bam	SM:MCF7_E2_WCE_rep1	LB:MCF7_E2_WCE_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/MCF7_E2_WCE_rep1/HI.1943.008.Index_13.MCF7_E2_WCE_ChIP8_GB_rep1.bam.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/MCF7_E2_WCE_rep1/HI.1943.008.Index_13.MCF7_E2_WCE_ChIP8_GB_rep1.bam/HI.1943.008.Index_13.MCF7_E2_WCE_ChIP8_GB_rep1.bam.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.HI.1943.008.Index_13.MCF7_E2_WCE_ChIP8_GB_rep1.bam.ed6c9e3e0ddd9935ff7d0f50d4db1322.mugqic.done
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
# JOB: bwa_mem_picard_sort_sam_8_JOB_ID: bwa_mem_picard_sort_sam.HI.1997.004.Index_3.MCF7_E2_NIPBL_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.1997.004.Index_3.MCF7_E2_NIPBL_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_8_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.1997.004.Index_3.MCF7_E2_NIPBL_ChIP8_GB_rep1.bam.3a60f4884e0f749c5b73af967ed678df.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.1997.004.Index_3.MCF7_E2_NIPBL_ChIP8_GB_rep1.bam.3a60f4884e0f749c5b73af967ed678df.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p alignment/MCF7_E2_NIPBL_rep1/HI.1997.004.Index_3.MCF7_E2_NIPBL_ChIP8_GB_rep1.bam && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:HI.1997.004.Index_3.MCF7_E2_NIPBL_ChIP8_GB_rep1.bam	SM:MCF7_E2_NIPBL_rep1	LB:MCF7_E2_NIPBL_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/MCF7_E2_NIPBL_rep1/HI.1997.004.Index_3.MCF7_E2_NIPBL_ChIP8_GB_rep1.bam.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/MCF7_E2_NIPBL_rep1/HI.1997.004.Index_3.MCF7_E2_NIPBL_ChIP8_GB_rep1.bam/HI.1997.004.Index_3.MCF7_E2_NIPBL_ChIP8_GB_rep1.bam.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.HI.1997.004.Index_3.MCF7_E2_NIPBL_ChIP8_GB_rep1.bam.3a60f4884e0f749c5b73af967ed678df.mugqic.done
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
# JOB: bwa_mem_picard_sort_sam_9_JOB_ID: bwa_mem_picard_sort_sam.HI.1997.004.Index_5.MCF7_E2_SMC1_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.1997.004.Index_5.MCF7_E2_SMC1_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_9_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.1997.004.Index_5.MCF7_E2_SMC1_ChIP8_GB_rep1.bam.5a6afc37bf07ec91acbc4c9ef496da18.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.1997.004.Index_5.MCF7_E2_SMC1_ChIP8_GB_rep1.bam.5a6afc37bf07ec91acbc4c9ef496da18.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p alignment/MCF7_E2_SMC1A_rep1/HI.1997.004.Index_5.MCF7_E2_SMC1_ChIP8_GB_rep1.bam && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:HI.1997.004.Index_5.MCF7_E2_SMC1_ChIP8_GB_rep1.bam	SM:MCF7_E2_SMC1A_rep1	LB:MCF7_E2_SMC1A_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/MCF7_E2_SMC1A_rep1/HI.1997.004.Index_5.MCF7_E2_SMC1_ChIP8_GB_rep1.bam.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/MCF7_E2_SMC1A_rep1/HI.1997.004.Index_5.MCF7_E2_SMC1_ChIP8_GB_rep1.bam/HI.1997.004.Index_5.MCF7_E2_SMC1_ChIP8_GB_rep1.bam.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.HI.1997.004.Index_5.MCF7_E2_SMC1_ChIP8_GB_rep1.bam.5a6afc37bf07ec91acbc4c9ef496da18.mugqic.done
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
# JOB: bwa_mem_picard_sort_sam_10_JOB_ID: bwa_mem_picard_sort_sam.HI.1943.008.Index_18.MCF7_E2_MED1_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.1943.008.Index_18.MCF7_E2_MED1_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_10_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.1943.008.Index_18.MCF7_E2_MED1_ChIP8_GB_rep1.bam.d26e48cbf19d3daf4ff94679f07c2beb.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.1943.008.Index_18.MCF7_E2_MED1_ChIP8_GB_rep1.bam.d26e48cbf19d3daf4ff94679f07c2beb.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p alignment/MCF7_E2_MED1_rep1/HI.1943.008.Index_18.MCF7_E2_MED1_ChIP8_GB_rep1.bam && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:HI.1943.008.Index_18.MCF7_E2_MED1_ChIP8_GB_rep1.bam	SM:MCF7_E2_MED1_rep1	LB:MCF7_E2_MED1_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/MCF7_E2_MED1_rep1/HI.1943.008.Index_18.MCF7_E2_MED1_ChIP8_GB_rep1.bam.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/MCF7_E2_MED1_rep1/HI.1943.008.Index_18.MCF7_E2_MED1_ChIP8_GB_rep1.bam/HI.1943.008.Index_18.MCF7_E2_MED1_ChIP8_GB_rep1.bam.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.HI.1943.008.Index_18.MCF7_E2_MED1_ChIP8_GB_rep1.bam.d26e48cbf19d3daf4ff94679f07c2beb.mugqic.done
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
# JOB: bwa_mem_picard_sort_sam_11_JOB_ID: bwa_mem_picard_sort_sam.HI.1997.005.Index_11.MCF7_E2_POL2_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.1997.005.Index_11.MCF7_E2_POL2_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_11_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.1997.005.Index_11.MCF7_E2_POL2_ChIP8_GB_rep1.bam.a04b1ea2cf6380f70b05dc3b61ecff3b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.1997.005.Index_11.MCF7_E2_POL2_ChIP8_GB_rep1.bam.a04b1ea2cf6380f70b05dc3b61ecff3b.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p alignment/MCF7_E2_POL2_rep1/HI.1997.005.Index_11.MCF7_E2_POL2_ChIP8_GB_rep1.bam && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:HI.1997.005.Index_11.MCF7_E2_POL2_ChIP8_GB_rep1.bam	SM:MCF7_E2_POL2_rep1	LB:MCF7_E2_POL2_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/MCF7_E2_POL2_rep1/HI.1997.005.Index_11.MCF7_E2_POL2_ChIP8_GB_rep1.bam.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/MCF7_E2_POL2_rep1/HI.1997.005.Index_11.MCF7_E2_POL2_ChIP8_GB_rep1.bam/HI.1997.005.Index_11.MCF7_E2_POL2_ChIP8_GB_rep1.bam.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.HI.1997.005.Index_11.MCF7_E2_POL2_ChIP8_GB_rep1.bam.a04b1ea2cf6380f70b05dc3b61ecff3b.mugqic.done
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
# JOB: bwa_mem_picard_sort_sam_12_JOB_ID: bwa_mem_picard_sort_sam.HI.1997.005.Index_21.MCF7_E2_ERA_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.1997.005.Index_21.MCF7_E2_ERA_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_12_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.1997.005.Index_21.MCF7_E2_ERA_ChIP8_GB_rep1.bam.0bcef63195d3dcfce3e0f81587d71a17.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.1997.005.Index_21.MCF7_E2_ERA_ChIP8_GB_rep1.bam.0bcef63195d3dcfce3e0f81587d71a17.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p alignment/MCF7_E2_ERA_rep1/HI.1997.005.Index_21.MCF7_E2_ERA_ChIP8_GB_rep1.bam && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:HI.1997.005.Index_21.MCF7_E2_ERA_ChIP8_GB_rep1.bam	SM:MCF7_E2_ERA_rep1	LB:MCF7_E2_ERA_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/MCF7_E2_ERA_rep1/HI.1997.005.Index_21.MCF7_E2_ERA_ChIP8_GB_rep1.bam.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/MCF7_E2_ERA_rep1/HI.1997.005.Index_21.MCF7_E2_ERA_ChIP8_GB_rep1.bam/HI.1997.005.Index_21.MCF7_E2_ERA_ChIP8_GB_rep1.bam.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.HI.1997.005.Index_21.MCF7_E2_ERA_ChIP8_GB_rep1.bam.0bcef63195d3dcfce3e0f81587d71a17.mugqic.done
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
# JOB: bwa_mem_picard_sort_sam_13_JOB_ID: bwa_mem_picard_sort_sam.HI.2791.003.Index_5.MCF7_WCE_C_rep2.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.2791.003.Index_5.MCF7_WCE_C_rep2.bam
JOB_DEPENDENCIES=$trimmomatic_13_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.2791.003.Index_5.MCF7_WCE_C_rep2.bam.acb320d6006f526b46b4e4c5d864cab4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.2791.003.Index_5.MCF7_WCE_C_rep2.bam.acb320d6006f526b46b4e4c5d864cab4.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p alignment/MCF7_CTRL_WCE_rep2/HI.2791.003.Index_5.MCF7_WCE_C_rep2.bam && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:HI.2791.003.Index_5.MCF7_WCE_C_rep2.bam	SM:MCF7_CTRL_WCE_rep2	LB:MCF7_CTRL_WCE_rep2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/MCF7_CTRL_WCE_rep2/HI.2791.003.Index_5.MCF7_WCE_C_rep2.bam.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/MCF7_CTRL_WCE_rep2/HI.2791.003.Index_5.MCF7_WCE_C_rep2.bam/HI.2791.003.Index_5.MCF7_WCE_C_rep2.bam.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.HI.2791.003.Index_5.MCF7_WCE_C_rep2.bam.acb320d6006f526b46b4e4c5d864cab4.mugqic.done
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
# JOB: bwa_mem_picard_sort_sam_14_JOB_ID: bwa_mem_picard_sort_sam.HI.2791.003.Index_9.MCF7_1_NIPBL_C_rep2.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.2791.003.Index_9.MCF7_1_NIPBL_C_rep2.bam
JOB_DEPENDENCIES=$trimmomatic_14_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.2791.003.Index_9.MCF7_1_NIPBL_C_rep2.bam.15bb8474dd42777b8b50297a98f63d17.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.2791.003.Index_9.MCF7_1_NIPBL_C_rep2.bam.15bb8474dd42777b8b50297a98f63d17.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p alignment/MCF7_CTRL_NIPBL_rep2/HI.2791.003.Index_9.MCF7_1_NIPBL_C_rep2.bam && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:HI.2791.003.Index_9.MCF7_1_NIPBL_C_rep2.bam	SM:MCF7_CTRL_NIPBL_rep2	LB:MCF7_CTRL_NIPBL_rep2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/MCF7_CTRL_NIPBL_rep2/HI.2791.003.Index_9.MCF7_1_NIPBL_C_rep2.bam.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/MCF7_CTRL_NIPBL_rep2/HI.2791.003.Index_9.MCF7_1_NIPBL_C_rep2.bam/HI.2791.003.Index_9.MCF7_1_NIPBL_C_rep2.bam.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.HI.2791.003.Index_9.MCF7_1_NIPBL_C_rep2.bam.15bb8474dd42777b8b50297a98f63d17.mugqic.done
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=64G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$bwa_mem_picard_sort_sam_14_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_15_JOB_ID: bwa_mem_picard_sort_sam.HI.2791.003.Index_11.MCF7_4_Smc1_C_rep2.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.2791.003.Index_11.MCF7_4_Smc1_C_rep2.bam
JOB_DEPENDENCIES=$trimmomatic_15_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.2791.003.Index_11.MCF7_4_Smc1_C_rep2.bam.3d183fa56a982f80862c7226a28949ec.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.2791.003.Index_11.MCF7_4_Smc1_C_rep2.bam.3d183fa56a982f80862c7226a28949ec.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p alignment/MCF7_CTRL_SMC1A_rep2/HI.2791.003.Index_11.MCF7_4_Smc1_C_rep2.bam && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:HI.2791.003.Index_11.MCF7_4_Smc1_C_rep2.bam	SM:MCF7_CTRL_SMC1A_rep2	LB:MCF7_CTRL_SMC1A_rep2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/MCF7_CTRL_SMC1A_rep2/HI.2791.003.Index_11.MCF7_4_Smc1_C_rep2.bam.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/MCF7_CTRL_SMC1A_rep2/HI.2791.003.Index_11.MCF7_4_Smc1_C_rep2.bam/HI.2791.003.Index_11.MCF7_4_Smc1_C_rep2.bam.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.HI.2791.003.Index_11.MCF7_4_Smc1_C_rep2.bam.3d183fa56a982f80862c7226a28949ec.mugqic.done
)
bwa_mem_picard_sort_sam_15_JOB_ID=$(echo "#! /bin/bash
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
echo "$bwa_mem_picard_sort_sam_15_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_16_JOB_ID: bwa_mem_picard_sort_sam.HI.2791.003.Index_8.MCF7_2_MED1_C_rep2.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.2791.003.Index_8.MCF7_2_MED1_C_rep2.bam
JOB_DEPENDENCIES=$trimmomatic_16_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.2791.003.Index_8.MCF7_2_MED1_C_rep2.bam.a0d659649d749993a4b620b9d87d76cb.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.2791.003.Index_8.MCF7_2_MED1_C_rep2.bam.a0d659649d749993a4b620b9d87d76cb.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p alignment/MCF7_CTRL_MED1_rep2/HI.2791.003.Index_8.MCF7_2_MED1_C_rep2.bam && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:HI.2791.003.Index_8.MCF7_2_MED1_C_rep2.bam	SM:MCF7_CTRL_MED1_rep2	LB:MCF7_CTRL_MED1_rep2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/MCF7_CTRL_MED1_rep2/HI.2791.003.Index_8.MCF7_2_MED1_C_rep2.bam.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/MCF7_CTRL_MED1_rep2/HI.2791.003.Index_8.MCF7_2_MED1_C_rep2.bam/HI.2791.003.Index_8.MCF7_2_MED1_C_rep2.bam.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.HI.2791.003.Index_8.MCF7_2_MED1_C_rep2.bam.a0d659649d749993a4b620b9d87d76cb.mugqic.done
)
bwa_mem_picard_sort_sam_16_JOB_ID=$(echo "#! /bin/bash
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
echo "$bwa_mem_picard_sort_sam_16_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_17_JOB_ID: bwa_mem_picard_sort_sam.HI.2791.003.Index_10.MCF7_3_Era_C_rep2.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.2791.003.Index_10.MCF7_3_Era_C_rep2.bam
JOB_DEPENDENCIES=$trimmomatic_17_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.2791.003.Index_10.MCF7_3_Era_C_rep2.bam.396b9f6ad0b3ac347b7c97df88d55ae9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.2791.003.Index_10.MCF7_3_Era_C_rep2.bam.396b9f6ad0b3ac347b7c97df88d55ae9.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p alignment/MCF7_CTRL_ERA_rep2/HI.2791.003.Index_10.MCF7_3_Era_C_rep2.bam && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:HI.2791.003.Index_10.MCF7_3_Era_C_rep2.bam	SM:MCF7_CTRL_ERA_rep2	LB:MCF7_CTRL_ERA_rep2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/MCF7_CTRL_ERA_rep2/HI.2791.003.Index_10.MCF7_3_Era_C_rep2.bam.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/MCF7_CTRL_ERA_rep2/HI.2791.003.Index_10.MCF7_3_Era_C_rep2.bam/HI.2791.003.Index_10.MCF7_3_Era_C_rep2.bam.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.HI.2791.003.Index_10.MCF7_3_Era_C_rep2.bam.396b9f6ad0b3ac347b7c97df88d55ae9.mugqic.done
)
bwa_mem_picard_sort_sam_17_JOB_ID=$(echo "#! /bin/bash
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
echo "$bwa_mem_picard_sort_sam_17_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_18_JOB_ID: bwa_mem_picard_sort_sam.HI.2791.003.Index_4.MCF7_WCE_E_rep2.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.2791.003.Index_4.MCF7_WCE_E_rep2.bam
JOB_DEPENDENCIES=$trimmomatic_18_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.2791.003.Index_4.MCF7_WCE_E_rep2.bam.bf675017e80aebb833176a00bda4603b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.2791.003.Index_4.MCF7_WCE_E_rep2.bam.bf675017e80aebb833176a00bda4603b.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p alignment/MCF7_E2_WCE_rep2/HI.2791.003.Index_4.MCF7_WCE_E_rep2.bam && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:HI.2791.003.Index_4.MCF7_WCE_E_rep2.bam	SM:MCF7_E2_WCE_rep2	LB:MCF7_E2_WCE_rep2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/MCF7_E2_WCE_rep2/HI.2791.003.Index_4.MCF7_WCE_E_rep2.bam.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/MCF7_E2_WCE_rep2/HI.2791.003.Index_4.MCF7_WCE_E_rep2.bam/HI.2791.003.Index_4.MCF7_WCE_E_rep2.bam.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.HI.2791.003.Index_4.MCF7_WCE_E_rep2.bam.bf675017e80aebb833176a00bda4603b.mugqic.done
)
bwa_mem_picard_sort_sam_18_JOB_ID=$(echo "#! /bin/bash
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
echo "$bwa_mem_picard_sort_sam_18_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_19_JOB_ID: bwa_mem_picard_sort_sam.HI.2791.003.Index_12.MCF7_1_NIPBL_E_rep2.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.2791.003.Index_12.MCF7_1_NIPBL_E_rep2.bam
JOB_DEPENDENCIES=$trimmomatic_19_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.2791.003.Index_12.MCF7_1_NIPBL_E_rep2.bam.fefc3b7256a25438531b6680b146b73e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.2791.003.Index_12.MCF7_1_NIPBL_E_rep2.bam.fefc3b7256a25438531b6680b146b73e.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p alignment/MCF7_E2_NIPBL_rep2/HI.2791.003.Index_12.MCF7_1_NIPBL_E_rep2.bam && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:HI.2791.003.Index_12.MCF7_1_NIPBL_E_rep2.bam	SM:MCF7_E2_NIPBL_rep2	LB:MCF7_E2_NIPBL_rep2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/MCF7_E2_NIPBL_rep2/HI.2791.003.Index_12.MCF7_1_NIPBL_E_rep2.bam.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/MCF7_E2_NIPBL_rep2/HI.2791.003.Index_12.MCF7_1_NIPBL_E_rep2.bam/HI.2791.003.Index_12.MCF7_1_NIPBL_E_rep2.bam.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.HI.2791.003.Index_12.MCF7_1_NIPBL_E_rep2.bam.fefc3b7256a25438531b6680b146b73e.mugqic.done
)
bwa_mem_picard_sort_sam_19_JOB_ID=$(echo "#! /bin/bash
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
echo "$bwa_mem_picard_sort_sam_19_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_20_JOB_ID: bwa_mem_picard_sort_sam.HI.2791.003.Index_3.MCF7_4_Smc1_E_rep2.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.2791.003.Index_3.MCF7_4_Smc1_E_rep2.bam
JOB_DEPENDENCIES=$trimmomatic_20_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.2791.003.Index_3.MCF7_4_Smc1_E_rep2.bam.2eaa6d2e429585bded15aa3db9953b4a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.2791.003.Index_3.MCF7_4_Smc1_E_rep2.bam.2eaa6d2e429585bded15aa3db9953b4a.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p alignment/MCF7_E2_SMC1A_rep2/HI.2791.003.Index_3.MCF7_4_Smc1_E_rep2.bam && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:HI.2791.003.Index_3.MCF7_4_Smc1_E_rep2.bam	SM:MCF7_E2_SMC1A_rep2	LB:MCF7_E2_SMC1A_rep2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/MCF7_E2_SMC1A_rep2/HI.2791.003.Index_3.MCF7_4_Smc1_E_rep2.bam.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/MCF7_E2_SMC1A_rep2/HI.2791.003.Index_3.MCF7_4_Smc1_E_rep2.bam/HI.2791.003.Index_3.MCF7_4_Smc1_E_rep2.bam.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.HI.2791.003.Index_3.MCF7_4_Smc1_E_rep2.bam.2eaa6d2e429585bded15aa3db9953b4a.mugqic.done
)
bwa_mem_picard_sort_sam_20_JOB_ID=$(echo "#! /bin/bash
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
echo "$bwa_mem_picard_sort_sam_20_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_21_JOB_ID: bwa_mem_picard_sort_sam.HI.2791.003.Index_19.MCF7_2_MED1_E_rep2.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.2791.003.Index_19.MCF7_2_MED1_E_rep2.bam
JOB_DEPENDENCIES=$trimmomatic_21_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.2791.003.Index_19.MCF7_2_MED1_E_rep2.bam.a7acc8fbb3dd993b205d6b7dcad7833a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.2791.003.Index_19.MCF7_2_MED1_E_rep2.bam.a7acc8fbb3dd993b205d6b7dcad7833a.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p alignment/MCF7_E2_MED1_rep2/HI.2791.003.Index_19.MCF7_2_MED1_E_rep2.bam && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:HI.2791.003.Index_19.MCF7_2_MED1_E_rep2.bam	SM:MCF7_E2_MED1_rep2	LB:MCF7_E2_MED1_rep2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/MCF7_E2_MED1_rep2/HI.2791.003.Index_19.MCF7_2_MED1_E_rep2.bam.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/MCF7_E2_MED1_rep2/HI.2791.003.Index_19.MCF7_2_MED1_E_rep2.bam/HI.2791.003.Index_19.MCF7_2_MED1_E_rep2.bam.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.HI.2791.003.Index_19.MCF7_2_MED1_E_rep2.bam.a7acc8fbb3dd993b205d6b7dcad7833a.mugqic.done
)
bwa_mem_picard_sort_sam_21_JOB_ID=$(echo "#! /bin/bash
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
echo "$bwa_mem_picard_sort_sam_21_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_22_JOB_ID: bwa_mem_picard_sort_sam.HI.2791.003.Index_1.MCF7_3_Era_E_rep2.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.2791.003.Index_1.MCF7_3_Era_E_rep2.bam
JOB_DEPENDENCIES=$trimmomatic_22_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.2791.003.Index_1.MCF7_3_Era_E_rep2.bam.2dbf4594c93b1b1ced17e6d3168f4363.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.2791.003.Index_1.MCF7_3_Era_E_rep2.bam.2dbf4594c93b1b1ced17e6d3168f4363.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p alignment/MCF7_E2_ERA_rep2/HI.2791.003.Index_1.MCF7_3_Era_E_rep2.bam && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:HI.2791.003.Index_1.MCF7_3_Era_E_rep2.bam	SM:MCF7_E2_ERA_rep2	LB:MCF7_E2_ERA_rep2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/MCF7_E2_ERA_rep2/HI.2791.003.Index_1.MCF7_3_Era_E_rep2.bam.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/MCF7_E2_ERA_rep2/HI.2791.003.Index_1.MCF7_3_Era_E_rep2.bam/HI.2791.003.Index_1.MCF7_3_Era_E_rep2.bam.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.HI.2791.003.Index_1.MCF7_3_Era_E_rep2.bam.2dbf4594c93b1b1ced17e6d3168f4363.mugqic.done
)
bwa_mem_picard_sort_sam_22_JOB_ID=$(echo "#! /bin/bash
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
echo "$bwa_mem_picard_sort_sam_22_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_23_JOB_ID: bwa_mem_picard_sort_sam.SRR1193526.fastq
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.SRR1193526.fastq
JOB_DEPENDENCIES=$trimmomatic_23_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.SRR1193526.fastq.5c86b630fd5a631c377c441abdffb309.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.SRR1193526.fastq.5c86b630fd5a631c377c441abdffb309.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p alignment/MCF7_CTRL_BRD4_rep1/SRR1193526.fastq && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:SRR1193526.fastq	SM:MCF7_CTRL_BRD4_rep1	LB:MCF7_CTRL_BRD4_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/MCF7_CTRL_BRD4_rep1/SRR1193526.fastq.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/MCF7_CTRL_BRD4_rep1/SRR1193526.fastq/SRR1193526.fastq.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.SRR1193526.fastq.5c86b630fd5a631c377c441abdffb309.mugqic.done
)
bwa_mem_picard_sort_sam_23_JOB_ID=$(echo "#! /bin/bash
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
echo "$bwa_mem_picard_sort_sam_23_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_24_JOB_ID: bwa_mem_picard_sort_sam.SRR1193527.fastq
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.SRR1193527.fastq
JOB_DEPENDENCIES=$trimmomatic_24_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.SRR1193527.fastq.4f9f7cba1bc2bdbd3343344d3a43c137.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.SRR1193527.fastq.4f9f7cba1bc2bdbd3343344d3a43c137.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p alignment/MCF7_CTRL_BRD4_rep2/SRR1193527.fastq && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:SRR1193527.fastq	SM:MCF7_CTRL_BRD4_rep2	LB:MCF7_CTRL_BRD4_rep2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/MCF7_CTRL_BRD4_rep2/SRR1193527.fastq.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/MCF7_CTRL_BRD4_rep2/SRR1193527.fastq/SRR1193527.fastq.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.SRR1193527.fastq.4f9f7cba1bc2bdbd3343344d3a43c137.mugqic.done
)
bwa_mem_picard_sort_sam_24_JOB_ID=$(echo "#! /bin/bash
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
echo "$bwa_mem_picard_sort_sam_24_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_25_JOB_ID: bwa_mem_picard_sort_sam.SRR1193528.fastq
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.SRR1193528.fastq
JOB_DEPENDENCIES=$trimmomatic_25_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.SRR1193528.fastq.7e3061ad08589ea6025f61a654efd4a7.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.SRR1193528.fastq.7e3061ad08589ea6025f61a654efd4a7.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p alignment/MCF7_CTRL_BRD4_rep3/SRR1193528.fastq && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:SRR1193528.fastq	SM:MCF7_CTRL_BRD4_rep3	LB:MCF7_CTRL_BRD4_rep3	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/MCF7_CTRL_BRD4_rep3/SRR1193528.fastq.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/MCF7_CTRL_BRD4_rep3/SRR1193528.fastq/SRR1193528.fastq.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.SRR1193528.fastq.7e3061ad08589ea6025f61a654efd4a7.mugqic.done
)
bwa_mem_picard_sort_sam_25_JOB_ID=$(echo "#! /bin/bash
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
echo "$bwa_mem_picard_sort_sam_25_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_26_JOB_ID: bwa_mem_picard_sort_sam.SRR1193529.fastq
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.SRR1193529.fastq
JOB_DEPENDENCIES=$trimmomatic_26_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.SRR1193529.fastq.395e023ad587dcb21b6c76e86ee7549e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.SRR1193529.fastq.395e023ad587dcb21b6c76e86ee7549e.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p alignment/MCF7_E2_BRD4_rep1/SRR1193529.fastq && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:SRR1193529.fastq	SM:MCF7_E2_BRD4_rep1	LB:MCF7_E2_BRD4_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/MCF7_E2_BRD4_rep1/SRR1193529.fastq.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/MCF7_E2_BRD4_rep1/SRR1193529.fastq/SRR1193529.fastq.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.SRR1193529.fastq.395e023ad587dcb21b6c76e86ee7549e.mugqic.done
)
bwa_mem_picard_sort_sam_26_JOB_ID=$(echo "#! /bin/bash
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
echo "$bwa_mem_picard_sort_sam_26_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_27_JOB_ID: bwa_mem_picard_sort_sam.SRR1193530.fastq
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.SRR1193530.fastq
JOB_DEPENDENCIES=$trimmomatic_27_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.SRR1193530.fastq.edf375cd91721168a7e87af4b146bf77.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.SRR1193530.fastq.edf375cd91721168a7e87af4b146bf77.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p alignment/MCF7_E2_BRD4_rep2/SRR1193530.fastq && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:SRR1193530.fastq	SM:MCF7_E2_BRD4_rep2	LB:MCF7_E2_BRD4_rep2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/MCF7_E2_BRD4_rep2/SRR1193530.fastq.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/MCF7_E2_BRD4_rep2/SRR1193530.fastq/SRR1193530.fastq.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.SRR1193530.fastq.edf375cd91721168a7e87af4b146bf77.mugqic.done
)
bwa_mem_picard_sort_sam_27_JOB_ID=$(echo "#! /bin/bash
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
echo "$bwa_mem_picard_sort_sam_27_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_28_JOB_ID: bwa_mem_picard_sort_sam.SRR1193531.fastq
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.SRR1193531.fastq
JOB_DEPENDENCIES=$trimmomatic_28_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.SRR1193531.fastq.acee69de539c6a721b3d006075a80072.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.SRR1193531.fastq.acee69de539c6a721b3d006075a80072.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p alignment/MCF7_E2_BRD4_rep3/SRR1193531.fastq && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:SRR1193531.fastq	SM:MCF7_E2_BRD4_rep3	LB:MCF7_E2_BRD4_rep3	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/MCF7_E2_BRD4_rep3/SRR1193531.fastq.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/MCF7_E2_BRD4_rep3/SRR1193531.fastq/SRR1193531.fastq.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.SRR1193531.fastq.acee69de539c6a721b3d006075a80072.mugqic.done
)
bwa_mem_picard_sort_sam_28_JOB_ID=$(echo "#! /bin/bash
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
echo "$bwa_mem_picard_sort_sam_28_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_29_JOB_ID: bwa_mem_picard_sort_sam.SRR1193562.fastq
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.SRR1193562.fastq
JOB_DEPENDENCIES=$trimmomatic_29_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.SRR1193562.fastq.746a5e4e8aff132851e508d3fbe3ce13.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.SRR1193562.fastq.746a5e4e8aff132851e508d3fbe3ce13.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p alignment/MCF7_CTRL_WCE/SRR1193562.fastq && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:SRR1193562.fastq	SM:MCF7_CTRL_WCE	LB:MCF7_CTRL_WCE	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/MCF7_CTRL_WCE/SRR1193562.fastq.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/MCF7_CTRL_WCE/SRR1193562.fastq/SRR1193562.fastq.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.SRR1193562.fastq.746a5e4e8aff132851e508d3fbe3ce13.mugqic.done
)
bwa_mem_picard_sort_sam_29_JOB_ID=$(echo "#! /bin/bash
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
echo "$bwa_mem_picard_sort_sam_29_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_30_JOB_ID: bwa_mem_picard_sort_sam.SRR1193563.fastq
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.SRR1193563.fastq
JOB_DEPENDENCIES=$trimmomatic_30_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.SRR1193563.fastq.3f478eaf05003fd262d123b493772135.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.SRR1193563.fastq.3f478eaf05003fd262d123b493772135.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p alignment/MCF7_E2_WCE/SRR1193563.fastq && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:SRR1193563.fastq	SM:MCF7_E2_WCE	LB:MCF7_E2_WCE	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/MCF7_E2_WCE/SRR1193563.fastq.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/MCF7_E2_WCE/SRR1193563.fastq/SRR1193563.fastq.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.SRR1193563.fastq.3f478eaf05003fd262d123b493772135.mugqic.done
)
bwa_mem_picard_sort_sam_30_JOB_ID=$(echo "#! /bin/bash
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
echo "$bwa_mem_picard_sort_sam_30_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_31_JOB_ID: bwa_mem_picard_sort_sam_report
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam_report
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_1_JOB_ID:$bwa_mem_picard_sort_sam_2_JOB_ID:$bwa_mem_picard_sort_sam_3_JOB_ID:$bwa_mem_picard_sort_sam_4_JOB_ID:$bwa_mem_picard_sort_sam_5_JOB_ID:$bwa_mem_picard_sort_sam_6_JOB_ID:$bwa_mem_picard_sort_sam_7_JOB_ID:$bwa_mem_picard_sort_sam_8_JOB_ID:$bwa_mem_picard_sort_sam_9_JOB_ID:$bwa_mem_picard_sort_sam_10_JOB_ID:$bwa_mem_picard_sort_sam_11_JOB_ID:$bwa_mem_picard_sort_sam_12_JOB_ID:$bwa_mem_picard_sort_sam_13_JOB_ID:$bwa_mem_picard_sort_sam_14_JOB_ID:$bwa_mem_picard_sort_sam_15_JOB_ID:$bwa_mem_picard_sort_sam_16_JOB_ID:$bwa_mem_picard_sort_sam_17_JOB_ID:$bwa_mem_picard_sort_sam_18_JOB_ID:$bwa_mem_picard_sort_sam_19_JOB_ID:$bwa_mem_picard_sort_sam_20_JOB_ID:$bwa_mem_picard_sort_sam_21_JOB_ID:$bwa_mem_picard_sort_sam_22_JOB_ID:$bwa_mem_picard_sort_sam_23_JOB_ID:$bwa_mem_picard_sort_sam_24_JOB_ID:$bwa_mem_picard_sort_sam_25_JOB_ID:$bwa_mem_picard_sort_sam_26_JOB_ID:$bwa_mem_picard_sort_sam_27_JOB_ID:$bwa_mem_picard_sort_sam_28_JOB_ID:$bwa_mem_picard_sort_sam_29_JOB_ID:$bwa_mem_picard_sort_sam_30_JOB_ID
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
bwa_mem_picard_sort_sam_31_JOB_ID=$(echo "#! /bin/bash
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
echo "$bwa_mem_picard_sort_sam_31_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# STEP: samtools_view_filter
#-------------------------------------------------------------------------------
STEP=samtools_view_filter
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_1_JOB_ID: samtools_view_filter.HI.1943.008.Index_2.MCF7_CTRL_WCE_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.1943.008.Index_2.MCF7_CTRL_WCE_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_1_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.1943.008.Index_2.MCF7_CTRL_WCE_ChIP8_GB_rep1.bam.08b09694537291c3493422a2612e1d5c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.1943.008.Index_2.MCF7_CTRL_WCE_ChIP8_GB_rep1.bam.08b09694537291c3493422a2612e1d5c.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/MCF7_CTRL_WCE_rep1/HI.1943.008.Index_2.MCF7_CTRL_WCE_ChIP8_GB_rep1.bam/HI.1943.008.Index_2.MCF7_CTRL_WCE_ChIP8_GB_rep1.bam.sorted.bam \
  > alignment/MCF7_CTRL_WCE_rep1/HI.1943.008.Index_2.MCF7_CTRL_WCE_ChIP8_GB_rep1.bam/HI.1943.008.Index_2.MCF7_CTRL_WCE_ChIP8_GB_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.1943.008.Index_2.MCF7_CTRL_WCE_ChIP8_GB_rep1.bam.08b09694537291c3493422a2612e1d5c.mugqic.done
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
# JOB: samtools_view_filter_2_JOB_ID: samtools_view_filter.HI.1997.004.Index_1.MCF7_CTRL_NIPBL_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.1997.004.Index_1.MCF7_CTRL_NIPBL_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_2_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.1997.004.Index_1.MCF7_CTRL_NIPBL_ChIP8_GB_rep1.bam.09271f42b3976ff5d65552272c1038e1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.1997.004.Index_1.MCF7_CTRL_NIPBL_ChIP8_GB_rep1.bam.09271f42b3976ff5d65552272c1038e1.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/MCF7_CTRL_NIPBL_rep1/HI.1997.004.Index_1.MCF7_CTRL_NIPBL_ChIP8_GB_rep1.bam/HI.1997.004.Index_1.MCF7_CTRL_NIPBL_ChIP8_GB_rep1.bam.sorted.bam \
  > alignment/MCF7_CTRL_NIPBL_rep1/HI.1997.004.Index_1.MCF7_CTRL_NIPBL_ChIP8_GB_rep1.bam/HI.1997.004.Index_1.MCF7_CTRL_NIPBL_ChIP8_GB_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.1997.004.Index_1.MCF7_CTRL_NIPBL_ChIP8_GB_rep1.bam.09271f42b3976ff5d65552272c1038e1.mugqic.done
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
# JOB: samtools_view_filter_3_JOB_ID: samtools_view_filter.HI.1997.004.Index_4.MCF7_CTRL_SMC1_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.1997.004.Index_4.MCF7_CTRL_SMC1_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_3_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.1997.004.Index_4.MCF7_CTRL_SMC1_ChIP8_GB_rep1.bam.0055033faa228aecabe9676d9a31894d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.1997.004.Index_4.MCF7_CTRL_SMC1_ChIP8_GB_rep1.bam.0055033faa228aecabe9676d9a31894d.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/MCF7_CTRL_SMC1A_rep1/HI.1997.004.Index_4.MCF7_CTRL_SMC1_ChIP8_GB_rep1.bam/HI.1997.004.Index_4.MCF7_CTRL_SMC1_ChIP8_GB_rep1.bam.sorted.bam \
  > alignment/MCF7_CTRL_SMC1A_rep1/HI.1997.004.Index_4.MCF7_CTRL_SMC1_ChIP8_GB_rep1.bam/HI.1997.004.Index_4.MCF7_CTRL_SMC1_ChIP8_GB_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.1997.004.Index_4.MCF7_CTRL_SMC1_ChIP8_GB_rep1.bam.0055033faa228aecabe9676d9a31894d.mugqic.done
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
# JOB: samtools_view_filter_4_JOB_ID: samtools_view_filter.HI.1943.008.Index_7.MCF7_CTRL_MED1_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.1943.008.Index_7.MCF7_CTRL_MED1_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_4_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.1943.008.Index_7.MCF7_CTRL_MED1_ChIP8_GB_rep1.bam.02ec6f9f7da6a465bbbcf25c9590c9eb.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.1943.008.Index_7.MCF7_CTRL_MED1_ChIP8_GB_rep1.bam.02ec6f9f7da6a465bbbcf25c9590c9eb.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/MCF7_CTRL_MED1_rep1/HI.1943.008.Index_7.MCF7_CTRL_MED1_ChIP8_GB_rep1.bam/HI.1943.008.Index_7.MCF7_CTRL_MED1_ChIP8_GB_rep1.bam.sorted.bam \
  > alignment/MCF7_CTRL_MED1_rep1/HI.1943.008.Index_7.MCF7_CTRL_MED1_ChIP8_GB_rep1.bam/HI.1943.008.Index_7.MCF7_CTRL_MED1_ChIP8_GB_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.1943.008.Index_7.MCF7_CTRL_MED1_ChIP8_GB_rep1.bam.02ec6f9f7da6a465bbbcf25c9590c9eb.mugqic.done
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
# JOB: samtools_view_filter_5_JOB_ID: samtools_view_filter.HI.1997.005.Index_10.MCF7_CTRL_POL2_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.1997.005.Index_10.MCF7_CTRL_POL2_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_5_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.1997.005.Index_10.MCF7_CTRL_POL2_ChIP8_GB_rep1.bam.af8cd0cd78d80a542be03bb46c439428.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.1997.005.Index_10.MCF7_CTRL_POL2_ChIP8_GB_rep1.bam.af8cd0cd78d80a542be03bb46c439428.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/MCF7_CTRL_POL2_rep1/HI.1997.005.Index_10.MCF7_CTRL_POL2_ChIP8_GB_rep1.bam/HI.1997.005.Index_10.MCF7_CTRL_POL2_ChIP8_GB_rep1.bam.sorted.bam \
  > alignment/MCF7_CTRL_POL2_rep1/HI.1997.005.Index_10.MCF7_CTRL_POL2_ChIP8_GB_rep1.bam/HI.1997.005.Index_10.MCF7_CTRL_POL2_ChIP8_GB_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.1997.005.Index_10.MCF7_CTRL_POL2_ChIP8_GB_rep1.bam.af8cd0cd78d80a542be03bb46c439428.mugqic.done
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
# JOB: samtools_view_filter_6_JOB_ID: samtools_view_filter.HI.1997.005.Index_25.MCF7_CTRL_ERA_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.1997.005.Index_25.MCF7_CTRL_ERA_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_6_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.1997.005.Index_25.MCF7_CTRL_ERA_ChIP8_GB_rep1.bam.92d40c95c17c4a546241ed7275837005.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.1997.005.Index_25.MCF7_CTRL_ERA_ChIP8_GB_rep1.bam.92d40c95c17c4a546241ed7275837005.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/MCF7_CTRL_ERA_rep1/HI.1997.005.Index_25.MCF7_CTRL_ERA_ChIP8_GB_rep1.bam/HI.1997.005.Index_25.MCF7_CTRL_ERA_ChIP8_GB_rep1.bam.sorted.bam \
  > alignment/MCF7_CTRL_ERA_rep1/HI.1997.005.Index_25.MCF7_CTRL_ERA_ChIP8_GB_rep1.bam/HI.1997.005.Index_25.MCF7_CTRL_ERA_ChIP8_GB_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.1997.005.Index_25.MCF7_CTRL_ERA_ChIP8_GB_rep1.bam.92d40c95c17c4a546241ed7275837005.mugqic.done
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
# JOB: samtools_view_filter_7_JOB_ID: samtools_view_filter.HI.1943.008.Index_13.MCF7_E2_WCE_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.1943.008.Index_13.MCF7_E2_WCE_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_7_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.1943.008.Index_13.MCF7_E2_WCE_ChIP8_GB_rep1.bam.5855af65b663672028bc22468bbfe15d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.1943.008.Index_13.MCF7_E2_WCE_ChIP8_GB_rep1.bam.5855af65b663672028bc22468bbfe15d.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/MCF7_E2_WCE_rep1/HI.1943.008.Index_13.MCF7_E2_WCE_ChIP8_GB_rep1.bam/HI.1943.008.Index_13.MCF7_E2_WCE_ChIP8_GB_rep1.bam.sorted.bam \
  > alignment/MCF7_E2_WCE_rep1/HI.1943.008.Index_13.MCF7_E2_WCE_ChIP8_GB_rep1.bam/HI.1943.008.Index_13.MCF7_E2_WCE_ChIP8_GB_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.1943.008.Index_13.MCF7_E2_WCE_ChIP8_GB_rep1.bam.5855af65b663672028bc22468bbfe15d.mugqic.done
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
# JOB: samtools_view_filter_8_JOB_ID: samtools_view_filter.HI.1997.004.Index_3.MCF7_E2_NIPBL_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.1997.004.Index_3.MCF7_E2_NIPBL_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_8_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.1997.004.Index_3.MCF7_E2_NIPBL_ChIP8_GB_rep1.bam.6ff85ffbc0659b1bd722d0d9a4278f6c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.1997.004.Index_3.MCF7_E2_NIPBL_ChIP8_GB_rep1.bam.6ff85ffbc0659b1bd722d0d9a4278f6c.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/MCF7_E2_NIPBL_rep1/HI.1997.004.Index_3.MCF7_E2_NIPBL_ChIP8_GB_rep1.bam/HI.1997.004.Index_3.MCF7_E2_NIPBL_ChIP8_GB_rep1.bam.sorted.bam \
  > alignment/MCF7_E2_NIPBL_rep1/HI.1997.004.Index_3.MCF7_E2_NIPBL_ChIP8_GB_rep1.bam/HI.1997.004.Index_3.MCF7_E2_NIPBL_ChIP8_GB_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.1997.004.Index_3.MCF7_E2_NIPBL_ChIP8_GB_rep1.bam.6ff85ffbc0659b1bd722d0d9a4278f6c.mugqic.done
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
# JOB: samtools_view_filter_9_JOB_ID: samtools_view_filter.HI.1997.004.Index_5.MCF7_E2_SMC1_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.1997.004.Index_5.MCF7_E2_SMC1_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_9_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.1997.004.Index_5.MCF7_E2_SMC1_ChIP8_GB_rep1.bam.08391c5ccba36e361849ea5f9bd77342.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.1997.004.Index_5.MCF7_E2_SMC1_ChIP8_GB_rep1.bam.08391c5ccba36e361849ea5f9bd77342.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/MCF7_E2_SMC1A_rep1/HI.1997.004.Index_5.MCF7_E2_SMC1_ChIP8_GB_rep1.bam/HI.1997.004.Index_5.MCF7_E2_SMC1_ChIP8_GB_rep1.bam.sorted.bam \
  > alignment/MCF7_E2_SMC1A_rep1/HI.1997.004.Index_5.MCF7_E2_SMC1_ChIP8_GB_rep1.bam/HI.1997.004.Index_5.MCF7_E2_SMC1_ChIP8_GB_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.1997.004.Index_5.MCF7_E2_SMC1_ChIP8_GB_rep1.bam.08391c5ccba36e361849ea5f9bd77342.mugqic.done
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
# JOB: samtools_view_filter_10_JOB_ID: samtools_view_filter.HI.1943.008.Index_18.MCF7_E2_MED1_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.1943.008.Index_18.MCF7_E2_MED1_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_10_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.1943.008.Index_18.MCF7_E2_MED1_ChIP8_GB_rep1.bam.0c0d358cc58358d0945f51133fbc5ab8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.1943.008.Index_18.MCF7_E2_MED1_ChIP8_GB_rep1.bam.0c0d358cc58358d0945f51133fbc5ab8.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/MCF7_E2_MED1_rep1/HI.1943.008.Index_18.MCF7_E2_MED1_ChIP8_GB_rep1.bam/HI.1943.008.Index_18.MCF7_E2_MED1_ChIP8_GB_rep1.bam.sorted.bam \
  > alignment/MCF7_E2_MED1_rep1/HI.1943.008.Index_18.MCF7_E2_MED1_ChIP8_GB_rep1.bam/HI.1943.008.Index_18.MCF7_E2_MED1_ChIP8_GB_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.1943.008.Index_18.MCF7_E2_MED1_ChIP8_GB_rep1.bam.0c0d358cc58358d0945f51133fbc5ab8.mugqic.done
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
# JOB: samtools_view_filter_11_JOB_ID: samtools_view_filter.HI.1997.005.Index_11.MCF7_E2_POL2_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.1997.005.Index_11.MCF7_E2_POL2_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_11_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.1997.005.Index_11.MCF7_E2_POL2_ChIP8_GB_rep1.bam.7ec514240ec217ab2be29d4fe5dc80fe.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.1997.005.Index_11.MCF7_E2_POL2_ChIP8_GB_rep1.bam.7ec514240ec217ab2be29d4fe5dc80fe.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/MCF7_E2_POL2_rep1/HI.1997.005.Index_11.MCF7_E2_POL2_ChIP8_GB_rep1.bam/HI.1997.005.Index_11.MCF7_E2_POL2_ChIP8_GB_rep1.bam.sorted.bam \
  > alignment/MCF7_E2_POL2_rep1/HI.1997.005.Index_11.MCF7_E2_POL2_ChIP8_GB_rep1.bam/HI.1997.005.Index_11.MCF7_E2_POL2_ChIP8_GB_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.1997.005.Index_11.MCF7_E2_POL2_ChIP8_GB_rep1.bam.7ec514240ec217ab2be29d4fe5dc80fe.mugqic.done
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
# JOB: samtools_view_filter_12_JOB_ID: samtools_view_filter.HI.1997.005.Index_21.MCF7_E2_ERA_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.1997.005.Index_21.MCF7_E2_ERA_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_12_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.1997.005.Index_21.MCF7_E2_ERA_ChIP8_GB_rep1.bam.017c84148833bbc39002ffe9fc08f8b6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.1997.005.Index_21.MCF7_E2_ERA_ChIP8_GB_rep1.bam.017c84148833bbc39002ffe9fc08f8b6.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/MCF7_E2_ERA_rep1/HI.1997.005.Index_21.MCF7_E2_ERA_ChIP8_GB_rep1.bam/HI.1997.005.Index_21.MCF7_E2_ERA_ChIP8_GB_rep1.bam.sorted.bam \
  > alignment/MCF7_E2_ERA_rep1/HI.1997.005.Index_21.MCF7_E2_ERA_ChIP8_GB_rep1.bam/HI.1997.005.Index_21.MCF7_E2_ERA_ChIP8_GB_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.1997.005.Index_21.MCF7_E2_ERA_ChIP8_GB_rep1.bam.017c84148833bbc39002ffe9fc08f8b6.mugqic.done
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
# JOB: samtools_view_filter_13_JOB_ID: samtools_view_filter.HI.2791.003.Index_5.MCF7_WCE_C_rep2.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.2791.003.Index_5.MCF7_WCE_C_rep2.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_13_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.2791.003.Index_5.MCF7_WCE_C_rep2.bam.285109b581d01b992304b770cf486110.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.2791.003.Index_5.MCF7_WCE_C_rep2.bam.285109b581d01b992304b770cf486110.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/MCF7_CTRL_WCE_rep2/HI.2791.003.Index_5.MCF7_WCE_C_rep2.bam/HI.2791.003.Index_5.MCF7_WCE_C_rep2.bam.sorted.bam \
  > alignment/MCF7_CTRL_WCE_rep2/HI.2791.003.Index_5.MCF7_WCE_C_rep2.bam/HI.2791.003.Index_5.MCF7_WCE_C_rep2.bam.sorted.filtered.bam
samtools_view_filter.HI.2791.003.Index_5.MCF7_WCE_C_rep2.bam.285109b581d01b992304b770cf486110.mugqic.done
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
# JOB: samtools_view_filter_14_JOB_ID: samtools_view_filter.HI.2791.003.Index_9.MCF7_1_NIPBL_C_rep2.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.2791.003.Index_9.MCF7_1_NIPBL_C_rep2.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_14_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.2791.003.Index_9.MCF7_1_NIPBL_C_rep2.bam.1855dccd57aa5f85829c354165e9d47a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.2791.003.Index_9.MCF7_1_NIPBL_C_rep2.bam.1855dccd57aa5f85829c354165e9d47a.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/MCF7_CTRL_NIPBL_rep2/HI.2791.003.Index_9.MCF7_1_NIPBL_C_rep2.bam/HI.2791.003.Index_9.MCF7_1_NIPBL_C_rep2.bam.sorted.bam \
  > alignment/MCF7_CTRL_NIPBL_rep2/HI.2791.003.Index_9.MCF7_1_NIPBL_C_rep2.bam/HI.2791.003.Index_9.MCF7_1_NIPBL_C_rep2.bam.sorted.filtered.bam
samtools_view_filter.HI.2791.003.Index_9.MCF7_1_NIPBL_C_rep2.bam.1855dccd57aa5f85829c354165e9d47a.mugqic.done
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
# JOB: samtools_view_filter_15_JOB_ID: samtools_view_filter.HI.2791.003.Index_11.MCF7_4_Smc1_C_rep2.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.2791.003.Index_11.MCF7_4_Smc1_C_rep2.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_15_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.2791.003.Index_11.MCF7_4_Smc1_C_rep2.bam.bc6ba516fc322d456e3824ee9786e789.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.2791.003.Index_11.MCF7_4_Smc1_C_rep2.bam.bc6ba516fc322d456e3824ee9786e789.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/MCF7_CTRL_SMC1A_rep2/HI.2791.003.Index_11.MCF7_4_Smc1_C_rep2.bam/HI.2791.003.Index_11.MCF7_4_Smc1_C_rep2.bam.sorted.bam \
  > alignment/MCF7_CTRL_SMC1A_rep2/HI.2791.003.Index_11.MCF7_4_Smc1_C_rep2.bam/HI.2791.003.Index_11.MCF7_4_Smc1_C_rep2.bam.sorted.filtered.bam
samtools_view_filter.HI.2791.003.Index_11.MCF7_4_Smc1_C_rep2.bam.bc6ba516fc322d456e3824ee9786e789.mugqic.done
)
samtools_view_filter_15_JOB_ID=$(echo "#! /bin/bash
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
echo "$samtools_view_filter_15_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_16_JOB_ID: samtools_view_filter.HI.2791.003.Index_8.MCF7_2_MED1_C_rep2.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.2791.003.Index_8.MCF7_2_MED1_C_rep2.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_16_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.2791.003.Index_8.MCF7_2_MED1_C_rep2.bam.e7287d1a03169d6465df22648020f1f4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.2791.003.Index_8.MCF7_2_MED1_C_rep2.bam.e7287d1a03169d6465df22648020f1f4.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/MCF7_CTRL_MED1_rep2/HI.2791.003.Index_8.MCF7_2_MED1_C_rep2.bam/HI.2791.003.Index_8.MCF7_2_MED1_C_rep2.bam.sorted.bam \
  > alignment/MCF7_CTRL_MED1_rep2/HI.2791.003.Index_8.MCF7_2_MED1_C_rep2.bam/HI.2791.003.Index_8.MCF7_2_MED1_C_rep2.bam.sorted.filtered.bam
samtools_view_filter.HI.2791.003.Index_8.MCF7_2_MED1_C_rep2.bam.e7287d1a03169d6465df22648020f1f4.mugqic.done
)
samtools_view_filter_16_JOB_ID=$(echo "#! /bin/bash
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
echo "$samtools_view_filter_16_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_17_JOB_ID: samtools_view_filter.HI.2791.003.Index_10.MCF7_3_Era_C_rep2.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.2791.003.Index_10.MCF7_3_Era_C_rep2.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_17_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.2791.003.Index_10.MCF7_3_Era_C_rep2.bam.7ff572a0cd8b52283d0c4830efa36cca.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.2791.003.Index_10.MCF7_3_Era_C_rep2.bam.7ff572a0cd8b52283d0c4830efa36cca.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/MCF7_CTRL_ERA_rep2/HI.2791.003.Index_10.MCF7_3_Era_C_rep2.bam/HI.2791.003.Index_10.MCF7_3_Era_C_rep2.bam.sorted.bam \
  > alignment/MCF7_CTRL_ERA_rep2/HI.2791.003.Index_10.MCF7_3_Era_C_rep2.bam/HI.2791.003.Index_10.MCF7_3_Era_C_rep2.bam.sorted.filtered.bam
samtools_view_filter.HI.2791.003.Index_10.MCF7_3_Era_C_rep2.bam.7ff572a0cd8b52283d0c4830efa36cca.mugqic.done
)
samtools_view_filter_17_JOB_ID=$(echo "#! /bin/bash
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
echo "$samtools_view_filter_17_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_18_JOB_ID: samtools_view_filter.HI.2791.003.Index_4.MCF7_WCE_E_rep2.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.2791.003.Index_4.MCF7_WCE_E_rep2.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_18_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.2791.003.Index_4.MCF7_WCE_E_rep2.bam.9c6aa82dc6ceb299d7546e48587f3e6c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.2791.003.Index_4.MCF7_WCE_E_rep2.bam.9c6aa82dc6ceb299d7546e48587f3e6c.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/MCF7_E2_WCE_rep2/HI.2791.003.Index_4.MCF7_WCE_E_rep2.bam/HI.2791.003.Index_4.MCF7_WCE_E_rep2.bam.sorted.bam \
  > alignment/MCF7_E2_WCE_rep2/HI.2791.003.Index_4.MCF7_WCE_E_rep2.bam/HI.2791.003.Index_4.MCF7_WCE_E_rep2.bam.sorted.filtered.bam
samtools_view_filter.HI.2791.003.Index_4.MCF7_WCE_E_rep2.bam.9c6aa82dc6ceb299d7546e48587f3e6c.mugqic.done
)
samtools_view_filter_18_JOB_ID=$(echo "#! /bin/bash
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
echo "$samtools_view_filter_18_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_19_JOB_ID: samtools_view_filter.HI.2791.003.Index_12.MCF7_1_NIPBL_E_rep2.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.2791.003.Index_12.MCF7_1_NIPBL_E_rep2.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_19_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.2791.003.Index_12.MCF7_1_NIPBL_E_rep2.bam.02326a1cf513b1be8da5f87b15956d1a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.2791.003.Index_12.MCF7_1_NIPBL_E_rep2.bam.02326a1cf513b1be8da5f87b15956d1a.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/MCF7_E2_NIPBL_rep2/HI.2791.003.Index_12.MCF7_1_NIPBL_E_rep2.bam/HI.2791.003.Index_12.MCF7_1_NIPBL_E_rep2.bam.sorted.bam \
  > alignment/MCF7_E2_NIPBL_rep2/HI.2791.003.Index_12.MCF7_1_NIPBL_E_rep2.bam/HI.2791.003.Index_12.MCF7_1_NIPBL_E_rep2.bam.sorted.filtered.bam
samtools_view_filter.HI.2791.003.Index_12.MCF7_1_NIPBL_E_rep2.bam.02326a1cf513b1be8da5f87b15956d1a.mugqic.done
)
samtools_view_filter_19_JOB_ID=$(echo "#! /bin/bash
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
echo "$samtools_view_filter_19_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_20_JOB_ID: samtools_view_filter.HI.2791.003.Index_3.MCF7_4_Smc1_E_rep2.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.2791.003.Index_3.MCF7_4_Smc1_E_rep2.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_20_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.2791.003.Index_3.MCF7_4_Smc1_E_rep2.bam.6255177b6fa67002c46bd2638306c72e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.2791.003.Index_3.MCF7_4_Smc1_E_rep2.bam.6255177b6fa67002c46bd2638306c72e.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/MCF7_E2_SMC1A_rep2/HI.2791.003.Index_3.MCF7_4_Smc1_E_rep2.bam/HI.2791.003.Index_3.MCF7_4_Smc1_E_rep2.bam.sorted.bam \
  > alignment/MCF7_E2_SMC1A_rep2/HI.2791.003.Index_3.MCF7_4_Smc1_E_rep2.bam/HI.2791.003.Index_3.MCF7_4_Smc1_E_rep2.bam.sorted.filtered.bam
samtools_view_filter.HI.2791.003.Index_3.MCF7_4_Smc1_E_rep2.bam.6255177b6fa67002c46bd2638306c72e.mugqic.done
)
samtools_view_filter_20_JOB_ID=$(echo "#! /bin/bash
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
echo "$samtools_view_filter_20_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_21_JOB_ID: samtools_view_filter.HI.2791.003.Index_19.MCF7_2_MED1_E_rep2.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.2791.003.Index_19.MCF7_2_MED1_E_rep2.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_21_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.2791.003.Index_19.MCF7_2_MED1_E_rep2.bam.29c7b2ad1388257a304f1b03380d0d3d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.2791.003.Index_19.MCF7_2_MED1_E_rep2.bam.29c7b2ad1388257a304f1b03380d0d3d.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/MCF7_E2_MED1_rep2/HI.2791.003.Index_19.MCF7_2_MED1_E_rep2.bam/HI.2791.003.Index_19.MCF7_2_MED1_E_rep2.bam.sorted.bam \
  > alignment/MCF7_E2_MED1_rep2/HI.2791.003.Index_19.MCF7_2_MED1_E_rep2.bam/HI.2791.003.Index_19.MCF7_2_MED1_E_rep2.bam.sorted.filtered.bam
samtools_view_filter.HI.2791.003.Index_19.MCF7_2_MED1_E_rep2.bam.29c7b2ad1388257a304f1b03380d0d3d.mugqic.done
)
samtools_view_filter_21_JOB_ID=$(echo "#! /bin/bash
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
echo "$samtools_view_filter_21_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_22_JOB_ID: samtools_view_filter.HI.2791.003.Index_1.MCF7_3_Era_E_rep2.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.2791.003.Index_1.MCF7_3_Era_E_rep2.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_22_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.2791.003.Index_1.MCF7_3_Era_E_rep2.bam.368edb2f7b272a7488df02d88b27473a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.2791.003.Index_1.MCF7_3_Era_E_rep2.bam.368edb2f7b272a7488df02d88b27473a.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/MCF7_E2_ERA_rep2/HI.2791.003.Index_1.MCF7_3_Era_E_rep2.bam/HI.2791.003.Index_1.MCF7_3_Era_E_rep2.bam.sorted.bam \
  > alignment/MCF7_E2_ERA_rep2/HI.2791.003.Index_1.MCF7_3_Era_E_rep2.bam/HI.2791.003.Index_1.MCF7_3_Era_E_rep2.bam.sorted.filtered.bam
samtools_view_filter.HI.2791.003.Index_1.MCF7_3_Era_E_rep2.bam.368edb2f7b272a7488df02d88b27473a.mugqic.done
)
samtools_view_filter_22_JOB_ID=$(echo "#! /bin/bash
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
echo "$samtools_view_filter_22_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_23_JOB_ID: samtools_view_filter.SRR1193526.fastq
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.SRR1193526.fastq
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_23_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.SRR1193526.fastq.f1365f26b11cc1e5090fd3b09351baf8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.SRR1193526.fastq.f1365f26b11cc1e5090fd3b09351baf8.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/MCF7_CTRL_BRD4_rep1/SRR1193526.fastq/SRR1193526.fastq.sorted.bam \
  > alignment/MCF7_CTRL_BRD4_rep1/SRR1193526.fastq/SRR1193526.fastq.sorted.filtered.bam
samtools_view_filter.SRR1193526.fastq.f1365f26b11cc1e5090fd3b09351baf8.mugqic.done
)
samtools_view_filter_23_JOB_ID=$(echo "#! /bin/bash
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
echo "$samtools_view_filter_23_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_24_JOB_ID: samtools_view_filter.SRR1193527.fastq
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.SRR1193527.fastq
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_24_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.SRR1193527.fastq.68f5925cf15f5e8d1e7059295d9f3a94.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.SRR1193527.fastq.68f5925cf15f5e8d1e7059295d9f3a94.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/MCF7_CTRL_BRD4_rep2/SRR1193527.fastq/SRR1193527.fastq.sorted.bam \
  > alignment/MCF7_CTRL_BRD4_rep2/SRR1193527.fastq/SRR1193527.fastq.sorted.filtered.bam
samtools_view_filter.SRR1193527.fastq.68f5925cf15f5e8d1e7059295d9f3a94.mugqic.done
)
samtools_view_filter_24_JOB_ID=$(echo "#! /bin/bash
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
echo "$samtools_view_filter_24_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_25_JOB_ID: samtools_view_filter.SRR1193528.fastq
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.SRR1193528.fastq
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_25_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.SRR1193528.fastq.0ba7405af4ecd57c7fbbf58e842e031d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.SRR1193528.fastq.0ba7405af4ecd57c7fbbf58e842e031d.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/MCF7_CTRL_BRD4_rep3/SRR1193528.fastq/SRR1193528.fastq.sorted.bam \
  > alignment/MCF7_CTRL_BRD4_rep3/SRR1193528.fastq/SRR1193528.fastq.sorted.filtered.bam
samtools_view_filter.SRR1193528.fastq.0ba7405af4ecd57c7fbbf58e842e031d.mugqic.done
)
samtools_view_filter_25_JOB_ID=$(echo "#! /bin/bash
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
echo "$samtools_view_filter_25_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_26_JOB_ID: samtools_view_filter.SRR1193529.fastq
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.SRR1193529.fastq
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_26_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.SRR1193529.fastq.c9b05818d6bf8f72564f0dcee4327eec.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.SRR1193529.fastq.c9b05818d6bf8f72564f0dcee4327eec.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/MCF7_E2_BRD4_rep1/SRR1193529.fastq/SRR1193529.fastq.sorted.bam \
  > alignment/MCF7_E2_BRD4_rep1/SRR1193529.fastq/SRR1193529.fastq.sorted.filtered.bam
samtools_view_filter.SRR1193529.fastq.c9b05818d6bf8f72564f0dcee4327eec.mugqic.done
)
samtools_view_filter_26_JOB_ID=$(echo "#! /bin/bash
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
echo "$samtools_view_filter_26_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_27_JOB_ID: samtools_view_filter.SRR1193530.fastq
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.SRR1193530.fastq
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_27_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.SRR1193530.fastq.bb06eb38d1a3d0a5ec99f4f44a56b3f9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.SRR1193530.fastq.bb06eb38d1a3d0a5ec99f4f44a56b3f9.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/MCF7_E2_BRD4_rep2/SRR1193530.fastq/SRR1193530.fastq.sorted.bam \
  > alignment/MCF7_E2_BRD4_rep2/SRR1193530.fastq/SRR1193530.fastq.sorted.filtered.bam
samtools_view_filter.SRR1193530.fastq.bb06eb38d1a3d0a5ec99f4f44a56b3f9.mugqic.done
)
samtools_view_filter_27_JOB_ID=$(echo "#! /bin/bash
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
echo "$samtools_view_filter_27_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_28_JOB_ID: samtools_view_filter.SRR1193531.fastq
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.SRR1193531.fastq
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_28_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.SRR1193531.fastq.15353734fef5f71f777be32522e51060.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.SRR1193531.fastq.15353734fef5f71f777be32522e51060.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/MCF7_E2_BRD4_rep3/SRR1193531.fastq/SRR1193531.fastq.sorted.bam \
  > alignment/MCF7_E2_BRD4_rep3/SRR1193531.fastq/SRR1193531.fastq.sorted.filtered.bam
samtools_view_filter.SRR1193531.fastq.15353734fef5f71f777be32522e51060.mugqic.done
)
samtools_view_filter_28_JOB_ID=$(echo "#! /bin/bash
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
echo "$samtools_view_filter_28_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_29_JOB_ID: samtools_view_filter.SRR1193562.fastq
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.SRR1193562.fastq
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_29_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.SRR1193562.fastq.5f5840048204da76c20412a012dc1505.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.SRR1193562.fastq.5f5840048204da76c20412a012dc1505.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/MCF7_CTRL_WCE/SRR1193562.fastq/SRR1193562.fastq.sorted.bam \
  > alignment/MCF7_CTRL_WCE/SRR1193562.fastq/SRR1193562.fastq.sorted.filtered.bam
samtools_view_filter.SRR1193562.fastq.5f5840048204da76c20412a012dc1505.mugqic.done
)
samtools_view_filter_29_JOB_ID=$(echo "#! /bin/bash
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
echo "$samtools_view_filter_29_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_30_JOB_ID: samtools_view_filter.SRR1193563.fastq
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.SRR1193563.fastq
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_30_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.SRR1193563.fastq.36bfd398e2165ea0fa9b14c8de44a661.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.SRR1193563.fastq.36bfd398e2165ea0fa9b14c8de44a661.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/MCF7_E2_WCE/SRR1193563.fastq/SRR1193563.fastq.sorted.bam \
  > alignment/MCF7_E2_WCE/SRR1193563.fastq/SRR1193563.fastq.sorted.filtered.bam
samtools_view_filter.SRR1193563.fastq.36bfd398e2165ea0fa9b14c8de44a661.mugqic.done
)
samtools_view_filter_30_JOB_ID=$(echo "#! /bin/bash
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
echo "$samtools_view_filter_30_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_31_JOB_ID: samtools_view_filter_report
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter_report
JOB_DEPENDENCIES=$samtools_view_filter_1_JOB_ID:$samtools_view_filter_2_JOB_ID:$samtools_view_filter_3_JOB_ID:$samtools_view_filter_4_JOB_ID:$samtools_view_filter_5_JOB_ID:$samtools_view_filter_6_JOB_ID:$samtools_view_filter_7_JOB_ID:$samtools_view_filter_8_JOB_ID:$samtools_view_filter_9_JOB_ID:$samtools_view_filter_10_JOB_ID:$samtools_view_filter_11_JOB_ID:$samtools_view_filter_12_JOB_ID:$samtools_view_filter_13_JOB_ID:$samtools_view_filter_14_JOB_ID:$samtools_view_filter_15_JOB_ID:$samtools_view_filter_16_JOB_ID:$samtools_view_filter_17_JOB_ID:$samtools_view_filter_18_JOB_ID:$samtools_view_filter_19_JOB_ID:$samtools_view_filter_20_JOB_ID:$samtools_view_filter_21_JOB_ID:$samtools_view_filter_22_JOB_ID:$samtools_view_filter_23_JOB_ID:$samtools_view_filter_24_JOB_ID:$samtools_view_filter_25_JOB_ID:$samtools_view_filter_26_JOB_ID:$samtools_view_filter_27_JOB_ID:$samtools_view_filter_28_JOB_ID:$samtools_view_filter_29_JOB_ID:$samtools_view_filter_30_JOB_ID
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
samtools_view_filter_31_JOB_ID=$(echo "#! /bin/bash
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
echo "$samtools_view_filter_31_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# STEP: picard_merge_sam_files
#-------------------------------------------------------------------------------
STEP=picard_merge_sam_files
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_1_JOB_ID: symlink_readset_sample_bam.MCF7_CTRL_WCE_rep1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.MCF7_CTRL_WCE_rep1
JOB_DEPENDENCIES=$samtools_view_filter_1_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.MCF7_CTRL_WCE_rep1.e331ce51cbb0291cd543e7d2ef26faf5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.MCF7_CTRL_WCE_rep1.e331ce51cbb0291cd543e7d2ef26faf5.mugqic.done'
mkdir -p alignment/MCF7_CTRL_WCE_rep1 && \
ln -s -f HI.1943.008.Index_2.MCF7_CTRL_WCE_ChIP8_GB_rep1.bam/HI.1943.008.Index_2.MCF7_CTRL_WCE_ChIP8_GB_rep1.bam.sorted.filtered.bam alignment/MCF7_CTRL_WCE_rep1/MCF7_CTRL_WCE_rep1.merged.bam
symlink_readset_sample_bam.MCF7_CTRL_WCE_rep1.e331ce51cbb0291cd543e7d2ef26faf5.mugqic.done
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
# JOB: picard_merge_sam_files_2_JOB_ID: symlink_readset_sample_bam.MCF7_CTRL_NIPBL_rep1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.MCF7_CTRL_NIPBL_rep1
JOB_DEPENDENCIES=$samtools_view_filter_2_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.MCF7_CTRL_NIPBL_rep1.9ae40456c62e70157613d63271aa8d1b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.MCF7_CTRL_NIPBL_rep1.9ae40456c62e70157613d63271aa8d1b.mugqic.done'
mkdir -p alignment/MCF7_CTRL_NIPBL_rep1 && \
ln -s -f HI.1997.004.Index_1.MCF7_CTRL_NIPBL_ChIP8_GB_rep1.bam/HI.1997.004.Index_1.MCF7_CTRL_NIPBL_ChIP8_GB_rep1.bam.sorted.filtered.bam alignment/MCF7_CTRL_NIPBL_rep1/MCF7_CTRL_NIPBL_rep1.merged.bam
symlink_readset_sample_bam.MCF7_CTRL_NIPBL_rep1.9ae40456c62e70157613d63271aa8d1b.mugqic.done
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
# JOB: picard_merge_sam_files_3_JOB_ID: symlink_readset_sample_bam.MCF7_CTRL_SMC1A_rep1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.MCF7_CTRL_SMC1A_rep1
JOB_DEPENDENCIES=$samtools_view_filter_3_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.MCF7_CTRL_SMC1A_rep1.5f9f1dc7494f54579ef5589ffbbe4fe3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.MCF7_CTRL_SMC1A_rep1.5f9f1dc7494f54579ef5589ffbbe4fe3.mugqic.done'
mkdir -p alignment/MCF7_CTRL_SMC1A_rep1 && \
ln -s -f HI.1997.004.Index_4.MCF7_CTRL_SMC1_ChIP8_GB_rep1.bam/HI.1997.004.Index_4.MCF7_CTRL_SMC1_ChIP8_GB_rep1.bam.sorted.filtered.bam alignment/MCF7_CTRL_SMC1A_rep1/MCF7_CTRL_SMC1A_rep1.merged.bam
symlink_readset_sample_bam.MCF7_CTRL_SMC1A_rep1.5f9f1dc7494f54579ef5589ffbbe4fe3.mugqic.done
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
# JOB: picard_merge_sam_files_4_JOB_ID: symlink_readset_sample_bam.MCF7_CTRL_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.MCF7_CTRL_MED1_rep1
JOB_DEPENDENCIES=$samtools_view_filter_4_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.MCF7_CTRL_MED1_rep1.fe0f51d58a8c787646ad53ab5fcc6030.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.MCF7_CTRL_MED1_rep1.fe0f51d58a8c787646ad53ab5fcc6030.mugqic.done'
mkdir -p alignment/MCF7_CTRL_MED1_rep1 && \
ln -s -f HI.1943.008.Index_7.MCF7_CTRL_MED1_ChIP8_GB_rep1.bam/HI.1943.008.Index_7.MCF7_CTRL_MED1_ChIP8_GB_rep1.bam.sorted.filtered.bam alignment/MCF7_CTRL_MED1_rep1/MCF7_CTRL_MED1_rep1.merged.bam
symlink_readset_sample_bam.MCF7_CTRL_MED1_rep1.fe0f51d58a8c787646ad53ab5fcc6030.mugqic.done
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
# JOB: picard_merge_sam_files_5_JOB_ID: symlink_readset_sample_bam.MCF7_CTRL_POL2_rep1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.MCF7_CTRL_POL2_rep1
JOB_DEPENDENCIES=$samtools_view_filter_5_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.MCF7_CTRL_POL2_rep1.cb14f2f38313aa614dde9c24c732cd2e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.MCF7_CTRL_POL2_rep1.cb14f2f38313aa614dde9c24c732cd2e.mugqic.done'
mkdir -p alignment/MCF7_CTRL_POL2_rep1 && \
ln -s -f HI.1997.005.Index_10.MCF7_CTRL_POL2_ChIP8_GB_rep1.bam/HI.1997.005.Index_10.MCF7_CTRL_POL2_ChIP8_GB_rep1.bam.sorted.filtered.bam alignment/MCF7_CTRL_POL2_rep1/MCF7_CTRL_POL2_rep1.merged.bam
symlink_readset_sample_bam.MCF7_CTRL_POL2_rep1.cb14f2f38313aa614dde9c24c732cd2e.mugqic.done
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
# JOB: picard_merge_sam_files_6_JOB_ID: symlink_readset_sample_bam.MCF7_CTRL_ERA_rep1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.MCF7_CTRL_ERA_rep1
JOB_DEPENDENCIES=$samtools_view_filter_6_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.MCF7_CTRL_ERA_rep1.a70d42844147932ed3ffbc426fff8011.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.MCF7_CTRL_ERA_rep1.a70d42844147932ed3ffbc426fff8011.mugqic.done'
mkdir -p alignment/MCF7_CTRL_ERA_rep1 && \
ln -s -f HI.1997.005.Index_25.MCF7_CTRL_ERA_ChIP8_GB_rep1.bam/HI.1997.005.Index_25.MCF7_CTRL_ERA_ChIP8_GB_rep1.bam.sorted.filtered.bam alignment/MCF7_CTRL_ERA_rep1/MCF7_CTRL_ERA_rep1.merged.bam
symlink_readset_sample_bam.MCF7_CTRL_ERA_rep1.a70d42844147932ed3ffbc426fff8011.mugqic.done
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
# JOB: picard_merge_sam_files_7_JOB_ID: symlink_readset_sample_bam.MCF7_E2_WCE_rep1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.MCF7_E2_WCE_rep1
JOB_DEPENDENCIES=$samtools_view_filter_7_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.MCF7_E2_WCE_rep1.a8a94ae0de907276cc8a891e27ff0487.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.MCF7_E2_WCE_rep1.a8a94ae0de907276cc8a891e27ff0487.mugqic.done'
mkdir -p alignment/MCF7_E2_WCE_rep1 && \
ln -s -f HI.1943.008.Index_13.MCF7_E2_WCE_ChIP8_GB_rep1.bam/HI.1943.008.Index_13.MCF7_E2_WCE_ChIP8_GB_rep1.bam.sorted.filtered.bam alignment/MCF7_E2_WCE_rep1/MCF7_E2_WCE_rep1.merged.bam
symlink_readset_sample_bam.MCF7_E2_WCE_rep1.a8a94ae0de907276cc8a891e27ff0487.mugqic.done
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
# JOB: picard_merge_sam_files_8_JOB_ID: symlink_readset_sample_bam.MCF7_E2_NIPBL_rep1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.MCF7_E2_NIPBL_rep1
JOB_DEPENDENCIES=$samtools_view_filter_8_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.MCF7_E2_NIPBL_rep1.5648f253e2ed815000a9f97b33d19634.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.MCF7_E2_NIPBL_rep1.5648f253e2ed815000a9f97b33d19634.mugqic.done'
mkdir -p alignment/MCF7_E2_NIPBL_rep1 && \
ln -s -f HI.1997.004.Index_3.MCF7_E2_NIPBL_ChIP8_GB_rep1.bam/HI.1997.004.Index_3.MCF7_E2_NIPBL_ChIP8_GB_rep1.bam.sorted.filtered.bam alignment/MCF7_E2_NIPBL_rep1/MCF7_E2_NIPBL_rep1.merged.bam
symlink_readset_sample_bam.MCF7_E2_NIPBL_rep1.5648f253e2ed815000a9f97b33d19634.mugqic.done
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
# JOB: picard_merge_sam_files_9_JOB_ID: symlink_readset_sample_bam.MCF7_E2_SMC1A_rep1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.MCF7_E2_SMC1A_rep1
JOB_DEPENDENCIES=$samtools_view_filter_9_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.MCF7_E2_SMC1A_rep1.fa746f0a20a89170cc48717d710cf611.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.MCF7_E2_SMC1A_rep1.fa746f0a20a89170cc48717d710cf611.mugqic.done'
mkdir -p alignment/MCF7_E2_SMC1A_rep1 && \
ln -s -f HI.1997.004.Index_5.MCF7_E2_SMC1_ChIP8_GB_rep1.bam/HI.1997.004.Index_5.MCF7_E2_SMC1_ChIP8_GB_rep1.bam.sorted.filtered.bam alignment/MCF7_E2_SMC1A_rep1/MCF7_E2_SMC1A_rep1.merged.bam
symlink_readset_sample_bam.MCF7_E2_SMC1A_rep1.fa746f0a20a89170cc48717d710cf611.mugqic.done
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
# JOB: picard_merge_sam_files_10_JOB_ID: symlink_readset_sample_bam.MCF7_E2_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.MCF7_E2_MED1_rep1
JOB_DEPENDENCIES=$samtools_view_filter_10_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.MCF7_E2_MED1_rep1.4a33e82b1e2b2f2fb858406f31af7102.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.MCF7_E2_MED1_rep1.4a33e82b1e2b2f2fb858406f31af7102.mugqic.done'
mkdir -p alignment/MCF7_E2_MED1_rep1 && \
ln -s -f HI.1943.008.Index_18.MCF7_E2_MED1_ChIP8_GB_rep1.bam/HI.1943.008.Index_18.MCF7_E2_MED1_ChIP8_GB_rep1.bam.sorted.filtered.bam alignment/MCF7_E2_MED1_rep1/MCF7_E2_MED1_rep1.merged.bam
symlink_readset_sample_bam.MCF7_E2_MED1_rep1.4a33e82b1e2b2f2fb858406f31af7102.mugqic.done
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
# JOB: picard_merge_sam_files_11_JOB_ID: symlink_readset_sample_bam.MCF7_E2_POL2_rep1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.MCF7_E2_POL2_rep1
JOB_DEPENDENCIES=$samtools_view_filter_11_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.MCF7_E2_POL2_rep1.56b4c313d4b10773cd958ca6a6872f98.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.MCF7_E2_POL2_rep1.56b4c313d4b10773cd958ca6a6872f98.mugqic.done'
mkdir -p alignment/MCF7_E2_POL2_rep1 && \
ln -s -f HI.1997.005.Index_11.MCF7_E2_POL2_ChIP8_GB_rep1.bam/HI.1997.005.Index_11.MCF7_E2_POL2_ChIP8_GB_rep1.bam.sorted.filtered.bam alignment/MCF7_E2_POL2_rep1/MCF7_E2_POL2_rep1.merged.bam
symlink_readset_sample_bam.MCF7_E2_POL2_rep1.56b4c313d4b10773cd958ca6a6872f98.mugqic.done
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
# JOB: picard_merge_sam_files_12_JOB_ID: symlink_readset_sample_bam.MCF7_E2_ERA_rep1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.MCF7_E2_ERA_rep1
JOB_DEPENDENCIES=$samtools_view_filter_12_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.MCF7_E2_ERA_rep1.1f80ccbef062b19a2edafeee8e25ac79.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.MCF7_E2_ERA_rep1.1f80ccbef062b19a2edafeee8e25ac79.mugqic.done'
mkdir -p alignment/MCF7_E2_ERA_rep1 && \
ln -s -f HI.1997.005.Index_21.MCF7_E2_ERA_ChIP8_GB_rep1.bam/HI.1997.005.Index_21.MCF7_E2_ERA_ChIP8_GB_rep1.bam.sorted.filtered.bam alignment/MCF7_E2_ERA_rep1/MCF7_E2_ERA_rep1.merged.bam
symlink_readset_sample_bam.MCF7_E2_ERA_rep1.1f80ccbef062b19a2edafeee8e25ac79.mugqic.done
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
# JOB: picard_merge_sam_files_13_JOB_ID: symlink_readset_sample_bam.MCF7_CTRL_WCE_rep2
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.MCF7_CTRL_WCE_rep2
JOB_DEPENDENCIES=$samtools_view_filter_13_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.MCF7_CTRL_WCE_rep2.681cc6f763b0b14bced568afc0924d28.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.MCF7_CTRL_WCE_rep2.681cc6f763b0b14bced568afc0924d28.mugqic.done'
mkdir -p alignment/MCF7_CTRL_WCE_rep2 && \
ln -s -f HI.2791.003.Index_5.MCF7_WCE_C_rep2.bam/HI.2791.003.Index_5.MCF7_WCE_C_rep2.bam.sorted.filtered.bam alignment/MCF7_CTRL_WCE_rep2/MCF7_CTRL_WCE_rep2.merged.bam
symlink_readset_sample_bam.MCF7_CTRL_WCE_rep2.681cc6f763b0b14bced568afc0924d28.mugqic.done
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
# JOB: picard_merge_sam_files_14_JOB_ID: symlink_readset_sample_bam.MCF7_CTRL_NIPBL_rep2
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.MCF7_CTRL_NIPBL_rep2
JOB_DEPENDENCIES=$samtools_view_filter_14_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.MCF7_CTRL_NIPBL_rep2.076d1f69d4b1bdd6a7d1490e9a8cb708.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.MCF7_CTRL_NIPBL_rep2.076d1f69d4b1bdd6a7d1490e9a8cb708.mugqic.done'
mkdir -p alignment/MCF7_CTRL_NIPBL_rep2 && \
ln -s -f HI.2791.003.Index_9.MCF7_1_NIPBL_C_rep2.bam/HI.2791.003.Index_9.MCF7_1_NIPBL_C_rep2.bam.sorted.filtered.bam alignment/MCF7_CTRL_NIPBL_rep2/MCF7_CTRL_NIPBL_rep2.merged.bam
symlink_readset_sample_bam.MCF7_CTRL_NIPBL_rep2.076d1f69d4b1bdd6a7d1490e9a8cb708.mugqic.done
)
picard_merge_sam_files_14_JOB_ID=$(echo "#! /bin/bash
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
echo "$picard_merge_sam_files_14_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_15_JOB_ID: symlink_readset_sample_bam.MCF7_CTRL_SMC1A_rep2
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.MCF7_CTRL_SMC1A_rep2
JOB_DEPENDENCIES=$samtools_view_filter_15_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.MCF7_CTRL_SMC1A_rep2.6fb9b76a9365d315eb04aa025097883a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.MCF7_CTRL_SMC1A_rep2.6fb9b76a9365d315eb04aa025097883a.mugqic.done'
mkdir -p alignment/MCF7_CTRL_SMC1A_rep2 && \
ln -s -f HI.2791.003.Index_11.MCF7_4_Smc1_C_rep2.bam/HI.2791.003.Index_11.MCF7_4_Smc1_C_rep2.bam.sorted.filtered.bam alignment/MCF7_CTRL_SMC1A_rep2/MCF7_CTRL_SMC1A_rep2.merged.bam
symlink_readset_sample_bam.MCF7_CTRL_SMC1A_rep2.6fb9b76a9365d315eb04aa025097883a.mugqic.done
)
picard_merge_sam_files_15_JOB_ID=$(echo "#! /bin/bash
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
echo "$picard_merge_sam_files_15_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_16_JOB_ID: symlink_readset_sample_bam.MCF7_CTRL_MED1_rep2
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.MCF7_CTRL_MED1_rep2
JOB_DEPENDENCIES=$samtools_view_filter_16_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.MCF7_CTRL_MED1_rep2.caaf6229d0b00570a1c1f7783adaa929.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.MCF7_CTRL_MED1_rep2.caaf6229d0b00570a1c1f7783adaa929.mugqic.done'
mkdir -p alignment/MCF7_CTRL_MED1_rep2 && \
ln -s -f HI.2791.003.Index_8.MCF7_2_MED1_C_rep2.bam/HI.2791.003.Index_8.MCF7_2_MED1_C_rep2.bam.sorted.filtered.bam alignment/MCF7_CTRL_MED1_rep2/MCF7_CTRL_MED1_rep2.merged.bam
symlink_readset_sample_bam.MCF7_CTRL_MED1_rep2.caaf6229d0b00570a1c1f7783adaa929.mugqic.done
)
picard_merge_sam_files_16_JOB_ID=$(echo "#! /bin/bash
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
echo "$picard_merge_sam_files_16_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_17_JOB_ID: symlink_readset_sample_bam.MCF7_CTRL_ERA_rep2
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.MCF7_CTRL_ERA_rep2
JOB_DEPENDENCIES=$samtools_view_filter_17_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.MCF7_CTRL_ERA_rep2.c1eb5070471f3e35e6cb4724624bdbca.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.MCF7_CTRL_ERA_rep2.c1eb5070471f3e35e6cb4724624bdbca.mugqic.done'
mkdir -p alignment/MCF7_CTRL_ERA_rep2 && \
ln -s -f HI.2791.003.Index_10.MCF7_3_Era_C_rep2.bam/HI.2791.003.Index_10.MCF7_3_Era_C_rep2.bam.sorted.filtered.bam alignment/MCF7_CTRL_ERA_rep2/MCF7_CTRL_ERA_rep2.merged.bam
symlink_readset_sample_bam.MCF7_CTRL_ERA_rep2.c1eb5070471f3e35e6cb4724624bdbca.mugqic.done
)
picard_merge_sam_files_17_JOB_ID=$(echo "#! /bin/bash
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
echo "$picard_merge_sam_files_17_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_18_JOB_ID: symlink_readset_sample_bam.MCF7_E2_WCE_rep2
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.MCF7_E2_WCE_rep2
JOB_DEPENDENCIES=$samtools_view_filter_18_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.MCF7_E2_WCE_rep2.be4e13e0cf6ca8c72be8b477a98da6ab.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.MCF7_E2_WCE_rep2.be4e13e0cf6ca8c72be8b477a98da6ab.mugqic.done'
mkdir -p alignment/MCF7_E2_WCE_rep2 && \
ln -s -f HI.2791.003.Index_4.MCF7_WCE_E_rep2.bam/HI.2791.003.Index_4.MCF7_WCE_E_rep2.bam.sorted.filtered.bam alignment/MCF7_E2_WCE_rep2/MCF7_E2_WCE_rep2.merged.bam
symlink_readset_sample_bam.MCF7_E2_WCE_rep2.be4e13e0cf6ca8c72be8b477a98da6ab.mugqic.done
)
picard_merge_sam_files_18_JOB_ID=$(echo "#! /bin/bash
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
echo "$picard_merge_sam_files_18_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_19_JOB_ID: symlink_readset_sample_bam.MCF7_E2_NIPBL_rep2
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.MCF7_E2_NIPBL_rep2
JOB_DEPENDENCIES=$samtools_view_filter_19_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.MCF7_E2_NIPBL_rep2.208276d4d5195be4abfb76f426458305.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.MCF7_E2_NIPBL_rep2.208276d4d5195be4abfb76f426458305.mugqic.done'
mkdir -p alignment/MCF7_E2_NIPBL_rep2 && \
ln -s -f HI.2791.003.Index_12.MCF7_1_NIPBL_E_rep2.bam/HI.2791.003.Index_12.MCF7_1_NIPBL_E_rep2.bam.sorted.filtered.bam alignment/MCF7_E2_NIPBL_rep2/MCF7_E2_NIPBL_rep2.merged.bam
symlink_readset_sample_bam.MCF7_E2_NIPBL_rep2.208276d4d5195be4abfb76f426458305.mugqic.done
)
picard_merge_sam_files_19_JOB_ID=$(echo "#! /bin/bash
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
echo "$picard_merge_sam_files_19_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_20_JOB_ID: symlink_readset_sample_bam.MCF7_E2_SMC1A_rep2
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.MCF7_E2_SMC1A_rep2
JOB_DEPENDENCIES=$samtools_view_filter_20_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.MCF7_E2_SMC1A_rep2.5641c367562728968531e7239b08a586.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.MCF7_E2_SMC1A_rep2.5641c367562728968531e7239b08a586.mugqic.done'
mkdir -p alignment/MCF7_E2_SMC1A_rep2 && \
ln -s -f HI.2791.003.Index_3.MCF7_4_Smc1_E_rep2.bam/HI.2791.003.Index_3.MCF7_4_Smc1_E_rep2.bam.sorted.filtered.bam alignment/MCF7_E2_SMC1A_rep2/MCF7_E2_SMC1A_rep2.merged.bam
symlink_readset_sample_bam.MCF7_E2_SMC1A_rep2.5641c367562728968531e7239b08a586.mugqic.done
)
picard_merge_sam_files_20_JOB_ID=$(echo "#! /bin/bash
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
echo "$picard_merge_sam_files_20_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_21_JOB_ID: symlink_readset_sample_bam.MCF7_E2_MED1_rep2
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.MCF7_E2_MED1_rep2
JOB_DEPENDENCIES=$samtools_view_filter_21_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.MCF7_E2_MED1_rep2.b53331a11e989e9c335440707dc54209.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.MCF7_E2_MED1_rep2.b53331a11e989e9c335440707dc54209.mugqic.done'
mkdir -p alignment/MCF7_E2_MED1_rep2 && \
ln -s -f HI.2791.003.Index_19.MCF7_2_MED1_E_rep2.bam/HI.2791.003.Index_19.MCF7_2_MED1_E_rep2.bam.sorted.filtered.bam alignment/MCF7_E2_MED1_rep2/MCF7_E2_MED1_rep2.merged.bam
symlink_readset_sample_bam.MCF7_E2_MED1_rep2.b53331a11e989e9c335440707dc54209.mugqic.done
)
picard_merge_sam_files_21_JOB_ID=$(echo "#! /bin/bash
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
echo "$picard_merge_sam_files_21_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_22_JOB_ID: symlink_readset_sample_bam.MCF7_E2_ERA_rep2
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.MCF7_E2_ERA_rep2
JOB_DEPENDENCIES=$samtools_view_filter_22_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.MCF7_E2_ERA_rep2.bad29ac6c727bb5baa1c67c1d4dee333.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.MCF7_E2_ERA_rep2.bad29ac6c727bb5baa1c67c1d4dee333.mugqic.done'
mkdir -p alignment/MCF7_E2_ERA_rep2 && \
ln -s -f HI.2791.003.Index_1.MCF7_3_Era_E_rep2.bam/HI.2791.003.Index_1.MCF7_3_Era_E_rep2.bam.sorted.filtered.bam alignment/MCF7_E2_ERA_rep2/MCF7_E2_ERA_rep2.merged.bam
symlink_readset_sample_bam.MCF7_E2_ERA_rep2.bad29ac6c727bb5baa1c67c1d4dee333.mugqic.done
)
picard_merge_sam_files_22_JOB_ID=$(echo "#! /bin/bash
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
echo "$picard_merge_sam_files_22_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_23_JOB_ID: symlink_readset_sample_bam.MCF7_CTRL_BRD4_rep1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.MCF7_CTRL_BRD4_rep1
JOB_DEPENDENCIES=$samtools_view_filter_23_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.MCF7_CTRL_BRD4_rep1.53368124a104a1c2dba50e624cc7b02b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.MCF7_CTRL_BRD4_rep1.53368124a104a1c2dba50e624cc7b02b.mugqic.done'
mkdir -p alignment/MCF7_CTRL_BRD4_rep1 && \
ln -s -f SRR1193526.fastq/SRR1193526.fastq.sorted.filtered.bam alignment/MCF7_CTRL_BRD4_rep1/MCF7_CTRL_BRD4_rep1.merged.bam
symlink_readset_sample_bam.MCF7_CTRL_BRD4_rep1.53368124a104a1c2dba50e624cc7b02b.mugqic.done
)
picard_merge_sam_files_23_JOB_ID=$(echo "#! /bin/bash
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
echo "$picard_merge_sam_files_23_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_24_JOB_ID: symlink_readset_sample_bam.MCF7_CTRL_BRD4_rep2
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.MCF7_CTRL_BRD4_rep2
JOB_DEPENDENCIES=$samtools_view_filter_24_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.MCF7_CTRL_BRD4_rep2.cfed391fa534fabd1c25afa6e8afbcba.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.MCF7_CTRL_BRD4_rep2.cfed391fa534fabd1c25afa6e8afbcba.mugqic.done'
mkdir -p alignment/MCF7_CTRL_BRD4_rep2 && \
ln -s -f SRR1193527.fastq/SRR1193527.fastq.sorted.filtered.bam alignment/MCF7_CTRL_BRD4_rep2/MCF7_CTRL_BRD4_rep2.merged.bam
symlink_readset_sample_bam.MCF7_CTRL_BRD4_rep2.cfed391fa534fabd1c25afa6e8afbcba.mugqic.done
)
picard_merge_sam_files_24_JOB_ID=$(echo "#! /bin/bash
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
echo "$picard_merge_sam_files_24_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_25_JOB_ID: symlink_readset_sample_bam.MCF7_CTRL_BRD4_rep3
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.MCF7_CTRL_BRD4_rep3
JOB_DEPENDENCIES=$samtools_view_filter_25_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.MCF7_CTRL_BRD4_rep3.17f9965fe3428502dcc019dba5efb011.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.MCF7_CTRL_BRD4_rep3.17f9965fe3428502dcc019dba5efb011.mugqic.done'
mkdir -p alignment/MCF7_CTRL_BRD4_rep3 && \
ln -s -f SRR1193528.fastq/SRR1193528.fastq.sorted.filtered.bam alignment/MCF7_CTRL_BRD4_rep3/MCF7_CTRL_BRD4_rep3.merged.bam
symlink_readset_sample_bam.MCF7_CTRL_BRD4_rep3.17f9965fe3428502dcc019dba5efb011.mugqic.done
)
picard_merge_sam_files_25_JOB_ID=$(echo "#! /bin/bash
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
echo "$picard_merge_sam_files_25_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_26_JOB_ID: symlink_readset_sample_bam.MCF7_E2_BRD4_rep1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.MCF7_E2_BRD4_rep1
JOB_DEPENDENCIES=$samtools_view_filter_26_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.MCF7_E2_BRD4_rep1.c9e1d92fdbab560cb138211853b6c3e6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.MCF7_E2_BRD4_rep1.c9e1d92fdbab560cb138211853b6c3e6.mugqic.done'
mkdir -p alignment/MCF7_E2_BRD4_rep1 && \
ln -s -f SRR1193529.fastq/SRR1193529.fastq.sorted.filtered.bam alignment/MCF7_E2_BRD4_rep1/MCF7_E2_BRD4_rep1.merged.bam
symlink_readset_sample_bam.MCF7_E2_BRD4_rep1.c9e1d92fdbab560cb138211853b6c3e6.mugqic.done
)
picard_merge_sam_files_26_JOB_ID=$(echo "#! /bin/bash
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
echo "$picard_merge_sam_files_26_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_27_JOB_ID: symlink_readset_sample_bam.MCF7_E2_BRD4_rep2
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.MCF7_E2_BRD4_rep2
JOB_DEPENDENCIES=$samtools_view_filter_27_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.MCF7_E2_BRD4_rep2.f683cde92abc2756dba0dd67eafbd2a1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.MCF7_E2_BRD4_rep2.f683cde92abc2756dba0dd67eafbd2a1.mugqic.done'
mkdir -p alignment/MCF7_E2_BRD4_rep2 && \
ln -s -f SRR1193530.fastq/SRR1193530.fastq.sorted.filtered.bam alignment/MCF7_E2_BRD4_rep2/MCF7_E2_BRD4_rep2.merged.bam
symlink_readset_sample_bam.MCF7_E2_BRD4_rep2.f683cde92abc2756dba0dd67eafbd2a1.mugqic.done
)
picard_merge_sam_files_27_JOB_ID=$(echo "#! /bin/bash
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
echo "$picard_merge_sam_files_27_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_28_JOB_ID: symlink_readset_sample_bam.MCF7_E2_BRD4_rep3
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.MCF7_E2_BRD4_rep3
JOB_DEPENDENCIES=$samtools_view_filter_28_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.MCF7_E2_BRD4_rep3.db8466648ce5a297004f6b4151a13a57.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.MCF7_E2_BRD4_rep3.db8466648ce5a297004f6b4151a13a57.mugqic.done'
mkdir -p alignment/MCF7_E2_BRD4_rep3 && \
ln -s -f SRR1193531.fastq/SRR1193531.fastq.sorted.filtered.bam alignment/MCF7_E2_BRD4_rep3/MCF7_E2_BRD4_rep3.merged.bam
symlink_readset_sample_bam.MCF7_E2_BRD4_rep3.db8466648ce5a297004f6b4151a13a57.mugqic.done
)
picard_merge_sam_files_28_JOB_ID=$(echo "#! /bin/bash
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
echo "$picard_merge_sam_files_28_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_29_JOB_ID: symlink_readset_sample_bam.MCF7_CTRL_WCE
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.MCF7_CTRL_WCE
JOB_DEPENDENCIES=$samtools_view_filter_29_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.MCF7_CTRL_WCE.e1019496d4ce1b6d577f562dfb96ce5f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.MCF7_CTRL_WCE.e1019496d4ce1b6d577f562dfb96ce5f.mugqic.done'
mkdir -p alignment/MCF7_CTRL_WCE && \
ln -s -f SRR1193562.fastq/SRR1193562.fastq.sorted.filtered.bam alignment/MCF7_CTRL_WCE/MCF7_CTRL_WCE.merged.bam
symlink_readset_sample_bam.MCF7_CTRL_WCE.e1019496d4ce1b6d577f562dfb96ce5f.mugqic.done
)
picard_merge_sam_files_29_JOB_ID=$(echo "#! /bin/bash
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
echo "$picard_merge_sam_files_29_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_30_JOB_ID: symlink_readset_sample_bam.MCF7_E2_WCE
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.MCF7_E2_WCE
JOB_DEPENDENCIES=$samtools_view_filter_30_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.MCF7_E2_WCE.c69bc69fecd26a4ea1f87b3c6f04e939.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.MCF7_E2_WCE.c69bc69fecd26a4ea1f87b3c6f04e939.mugqic.done'
mkdir -p alignment/MCF7_E2_WCE && \
ln -s -f SRR1193563.fastq/SRR1193563.fastq.sorted.filtered.bam alignment/MCF7_E2_WCE/MCF7_E2_WCE.merged.bam
symlink_readset_sample_bam.MCF7_E2_WCE.c69bc69fecd26a4ea1f87b3c6f04e939.mugqic.done
)
picard_merge_sam_files_30_JOB_ID=$(echo "#! /bin/bash
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
echo "$picard_merge_sam_files_30_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# STEP: picard_mark_duplicates
#-------------------------------------------------------------------------------
STEP=picard_mark_duplicates
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_1_JOB_ID: picard_mark_duplicates.MCF7_CTRL_WCE_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.MCF7_CTRL_WCE_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_1_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.MCF7_CTRL_WCE_rep1.acc294bf79b9e24d064cfbbc320563d8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.MCF7_CTRL_WCE_rep1.acc294bf79b9e24d064cfbbc320563d8.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/MCF7_CTRL_WCE_rep1/MCF7_CTRL_WCE_rep1.merged.bam \
 OUTPUT=alignment/MCF7_CTRL_WCE_rep1/MCF7_CTRL_WCE_rep1.sorted.dup.bam \
 METRICS_FILE=alignment/MCF7_CTRL_WCE_rep1/MCF7_CTRL_WCE_rep1.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.MCF7_CTRL_WCE_rep1.acc294bf79b9e24d064cfbbc320563d8.mugqic.done
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
# JOB: picard_mark_duplicates_2_JOB_ID: picard_mark_duplicates.MCF7_CTRL_NIPBL_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.MCF7_CTRL_NIPBL_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_2_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.MCF7_CTRL_NIPBL_rep1.d56bdb7e799cda8f60d96034547db7b4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.MCF7_CTRL_NIPBL_rep1.d56bdb7e799cda8f60d96034547db7b4.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/MCF7_CTRL_NIPBL_rep1/MCF7_CTRL_NIPBL_rep1.merged.bam \
 OUTPUT=alignment/MCF7_CTRL_NIPBL_rep1/MCF7_CTRL_NIPBL_rep1.sorted.dup.bam \
 METRICS_FILE=alignment/MCF7_CTRL_NIPBL_rep1/MCF7_CTRL_NIPBL_rep1.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.MCF7_CTRL_NIPBL_rep1.d56bdb7e799cda8f60d96034547db7b4.mugqic.done
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
# JOB: picard_mark_duplicates_3_JOB_ID: picard_mark_duplicates.MCF7_CTRL_SMC1A_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.MCF7_CTRL_SMC1A_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_3_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.MCF7_CTRL_SMC1A_rep1.a2a2acce881a6a702d0654a29f1f4d52.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.MCF7_CTRL_SMC1A_rep1.a2a2acce881a6a702d0654a29f1f4d52.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/MCF7_CTRL_SMC1A_rep1/MCF7_CTRL_SMC1A_rep1.merged.bam \
 OUTPUT=alignment/MCF7_CTRL_SMC1A_rep1/MCF7_CTRL_SMC1A_rep1.sorted.dup.bam \
 METRICS_FILE=alignment/MCF7_CTRL_SMC1A_rep1/MCF7_CTRL_SMC1A_rep1.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.MCF7_CTRL_SMC1A_rep1.a2a2acce881a6a702d0654a29f1f4d52.mugqic.done
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
# JOB: picard_mark_duplicates_4_JOB_ID: picard_mark_duplicates.MCF7_CTRL_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.MCF7_CTRL_MED1_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_4_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.MCF7_CTRL_MED1_rep1.446c2294be7fcbab0869e10be281e8e5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.MCF7_CTRL_MED1_rep1.446c2294be7fcbab0869e10be281e8e5.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/MCF7_CTRL_MED1_rep1/MCF7_CTRL_MED1_rep1.merged.bam \
 OUTPUT=alignment/MCF7_CTRL_MED1_rep1/MCF7_CTRL_MED1_rep1.sorted.dup.bam \
 METRICS_FILE=alignment/MCF7_CTRL_MED1_rep1/MCF7_CTRL_MED1_rep1.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.MCF7_CTRL_MED1_rep1.446c2294be7fcbab0869e10be281e8e5.mugqic.done
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
# JOB: picard_mark_duplicates_5_JOB_ID: picard_mark_duplicates.MCF7_CTRL_POL2_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.MCF7_CTRL_POL2_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_5_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.MCF7_CTRL_POL2_rep1.3ec25aa0695d442e64bd5163a5daaba4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.MCF7_CTRL_POL2_rep1.3ec25aa0695d442e64bd5163a5daaba4.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/MCF7_CTRL_POL2_rep1/MCF7_CTRL_POL2_rep1.merged.bam \
 OUTPUT=alignment/MCF7_CTRL_POL2_rep1/MCF7_CTRL_POL2_rep1.sorted.dup.bam \
 METRICS_FILE=alignment/MCF7_CTRL_POL2_rep1/MCF7_CTRL_POL2_rep1.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.MCF7_CTRL_POL2_rep1.3ec25aa0695d442e64bd5163a5daaba4.mugqic.done
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
# JOB: picard_mark_duplicates_6_JOB_ID: picard_mark_duplicates.MCF7_CTRL_ERA_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.MCF7_CTRL_ERA_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_6_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.MCF7_CTRL_ERA_rep1.16bcd077c3fd71a76a3813b77c80f0e6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.MCF7_CTRL_ERA_rep1.16bcd077c3fd71a76a3813b77c80f0e6.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/MCF7_CTRL_ERA_rep1/MCF7_CTRL_ERA_rep1.merged.bam \
 OUTPUT=alignment/MCF7_CTRL_ERA_rep1/MCF7_CTRL_ERA_rep1.sorted.dup.bam \
 METRICS_FILE=alignment/MCF7_CTRL_ERA_rep1/MCF7_CTRL_ERA_rep1.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.MCF7_CTRL_ERA_rep1.16bcd077c3fd71a76a3813b77c80f0e6.mugqic.done
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
# JOB: picard_mark_duplicates_7_JOB_ID: picard_mark_duplicates.MCF7_E2_WCE_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.MCF7_E2_WCE_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_7_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.MCF7_E2_WCE_rep1.6f610008fb3e974b868a4be2cfb0c8fb.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.MCF7_E2_WCE_rep1.6f610008fb3e974b868a4be2cfb0c8fb.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/MCF7_E2_WCE_rep1/MCF7_E2_WCE_rep1.merged.bam \
 OUTPUT=alignment/MCF7_E2_WCE_rep1/MCF7_E2_WCE_rep1.sorted.dup.bam \
 METRICS_FILE=alignment/MCF7_E2_WCE_rep1/MCF7_E2_WCE_rep1.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.MCF7_E2_WCE_rep1.6f610008fb3e974b868a4be2cfb0c8fb.mugqic.done
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
# JOB: picard_mark_duplicates_8_JOB_ID: picard_mark_duplicates.MCF7_E2_NIPBL_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.MCF7_E2_NIPBL_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_8_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.MCF7_E2_NIPBL_rep1.53658c83e0f066135fd370d3cf64d53d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.MCF7_E2_NIPBL_rep1.53658c83e0f066135fd370d3cf64d53d.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/MCF7_E2_NIPBL_rep1/MCF7_E2_NIPBL_rep1.merged.bam \
 OUTPUT=alignment/MCF7_E2_NIPBL_rep1/MCF7_E2_NIPBL_rep1.sorted.dup.bam \
 METRICS_FILE=alignment/MCF7_E2_NIPBL_rep1/MCF7_E2_NIPBL_rep1.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.MCF7_E2_NIPBL_rep1.53658c83e0f066135fd370d3cf64d53d.mugqic.done
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
# JOB: picard_mark_duplicates_9_JOB_ID: picard_mark_duplicates.MCF7_E2_SMC1A_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.MCF7_E2_SMC1A_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_9_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.MCF7_E2_SMC1A_rep1.07f6b7e522ae92ae1d08560625a31d9b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.MCF7_E2_SMC1A_rep1.07f6b7e522ae92ae1d08560625a31d9b.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/MCF7_E2_SMC1A_rep1/MCF7_E2_SMC1A_rep1.merged.bam \
 OUTPUT=alignment/MCF7_E2_SMC1A_rep1/MCF7_E2_SMC1A_rep1.sorted.dup.bam \
 METRICS_FILE=alignment/MCF7_E2_SMC1A_rep1/MCF7_E2_SMC1A_rep1.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.MCF7_E2_SMC1A_rep1.07f6b7e522ae92ae1d08560625a31d9b.mugqic.done
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
# JOB: picard_mark_duplicates_10_JOB_ID: picard_mark_duplicates.MCF7_E2_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.MCF7_E2_MED1_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_10_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.MCF7_E2_MED1_rep1.d2dc42ce983f005377253fdc85b3ef31.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.MCF7_E2_MED1_rep1.d2dc42ce983f005377253fdc85b3ef31.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/MCF7_E2_MED1_rep1/MCF7_E2_MED1_rep1.merged.bam \
 OUTPUT=alignment/MCF7_E2_MED1_rep1/MCF7_E2_MED1_rep1.sorted.dup.bam \
 METRICS_FILE=alignment/MCF7_E2_MED1_rep1/MCF7_E2_MED1_rep1.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.MCF7_E2_MED1_rep1.d2dc42ce983f005377253fdc85b3ef31.mugqic.done
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
# JOB: picard_mark_duplicates_11_JOB_ID: picard_mark_duplicates.MCF7_E2_POL2_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.MCF7_E2_POL2_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_11_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.MCF7_E2_POL2_rep1.d5404d463fc7b3a5a2f5347d3ce32ae6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.MCF7_E2_POL2_rep1.d5404d463fc7b3a5a2f5347d3ce32ae6.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/MCF7_E2_POL2_rep1/MCF7_E2_POL2_rep1.merged.bam \
 OUTPUT=alignment/MCF7_E2_POL2_rep1/MCF7_E2_POL2_rep1.sorted.dup.bam \
 METRICS_FILE=alignment/MCF7_E2_POL2_rep1/MCF7_E2_POL2_rep1.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.MCF7_E2_POL2_rep1.d5404d463fc7b3a5a2f5347d3ce32ae6.mugqic.done
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
# JOB: picard_mark_duplicates_12_JOB_ID: picard_mark_duplicates.MCF7_E2_ERA_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.MCF7_E2_ERA_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_12_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.MCF7_E2_ERA_rep1.7e5d92100d2fdda51295473fdc364f8d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.MCF7_E2_ERA_rep1.7e5d92100d2fdda51295473fdc364f8d.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/MCF7_E2_ERA_rep1/MCF7_E2_ERA_rep1.merged.bam \
 OUTPUT=alignment/MCF7_E2_ERA_rep1/MCF7_E2_ERA_rep1.sorted.dup.bam \
 METRICS_FILE=alignment/MCF7_E2_ERA_rep1/MCF7_E2_ERA_rep1.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.MCF7_E2_ERA_rep1.7e5d92100d2fdda51295473fdc364f8d.mugqic.done
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
# JOB: picard_mark_duplicates_13_JOB_ID: picard_mark_duplicates.MCF7_CTRL_WCE_rep2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.MCF7_CTRL_WCE_rep2
JOB_DEPENDENCIES=$picard_merge_sam_files_13_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.MCF7_CTRL_WCE_rep2.3f72cc3863cddaef5db091eb1cd440cc.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.MCF7_CTRL_WCE_rep2.3f72cc3863cddaef5db091eb1cd440cc.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/MCF7_CTRL_WCE_rep2/MCF7_CTRL_WCE_rep2.merged.bam \
 OUTPUT=alignment/MCF7_CTRL_WCE_rep2/MCF7_CTRL_WCE_rep2.sorted.dup.bam \
 METRICS_FILE=alignment/MCF7_CTRL_WCE_rep2/MCF7_CTRL_WCE_rep2.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.MCF7_CTRL_WCE_rep2.3f72cc3863cddaef5db091eb1cd440cc.mugqic.done
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
# JOB: picard_mark_duplicates_14_JOB_ID: picard_mark_duplicates.MCF7_CTRL_NIPBL_rep2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.MCF7_CTRL_NIPBL_rep2
JOB_DEPENDENCIES=$picard_merge_sam_files_14_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.MCF7_CTRL_NIPBL_rep2.09ceb0df5e9572a3f50ad0e297bfcf8d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.MCF7_CTRL_NIPBL_rep2.09ceb0df5e9572a3f50ad0e297bfcf8d.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/MCF7_CTRL_NIPBL_rep2/MCF7_CTRL_NIPBL_rep2.merged.bam \
 OUTPUT=alignment/MCF7_CTRL_NIPBL_rep2/MCF7_CTRL_NIPBL_rep2.sorted.dup.bam \
 METRICS_FILE=alignment/MCF7_CTRL_NIPBL_rep2/MCF7_CTRL_NIPBL_rep2.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.MCF7_CTRL_NIPBL_rep2.09ceb0df5e9572a3f50ad0e297bfcf8d.mugqic.done
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=48:00:0 --mem=8G -N 1 -n 2 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$picard_mark_duplicates_14_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_15_JOB_ID: picard_mark_duplicates.MCF7_CTRL_SMC1A_rep2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.MCF7_CTRL_SMC1A_rep2
JOB_DEPENDENCIES=$picard_merge_sam_files_15_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.MCF7_CTRL_SMC1A_rep2.9224305084651aea63f8690b99d40da0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.MCF7_CTRL_SMC1A_rep2.9224305084651aea63f8690b99d40da0.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/MCF7_CTRL_SMC1A_rep2/MCF7_CTRL_SMC1A_rep2.merged.bam \
 OUTPUT=alignment/MCF7_CTRL_SMC1A_rep2/MCF7_CTRL_SMC1A_rep2.sorted.dup.bam \
 METRICS_FILE=alignment/MCF7_CTRL_SMC1A_rep2/MCF7_CTRL_SMC1A_rep2.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.MCF7_CTRL_SMC1A_rep2.9224305084651aea63f8690b99d40da0.mugqic.done
)
picard_mark_duplicates_15_JOB_ID=$(echo "#! /bin/bash
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
echo "$picard_mark_duplicates_15_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_16_JOB_ID: picard_mark_duplicates.MCF7_CTRL_MED1_rep2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.MCF7_CTRL_MED1_rep2
JOB_DEPENDENCIES=$picard_merge_sam_files_16_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.MCF7_CTRL_MED1_rep2.00976321ee83e850d2aad7549fd412e9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.MCF7_CTRL_MED1_rep2.00976321ee83e850d2aad7549fd412e9.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/MCF7_CTRL_MED1_rep2/MCF7_CTRL_MED1_rep2.merged.bam \
 OUTPUT=alignment/MCF7_CTRL_MED1_rep2/MCF7_CTRL_MED1_rep2.sorted.dup.bam \
 METRICS_FILE=alignment/MCF7_CTRL_MED1_rep2/MCF7_CTRL_MED1_rep2.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.MCF7_CTRL_MED1_rep2.00976321ee83e850d2aad7549fd412e9.mugqic.done
)
picard_mark_duplicates_16_JOB_ID=$(echo "#! /bin/bash
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
echo "$picard_mark_duplicates_16_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_17_JOB_ID: picard_mark_duplicates.MCF7_CTRL_ERA_rep2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.MCF7_CTRL_ERA_rep2
JOB_DEPENDENCIES=$picard_merge_sam_files_17_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.MCF7_CTRL_ERA_rep2.5a285055f35e5e43df37e6ebe4c646d8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.MCF7_CTRL_ERA_rep2.5a285055f35e5e43df37e6ebe4c646d8.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/MCF7_CTRL_ERA_rep2/MCF7_CTRL_ERA_rep2.merged.bam \
 OUTPUT=alignment/MCF7_CTRL_ERA_rep2/MCF7_CTRL_ERA_rep2.sorted.dup.bam \
 METRICS_FILE=alignment/MCF7_CTRL_ERA_rep2/MCF7_CTRL_ERA_rep2.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.MCF7_CTRL_ERA_rep2.5a285055f35e5e43df37e6ebe4c646d8.mugqic.done
)
picard_mark_duplicates_17_JOB_ID=$(echo "#! /bin/bash
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
echo "$picard_mark_duplicates_17_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_18_JOB_ID: picard_mark_duplicates.MCF7_E2_WCE_rep2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.MCF7_E2_WCE_rep2
JOB_DEPENDENCIES=$picard_merge_sam_files_18_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.MCF7_E2_WCE_rep2.621f266e4e084b242bd7bf91a31a3cb3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.MCF7_E2_WCE_rep2.621f266e4e084b242bd7bf91a31a3cb3.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/MCF7_E2_WCE_rep2/MCF7_E2_WCE_rep2.merged.bam \
 OUTPUT=alignment/MCF7_E2_WCE_rep2/MCF7_E2_WCE_rep2.sorted.dup.bam \
 METRICS_FILE=alignment/MCF7_E2_WCE_rep2/MCF7_E2_WCE_rep2.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.MCF7_E2_WCE_rep2.621f266e4e084b242bd7bf91a31a3cb3.mugqic.done
)
picard_mark_duplicates_18_JOB_ID=$(echo "#! /bin/bash
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
echo "$picard_mark_duplicates_18_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_19_JOB_ID: picard_mark_duplicates.MCF7_E2_NIPBL_rep2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.MCF7_E2_NIPBL_rep2
JOB_DEPENDENCIES=$picard_merge_sam_files_19_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.MCF7_E2_NIPBL_rep2.a26e13db541b82e923fe751f9dbd87a7.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.MCF7_E2_NIPBL_rep2.a26e13db541b82e923fe751f9dbd87a7.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/MCF7_E2_NIPBL_rep2/MCF7_E2_NIPBL_rep2.merged.bam \
 OUTPUT=alignment/MCF7_E2_NIPBL_rep2/MCF7_E2_NIPBL_rep2.sorted.dup.bam \
 METRICS_FILE=alignment/MCF7_E2_NIPBL_rep2/MCF7_E2_NIPBL_rep2.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.MCF7_E2_NIPBL_rep2.a26e13db541b82e923fe751f9dbd87a7.mugqic.done
)
picard_mark_duplicates_19_JOB_ID=$(echo "#! /bin/bash
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
echo "$picard_mark_duplicates_19_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_20_JOB_ID: picard_mark_duplicates.MCF7_E2_SMC1A_rep2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.MCF7_E2_SMC1A_rep2
JOB_DEPENDENCIES=$picard_merge_sam_files_20_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.MCF7_E2_SMC1A_rep2.6fef30d5cc1950e84dc9a29ecc65b99f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.MCF7_E2_SMC1A_rep2.6fef30d5cc1950e84dc9a29ecc65b99f.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/MCF7_E2_SMC1A_rep2/MCF7_E2_SMC1A_rep2.merged.bam \
 OUTPUT=alignment/MCF7_E2_SMC1A_rep2/MCF7_E2_SMC1A_rep2.sorted.dup.bam \
 METRICS_FILE=alignment/MCF7_E2_SMC1A_rep2/MCF7_E2_SMC1A_rep2.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.MCF7_E2_SMC1A_rep2.6fef30d5cc1950e84dc9a29ecc65b99f.mugqic.done
)
picard_mark_duplicates_20_JOB_ID=$(echo "#! /bin/bash
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
echo "$picard_mark_duplicates_20_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_21_JOB_ID: picard_mark_duplicates.MCF7_E2_MED1_rep2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.MCF7_E2_MED1_rep2
JOB_DEPENDENCIES=$picard_merge_sam_files_21_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.MCF7_E2_MED1_rep2.db5ab007d5d963d44f687793685010ff.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.MCF7_E2_MED1_rep2.db5ab007d5d963d44f687793685010ff.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/MCF7_E2_MED1_rep2/MCF7_E2_MED1_rep2.merged.bam \
 OUTPUT=alignment/MCF7_E2_MED1_rep2/MCF7_E2_MED1_rep2.sorted.dup.bam \
 METRICS_FILE=alignment/MCF7_E2_MED1_rep2/MCF7_E2_MED1_rep2.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.MCF7_E2_MED1_rep2.db5ab007d5d963d44f687793685010ff.mugqic.done
)
picard_mark_duplicates_21_JOB_ID=$(echo "#! /bin/bash
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
echo "$picard_mark_duplicates_21_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_22_JOB_ID: picard_mark_duplicates.MCF7_E2_ERA_rep2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.MCF7_E2_ERA_rep2
JOB_DEPENDENCIES=$picard_merge_sam_files_22_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.MCF7_E2_ERA_rep2.b30f587a8e1841a3d4e760e46e412b53.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.MCF7_E2_ERA_rep2.b30f587a8e1841a3d4e760e46e412b53.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/MCF7_E2_ERA_rep2/MCF7_E2_ERA_rep2.merged.bam \
 OUTPUT=alignment/MCF7_E2_ERA_rep2/MCF7_E2_ERA_rep2.sorted.dup.bam \
 METRICS_FILE=alignment/MCF7_E2_ERA_rep2/MCF7_E2_ERA_rep2.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.MCF7_E2_ERA_rep2.b30f587a8e1841a3d4e760e46e412b53.mugqic.done
)
picard_mark_duplicates_22_JOB_ID=$(echo "#! /bin/bash
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
echo "$picard_mark_duplicates_22_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_23_JOB_ID: picard_mark_duplicates.MCF7_CTRL_BRD4_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.MCF7_CTRL_BRD4_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_23_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.MCF7_CTRL_BRD4_rep1.91c5e2325ac6a2d18685685fc3e870ea.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.MCF7_CTRL_BRD4_rep1.91c5e2325ac6a2d18685685fc3e870ea.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/MCF7_CTRL_BRD4_rep1/MCF7_CTRL_BRD4_rep1.merged.bam \
 OUTPUT=alignment/MCF7_CTRL_BRD4_rep1/MCF7_CTRL_BRD4_rep1.sorted.dup.bam \
 METRICS_FILE=alignment/MCF7_CTRL_BRD4_rep1/MCF7_CTRL_BRD4_rep1.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.MCF7_CTRL_BRD4_rep1.91c5e2325ac6a2d18685685fc3e870ea.mugqic.done
)
picard_mark_duplicates_23_JOB_ID=$(echo "#! /bin/bash
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
echo "$picard_mark_duplicates_23_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_24_JOB_ID: picard_mark_duplicates.MCF7_CTRL_BRD4_rep2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.MCF7_CTRL_BRD4_rep2
JOB_DEPENDENCIES=$picard_merge_sam_files_24_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.MCF7_CTRL_BRD4_rep2.bb64301eb25015b8b6344a9787be9f3d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.MCF7_CTRL_BRD4_rep2.bb64301eb25015b8b6344a9787be9f3d.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/MCF7_CTRL_BRD4_rep2/MCF7_CTRL_BRD4_rep2.merged.bam \
 OUTPUT=alignment/MCF7_CTRL_BRD4_rep2/MCF7_CTRL_BRD4_rep2.sorted.dup.bam \
 METRICS_FILE=alignment/MCF7_CTRL_BRD4_rep2/MCF7_CTRL_BRD4_rep2.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.MCF7_CTRL_BRD4_rep2.bb64301eb25015b8b6344a9787be9f3d.mugqic.done
)
picard_mark_duplicates_24_JOB_ID=$(echo "#! /bin/bash
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
echo "$picard_mark_duplicates_24_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_25_JOB_ID: picard_mark_duplicates.MCF7_CTRL_BRD4_rep3
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.MCF7_CTRL_BRD4_rep3
JOB_DEPENDENCIES=$picard_merge_sam_files_25_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.MCF7_CTRL_BRD4_rep3.88119bb5684311a7c95326b9b9f3849b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.MCF7_CTRL_BRD4_rep3.88119bb5684311a7c95326b9b9f3849b.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/MCF7_CTRL_BRD4_rep3/MCF7_CTRL_BRD4_rep3.merged.bam \
 OUTPUT=alignment/MCF7_CTRL_BRD4_rep3/MCF7_CTRL_BRD4_rep3.sorted.dup.bam \
 METRICS_FILE=alignment/MCF7_CTRL_BRD4_rep3/MCF7_CTRL_BRD4_rep3.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.MCF7_CTRL_BRD4_rep3.88119bb5684311a7c95326b9b9f3849b.mugqic.done
)
picard_mark_duplicates_25_JOB_ID=$(echo "#! /bin/bash
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
echo "$picard_mark_duplicates_25_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_26_JOB_ID: picard_mark_duplicates.MCF7_E2_BRD4_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.MCF7_E2_BRD4_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_26_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.MCF7_E2_BRD4_rep1.42a4ad5aae7093ba449969ce9be8875e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.MCF7_E2_BRD4_rep1.42a4ad5aae7093ba449969ce9be8875e.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/MCF7_E2_BRD4_rep1/MCF7_E2_BRD4_rep1.merged.bam \
 OUTPUT=alignment/MCF7_E2_BRD4_rep1/MCF7_E2_BRD4_rep1.sorted.dup.bam \
 METRICS_FILE=alignment/MCF7_E2_BRD4_rep1/MCF7_E2_BRD4_rep1.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.MCF7_E2_BRD4_rep1.42a4ad5aae7093ba449969ce9be8875e.mugqic.done
)
picard_mark_duplicates_26_JOB_ID=$(echo "#! /bin/bash
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
echo "$picard_mark_duplicates_26_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_27_JOB_ID: picard_mark_duplicates.MCF7_E2_BRD4_rep2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.MCF7_E2_BRD4_rep2
JOB_DEPENDENCIES=$picard_merge_sam_files_27_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.MCF7_E2_BRD4_rep2.3c7346809be84c3c4687e00c9f657926.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.MCF7_E2_BRD4_rep2.3c7346809be84c3c4687e00c9f657926.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/MCF7_E2_BRD4_rep2/MCF7_E2_BRD4_rep2.merged.bam \
 OUTPUT=alignment/MCF7_E2_BRD4_rep2/MCF7_E2_BRD4_rep2.sorted.dup.bam \
 METRICS_FILE=alignment/MCF7_E2_BRD4_rep2/MCF7_E2_BRD4_rep2.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.MCF7_E2_BRD4_rep2.3c7346809be84c3c4687e00c9f657926.mugqic.done
)
picard_mark_duplicates_27_JOB_ID=$(echo "#! /bin/bash
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
echo "$picard_mark_duplicates_27_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_28_JOB_ID: picard_mark_duplicates.MCF7_E2_BRD4_rep3
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.MCF7_E2_BRD4_rep3
JOB_DEPENDENCIES=$picard_merge_sam_files_28_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.MCF7_E2_BRD4_rep3.55329576556a449eee60381e2e586c24.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.MCF7_E2_BRD4_rep3.55329576556a449eee60381e2e586c24.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/MCF7_E2_BRD4_rep3/MCF7_E2_BRD4_rep3.merged.bam \
 OUTPUT=alignment/MCF7_E2_BRD4_rep3/MCF7_E2_BRD4_rep3.sorted.dup.bam \
 METRICS_FILE=alignment/MCF7_E2_BRD4_rep3/MCF7_E2_BRD4_rep3.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.MCF7_E2_BRD4_rep3.55329576556a449eee60381e2e586c24.mugqic.done
)
picard_mark_duplicates_28_JOB_ID=$(echo "#! /bin/bash
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
echo "$picard_mark_duplicates_28_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_29_JOB_ID: picard_mark_duplicates.MCF7_CTRL_WCE
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.MCF7_CTRL_WCE
JOB_DEPENDENCIES=$picard_merge_sam_files_29_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.MCF7_CTRL_WCE.0daa994586ae717a919599e431d10819.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.MCF7_CTRL_WCE.0daa994586ae717a919599e431d10819.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/MCF7_CTRL_WCE/MCF7_CTRL_WCE.merged.bam \
 OUTPUT=alignment/MCF7_CTRL_WCE/MCF7_CTRL_WCE.sorted.dup.bam \
 METRICS_FILE=alignment/MCF7_CTRL_WCE/MCF7_CTRL_WCE.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.MCF7_CTRL_WCE.0daa994586ae717a919599e431d10819.mugqic.done
)
picard_mark_duplicates_29_JOB_ID=$(echo "#! /bin/bash
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
echo "$picard_mark_duplicates_29_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_30_JOB_ID: picard_mark_duplicates.MCF7_E2_WCE
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.MCF7_E2_WCE
JOB_DEPENDENCIES=$picard_merge_sam_files_30_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.MCF7_E2_WCE.790af8ca7ba72217d5cf59dca9e9d710.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.MCF7_E2_WCE.790af8ca7ba72217d5cf59dca9e9d710.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/MCF7_E2_WCE/MCF7_E2_WCE.merged.bam \
 OUTPUT=alignment/MCF7_E2_WCE/MCF7_E2_WCE.sorted.dup.bam \
 METRICS_FILE=alignment/MCF7_E2_WCE/MCF7_E2_WCE.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.MCF7_E2_WCE.790af8ca7ba72217d5cf59dca9e9d710.mugqic.done
)
picard_mark_duplicates_30_JOB_ID=$(echo "#! /bin/bash
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
echo "$picard_mark_duplicates_30_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_31_JOB_ID: picard_mark_duplicates_report
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates_report
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_2_JOB_ID:$picard_mark_duplicates_3_JOB_ID:$picard_mark_duplicates_4_JOB_ID:$picard_mark_duplicates_5_JOB_ID:$picard_mark_duplicates_6_JOB_ID:$picard_mark_duplicates_7_JOB_ID:$picard_mark_duplicates_8_JOB_ID:$picard_mark_duplicates_9_JOB_ID:$picard_mark_duplicates_10_JOB_ID:$picard_mark_duplicates_11_JOB_ID:$picard_mark_duplicates_12_JOB_ID:$picard_mark_duplicates_13_JOB_ID:$picard_mark_duplicates_14_JOB_ID:$picard_mark_duplicates_15_JOB_ID:$picard_mark_duplicates_16_JOB_ID:$picard_mark_duplicates_17_JOB_ID:$picard_mark_duplicates_18_JOB_ID:$picard_mark_duplicates_19_JOB_ID:$picard_mark_duplicates_20_JOB_ID:$picard_mark_duplicates_21_JOB_ID:$picard_mark_duplicates_22_JOB_ID:$picard_mark_duplicates_23_JOB_ID:$picard_mark_duplicates_24_JOB_ID:$picard_mark_duplicates_25_JOB_ID:$picard_mark_duplicates_26_JOB_ID:$picard_mark_duplicates_27_JOB_ID:$picard_mark_duplicates_28_JOB_ID:$picard_mark_duplicates_29_JOB_ID:$picard_mark_duplicates_30_JOB_ID
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
picard_mark_duplicates_31_JOB_ID=$(echo "#! /bin/bash
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
echo "$picard_mark_duplicates_31_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

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
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_2_JOB_ID:$picard_mark_duplicates_3_JOB_ID:$picard_mark_duplicates_4_JOB_ID:$picard_mark_duplicates_5_JOB_ID:$picard_mark_duplicates_6_JOB_ID:$picard_mark_duplicates_7_JOB_ID:$picard_mark_duplicates_8_JOB_ID:$picard_mark_duplicates_9_JOB_ID:$picard_mark_duplicates_10_JOB_ID:$picard_mark_duplicates_11_JOB_ID:$picard_mark_duplicates_12_JOB_ID:$picard_mark_duplicates_13_JOB_ID:$picard_mark_duplicates_14_JOB_ID:$picard_mark_duplicates_15_JOB_ID:$picard_mark_duplicates_16_JOB_ID:$picard_mark_duplicates_17_JOB_ID:$picard_mark_duplicates_18_JOB_ID:$picard_mark_duplicates_19_JOB_ID:$picard_mark_duplicates_20_JOB_ID:$picard_mark_duplicates_21_JOB_ID:$picard_mark_duplicates_22_JOB_ID:$picard_mark_duplicates_23_JOB_ID:$picard_mark_duplicates_24_JOB_ID:$picard_mark_duplicates_25_JOB_ID:$picard_mark_duplicates_26_JOB_ID:$picard_mark_duplicates_27_JOB_ID:$picard_mark_duplicates_28_JOB_ID:$picard_mark_duplicates_29_JOB_ID:$picard_mark_duplicates_30_JOB_ID
JOB_DONE=job_output/metrics/metrics.flagstat.94382b66febeb1852cb4c3bfb364fccb.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'metrics.flagstat.94382b66febeb1852cb4c3bfb364fccb.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools flagstat \
  alignment/MCF7_CTRL_WCE_rep1/MCF7_CTRL_WCE_rep1.sorted.dup.bam \
  > alignment/MCF7_CTRL_WCE_rep1/MCF7_CTRL_WCE_rep1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/MCF7_CTRL_NIPBL_rep1/MCF7_CTRL_NIPBL_rep1.sorted.dup.bam \
  > alignment/MCF7_CTRL_NIPBL_rep1/MCF7_CTRL_NIPBL_rep1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/MCF7_CTRL_SMC1A_rep1/MCF7_CTRL_SMC1A_rep1.sorted.dup.bam \
  > alignment/MCF7_CTRL_SMC1A_rep1/MCF7_CTRL_SMC1A_rep1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/MCF7_CTRL_MED1_rep1/MCF7_CTRL_MED1_rep1.sorted.dup.bam \
  > alignment/MCF7_CTRL_MED1_rep1/MCF7_CTRL_MED1_rep1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/MCF7_CTRL_POL2_rep1/MCF7_CTRL_POL2_rep1.sorted.dup.bam \
  > alignment/MCF7_CTRL_POL2_rep1/MCF7_CTRL_POL2_rep1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/MCF7_CTRL_ERA_rep1/MCF7_CTRL_ERA_rep1.sorted.dup.bam \
  > alignment/MCF7_CTRL_ERA_rep1/MCF7_CTRL_ERA_rep1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/MCF7_E2_WCE_rep1/MCF7_E2_WCE_rep1.sorted.dup.bam \
  > alignment/MCF7_E2_WCE_rep1/MCF7_E2_WCE_rep1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/MCF7_E2_NIPBL_rep1/MCF7_E2_NIPBL_rep1.sorted.dup.bam \
  > alignment/MCF7_E2_NIPBL_rep1/MCF7_E2_NIPBL_rep1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/MCF7_E2_SMC1A_rep1/MCF7_E2_SMC1A_rep1.sorted.dup.bam \
  > alignment/MCF7_E2_SMC1A_rep1/MCF7_E2_SMC1A_rep1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/MCF7_E2_MED1_rep1/MCF7_E2_MED1_rep1.sorted.dup.bam \
  > alignment/MCF7_E2_MED1_rep1/MCF7_E2_MED1_rep1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/MCF7_E2_POL2_rep1/MCF7_E2_POL2_rep1.sorted.dup.bam \
  > alignment/MCF7_E2_POL2_rep1/MCF7_E2_POL2_rep1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/MCF7_E2_ERA_rep1/MCF7_E2_ERA_rep1.sorted.dup.bam \
  > alignment/MCF7_E2_ERA_rep1/MCF7_E2_ERA_rep1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/MCF7_CTRL_WCE_rep2/MCF7_CTRL_WCE_rep2.sorted.dup.bam \
  > alignment/MCF7_CTRL_WCE_rep2/MCF7_CTRL_WCE_rep2.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/MCF7_CTRL_NIPBL_rep2/MCF7_CTRL_NIPBL_rep2.sorted.dup.bam \
  > alignment/MCF7_CTRL_NIPBL_rep2/MCF7_CTRL_NIPBL_rep2.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/MCF7_CTRL_SMC1A_rep2/MCF7_CTRL_SMC1A_rep2.sorted.dup.bam \
  > alignment/MCF7_CTRL_SMC1A_rep2/MCF7_CTRL_SMC1A_rep2.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/MCF7_CTRL_MED1_rep2/MCF7_CTRL_MED1_rep2.sorted.dup.bam \
  > alignment/MCF7_CTRL_MED1_rep2/MCF7_CTRL_MED1_rep2.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/MCF7_CTRL_ERA_rep2/MCF7_CTRL_ERA_rep2.sorted.dup.bam \
  > alignment/MCF7_CTRL_ERA_rep2/MCF7_CTRL_ERA_rep2.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/MCF7_E2_WCE_rep2/MCF7_E2_WCE_rep2.sorted.dup.bam \
  > alignment/MCF7_E2_WCE_rep2/MCF7_E2_WCE_rep2.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/MCF7_E2_NIPBL_rep2/MCF7_E2_NIPBL_rep2.sorted.dup.bam \
  > alignment/MCF7_E2_NIPBL_rep2/MCF7_E2_NIPBL_rep2.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/MCF7_E2_SMC1A_rep2/MCF7_E2_SMC1A_rep2.sorted.dup.bam \
  > alignment/MCF7_E2_SMC1A_rep2/MCF7_E2_SMC1A_rep2.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/MCF7_E2_MED1_rep2/MCF7_E2_MED1_rep2.sorted.dup.bam \
  > alignment/MCF7_E2_MED1_rep2/MCF7_E2_MED1_rep2.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/MCF7_E2_ERA_rep2/MCF7_E2_ERA_rep2.sorted.dup.bam \
  > alignment/MCF7_E2_ERA_rep2/MCF7_E2_ERA_rep2.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/MCF7_CTRL_BRD4_rep1/MCF7_CTRL_BRD4_rep1.sorted.dup.bam \
  > alignment/MCF7_CTRL_BRD4_rep1/MCF7_CTRL_BRD4_rep1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/MCF7_CTRL_BRD4_rep2/MCF7_CTRL_BRD4_rep2.sorted.dup.bam \
  > alignment/MCF7_CTRL_BRD4_rep2/MCF7_CTRL_BRD4_rep2.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/MCF7_CTRL_BRD4_rep3/MCF7_CTRL_BRD4_rep3.sorted.dup.bam \
  > alignment/MCF7_CTRL_BRD4_rep3/MCF7_CTRL_BRD4_rep3.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/MCF7_E2_BRD4_rep1/MCF7_E2_BRD4_rep1.sorted.dup.bam \
  > alignment/MCF7_E2_BRD4_rep1/MCF7_E2_BRD4_rep1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/MCF7_E2_BRD4_rep2/MCF7_E2_BRD4_rep2.sorted.dup.bam \
  > alignment/MCF7_E2_BRD4_rep2/MCF7_E2_BRD4_rep2.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/MCF7_E2_BRD4_rep3/MCF7_E2_BRD4_rep3.sorted.dup.bam \
  > alignment/MCF7_E2_BRD4_rep3/MCF7_E2_BRD4_rep3.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/MCF7_CTRL_WCE/MCF7_CTRL_WCE.sorted.dup.bam \
  > alignment/MCF7_CTRL_WCE/MCF7_CTRL_WCE.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/MCF7_E2_WCE/MCF7_E2_WCE.sorted.dup.bam \
  > alignment/MCF7_E2_WCE/MCF7_E2_WCE.sorted.dup.bam.flagstat
metrics.flagstat.94382b66febeb1852cb4c3bfb364fccb.mugqic.done
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
JOB_DONE=job_output/metrics/metrics_report.48a649426add707b86e5bccd32986028.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'metrics_report.48a649426add707b86e5bccd32986028.mugqic.done'
module load mugqic/pandoc/1.15.2 && \
for sample in MCF7_CTRL_WCE_rep1 MCF7_CTRL_NIPBL_rep1 MCF7_CTRL_SMC1A_rep1 MCF7_CTRL_MED1_rep1 MCF7_CTRL_POL2_rep1 MCF7_CTRL_ERA_rep1 MCF7_E2_WCE_rep1 MCF7_E2_NIPBL_rep1 MCF7_E2_SMC1A_rep1 MCF7_E2_MED1_rep1 MCF7_E2_POL2_rep1 MCF7_E2_ERA_rep1 MCF7_CTRL_WCE_rep2 MCF7_CTRL_NIPBL_rep2 MCF7_CTRL_SMC1A_rep2 MCF7_CTRL_MED1_rep2 MCF7_CTRL_ERA_rep2 MCF7_E2_WCE_rep2 MCF7_E2_NIPBL_rep2 MCF7_E2_SMC1A_rep2 MCF7_E2_MED1_rep2 MCF7_E2_ERA_rep2 MCF7_CTRL_BRD4_rep1 MCF7_CTRL_BRD4_rep2 MCF7_CTRL_BRD4_rep3 MCF7_E2_BRD4_rep1 MCF7_E2_BRD4_rep2 MCF7_E2_BRD4_rep3 MCF7_CTRL_WCE MCF7_E2_WCE
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

metrics_report.48a649426add707b86e5bccd32986028.mugqic.done
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
# JOB: homer_make_tag_directory_1_JOB_ID: homer_make_tag_directory.MCF7_CTRL_WCE_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.MCF7_CTRL_WCE_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.MCF7_CTRL_WCE_rep1.a5ddb2044067bca25d62a887efaced3a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.MCF7_CTRL_WCE_rep1.a5ddb2044067bca25d62a887efaced3a.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/MCF7_CTRL_WCE_rep1 \
            alignment/MCF7_CTRL_WCE_rep1/MCF7_CTRL_WCE_rep1.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.MCF7_CTRL_WCE_rep1.a5ddb2044067bca25d62a887efaced3a.mugqic.done
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
# JOB: homer_make_tag_directory_2_JOB_ID: homer_make_tag_directory.MCF7_CTRL_NIPBL_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.MCF7_CTRL_NIPBL_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_2_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.MCF7_CTRL_NIPBL_rep1.60b414baf903453df100488d7044854c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.MCF7_CTRL_NIPBL_rep1.60b414baf903453df100488d7044854c.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/MCF7_CTRL_NIPBL_rep1 \
            alignment/MCF7_CTRL_NIPBL_rep1/MCF7_CTRL_NIPBL_rep1.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.MCF7_CTRL_NIPBL_rep1.60b414baf903453df100488d7044854c.mugqic.done
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
# JOB: homer_make_tag_directory_3_JOB_ID: homer_make_tag_directory.MCF7_CTRL_SMC1A_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.MCF7_CTRL_SMC1A_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_3_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.MCF7_CTRL_SMC1A_rep1.9b89f441bf8ea073462cefd50561a034.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.MCF7_CTRL_SMC1A_rep1.9b89f441bf8ea073462cefd50561a034.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/MCF7_CTRL_SMC1A_rep1 \
            alignment/MCF7_CTRL_SMC1A_rep1/MCF7_CTRL_SMC1A_rep1.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.MCF7_CTRL_SMC1A_rep1.9b89f441bf8ea073462cefd50561a034.mugqic.done
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
# JOB: homer_make_tag_directory_4_JOB_ID: homer_make_tag_directory.MCF7_CTRL_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.MCF7_CTRL_MED1_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_4_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.MCF7_CTRL_MED1_rep1.0bc2e83fe1f0aaabaee4863791498040.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.MCF7_CTRL_MED1_rep1.0bc2e83fe1f0aaabaee4863791498040.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/MCF7_CTRL_MED1_rep1 \
            alignment/MCF7_CTRL_MED1_rep1/MCF7_CTRL_MED1_rep1.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.MCF7_CTRL_MED1_rep1.0bc2e83fe1f0aaabaee4863791498040.mugqic.done
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
# JOB: homer_make_tag_directory_5_JOB_ID: homer_make_tag_directory.MCF7_CTRL_POL2_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.MCF7_CTRL_POL2_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_5_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.MCF7_CTRL_POL2_rep1.34f9891787fac61f67a9d0bc30336e10.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.MCF7_CTRL_POL2_rep1.34f9891787fac61f67a9d0bc30336e10.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/MCF7_CTRL_POL2_rep1 \
            alignment/MCF7_CTRL_POL2_rep1/MCF7_CTRL_POL2_rep1.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.MCF7_CTRL_POL2_rep1.34f9891787fac61f67a9d0bc30336e10.mugqic.done
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
# JOB: homer_make_tag_directory_6_JOB_ID: homer_make_tag_directory.MCF7_CTRL_ERA_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.MCF7_CTRL_ERA_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_6_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.MCF7_CTRL_ERA_rep1.9a3961a9ff702c70a84e2fdeabc36ff2.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.MCF7_CTRL_ERA_rep1.9a3961a9ff702c70a84e2fdeabc36ff2.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/MCF7_CTRL_ERA_rep1 \
            alignment/MCF7_CTRL_ERA_rep1/MCF7_CTRL_ERA_rep1.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.MCF7_CTRL_ERA_rep1.9a3961a9ff702c70a84e2fdeabc36ff2.mugqic.done
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
# JOB: homer_make_tag_directory_7_JOB_ID: homer_make_tag_directory.MCF7_E2_WCE_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.MCF7_E2_WCE_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_7_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.MCF7_E2_WCE_rep1.b34f13567c281d4ee7dda2aafa692bf4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.MCF7_E2_WCE_rep1.b34f13567c281d4ee7dda2aafa692bf4.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/MCF7_E2_WCE_rep1 \
            alignment/MCF7_E2_WCE_rep1/MCF7_E2_WCE_rep1.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.MCF7_E2_WCE_rep1.b34f13567c281d4ee7dda2aafa692bf4.mugqic.done
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
# JOB: homer_make_tag_directory_8_JOB_ID: homer_make_tag_directory.MCF7_E2_NIPBL_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.MCF7_E2_NIPBL_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_8_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.MCF7_E2_NIPBL_rep1.23a52cc8464eb1b7f9355d1b9bf93d32.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.MCF7_E2_NIPBL_rep1.23a52cc8464eb1b7f9355d1b9bf93d32.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/MCF7_E2_NIPBL_rep1 \
            alignment/MCF7_E2_NIPBL_rep1/MCF7_E2_NIPBL_rep1.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.MCF7_E2_NIPBL_rep1.23a52cc8464eb1b7f9355d1b9bf93d32.mugqic.done
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
# JOB: homer_make_tag_directory_9_JOB_ID: homer_make_tag_directory.MCF7_E2_SMC1A_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.MCF7_E2_SMC1A_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_9_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.MCF7_E2_SMC1A_rep1.1ebac64a1fe4eaf67f5f06c565cb0455.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.MCF7_E2_SMC1A_rep1.1ebac64a1fe4eaf67f5f06c565cb0455.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/MCF7_E2_SMC1A_rep1 \
            alignment/MCF7_E2_SMC1A_rep1/MCF7_E2_SMC1A_rep1.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.MCF7_E2_SMC1A_rep1.1ebac64a1fe4eaf67f5f06c565cb0455.mugqic.done
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
# JOB: homer_make_tag_directory_10_JOB_ID: homer_make_tag_directory.MCF7_E2_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.MCF7_E2_MED1_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_10_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.MCF7_E2_MED1_rep1.129fa4702057491b3f3dce35bedc86ca.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.MCF7_E2_MED1_rep1.129fa4702057491b3f3dce35bedc86ca.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/MCF7_E2_MED1_rep1 \
            alignment/MCF7_E2_MED1_rep1/MCF7_E2_MED1_rep1.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.MCF7_E2_MED1_rep1.129fa4702057491b3f3dce35bedc86ca.mugqic.done
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
# JOB: homer_make_tag_directory_11_JOB_ID: homer_make_tag_directory.MCF7_E2_POL2_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.MCF7_E2_POL2_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_11_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.MCF7_E2_POL2_rep1.975a0ae0fb3014a7586096b738a30d51.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.MCF7_E2_POL2_rep1.975a0ae0fb3014a7586096b738a30d51.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/MCF7_E2_POL2_rep1 \
            alignment/MCF7_E2_POL2_rep1/MCF7_E2_POL2_rep1.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.MCF7_E2_POL2_rep1.975a0ae0fb3014a7586096b738a30d51.mugqic.done
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
# JOB: homer_make_tag_directory_12_JOB_ID: homer_make_tag_directory.MCF7_E2_ERA_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.MCF7_E2_ERA_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_12_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.MCF7_E2_ERA_rep1.01a49d8145003302f1f41ae26b62b97a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.MCF7_E2_ERA_rep1.01a49d8145003302f1f41ae26b62b97a.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/MCF7_E2_ERA_rep1 \
            alignment/MCF7_E2_ERA_rep1/MCF7_E2_ERA_rep1.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.MCF7_E2_ERA_rep1.01a49d8145003302f1f41ae26b62b97a.mugqic.done
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
# JOB: homer_make_tag_directory_13_JOB_ID: homer_make_tag_directory.MCF7_CTRL_WCE_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.MCF7_CTRL_WCE_rep2
JOB_DEPENDENCIES=$picard_mark_duplicates_13_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.MCF7_CTRL_WCE_rep2.86ba84053aa60541d76cedcbefccdbee.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.MCF7_CTRL_WCE_rep2.86ba84053aa60541d76cedcbefccdbee.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/MCF7_CTRL_WCE_rep2 \
            alignment/MCF7_CTRL_WCE_rep2/MCF7_CTRL_WCE_rep2.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.MCF7_CTRL_WCE_rep2.86ba84053aa60541d76cedcbefccdbee.mugqic.done
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
# JOB: homer_make_tag_directory_14_JOB_ID: homer_make_tag_directory.MCF7_CTRL_NIPBL_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.MCF7_CTRL_NIPBL_rep2
JOB_DEPENDENCIES=$picard_mark_duplicates_14_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.MCF7_CTRL_NIPBL_rep2.077b2a92c57336d1498f1142996ee876.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.MCF7_CTRL_NIPBL_rep2.077b2a92c57336d1498f1142996ee876.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/MCF7_CTRL_NIPBL_rep2 \
            alignment/MCF7_CTRL_NIPBL_rep2/MCF7_CTRL_NIPBL_rep2.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.MCF7_CTRL_NIPBL_rep2.077b2a92c57336d1498f1142996ee876.mugqic.done
)
homer_make_tag_directory_14_JOB_ID=$(echo "#! /bin/bash
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
echo "$homer_make_tag_directory_14_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_15_JOB_ID: homer_make_tag_directory.MCF7_CTRL_SMC1A_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.MCF7_CTRL_SMC1A_rep2
JOB_DEPENDENCIES=$picard_mark_duplicates_15_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.MCF7_CTRL_SMC1A_rep2.22dc7aab3e9cba814ef57b982b892299.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.MCF7_CTRL_SMC1A_rep2.22dc7aab3e9cba814ef57b982b892299.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/MCF7_CTRL_SMC1A_rep2 \
            alignment/MCF7_CTRL_SMC1A_rep2/MCF7_CTRL_SMC1A_rep2.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.MCF7_CTRL_SMC1A_rep2.22dc7aab3e9cba814ef57b982b892299.mugqic.done
)
homer_make_tag_directory_15_JOB_ID=$(echo "#! /bin/bash
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
echo "$homer_make_tag_directory_15_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_16_JOB_ID: homer_make_tag_directory.MCF7_CTRL_MED1_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.MCF7_CTRL_MED1_rep2
JOB_DEPENDENCIES=$picard_mark_duplicates_16_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.MCF7_CTRL_MED1_rep2.d93573147c40439bf574bd63e2026b8e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.MCF7_CTRL_MED1_rep2.d93573147c40439bf574bd63e2026b8e.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/MCF7_CTRL_MED1_rep2 \
            alignment/MCF7_CTRL_MED1_rep2/MCF7_CTRL_MED1_rep2.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.MCF7_CTRL_MED1_rep2.d93573147c40439bf574bd63e2026b8e.mugqic.done
)
homer_make_tag_directory_16_JOB_ID=$(echo "#! /bin/bash
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
echo "$homer_make_tag_directory_16_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_17_JOB_ID: homer_make_tag_directory.MCF7_CTRL_ERA_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.MCF7_CTRL_ERA_rep2
JOB_DEPENDENCIES=$picard_mark_duplicates_17_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.MCF7_CTRL_ERA_rep2.8c790fc17debf055416c85d9f0aaa351.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.MCF7_CTRL_ERA_rep2.8c790fc17debf055416c85d9f0aaa351.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/MCF7_CTRL_ERA_rep2 \
            alignment/MCF7_CTRL_ERA_rep2/MCF7_CTRL_ERA_rep2.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.MCF7_CTRL_ERA_rep2.8c790fc17debf055416c85d9f0aaa351.mugqic.done
)
homer_make_tag_directory_17_JOB_ID=$(echo "#! /bin/bash
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
echo "$homer_make_tag_directory_17_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_18_JOB_ID: homer_make_tag_directory.MCF7_E2_WCE_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.MCF7_E2_WCE_rep2
JOB_DEPENDENCIES=$picard_mark_duplicates_18_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.MCF7_E2_WCE_rep2.f52b5403245d1e2da9e8fa0aa6eed149.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.MCF7_E2_WCE_rep2.f52b5403245d1e2da9e8fa0aa6eed149.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/MCF7_E2_WCE_rep2 \
            alignment/MCF7_E2_WCE_rep2/MCF7_E2_WCE_rep2.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.MCF7_E2_WCE_rep2.f52b5403245d1e2da9e8fa0aa6eed149.mugqic.done
)
homer_make_tag_directory_18_JOB_ID=$(echo "#! /bin/bash
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
echo "$homer_make_tag_directory_18_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_19_JOB_ID: homer_make_tag_directory.MCF7_E2_NIPBL_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.MCF7_E2_NIPBL_rep2
JOB_DEPENDENCIES=$picard_mark_duplicates_19_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.MCF7_E2_NIPBL_rep2.fc7230f673dbc150003ba4fafb15b49e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.MCF7_E2_NIPBL_rep2.fc7230f673dbc150003ba4fafb15b49e.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/MCF7_E2_NIPBL_rep2 \
            alignment/MCF7_E2_NIPBL_rep2/MCF7_E2_NIPBL_rep2.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.MCF7_E2_NIPBL_rep2.fc7230f673dbc150003ba4fafb15b49e.mugqic.done
)
homer_make_tag_directory_19_JOB_ID=$(echo "#! /bin/bash
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
echo "$homer_make_tag_directory_19_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_20_JOB_ID: homer_make_tag_directory.MCF7_E2_SMC1A_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.MCF7_E2_SMC1A_rep2
JOB_DEPENDENCIES=$picard_mark_duplicates_20_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.MCF7_E2_SMC1A_rep2.7ea8a9856206486f17e73217a089fcbb.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.MCF7_E2_SMC1A_rep2.7ea8a9856206486f17e73217a089fcbb.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/MCF7_E2_SMC1A_rep2 \
            alignment/MCF7_E2_SMC1A_rep2/MCF7_E2_SMC1A_rep2.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.MCF7_E2_SMC1A_rep2.7ea8a9856206486f17e73217a089fcbb.mugqic.done
)
homer_make_tag_directory_20_JOB_ID=$(echo "#! /bin/bash
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
echo "$homer_make_tag_directory_20_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_21_JOB_ID: homer_make_tag_directory.MCF7_E2_MED1_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.MCF7_E2_MED1_rep2
JOB_DEPENDENCIES=$picard_mark_duplicates_21_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.MCF7_E2_MED1_rep2.37793a83a79a0e062e06898c0e5a1e7e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.MCF7_E2_MED1_rep2.37793a83a79a0e062e06898c0e5a1e7e.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/MCF7_E2_MED1_rep2 \
            alignment/MCF7_E2_MED1_rep2/MCF7_E2_MED1_rep2.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.MCF7_E2_MED1_rep2.37793a83a79a0e062e06898c0e5a1e7e.mugqic.done
)
homer_make_tag_directory_21_JOB_ID=$(echo "#! /bin/bash
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
echo "$homer_make_tag_directory_21_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_22_JOB_ID: homer_make_tag_directory.MCF7_E2_ERA_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.MCF7_E2_ERA_rep2
JOB_DEPENDENCIES=$picard_mark_duplicates_22_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.MCF7_E2_ERA_rep2.f9ef9cab65a222abcef5a5e56e3661c1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.MCF7_E2_ERA_rep2.f9ef9cab65a222abcef5a5e56e3661c1.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/MCF7_E2_ERA_rep2 \
            alignment/MCF7_E2_ERA_rep2/MCF7_E2_ERA_rep2.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.MCF7_E2_ERA_rep2.f9ef9cab65a222abcef5a5e56e3661c1.mugqic.done
)
homer_make_tag_directory_22_JOB_ID=$(echo "#! /bin/bash
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
echo "$homer_make_tag_directory_22_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_23_JOB_ID: homer_make_tag_directory.MCF7_CTRL_BRD4_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.MCF7_CTRL_BRD4_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_23_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.MCF7_CTRL_BRD4_rep1.8f3687ab5c92ec3b9a003f5d786b59ad.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.MCF7_CTRL_BRD4_rep1.8f3687ab5c92ec3b9a003f5d786b59ad.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/MCF7_CTRL_BRD4_rep1 \
            alignment/MCF7_CTRL_BRD4_rep1/MCF7_CTRL_BRD4_rep1.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.MCF7_CTRL_BRD4_rep1.8f3687ab5c92ec3b9a003f5d786b59ad.mugqic.done
)
homer_make_tag_directory_23_JOB_ID=$(echo "#! /bin/bash
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
echo "$homer_make_tag_directory_23_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_24_JOB_ID: homer_make_tag_directory.MCF7_CTRL_BRD4_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.MCF7_CTRL_BRD4_rep2
JOB_DEPENDENCIES=$picard_mark_duplicates_24_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.MCF7_CTRL_BRD4_rep2.2a4ac264e7855b5fd95436a56f385482.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.MCF7_CTRL_BRD4_rep2.2a4ac264e7855b5fd95436a56f385482.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/MCF7_CTRL_BRD4_rep2 \
            alignment/MCF7_CTRL_BRD4_rep2/MCF7_CTRL_BRD4_rep2.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.MCF7_CTRL_BRD4_rep2.2a4ac264e7855b5fd95436a56f385482.mugqic.done
)
homer_make_tag_directory_24_JOB_ID=$(echo "#! /bin/bash
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
echo "$homer_make_tag_directory_24_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_25_JOB_ID: homer_make_tag_directory.MCF7_CTRL_BRD4_rep3
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.MCF7_CTRL_BRD4_rep3
JOB_DEPENDENCIES=$picard_mark_duplicates_25_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.MCF7_CTRL_BRD4_rep3.320f8cf86a814484ddd1fcce4587c383.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.MCF7_CTRL_BRD4_rep3.320f8cf86a814484ddd1fcce4587c383.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/MCF7_CTRL_BRD4_rep3 \
            alignment/MCF7_CTRL_BRD4_rep3/MCF7_CTRL_BRD4_rep3.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.MCF7_CTRL_BRD4_rep3.320f8cf86a814484ddd1fcce4587c383.mugqic.done
)
homer_make_tag_directory_25_JOB_ID=$(echo "#! /bin/bash
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
echo "$homer_make_tag_directory_25_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_26_JOB_ID: homer_make_tag_directory.MCF7_E2_BRD4_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.MCF7_E2_BRD4_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_26_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.MCF7_E2_BRD4_rep1.3569051fa0102fe7660e577bb872ded1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.MCF7_E2_BRD4_rep1.3569051fa0102fe7660e577bb872ded1.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/MCF7_E2_BRD4_rep1 \
            alignment/MCF7_E2_BRD4_rep1/MCF7_E2_BRD4_rep1.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.MCF7_E2_BRD4_rep1.3569051fa0102fe7660e577bb872ded1.mugqic.done
)
homer_make_tag_directory_26_JOB_ID=$(echo "#! /bin/bash
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
echo "$homer_make_tag_directory_26_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_27_JOB_ID: homer_make_tag_directory.MCF7_E2_BRD4_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.MCF7_E2_BRD4_rep2
JOB_DEPENDENCIES=$picard_mark_duplicates_27_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.MCF7_E2_BRD4_rep2.c20f90443fec73637287e837f4f02ef2.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.MCF7_E2_BRD4_rep2.c20f90443fec73637287e837f4f02ef2.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/MCF7_E2_BRD4_rep2 \
            alignment/MCF7_E2_BRD4_rep2/MCF7_E2_BRD4_rep2.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.MCF7_E2_BRD4_rep2.c20f90443fec73637287e837f4f02ef2.mugqic.done
)
homer_make_tag_directory_27_JOB_ID=$(echo "#! /bin/bash
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
echo "$homer_make_tag_directory_27_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_28_JOB_ID: homer_make_tag_directory.MCF7_E2_BRD4_rep3
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.MCF7_E2_BRD4_rep3
JOB_DEPENDENCIES=$picard_mark_duplicates_28_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.MCF7_E2_BRD4_rep3.633869d583ab3a61d24ab10c8f21662f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.MCF7_E2_BRD4_rep3.633869d583ab3a61d24ab10c8f21662f.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/MCF7_E2_BRD4_rep3 \
            alignment/MCF7_E2_BRD4_rep3/MCF7_E2_BRD4_rep3.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.MCF7_E2_BRD4_rep3.633869d583ab3a61d24ab10c8f21662f.mugqic.done
)
homer_make_tag_directory_28_JOB_ID=$(echo "#! /bin/bash
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
echo "$homer_make_tag_directory_28_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_29_JOB_ID: homer_make_tag_directory.MCF7_CTRL_WCE
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.MCF7_CTRL_WCE
JOB_DEPENDENCIES=$picard_mark_duplicates_29_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.MCF7_CTRL_WCE.1fb83957ec87ae7ca5a40a9c6afdf94f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.MCF7_CTRL_WCE.1fb83957ec87ae7ca5a40a9c6afdf94f.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/MCF7_CTRL_WCE \
            alignment/MCF7_CTRL_WCE/MCF7_CTRL_WCE.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.MCF7_CTRL_WCE.1fb83957ec87ae7ca5a40a9c6afdf94f.mugqic.done
)
homer_make_tag_directory_29_JOB_ID=$(echo "#! /bin/bash
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
echo "$homer_make_tag_directory_29_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_30_JOB_ID: homer_make_tag_directory.MCF7_E2_WCE
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.MCF7_E2_WCE
JOB_DEPENDENCIES=$picard_mark_duplicates_30_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.MCF7_E2_WCE.100c36d16a814c1f9506e1f48f3bacdc.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.MCF7_E2_WCE.100c36d16a814c1f9506e1f48f3bacdc.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/MCF7_E2_WCE \
            alignment/MCF7_E2_WCE/MCF7_E2_WCE.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.MCF7_E2_WCE.100c36d16a814c1f9506e1f48f3bacdc.mugqic.done
)
homer_make_tag_directory_30_JOB_ID=$(echo "#! /bin/bash
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
echo "$homer_make_tag_directory_30_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

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
JOB_DEPENDENCIES=$homer_make_tag_directory_1_JOB_ID:$homer_make_tag_directory_2_JOB_ID:$homer_make_tag_directory_3_JOB_ID:$homer_make_tag_directory_4_JOB_ID:$homer_make_tag_directory_5_JOB_ID:$homer_make_tag_directory_6_JOB_ID:$homer_make_tag_directory_7_JOB_ID:$homer_make_tag_directory_8_JOB_ID:$homer_make_tag_directory_9_JOB_ID:$homer_make_tag_directory_10_JOB_ID:$homer_make_tag_directory_11_JOB_ID:$homer_make_tag_directory_12_JOB_ID:$homer_make_tag_directory_13_JOB_ID:$homer_make_tag_directory_14_JOB_ID:$homer_make_tag_directory_15_JOB_ID:$homer_make_tag_directory_16_JOB_ID:$homer_make_tag_directory_17_JOB_ID:$homer_make_tag_directory_18_JOB_ID:$homer_make_tag_directory_19_JOB_ID:$homer_make_tag_directory_20_JOB_ID:$homer_make_tag_directory_21_JOB_ID:$homer_make_tag_directory_22_JOB_ID:$homer_make_tag_directory_23_JOB_ID:$homer_make_tag_directory_24_JOB_ID:$homer_make_tag_directory_25_JOB_ID:$homer_make_tag_directory_26_JOB_ID:$homer_make_tag_directory_27_JOB_ID:$homer_make_tag_directory_28_JOB_ID:$homer_make_tag_directory_29_JOB_ID:$homer_make_tag_directory_30_JOB_ID
JOB_DONE=job_output/qc_metrics/qc_plots_R.54c5c6735089a46d435151523cb2a12f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'qc_plots_R.54c5c6735089a46d435151523cb2a12f.mugqic.done'
module load mugqic/mugqic_tools/2.1.9 mugqic/R_Bioconductor/3.5.0_3.7 && \
mkdir -p graphs && \
Rscript $R_TOOLS/chipSeqGenerateQCMetrics.R \
  ../../raw/chip-seq/design_2019-02-23.txt \
  /scratch/chris11/sb_cofactor/MCF7/output/chip-pipeline-GRCh38 && \
cp /cvmfs/soft.mugqic/CentOS6/software/genpipes/genpipes-3.1.2/bfx/report/ChipSeq.qc_metrics.md report/ChipSeq.qc_metrics.md && \
for sample in MCF7_CTRL_WCE_rep1 MCF7_CTRL_NIPBL_rep1 MCF7_CTRL_SMC1A_rep1 MCF7_CTRL_MED1_rep1 MCF7_CTRL_POL2_rep1 MCF7_CTRL_ERA_rep1 MCF7_E2_WCE_rep1 MCF7_E2_NIPBL_rep1 MCF7_E2_SMC1A_rep1 MCF7_E2_MED1_rep1 MCF7_E2_POL2_rep1 MCF7_E2_ERA_rep1 MCF7_CTRL_WCE_rep2 MCF7_CTRL_NIPBL_rep2 MCF7_CTRL_SMC1A_rep2 MCF7_CTRL_MED1_rep2 MCF7_CTRL_ERA_rep2 MCF7_E2_WCE_rep2 MCF7_E2_NIPBL_rep2 MCF7_E2_SMC1A_rep2 MCF7_E2_MED1_rep2 MCF7_E2_ERA_rep2 MCF7_CTRL_BRD4_rep1 MCF7_CTRL_BRD4_rep2 MCF7_CTRL_BRD4_rep3 MCF7_E2_BRD4_rep1 MCF7_E2_BRD4_rep2 MCF7_E2_BRD4_rep3 MCF7_CTRL_WCE MCF7_E2_WCE
do
  cp --parents graphs/${sample}_QC_Metrics.ps report/
  convert -rotate 90 graphs/${sample}_QC_Metrics.ps report/graphs/${sample}_QC_Metrics.png
  echo -e "----

![QC Metrics for Sample $sample ([download high-res image](graphs/${sample}_QC_Metrics.ps))](graphs/${sample}_QC_Metrics.png)
" \
  >> report/ChipSeq.qc_metrics.md
done
qc_plots_R.54c5c6735089a46d435151523cb2a12f.mugqic.done
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
# JOB: macs2_callpeak_1_JOB_ID: macs2_callpeak.MCF7_CTRL_NIPBL_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.MCF7_CTRL_NIPBL_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_2_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.MCF7_CTRL_NIPBL_rep1.984940451f2df63444901fd83acdfb1c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.MCF7_CTRL_NIPBL_rep1.984940451f2df63444901fd83acdfb1c.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/MCF7_CTRL_NIPBL_rep1 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/MCF7_CTRL_NIPBL_rep1/MCF7_CTRL_NIPBL_rep1.sorted.dup.bam \
  --control \
  alignment/MCF7_CTRL_WCE_rep1/MCF7_CTRL_WCE_rep1.sorted.dup.bam \
  --name peak_call/MCF7_CTRL_NIPBL_rep1/MCF7_CTRL_NIPBL_rep1 \
  >& peak_call/MCF7_CTRL_NIPBL_rep1/MCF7_CTRL_NIPBL_rep1.diag.macs.out
macs2_callpeak.MCF7_CTRL_NIPBL_rep1.984940451f2df63444901fd83acdfb1c.mugqic.done
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
# JOB: macs2_callpeak_2_JOB_ID: macs2_callpeak_bigBed.MCF7_CTRL_NIPBL_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.MCF7_CTRL_NIPBL_rep1
JOB_DEPENDENCIES=$macs2_callpeak_1_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.MCF7_CTRL_NIPBL_rep1.412265a47e56ac11995d0b1e11e20711.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.MCF7_CTRL_NIPBL_rep1.412265a47e56ac11995d0b1e11e20711.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/MCF7_CTRL_NIPBL_rep1/MCF7_CTRL_NIPBL_rep1_peaks.narrowPeak > peak_call/MCF7_CTRL_NIPBL_rep1/MCF7_CTRL_NIPBL_rep1_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/MCF7_CTRL_NIPBL_rep1/MCF7_CTRL_NIPBL_rep1_peaks.narrowPeak.bed \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/MCF7_CTRL_NIPBL_rep1/MCF7_CTRL_NIPBL_rep1_peaks.narrowPeak.bb
macs2_callpeak_bigBed.MCF7_CTRL_NIPBL_rep1.412265a47e56ac11995d0b1e11e20711.mugqic.done
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
# JOB: macs2_callpeak_3_JOB_ID: macs2_callpeak.MCF7_CTRL_SMC1A_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.MCF7_CTRL_SMC1A_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_3_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.MCF7_CTRL_SMC1A_rep1.98dbad5bc533843f5cc8f37cffc2ad6c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.MCF7_CTRL_SMC1A_rep1.98dbad5bc533843f5cc8f37cffc2ad6c.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/MCF7_CTRL_SMC1A_rep1 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/MCF7_CTRL_SMC1A_rep1/MCF7_CTRL_SMC1A_rep1.sorted.dup.bam \
  --control \
  alignment/MCF7_CTRL_WCE_rep1/MCF7_CTRL_WCE_rep1.sorted.dup.bam \
  --name peak_call/MCF7_CTRL_SMC1A_rep1/MCF7_CTRL_SMC1A_rep1 \
  >& peak_call/MCF7_CTRL_SMC1A_rep1/MCF7_CTRL_SMC1A_rep1.diag.macs.out
macs2_callpeak.MCF7_CTRL_SMC1A_rep1.98dbad5bc533843f5cc8f37cffc2ad6c.mugqic.done
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
# JOB: macs2_callpeak_4_JOB_ID: macs2_callpeak_bigBed.MCF7_CTRL_SMC1A_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.MCF7_CTRL_SMC1A_rep1
JOB_DEPENDENCIES=$macs2_callpeak_3_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.MCF7_CTRL_SMC1A_rep1.0f342db0869a1ca4c3d10e0996dfa1e0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.MCF7_CTRL_SMC1A_rep1.0f342db0869a1ca4c3d10e0996dfa1e0.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/MCF7_CTRL_SMC1A_rep1/MCF7_CTRL_SMC1A_rep1_peaks.narrowPeak > peak_call/MCF7_CTRL_SMC1A_rep1/MCF7_CTRL_SMC1A_rep1_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/MCF7_CTRL_SMC1A_rep1/MCF7_CTRL_SMC1A_rep1_peaks.narrowPeak.bed \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/MCF7_CTRL_SMC1A_rep1/MCF7_CTRL_SMC1A_rep1_peaks.narrowPeak.bb
macs2_callpeak_bigBed.MCF7_CTRL_SMC1A_rep1.0f342db0869a1ca4c3d10e0996dfa1e0.mugqic.done
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
# JOB: macs2_callpeak_5_JOB_ID: macs2_callpeak.MCF7_CTRL_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.MCF7_CTRL_MED1_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_4_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.MCF7_CTRL_MED1_rep1.f3177c8b2f97ad6cb7881d9aa516fba1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.MCF7_CTRL_MED1_rep1.f3177c8b2f97ad6cb7881d9aa516fba1.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/MCF7_CTRL_MED1_rep1 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/MCF7_CTRL_MED1_rep1/MCF7_CTRL_MED1_rep1.sorted.dup.bam \
  --control \
  alignment/MCF7_CTRL_WCE_rep1/MCF7_CTRL_WCE_rep1.sorted.dup.bam \
  --name peak_call/MCF7_CTRL_MED1_rep1/MCF7_CTRL_MED1_rep1 \
  >& peak_call/MCF7_CTRL_MED1_rep1/MCF7_CTRL_MED1_rep1.diag.macs.out
macs2_callpeak.MCF7_CTRL_MED1_rep1.f3177c8b2f97ad6cb7881d9aa516fba1.mugqic.done
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
# JOB: macs2_callpeak_6_JOB_ID: macs2_callpeak_bigBed.MCF7_CTRL_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.MCF7_CTRL_MED1_rep1
JOB_DEPENDENCIES=$macs2_callpeak_5_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.MCF7_CTRL_MED1_rep1.62e71aca56657f3d12fb905878cbbcee.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.MCF7_CTRL_MED1_rep1.62e71aca56657f3d12fb905878cbbcee.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/MCF7_CTRL_MED1_rep1/MCF7_CTRL_MED1_rep1_peaks.narrowPeak > peak_call/MCF7_CTRL_MED1_rep1/MCF7_CTRL_MED1_rep1_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/MCF7_CTRL_MED1_rep1/MCF7_CTRL_MED1_rep1_peaks.narrowPeak.bed \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/MCF7_CTRL_MED1_rep1/MCF7_CTRL_MED1_rep1_peaks.narrowPeak.bb
macs2_callpeak_bigBed.MCF7_CTRL_MED1_rep1.62e71aca56657f3d12fb905878cbbcee.mugqic.done
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
# JOB: macs2_callpeak_7_JOB_ID: macs2_callpeak.MCF7_CTRL_POL2_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.MCF7_CTRL_POL2_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_5_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.MCF7_CTRL_POL2_rep1.0400a982b79fb5a3f8a5d2f2c1d51ae4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.MCF7_CTRL_POL2_rep1.0400a982b79fb5a3f8a5d2f2c1d51ae4.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/MCF7_CTRL_POL2_rep1 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/MCF7_CTRL_POL2_rep1/MCF7_CTRL_POL2_rep1.sorted.dup.bam \
  --control \
  alignment/MCF7_CTRL_WCE_rep1/MCF7_CTRL_WCE_rep1.sorted.dup.bam \
  --name peak_call/MCF7_CTRL_POL2_rep1/MCF7_CTRL_POL2_rep1 \
  >& peak_call/MCF7_CTRL_POL2_rep1/MCF7_CTRL_POL2_rep1.diag.macs.out
macs2_callpeak.MCF7_CTRL_POL2_rep1.0400a982b79fb5a3f8a5d2f2c1d51ae4.mugqic.done
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
# JOB: macs2_callpeak_8_JOB_ID: macs2_callpeak_bigBed.MCF7_CTRL_POL2_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.MCF7_CTRL_POL2_rep1
JOB_DEPENDENCIES=$macs2_callpeak_7_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.MCF7_CTRL_POL2_rep1.4b9a1bb93a4fc948cb1ed045d42bd02c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.MCF7_CTRL_POL2_rep1.4b9a1bb93a4fc948cb1ed045d42bd02c.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/MCF7_CTRL_POL2_rep1/MCF7_CTRL_POL2_rep1_peaks.narrowPeak > peak_call/MCF7_CTRL_POL2_rep1/MCF7_CTRL_POL2_rep1_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/MCF7_CTRL_POL2_rep1/MCF7_CTRL_POL2_rep1_peaks.narrowPeak.bed \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/MCF7_CTRL_POL2_rep1/MCF7_CTRL_POL2_rep1_peaks.narrowPeak.bb
macs2_callpeak_bigBed.MCF7_CTRL_POL2_rep1.4b9a1bb93a4fc948cb1ed045d42bd02c.mugqic.done
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
# JOB: macs2_callpeak_9_JOB_ID: macs2_callpeak.MCF7_CTRL_ERA_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.MCF7_CTRL_ERA_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_6_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.MCF7_CTRL_ERA_rep1.f3034b47a546197ce43b1f3bc44d3efb.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.MCF7_CTRL_ERA_rep1.f3034b47a546197ce43b1f3bc44d3efb.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/MCF7_CTRL_ERA_rep1 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/MCF7_CTRL_ERA_rep1/MCF7_CTRL_ERA_rep1.sorted.dup.bam \
  --control \
  alignment/MCF7_CTRL_WCE_rep1/MCF7_CTRL_WCE_rep1.sorted.dup.bam \
  --name peak_call/MCF7_CTRL_ERA_rep1/MCF7_CTRL_ERA_rep1 \
  >& peak_call/MCF7_CTRL_ERA_rep1/MCF7_CTRL_ERA_rep1.diag.macs.out
macs2_callpeak.MCF7_CTRL_ERA_rep1.f3034b47a546197ce43b1f3bc44d3efb.mugqic.done
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
# JOB: macs2_callpeak_10_JOB_ID: macs2_callpeak_bigBed.MCF7_CTRL_ERA_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.MCF7_CTRL_ERA_rep1
JOB_DEPENDENCIES=$macs2_callpeak_9_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.MCF7_CTRL_ERA_rep1.bff81ff9e4eaa6a8ab8506b459683eea.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.MCF7_CTRL_ERA_rep1.bff81ff9e4eaa6a8ab8506b459683eea.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/MCF7_CTRL_ERA_rep1/MCF7_CTRL_ERA_rep1_peaks.narrowPeak > peak_call/MCF7_CTRL_ERA_rep1/MCF7_CTRL_ERA_rep1_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/MCF7_CTRL_ERA_rep1/MCF7_CTRL_ERA_rep1_peaks.narrowPeak.bed \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/MCF7_CTRL_ERA_rep1/MCF7_CTRL_ERA_rep1_peaks.narrowPeak.bb
macs2_callpeak_bigBed.MCF7_CTRL_ERA_rep1.bff81ff9e4eaa6a8ab8506b459683eea.mugqic.done
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
# JOB: macs2_callpeak_11_JOB_ID: macs2_callpeak.MCF7_E2_NIPBL_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.MCF7_E2_NIPBL_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_7_JOB_ID:$picard_mark_duplicates_8_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.MCF7_E2_NIPBL_rep1.9371d05d6352c1d6bd9f210687dac453.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.MCF7_E2_NIPBL_rep1.9371d05d6352c1d6bd9f210687dac453.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/MCF7_E2_NIPBL_rep1 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/MCF7_E2_NIPBL_rep1/MCF7_E2_NIPBL_rep1.sorted.dup.bam \
  --control \
  alignment/MCF7_E2_WCE_rep1/MCF7_E2_WCE_rep1.sorted.dup.bam \
  --name peak_call/MCF7_E2_NIPBL_rep1/MCF7_E2_NIPBL_rep1 \
  >& peak_call/MCF7_E2_NIPBL_rep1/MCF7_E2_NIPBL_rep1.diag.macs.out
macs2_callpeak.MCF7_E2_NIPBL_rep1.9371d05d6352c1d6bd9f210687dac453.mugqic.done
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
# JOB: macs2_callpeak_12_JOB_ID: macs2_callpeak_bigBed.MCF7_E2_NIPBL_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.MCF7_E2_NIPBL_rep1
JOB_DEPENDENCIES=$macs2_callpeak_11_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.MCF7_E2_NIPBL_rep1.e73a260d756713683a9ae8f861843b24.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.MCF7_E2_NIPBL_rep1.e73a260d756713683a9ae8f861843b24.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/MCF7_E2_NIPBL_rep1/MCF7_E2_NIPBL_rep1_peaks.narrowPeak > peak_call/MCF7_E2_NIPBL_rep1/MCF7_E2_NIPBL_rep1_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/MCF7_E2_NIPBL_rep1/MCF7_E2_NIPBL_rep1_peaks.narrowPeak.bed \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/MCF7_E2_NIPBL_rep1/MCF7_E2_NIPBL_rep1_peaks.narrowPeak.bb
macs2_callpeak_bigBed.MCF7_E2_NIPBL_rep1.e73a260d756713683a9ae8f861843b24.mugqic.done
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
# JOB: macs2_callpeak_13_JOB_ID: macs2_callpeak.MCF7_E2_SMC1A_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.MCF7_E2_SMC1A_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_7_JOB_ID:$picard_mark_duplicates_9_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.MCF7_E2_SMC1A_rep1.dc5f9b6bc054daa9c7e3f4b9d144a99d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.MCF7_E2_SMC1A_rep1.dc5f9b6bc054daa9c7e3f4b9d144a99d.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/MCF7_E2_SMC1A_rep1 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/MCF7_E2_SMC1A_rep1/MCF7_E2_SMC1A_rep1.sorted.dup.bam \
  --control \
  alignment/MCF7_E2_WCE_rep1/MCF7_E2_WCE_rep1.sorted.dup.bam \
  --name peak_call/MCF7_E2_SMC1A_rep1/MCF7_E2_SMC1A_rep1 \
  >& peak_call/MCF7_E2_SMC1A_rep1/MCF7_E2_SMC1A_rep1.diag.macs.out
macs2_callpeak.MCF7_E2_SMC1A_rep1.dc5f9b6bc054daa9c7e3f4b9d144a99d.mugqic.done
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
# JOB: macs2_callpeak_14_JOB_ID: macs2_callpeak_bigBed.MCF7_E2_SMC1A_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.MCF7_E2_SMC1A_rep1
JOB_DEPENDENCIES=$macs2_callpeak_13_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.MCF7_E2_SMC1A_rep1.4a0c2b0e007846aef1fda10ca982f82a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.MCF7_E2_SMC1A_rep1.4a0c2b0e007846aef1fda10ca982f82a.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/MCF7_E2_SMC1A_rep1/MCF7_E2_SMC1A_rep1_peaks.narrowPeak > peak_call/MCF7_E2_SMC1A_rep1/MCF7_E2_SMC1A_rep1_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/MCF7_E2_SMC1A_rep1/MCF7_E2_SMC1A_rep1_peaks.narrowPeak.bed \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/MCF7_E2_SMC1A_rep1/MCF7_E2_SMC1A_rep1_peaks.narrowPeak.bb
macs2_callpeak_bigBed.MCF7_E2_SMC1A_rep1.4a0c2b0e007846aef1fda10ca982f82a.mugqic.done
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
# JOB: macs2_callpeak_15_JOB_ID: macs2_callpeak.MCF7_E2_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.MCF7_E2_MED1_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_7_JOB_ID:$picard_mark_duplicates_10_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.MCF7_E2_MED1_rep1.0174edee9c106677c9131a88bd72006f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.MCF7_E2_MED1_rep1.0174edee9c106677c9131a88bd72006f.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/MCF7_E2_MED1_rep1 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/MCF7_E2_MED1_rep1/MCF7_E2_MED1_rep1.sorted.dup.bam \
  --control \
  alignment/MCF7_E2_WCE_rep1/MCF7_E2_WCE_rep1.sorted.dup.bam \
  --name peak_call/MCF7_E2_MED1_rep1/MCF7_E2_MED1_rep1 \
  >& peak_call/MCF7_E2_MED1_rep1/MCF7_E2_MED1_rep1.diag.macs.out
macs2_callpeak.MCF7_E2_MED1_rep1.0174edee9c106677c9131a88bd72006f.mugqic.done
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
# JOB: macs2_callpeak_16_JOB_ID: macs2_callpeak_bigBed.MCF7_E2_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.MCF7_E2_MED1_rep1
JOB_DEPENDENCIES=$macs2_callpeak_15_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.MCF7_E2_MED1_rep1.a3488db2d66e5cedf01f8e02cf82d08e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.MCF7_E2_MED1_rep1.a3488db2d66e5cedf01f8e02cf82d08e.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/MCF7_E2_MED1_rep1/MCF7_E2_MED1_rep1_peaks.narrowPeak > peak_call/MCF7_E2_MED1_rep1/MCF7_E2_MED1_rep1_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/MCF7_E2_MED1_rep1/MCF7_E2_MED1_rep1_peaks.narrowPeak.bed \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/MCF7_E2_MED1_rep1/MCF7_E2_MED1_rep1_peaks.narrowPeak.bb
macs2_callpeak_bigBed.MCF7_E2_MED1_rep1.a3488db2d66e5cedf01f8e02cf82d08e.mugqic.done
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
# JOB: macs2_callpeak_17_JOB_ID: macs2_callpeak.MCF7_E2_POL2_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.MCF7_E2_POL2_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_7_JOB_ID:$picard_mark_duplicates_11_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.MCF7_E2_POL2_rep1.acbbdb1d1405c3a6349c517a1151d7d3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.MCF7_E2_POL2_rep1.acbbdb1d1405c3a6349c517a1151d7d3.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/MCF7_E2_POL2_rep1 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/MCF7_E2_POL2_rep1/MCF7_E2_POL2_rep1.sorted.dup.bam \
  --control \
  alignment/MCF7_E2_WCE_rep1/MCF7_E2_WCE_rep1.sorted.dup.bam \
  --name peak_call/MCF7_E2_POL2_rep1/MCF7_E2_POL2_rep1 \
  >& peak_call/MCF7_E2_POL2_rep1/MCF7_E2_POL2_rep1.diag.macs.out
macs2_callpeak.MCF7_E2_POL2_rep1.acbbdb1d1405c3a6349c517a1151d7d3.mugqic.done
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
# JOB: macs2_callpeak_18_JOB_ID: macs2_callpeak_bigBed.MCF7_E2_POL2_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.MCF7_E2_POL2_rep1
JOB_DEPENDENCIES=$macs2_callpeak_17_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.MCF7_E2_POL2_rep1.92259c2200da09ac7502bb1cb0362d07.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.MCF7_E2_POL2_rep1.92259c2200da09ac7502bb1cb0362d07.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/MCF7_E2_POL2_rep1/MCF7_E2_POL2_rep1_peaks.narrowPeak > peak_call/MCF7_E2_POL2_rep1/MCF7_E2_POL2_rep1_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/MCF7_E2_POL2_rep1/MCF7_E2_POL2_rep1_peaks.narrowPeak.bed \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/MCF7_E2_POL2_rep1/MCF7_E2_POL2_rep1_peaks.narrowPeak.bb
macs2_callpeak_bigBed.MCF7_E2_POL2_rep1.92259c2200da09ac7502bb1cb0362d07.mugqic.done
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
# JOB: macs2_callpeak_19_JOB_ID: macs2_callpeak.MCF7_E2_ERA_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.MCF7_E2_ERA_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_7_JOB_ID:$picard_mark_duplicates_12_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.MCF7_E2_ERA_rep1.cf83b6d8e9b05353b3f3724db66264f7.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.MCF7_E2_ERA_rep1.cf83b6d8e9b05353b3f3724db66264f7.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/MCF7_E2_ERA_rep1 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/MCF7_E2_ERA_rep1/MCF7_E2_ERA_rep1.sorted.dup.bam \
  --control \
  alignment/MCF7_E2_WCE_rep1/MCF7_E2_WCE_rep1.sorted.dup.bam \
  --name peak_call/MCF7_E2_ERA_rep1/MCF7_E2_ERA_rep1 \
  >& peak_call/MCF7_E2_ERA_rep1/MCF7_E2_ERA_rep1.diag.macs.out
macs2_callpeak.MCF7_E2_ERA_rep1.cf83b6d8e9b05353b3f3724db66264f7.mugqic.done
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
# JOB: macs2_callpeak_20_JOB_ID: macs2_callpeak_bigBed.MCF7_E2_ERA_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.MCF7_E2_ERA_rep1
JOB_DEPENDENCIES=$macs2_callpeak_19_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.MCF7_E2_ERA_rep1.0fd51d8f00be7d2633a3ea5e5ed70bbc.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.MCF7_E2_ERA_rep1.0fd51d8f00be7d2633a3ea5e5ed70bbc.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/MCF7_E2_ERA_rep1/MCF7_E2_ERA_rep1_peaks.narrowPeak > peak_call/MCF7_E2_ERA_rep1/MCF7_E2_ERA_rep1_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/MCF7_E2_ERA_rep1/MCF7_E2_ERA_rep1_peaks.narrowPeak.bed \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/MCF7_E2_ERA_rep1/MCF7_E2_ERA_rep1_peaks.narrowPeak.bb
macs2_callpeak_bigBed.MCF7_E2_ERA_rep1.0fd51d8f00be7d2633a3ea5e5ed70bbc.mugqic.done
)
macs2_callpeak_20_JOB_ID=$(echo "#! /bin/bash
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
echo "$macs2_callpeak_20_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_21_JOB_ID: macs2_callpeak.MCF7_CTRL_NIPBL_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.MCF7_CTRL_NIPBL_rep2
JOB_DEPENDENCIES=$picard_mark_duplicates_13_JOB_ID:$picard_mark_duplicates_14_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.MCF7_CTRL_NIPBL_rep2.9cbea0a7c123522a67c1fa93ed03b079.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.MCF7_CTRL_NIPBL_rep2.9cbea0a7c123522a67c1fa93ed03b079.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/MCF7_CTRL_NIPBL_rep2 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/MCF7_CTRL_NIPBL_rep2/MCF7_CTRL_NIPBL_rep2.sorted.dup.bam \
  --control \
  alignment/MCF7_CTRL_WCE_rep2/MCF7_CTRL_WCE_rep2.sorted.dup.bam \
  --name peak_call/MCF7_CTRL_NIPBL_rep2/MCF7_CTRL_NIPBL_rep2 \
  >& peak_call/MCF7_CTRL_NIPBL_rep2/MCF7_CTRL_NIPBL_rep2.diag.macs.out
macs2_callpeak.MCF7_CTRL_NIPBL_rep2.9cbea0a7c123522a67c1fa93ed03b079.mugqic.done
)
macs2_callpeak_21_JOB_ID=$(echo "#! /bin/bash
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
echo "$macs2_callpeak_21_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_22_JOB_ID: macs2_callpeak_bigBed.MCF7_CTRL_NIPBL_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.MCF7_CTRL_NIPBL_rep2
JOB_DEPENDENCIES=$macs2_callpeak_21_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.MCF7_CTRL_NIPBL_rep2.2b0b5992276f211097f220609d17a985.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.MCF7_CTRL_NIPBL_rep2.2b0b5992276f211097f220609d17a985.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/MCF7_CTRL_NIPBL_rep2/MCF7_CTRL_NIPBL_rep2_peaks.narrowPeak > peak_call/MCF7_CTRL_NIPBL_rep2/MCF7_CTRL_NIPBL_rep2_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/MCF7_CTRL_NIPBL_rep2/MCF7_CTRL_NIPBL_rep2_peaks.narrowPeak.bed \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/MCF7_CTRL_NIPBL_rep2/MCF7_CTRL_NIPBL_rep2_peaks.narrowPeak.bb
macs2_callpeak_bigBed.MCF7_CTRL_NIPBL_rep2.2b0b5992276f211097f220609d17a985.mugqic.done
)
macs2_callpeak_22_JOB_ID=$(echo "#! /bin/bash
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
echo "$macs2_callpeak_22_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_23_JOB_ID: macs2_callpeak.MCF7_CTRL_SMC1A_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.MCF7_CTRL_SMC1A_rep2
JOB_DEPENDENCIES=$picard_mark_duplicates_13_JOB_ID:$picard_mark_duplicates_15_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.MCF7_CTRL_SMC1A_rep2.384ebd58a3743ea46004350f801f8b5b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.MCF7_CTRL_SMC1A_rep2.384ebd58a3743ea46004350f801f8b5b.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/MCF7_CTRL_SMC1A_rep2 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/MCF7_CTRL_SMC1A_rep2/MCF7_CTRL_SMC1A_rep2.sorted.dup.bam \
  --control \
  alignment/MCF7_CTRL_WCE_rep2/MCF7_CTRL_WCE_rep2.sorted.dup.bam \
  --name peak_call/MCF7_CTRL_SMC1A_rep2/MCF7_CTRL_SMC1A_rep2 \
  >& peak_call/MCF7_CTRL_SMC1A_rep2/MCF7_CTRL_SMC1A_rep2.diag.macs.out
macs2_callpeak.MCF7_CTRL_SMC1A_rep2.384ebd58a3743ea46004350f801f8b5b.mugqic.done
)
macs2_callpeak_23_JOB_ID=$(echo "#! /bin/bash
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
echo "$macs2_callpeak_23_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_24_JOB_ID: macs2_callpeak_bigBed.MCF7_CTRL_SMC1A_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.MCF7_CTRL_SMC1A_rep2
JOB_DEPENDENCIES=$macs2_callpeak_23_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.MCF7_CTRL_SMC1A_rep2.b64eb38ffaa12d2dfba6da1e82019586.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.MCF7_CTRL_SMC1A_rep2.b64eb38ffaa12d2dfba6da1e82019586.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/MCF7_CTRL_SMC1A_rep2/MCF7_CTRL_SMC1A_rep2_peaks.narrowPeak > peak_call/MCF7_CTRL_SMC1A_rep2/MCF7_CTRL_SMC1A_rep2_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/MCF7_CTRL_SMC1A_rep2/MCF7_CTRL_SMC1A_rep2_peaks.narrowPeak.bed \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/MCF7_CTRL_SMC1A_rep2/MCF7_CTRL_SMC1A_rep2_peaks.narrowPeak.bb
macs2_callpeak_bigBed.MCF7_CTRL_SMC1A_rep2.b64eb38ffaa12d2dfba6da1e82019586.mugqic.done
)
macs2_callpeak_24_JOB_ID=$(echo "#! /bin/bash
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
echo "$macs2_callpeak_24_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_25_JOB_ID: macs2_callpeak.MCF7_CTRL_MED1_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.MCF7_CTRL_MED1_rep2
JOB_DEPENDENCIES=$picard_mark_duplicates_13_JOB_ID:$picard_mark_duplicates_16_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.MCF7_CTRL_MED1_rep2.00849fdf28481ebd82fbe4048ceceda5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.MCF7_CTRL_MED1_rep2.00849fdf28481ebd82fbe4048ceceda5.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/MCF7_CTRL_MED1_rep2 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/MCF7_CTRL_MED1_rep2/MCF7_CTRL_MED1_rep2.sorted.dup.bam \
  --control \
  alignment/MCF7_CTRL_WCE_rep2/MCF7_CTRL_WCE_rep2.sorted.dup.bam \
  --name peak_call/MCF7_CTRL_MED1_rep2/MCF7_CTRL_MED1_rep2 \
  >& peak_call/MCF7_CTRL_MED1_rep2/MCF7_CTRL_MED1_rep2.diag.macs.out
macs2_callpeak.MCF7_CTRL_MED1_rep2.00849fdf28481ebd82fbe4048ceceda5.mugqic.done
)
macs2_callpeak_25_JOB_ID=$(echo "#! /bin/bash
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
echo "$macs2_callpeak_25_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_26_JOB_ID: macs2_callpeak_bigBed.MCF7_CTRL_MED1_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.MCF7_CTRL_MED1_rep2
JOB_DEPENDENCIES=$macs2_callpeak_25_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.MCF7_CTRL_MED1_rep2.50d7eabbc3582d2132a9a26b0eedd63f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.MCF7_CTRL_MED1_rep2.50d7eabbc3582d2132a9a26b0eedd63f.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/MCF7_CTRL_MED1_rep2/MCF7_CTRL_MED1_rep2_peaks.narrowPeak > peak_call/MCF7_CTRL_MED1_rep2/MCF7_CTRL_MED1_rep2_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/MCF7_CTRL_MED1_rep2/MCF7_CTRL_MED1_rep2_peaks.narrowPeak.bed \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/MCF7_CTRL_MED1_rep2/MCF7_CTRL_MED1_rep2_peaks.narrowPeak.bb
macs2_callpeak_bigBed.MCF7_CTRL_MED1_rep2.50d7eabbc3582d2132a9a26b0eedd63f.mugqic.done
)
macs2_callpeak_26_JOB_ID=$(echo "#! /bin/bash
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
echo "$macs2_callpeak_26_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_27_JOB_ID: macs2_callpeak.MCF7_CTRL_ERA_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.MCF7_CTRL_ERA_rep2
JOB_DEPENDENCIES=$picard_mark_duplicates_13_JOB_ID:$picard_mark_duplicates_17_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.MCF7_CTRL_ERA_rep2.882d40cdda731398e407b7e6f8e5dd8f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.MCF7_CTRL_ERA_rep2.882d40cdda731398e407b7e6f8e5dd8f.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/MCF7_CTRL_ERA_rep2 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/MCF7_CTRL_ERA_rep2/MCF7_CTRL_ERA_rep2.sorted.dup.bam \
  --control \
  alignment/MCF7_CTRL_WCE_rep2/MCF7_CTRL_WCE_rep2.sorted.dup.bam \
  --name peak_call/MCF7_CTRL_ERA_rep2/MCF7_CTRL_ERA_rep2 \
  >& peak_call/MCF7_CTRL_ERA_rep2/MCF7_CTRL_ERA_rep2.diag.macs.out
macs2_callpeak.MCF7_CTRL_ERA_rep2.882d40cdda731398e407b7e6f8e5dd8f.mugqic.done
)
macs2_callpeak_27_JOB_ID=$(echo "#! /bin/bash
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
echo "$macs2_callpeak_27_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_28_JOB_ID: macs2_callpeak_bigBed.MCF7_CTRL_ERA_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.MCF7_CTRL_ERA_rep2
JOB_DEPENDENCIES=$macs2_callpeak_27_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.MCF7_CTRL_ERA_rep2.5d9075bbc979e6d7d71dca4a7cfba15b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.MCF7_CTRL_ERA_rep2.5d9075bbc979e6d7d71dca4a7cfba15b.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/MCF7_CTRL_ERA_rep2/MCF7_CTRL_ERA_rep2_peaks.narrowPeak > peak_call/MCF7_CTRL_ERA_rep2/MCF7_CTRL_ERA_rep2_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/MCF7_CTRL_ERA_rep2/MCF7_CTRL_ERA_rep2_peaks.narrowPeak.bed \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/MCF7_CTRL_ERA_rep2/MCF7_CTRL_ERA_rep2_peaks.narrowPeak.bb
macs2_callpeak_bigBed.MCF7_CTRL_ERA_rep2.5d9075bbc979e6d7d71dca4a7cfba15b.mugqic.done
)
macs2_callpeak_28_JOB_ID=$(echo "#! /bin/bash
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
echo "$macs2_callpeak_28_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_29_JOB_ID: macs2_callpeak.MCF7_E2_NIPBL_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.MCF7_E2_NIPBL_rep2
JOB_DEPENDENCIES=$picard_mark_duplicates_18_JOB_ID:$picard_mark_duplicates_19_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.MCF7_E2_NIPBL_rep2.0a964b4c8c05af298aa7b497d164ba76.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.MCF7_E2_NIPBL_rep2.0a964b4c8c05af298aa7b497d164ba76.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/MCF7_E2_NIPBL_rep2 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/MCF7_E2_NIPBL_rep2/MCF7_E2_NIPBL_rep2.sorted.dup.bam \
  --control \
  alignment/MCF7_E2_WCE_rep2/MCF7_E2_WCE_rep2.sorted.dup.bam \
  --name peak_call/MCF7_E2_NIPBL_rep2/MCF7_E2_NIPBL_rep2 \
  >& peak_call/MCF7_E2_NIPBL_rep2/MCF7_E2_NIPBL_rep2.diag.macs.out
macs2_callpeak.MCF7_E2_NIPBL_rep2.0a964b4c8c05af298aa7b497d164ba76.mugqic.done
)
macs2_callpeak_29_JOB_ID=$(echo "#! /bin/bash
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
echo "$macs2_callpeak_29_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_30_JOB_ID: macs2_callpeak_bigBed.MCF7_E2_NIPBL_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.MCF7_E2_NIPBL_rep2
JOB_DEPENDENCIES=$macs2_callpeak_29_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.MCF7_E2_NIPBL_rep2.ae9181ef7cecf438142cd241a70fa0d3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.MCF7_E2_NIPBL_rep2.ae9181ef7cecf438142cd241a70fa0d3.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/MCF7_E2_NIPBL_rep2/MCF7_E2_NIPBL_rep2_peaks.narrowPeak > peak_call/MCF7_E2_NIPBL_rep2/MCF7_E2_NIPBL_rep2_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/MCF7_E2_NIPBL_rep2/MCF7_E2_NIPBL_rep2_peaks.narrowPeak.bed \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/MCF7_E2_NIPBL_rep2/MCF7_E2_NIPBL_rep2_peaks.narrowPeak.bb
macs2_callpeak_bigBed.MCF7_E2_NIPBL_rep2.ae9181ef7cecf438142cd241a70fa0d3.mugqic.done
)
macs2_callpeak_30_JOB_ID=$(echo "#! /bin/bash
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
echo "$macs2_callpeak_30_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_31_JOB_ID: macs2_callpeak.MCF7_E2_SMC1A_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.MCF7_E2_SMC1A_rep2
JOB_DEPENDENCIES=$picard_mark_duplicates_18_JOB_ID:$picard_mark_duplicates_20_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.MCF7_E2_SMC1A_rep2.53f9d74e9a241584a889e3d85c1a37eb.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.MCF7_E2_SMC1A_rep2.53f9d74e9a241584a889e3d85c1a37eb.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/MCF7_E2_SMC1A_rep2 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/MCF7_E2_SMC1A_rep2/MCF7_E2_SMC1A_rep2.sorted.dup.bam \
  --control \
  alignment/MCF7_E2_WCE_rep2/MCF7_E2_WCE_rep2.sorted.dup.bam \
  --name peak_call/MCF7_E2_SMC1A_rep2/MCF7_E2_SMC1A_rep2 \
  >& peak_call/MCF7_E2_SMC1A_rep2/MCF7_E2_SMC1A_rep2.diag.macs.out
macs2_callpeak.MCF7_E2_SMC1A_rep2.53f9d74e9a241584a889e3d85c1a37eb.mugqic.done
)
macs2_callpeak_31_JOB_ID=$(echo "#! /bin/bash
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
echo "$macs2_callpeak_31_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_32_JOB_ID: macs2_callpeak_bigBed.MCF7_E2_SMC1A_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.MCF7_E2_SMC1A_rep2
JOB_DEPENDENCIES=$macs2_callpeak_31_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.MCF7_E2_SMC1A_rep2.d442cb9cb1257816a9e13950ed4cf156.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.MCF7_E2_SMC1A_rep2.d442cb9cb1257816a9e13950ed4cf156.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/MCF7_E2_SMC1A_rep2/MCF7_E2_SMC1A_rep2_peaks.narrowPeak > peak_call/MCF7_E2_SMC1A_rep2/MCF7_E2_SMC1A_rep2_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/MCF7_E2_SMC1A_rep2/MCF7_E2_SMC1A_rep2_peaks.narrowPeak.bed \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/MCF7_E2_SMC1A_rep2/MCF7_E2_SMC1A_rep2_peaks.narrowPeak.bb
macs2_callpeak_bigBed.MCF7_E2_SMC1A_rep2.d442cb9cb1257816a9e13950ed4cf156.mugqic.done
)
macs2_callpeak_32_JOB_ID=$(echo "#! /bin/bash
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
echo "$macs2_callpeak_32_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_33_JOB_ID: macs2_callpeak.MCF7_E2_MED1_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.MCF7_E2_MED1_rep2
JOB_DEPENDENCIES=$picard_mark_duplicates_18_JOB_ID:$picard_mark_duplicates_21_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.MCF7_E2_MED1_rep2.abab655e225f7e81926c4ff5bbb46897.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.MCF7_E2_MED1_rep2.abab655e225f7e81926c4ff5bbb46897.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/MCF7_E2_MED1_rep2 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/MCF7_E2_MED1_rep2/MCF7_E2_MED1_rep2.sorted.dup.bam \
  --control \
  alignment/MCF7_E2_WCE_rep2/MCF7_E2_WCE_rep2.sorted.dup.bam \
  --name peak_call/MCF7_E2_MED1_rep2/MCF7_E2_MED1_rep2 \
  >& peak_call/MCF7_E2_MED1_rep2/MCF7_E2_MED1_rep2.diag.macs.out
macs2_callpeak.MCF7_E2_MED1_rep2.abab655e225f7e81926c4ff5bbb46897.mugqic.done
)
macs2_callpeak_33_JOB_ID=$(echo "#! /bin/bash
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
echo "$macs2_callpeak_33_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_34_JOB_ID: macs2_callpeak_bigBed.MCF7_E2_MED1_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.MCF7_E2_MED1_rep2
JOB_DEPENDENCIES=$macs2_callpeak_33_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.MCF7_E2_MED1_rep2.e14312ac41e7531ef780563ba06072d8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.MCF7_E2_MED1_rep2.e14312ac41e7531ef780563ba06072d8.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/MCF7_E2_MED1_rep2/MCF7_E2_MED1_rep2_peaks.narrowPeak > peak_call/MCF7_E2_MED1_rep2/MCF7_E2_MED1_rep2_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/MCF7_E2_MED1_rep2/MCF7_E2_MED1_rep2_peaks.narrowPeak.bed \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/MCF7_E2_MED1_rep2/MCF7_E2_MED1_rep2_peaks.narrowPeak.bb
macs2_callpeak_bigBed.MCF7_E2_MED1_rep2.e14312ac41e7531ef780563ba06072d8.mugqic.done
)
macs2_callpeak_34_JOB_ID=$(echo "#! /bin/bash
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
echo "$macs2_callpeak_34_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_35_JOB_ID: macs2_callpeak.MCF7_E2_ERA_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.MCF7_E2_ERA_rep2
JOB_DEPENDENCIES=$picard_mark_duplicates_18_JOB_ID:$picard_mark_duplicates_22_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.MCF7_E2_ERA_rep2.a89b1777747bcf8c2da1d5477c14ec91.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.MCF7_E2_ERA_rep2.a89b1777747bcf8c2da1d5477c14ec91.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/MCF7_E2_ERA_rep2 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/MCF7_E2_ERA_rep2/MCF7_E2_ERA_rep2.sorted.dup.bam \
  --control \
  alignment/MCF7_E2_WCE_rep2/MCF7_E2_WCE_rep2.sorted.dup.bam \
  --name peak_call/MCF7_E2_ERA_rep2/MCF7_E2_ERA_rep2 \
  >& peak_call/MCF7_E2_ERA_rep2/MCF7_E2_ERA_rep2.diag.macs.out
macs2_callpeak.MCF7_E2_ERA_rep2.a89b1777747bcf8c2da1d5477c14ec91.mugqic.done
)
macs2_callpeak_35_JOB_ID=$(echo "#! /bin/bash
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
echo "$macs2_callpeak_35_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_36_JOB_ID: macs2_callpeak_bigBed.MCF7_E2_ERA_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.MCF7_E2_ERA_rep2
JOB_DEPENDENCIES=$macs2_callpeak_35_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.MCF7_E2_ERA_rep2.729ad26b5d730fb3ab3666f0ca812a65.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.MCF7_E2_ERA_rep2.729ad26b5d730fb3ab3666f0ca812a65.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/MCF7_E2_ERA_rep2/MCF7_E2_ERA_rep2_peaks.narrowPeak > peak_call/MCF7_E2_ERA_rep2/MCF7_E2_ERA_rep2_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/MCF7_E2_ERA_rep2/MCF7_E2_ERA_rep2_peaks.narrowPeak.bed \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/MCF7_E2_ERA_rep2/MCF7_E2_ERA_rep2_peaks.narrowPeak.bb
macs2_callpeak_bigBed.MCF7_E2_ERA_rep2.729ad26b5d730fb3ab3666f0ca812a65.mugqic.done
)
macs2_callpeak_36_JOB_ID=$(echo "#! /bin/bash
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
echo "$macs2_callpeak_36_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_37_JOB_ID: macs2_callpeak.MCF7_CTRL_BRD4
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.MCF7_CTRL_BRD4
JOB_DEPENDENCIES=$picard_mark_duplicates_23_JOB_ID:$picard_mark_duplicates_24_JOB_ID:$picard_mark_duplicates_25_JOB_ID:$picard_mark_duplicates_29_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.MCF7_CTRL_BRD4.d460e42a8b1d8aa236077f3e84607a1e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.MCF7_CTRL_BRD4.d460e42a8b1d8aa236077f3e84607a1e.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/MCF7_CTRL_BRD4 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/MCF7_CTRL_BRD4_rep1/MCF7_CTRL_BRD4_rep1.sorted.dup.bam \
  alignment/MCF7_CTRL_BRD4_rep2/MCF7_CTRL_BRD4_rep2.sorted.dup.bam \
  alignment/MCF7_CTRL_BRD4_rep3/MCF7_CTRL_BRD4_rep3.sorted.dup.bam \
  --control \
  alignment/MCF7_CTRL_WCE/MCF7_CTRL_WCE.sorted.dup.bam \
  --name peak_call/MCF7_CTRL_BRD4/MCF7_CTRL_BRD4 \
  >& peak_call/MCF7_CTRL_BRD4/MCF7_CTRL_BRD4.diag.macs.out
macs2_callpeak.MCF7_CTRL_BRD4.d460e42a8b1d8aa236077f3e84607a1e.mugqic.done
)
macs2_callpeak_37_JOB_ID=$(echo "#! /bin/bash
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
echo "$macs2_callpeak_37_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_38_JOB_ID: macs2_callpeak_bigBed.MCF7_CTRL_BRD4
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.MCF7_CTRL_BRD4
JOB_DEPENDENCIES=$macs2_callpeak_37_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.MCF7_CTRL_BRD4.6988d971dfc079f0892812c1c913740d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.MCF7_CTRL_BRD4.6988d971dfc079f0892812c1c913740d.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/MCF7_CTRL_BRD4/MCF7_CTRL_BRD4_peaks.narrowPeak > peak_call/MCF7_CTRL_BRD4/MCF7_CTRL_BRD4_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/MCF7_CTRL_BRD4/MCF7_CTRL_BRD4_peaks.narrowPeak.bed \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/MCF7_CTRL_BRD4/MCF7_CTRL_BRD4_peaks.narrowPeak.bb
macs2_callpeak_bigBed.MCF7_CTRL_BRD4.6988d971dfc079f0892812c1c913740d.mugqic.done
)
macs2_callpeak_38_JOB_ID=$(echo "#! /bin/bash
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
echo "$macs2_callpeak_38_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_39_JOB_ID: macs2_callpeak.MCF7_E2_BRD4
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.MCF7_E2_BRD4
JOB_DEPENDENCIES=$picard_mark_duplicates_26_JOB_ID:$picard_mark_duplicates_27_JOB_ID:$picard_mark_duplicates_28_JOB_ID:$picard_mark_duplicates_30_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.MCF7_E2_BRD4.c2d43e9c248664e075b8e417a3d51627.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.MCF7_E2_BRD4.c2d43e9c248664e075b8e417a3d51627.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/MCF7_E2_BRD4 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/MCF7_E2_BRD4_rep1/MCF7_E2_BRD4_rep1.sorted.dup.bam \
  alignment/MCF7_E2_BRD4_rep2/MCF7_E2_BRD4_rep2.sorted.dup.bam \
  alignment/MCF7_E2_BRD4_rep3/MCF7_E2_BRD4_rep3.sorted.dup.bam \
  --control \
  alignment/MCF7_E2_WCE/MCF7_E2_WCE.sorted.dup.bam \
  --name peak_call/MCF7_E2_BRD4/MCF7_E2_BRD4 \
  >& peak_call/MCF7_E2_BRD4/MCF7_E2_BRD4.diag.macs.out
macs2_callpeak.MCF7_E2_BRD4.c2d43e9c248664e075b8e417a3d51627.mugqic.done
)
macs2_callpeak_39_JOB_ID=$(echo "#! /bin/bash
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
echo "$macs2_callpeak_39_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_40_JOB_ID: macs2_callpeak_bigBed.MCF7_E2_BRD4
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.MCF7_E2_BRD4
JOB_DEPENDENCIES=$macs2_callpeak_39_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.MCF7_E2_BRD4.dc96e81134f889ea4b58d1e23f838dbc.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.MCF7_E2_BRD4.dc96e81134f889ea4b58d1e23f838dbc.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/MCF7_E2_BRD4/MCF7_E2_BRD4_peaks.narrowPeak > peak_call/MCF7_E2_BRD4/MCF7_E2_BRD4_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/MCF7_E2_BRD4/MCF7_E2_BRD4_peaks.narrowPeak.bed \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/MCF7_E2_BRD4/MCF7_E2_BRD4_peaks.narrowPeak.bb
macs2_callpeak_bigBed.MCF7_E2_BRD4.dc96e81134f889ea4b58d1e23f838dbc.mugqic.done
)
macs2_callpeak_40_JOB_ID=$(echo "#! /bin/bash
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
echo "$macs2_callpeak_40_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_41_JOB_ID: macs2_callpeak_report
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_report
JOB_DEPENDENCIES=$macs2_callpeak_1_JOB_ID:$macs2_callpeak_3_JOB_ID:$macs2_callpeak_5_JOB_ID:$macs2_callpeak_7_JOB_ID:$macs2_callpeak_9_JOB_ID:$macs2_callpeak_11_JOB_ID:$macs2_callpeak_13_JOB_ID:$macs2_callpeak_15_JOB_ID:$macs2_callpeak_17_JOB_ID:$macs2_callpeak_19_JOB_ID:$macs2_callpeak_21_JOB_ID:$macs2_callpeak_23_JOB_ID:$macs2_callpeak_25_JOB_ID:$macs2_callpeak_27_JOB_ID:$macs2_callpeak_29_JOB_ID:$macs2_callpeak_31_JOB_ID:$macs2_callpeak_33_JOB_ID:$macs2_callpeak_35_JOB_ID:$macs2_callpeak_37_JOB_ID:$macs2_callpeak_39_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_report.28ca416ee5777b760a14b1c31d2db310.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_report.28ca416ee5777b760a14b1c31d2db310.mugqic.done'
mkdir -p report && \
cp /cvmfs/soft.mugqic/CentOS6/software/genpipes/genpipes-3.1.2/bfx/report/ChipSeq.macs2_callpeak.md report/ && \
for contrast in MCF7_CTRL_NIPBL_rep1 MCF7_CTRL_SMC1A_rep1 MCF7_CTRL_MED1_rep1 MCF7_CTRL_POL2_rep1 MCF7_CTRL_ERA_rep1 MCF7_E2_NIPBL_rep1 MCF7_E2_SMC1A_rep1 MCF7_E2_MED1_rep1 MCF7_E2_POL2_rep1 MCF7_E2_ERA_rep1 MCF7_CTRL_NIPBL_rep2 MCF7_CTRL_SMC1A_rep2 MCF7_CTRL_MED1_rep2 MCF7_CTRL_ERA_rep2 MCF7_E2_NIPBL_rep2 MCF7_E2_SMC1A_rep2 MCF7_E2_MED1_rep2 MCF7_E2_ERA_rep2 MCF7_CTRL_BRD4 MCF7_E2_BRD4
do
  cp -a --parents peak_call/$contrast/ report/ && \
  echo -e "* [Peak Calls File for Design $contrast](peak_call/$contrast/${contrast}_peaks.xls)" \
  >> report/ChipSeq.macs2_callpeak.md
done
macs2_callpeak_report.28ca416ee5777b760a14b1c31d2db310.mugqic.done
)
macs2_callpeak_41_JOB_ID=$(echo "#! /bin/bash
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
echo "$macs2_callpeak_41_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# Call home with pipeline statistics
#-------------------------------------------------------------------------------
LOG_MD5=$(echo $USER-'206.12.124.2-ChipSeq-MCF7_CTRL_ERA_rep2.HI.2791.003.Index_10.MCF7_3_Era_C_rep2.bam,MCF7_CTRL_ERA_rep1.HI.1997.005.Index_25.MCF7_CTRL_ERA_ChIP8_GB_rep1.bam,MCF7_CTRL_WCE_rep1.HI.1943.008.Index_2.MCF7_CTRL_WCE_ChIP8_GB_rep1.bam,MCF7_CTRL_WCE_rep2.HI.2791.003.Index_5.MCF7_WCE_C_rep2.bam,MCF7_E2_WCE.SRR1193563.fastq,MCF7_E2_NIPBL_rep2.HI.2791.003.Index_12.MCF7_1_NIPBL_E_rep2.bam,MCF7_E2_NIPBL_rep1.HI.1997.004.Index_3.MCF7_E2_NIPBL_ChIP8_GB_rep1.bam,MCF7_CTRL_WCE.SRR1193562.fastq,MCF7_CTRL_POL2_rep1.HI.1997.005.Index_10.MCF7_CTRL_POL2_ChIP8_GB_rep1.bam,MCF7_CTRL_BRD4_rep1.SRR1193526.fastq,MCF7_CTRL_BRD4_rep3.SRR1193528.fastq,MCF7_E2_POL2_rep1.HI.1997.005.Index_11.MCF7_E2_POL2_ChIP8_GB_rep1.bam,MCF7_E2_WCE_rep1.HI.1943.008.Index_13.MCF7_E2_WCE_ChIP8_GB_rep1.bam,MCF7_E2_MED1_rep1.HI.1943.008.Index_18.MCF7_E2_MED1_ChIP8_GB_rep1.bam,MCF7_E2_MED1_rep2.HI.2791.003.Index_19.MCF7_2_MED1_E_rep2.bam,MCF7_E2_WCE_rep2.HI.2791.003.Index_4.MCF7_WCE_E_rep2.bam,MCF7_E2_ERA_rep2.HI.2791.003.Index_1.MCF7_3_Era_E_rep2.bam,MCF7_E2_ERA_rep1.HI.1997.005.Index_21.MCF7_E2_ERA_ChIP8_GB_rep1.bam,MCF7_E2_BRD4_rep1.SRR1193529.fastq,MCF7_CTRL_BRD4_rep2.SRR1193527.fastq,MCF7_E2_BRD4_rep3.SRR1193531.fastq,MCF7_E2_BRD4_rep2.SRR1193530.fastq,MCF7_CTRL_SMC1A_rep1.HI.1997.004.Index_4.MCF7_CTRL_SMC1_ChIP8_GB_rep1.bam,MCF7_CTRL_SMC1A_rep2.HI.2791.003.Index_11.MCF7_4_Smc1_C_rep2.bam,MCF7_CTRL_NIPBL_rep2.HI.2791.003.Index_9.MCF7_1_NIPBL_C_rep2.bam,MCF7_CTRL_NIPBL_rep1.HI.1997.004.Index_1.MCF7_CTRL_NIPBL_ChIP8_GB_rep1.bam,MCF7_CTRL_MED1_rep1.HI.1943.008.Index_7.MCF7_CTRL_MED1_ChIP8_GB_rep1.bam,MCF7_E2_SMC1A_rep1.HI.1997.004.Index_5.MCF7_E2_SMC1_ChIP8_GB_rep1.bam,MCF7_CTRL_MED1_rep2.HI.2791.003.Index_8.MCF7_2_MED1_C_rep2.bam,MCF7_E2_SMC1A_rep2.HI.2791.003.Index_3.MCF7_4_Smc1_E_rep2.bam' | md5sum | awk '{ print $1 }')
wget "http://mugqic.hpc.mcgill.ca/cgi-bin/pipeline.cgi?hostname=cedar1.cedar.computecanada.ca&ip=206.12.124.2&pipeline=ChipSeq&steps=picard_sam_to_fastq,trimmomatic,merge_trimmomatic_stats,bwa_mem_picard_sort_sam,samtools_view_filter,picard_merge_sam_files,picard_mark_duplicates,metrics,homer_make_tag_directory,qc_metrics,macs2_callpeak&samples=30&md5=$LOG_MD5" --quiet --output-document=/dev/null

