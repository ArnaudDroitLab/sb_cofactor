#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

#-------------------------------------------------------------------------------
# ChipSeq SlurmScheduler Job Submission Bash script
# Version: 3.1.2-beta
# Created on: 2019-02-22T08:53:00
# Steps:
#   picard_sam_to_fastq: 5 jobs
#   trimmomatic: 9 jobs
#   merge_trimmomatic_stats: 1 job
#   bwa_mem_picard_sort_sam: 10 jobs
#   samtools_view_filter: 10 jobs
#   picard_merge_sam_files: 9 jobs
#   picard_mark_duplicates: 10 jobs
#   metrics: 2 jobs
#   homer_make_tag_directory: 9 jobs
#   qc_metrics: 1 job
#   macs2_callpeak: 13 jobs
#   TOTAL: 79 jobs
#-------------------------------------------------------------------------------

OUTPUT_DIR=/scratch/chris11/sb_cofactor/K562/output/chip-pipeline-GRCh38
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
# JOB: picard_sam_to_fastq_1_JOB_ID: picard_sam_to_fastq.HI.1385.007.Index_6.K562_WCE_MF_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.1385.007.Index_6.K562_WCE_MF_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.1385.007.Index_6.K562_WCE_MF_rep1.bam.813307903d3ba59216f9cca8d9f43c0e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.1385.007.Index_6.K562_WCE_MF_rep1.bam.813307903d3ba59216f9cca8d9f43c0e.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/scratch/chris11/sb_cofactor/K562/raw/chip-seq/HI.1385.007.Index_6.K562_WCE_MF_rep1.bam \
 FASTQ=/scratch/chris11/sb_cofactor/K562/raw/chip-seq/HI.1385.007.Index_6.K562_WCE_MF_rep1.single.fastq.gz
picard_sam_to_fastq.HI.1385.007.Index_6.K562_WCE_MF_rep1.bam.813307903d3ba59216f9cca8d9f43c0e.mugqic.done
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
# JOB: picard_sam_to_fastq_2_JOB_ID: picard_sam_to_fastq.HI.1385.007.Index_10.K562_CDK9_MF_ChIP2_4_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.1385.007.Index_10.K562_CDK9_MF_ChIP2_4_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.1385.007.Index_10.K562_CDK9_MF_ChIP2_4_rep1.bam.63347d9cdbf0a5c6e2db1ca35f83b098.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.1385.007.Index_10.K562_CDK9_MF_ChIP2_4_rep1.bam.63347d9cdbf0a5c6e2db1ca35f83b098.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/scratch/chris11/sb_cofactor/K562/raw/chip-seq/HI.1385.007.Index_10.K562_CDK9_MF_ChIP2_4_rep1.bam \
 FASTQ=/scratch/chris11/sb_cofactor/K562/raw/chip-seq/HI.1385.007.Index_10.K562_CDK9_MF_ChIP2_4_rep1.single.fastq.gz
picard_sam_to_fastq.HI.1385.007.Index_10.K562_CDK9_MF_ChIP2_4_rep1.bam.63347d9cdbf0a5c6e2db1ca35f83b098.mugqic.done
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
# JOB: picard_sam_to_fastq_3_JOB_ID: picard_sam_to_fastq.HI.1385.007.Index_9.K562_SMC1_MF_ChIP2_3_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.1385.007.Index_9.K562_SMC1_MF_ChIP2_3_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.1385.007.Index_9.K562_SMC1_MF_ChIP2_3_rep1.bam.56965aca5d9d608144f9587a3188bad3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.1385.007.Index_9.K562_SMC1_MF_ChIP2_3_rep1.bam.56965aca5d9d608144f9587a3188bad3.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/scratch/chris11/sb_cofactor/K562/raw/chip-seq/HI.1385.007.Index_9.K562_SMC1_MF_ChIP2_3_rep1.bam \
 FASTQ=/scratch/chris11/sb_cofactor/K562/raw/chip-seq/HI.1385.007.Index_9.K562_SMC1_MF_ChIP2_3_rep1.single.fastq.gz
picard_sam_to_fastq.HI.1385.007.Index_9.K562_SMC1_MF_ChIP2_3_rep1.bam.56965aca5d9d608144f9587a3188bad3.mugqic.done
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
# JOB: picard_sam_to_fastq_4_JOB_ID: picard_sam_to_fastq.HI.1385.007.Index_8.K562_MED1_MF_ChIP2_2_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.1385.007.Index_8.K562_MED1_MF_ChIP2_2_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.1385.007.Index_8.K562_MED1_MF_ChIP2_2_rep1.bam.e1b86aaeb7789b132ae479af8f31da6c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.1385.007.Index_8.K562_MED1_MF_ChIP2_2_rep1.bam.e1b86aaeb7789b132ae479af8f31da6c.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/scratch/chris11/sb_cofactor/K562/raw/chip-seq/HI.1385.007.Index_8.K562_MED1_MF_ChIP2_2_rep1.bam \
 FASTQ=/scratch/chris11/sb_cofactor/K562/raw/chip-seq/HI.1385.007.Index_8.K562_MED1_MF_ChIP2_2_rep1.single.fastq.gz
picard_sam_to_fastq.HI.1385.007.Index_8.K562_MED1_MF_ChIP2_2_rep1.bam.e1b86aaeb7789b132ae479af8f31da6c.mugqic.done
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
# JOB: picard_sam_to_fastq_5_JOB_ID: picard_sam_to_fastq.HI.1385.007.Index_7.K562_NIPBL_MF_ChIP2_1_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.1385.007.Index_7.K562_NIPBL_MF_ChIP2_1_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.1385.007.Index_7.K562_NIPBL_MF_ChIP2_1_rep1.bam.f9d7f9df58343cd99d291f59b98aa89f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.1385.007.Index_7.K562_NIPBL_MF_ChIP2_1_rep1.bam.f9d7f9df58343cd99d291f59b98aa89f.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/scratch/chris11/sb_cofactor/K562/raw/chip-seq/HI.1385.007.Index_7.K562_NIPBL_MF_ChIP2_1_rep1.bam \
 FASTQ=/scratch/chris11/sb_cofactor/K562/raw/chip-seq/HI.1385.007.Index_7.K562_NIPBL_MF_ChIP2_1_rep1.single.fastq.gz
picard_sam_to_fastq.HI.1385.007.Index_7.K562_NIPBL_MF_ChIP2_1_rep1.bam.f9d7f9df58343cd99d291f59b98aa89f.mugqic.done
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
# STEP: trimmomatic
#-------------------------------------------------------------------------------
STEP=trimmomatic
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: trimmomatic_1_JOB_ID: trimmomatic.HI.1385.007.Index_6.K562_WCE_MF_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.1385.007.Index_6.K562_WCE_MF_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_1_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.1385.007.Index_6.K562_WCE_MF_rep1.bam.fb3d22e043e787aeb055805a8b07ca3d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.1385.007.Index_6.K562_WCE_MF_rep1.bam.fb3d22e043e787aeb055805a8b07ca3d.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/K562_WCE_rep1 && \
`cat > trim/K562_WCE_rep1/HI.1385.007.Index_6.K562_WCE_MF_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /scratch/chris11/sb_cofactor/K562/raw/chip-seq/HI.1385.007.Index_6.K562_WCE_MF_rep1.single.fastq.gz \
  trim/K562_WCE_rep1/HI.1385.007.Index_6.K562_WCE_MF_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/K562_WCE_rep1/HI.1385.007.Index_6.K562_WCE_MF_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/K562_WCE_rep1/HI.1385.007.Index_6.K562_WCE_MF_rep1.bam.trim.log
trimmomatic.HI.1385.007.Index_6.K562_WCE_MF_rep1.bam.fb3d22e043e787aeb055805a8b07ca3d.mugqic.done
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
# JOB: trimmomatic_2_JOB_ID: trimmomatic.HI.1385.007.Index_10.K562_CDK9_MF_ChIP2_4_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.1385.007.Index_10.K562_CDK9_MF_ChIP2_4_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_2_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.1385.007.Index_10.K562_CDK9_MF_ChIP2_4_rep1.bam.329e2aa4fd9870603b7dadd5e181c7d9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.1385.007.Index_10.K562_CDK9_MF_ChIP2_4_rep1.bam.329e2aa4fd9870603b7dadd5e181c7d9.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/K562_CDK9_rep1 && \
`cat > trim/K562_CDK9_rep1/HI.1385.007.Index_10.K562_CDK9_MF_ChIP2_4_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /scratch/chris11/sb_cofactor/K562/raw/chip-seq/HI.1385.007.Index_10.K562_CDK9_MF_ChIP2_4_rep1.single.fastq.gz \
  trim/K562_CDK9_rep1/HI.1385.007.Index_10.K562_CDK9_MF_ChIP2_4_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/K562_CDK9_rep1/HI.1385.007.Index_10.K562_CDK9_MF_ChIP2_4_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/K562_CDK9_rep1/HI.1385.007.Index_10.K562_CDK9_MF_ChIP2_4_rep1.bam.trim.log
trimmomatic.HI.1385.007.Index_10.K562_CDK9_MF_ChIP2_4_rep1.bam.329e2aa4fd9870603b7dadd5e181c7d9.mugqic.done
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
# JOB: trimmomatic_3_JOB_ID: trimmomatic.HI.1385.007.Index_9.K562_SMC1_MF_ChIP2_3_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.1385.007.Index_9.K562_SMC1_MF_ChIP2_3_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_3_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.1385.007.Index_9.K562_SMC1_MF_ChIP2_3_rep1.bam.d34fa2b5f109c79273dd13db012f5b7e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.1385.007.Index_9.K562_SMC1_MF_ChIP2_3_rep1.bam.d34fa2b5f109c79273dd13db012f5b7e.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/K562_SMC1A_rep1 && \
`cat > trim/K562_SMC1A_rep1/HI.1385.007.Index_9.K562_SMC1_MF_ChIP2_3_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /scratch/chris11/sb_cofactor/K562/raw/chip-seq/HI.1385.007.Index_9.K562_SMC1_MF_ChIP2_3_rep1.single.fastq.gz \
  trim/K562_SMC1A_rep1/HI.1385.007.Index_9.K562_SMC1_MF_ChIP2_3_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/K562_SMC1A_rep1/HI.1385.007.Index_9.K562_SMC1_MF_ChIP2_3_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/K562_SMC1A_rep1/HI.1385.007.Index_9.K562_SMC1_MF_ChIP2_3_rep1.bam.trim.log
trimmomatic.HI.1385.007.Index_9.K562_SMC1_MF_ChIP2_3_rep1.bam.d34fa2b5f109c79273dd13db012f5b7e.mugqic.done
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
# JOB: trimmomatic_4_JOB_ID: trimmomatic.HI.1385.007.Index_8.K562_MED1_MF_ChIP2_2_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.1385.007.Index_8.K562_MED1_MF_ChIP2_2_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_4_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.1385.007.Index_8.K562_MED1_MF_ChIP2_2_rep1.bam.1ee76e3b2484ef55190e585fa5ca358f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.1385.007.Index_8.K562_MED1_MF_ChIP2_2_rep1.bam.1ee76e3b2484ef55190e585fa5ca358f.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/K562_MED1_rep1 && \
`cat > trim/K562_MED1_rep1/HI.1385.007.Index_8.K562_MED1_MF_ChIP2_2_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /scratch/chris11/sb_cofactor/K562/raw/chip-seq/HI.1385.007.Index_8.K562_MED1_MF_ChIP2_2_rep1.single.fastq.gz \
  trim/K562_MED1_rep1/HI.1385.007.Index_8.K562_MED1_MF_ChIP2_2_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/K562_MED1_rep1/HI.1385.007.Index_8.K562_MED1_MF_ChIP2_2_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/K562_MED1_rep1/HI.1385.007.Index_8.K562_MED1_MF_ChIP2_2_rep1.bam.trim.log
trimmomatic.HI.1385.007.Index_8.K562_MED1_MF_ChIP2_2_rep1.bam.1ee76e3b2484ef55190e585fa5ca358f.mugqic.done
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
# JOB: trimmomatic_5_JOB_ID: trimmomatic.HI.1385.007.Index_7.K562_NIPBL_MF_ChIP2_1_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.1385.007.Index_7.K562_NIPBL_MF_ChIP2_1_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_5_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.1385.007.Index_7.K562_NIPBL_MF_ChIP2_1_rep1.bam.2e47286e67b08dac6cff7fd9e058da05.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.1385.007.Index_7.K562_NIPBL_MF_ChIP2_1_rep1.bam.2e47286e67b08dac6cff7fd9e058da05.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/K562_NIPBL_rep1 && \
`cat > trim/K562_NIPBL_rep1/HI.1385.007.Index_7.K562_NIPBL_MF_ChIP2_1_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /scratch/chris11/sb_cofactor/K562/raw/chip-seq/HI.1385.007.Index_7.K562_NIPBL_MF_ChIP2_1_rep1.single.fastq.gz \
  trim/K562_NIPBL_rep1/HI.1385.007.Index_7.K562_NIPBL_MF_ChIP2_1_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/K562_NIPBL_rep1/HI.1385.007.Index_7.K562_NIPBL_MF_ChIP2_1_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/K562_NIPBL_rep1/HI.1385.007.Index_7.K562_NIPBL_MF_ChIP2_1_rep1.bam.trim.log
trimmomatic.HI.1385.007.Index_7.K562_NIPBL_MF_ChIP2_1_rep1.bam.2e47286e67b08dac6cff7fd9e058da05.mugqic.done
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
# JOB: trimmomatic_6_JOB_ID: trimmomatic.ENCFF335PHG.fastq
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.ENCFF335PHG.fastq
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.ENCFF335PHG.fastq.6e9fe717280ffdaec5eadce3b3a1909a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.ENCFF335PHG.fastq.6e9fe717280ffdaec5eadce3b3a1909a.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/K562_BRD4_rep1 && \
`cat > trim/K562_BRD4_rep1/ENCFF335PHG.fastq.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /scratch/chris11/sb_cofactor/K562/raw/chip-seq/ENCFF335PHG.fastq \
  trim/K562_BRD4_rep1/ENCFF335PHG.fastq.trim.single.fastq.gz \
  ILLUMINACLIP:trim/K562_BRD4_rep1/ENCFF335PHG.fastq.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/K562_BRD4_rep1/ENCFF335PHG.fastq.trim.log
trimmomatic.ENCFF335PHG.fastq.6e9fe717280ffdaec5eadce3b3a1909a.mugqic.done
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=24G -N 1 -n 6 | grep "[0-9]" | cut -d\  -f4)
echo "$trimmomatic_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: trimmomatic_7_JOB_ID: trimmomatic.ENCFF049ACQ.fastq
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.ENCFF049ACQ.fastq
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.ENCFF049ACQ.fastq.010ccb14fa7261eac73342cefe094780.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.ENCFF049ACQ.fastq.010ccb14fa7261eac73342cefe094780.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/K562_BRD4_rep2 && \
`cat > trim/K562_BRD4_rep2/ENCFF049ACQ.fastq.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /scratch/chris11/sb_cofactor/K562/raw/chip-seq/ENCFF049ACQ.fastq \
  trim/K562_BRD4_rep2/ENCFF049ACQ.fastq.trim.single.fastq.gz \
  ILLUMINACLIP:trim/K562_BRD4_rep2/ENCFF049ACQ.fastq.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/K562_BRD4_rep2/ENCFF049ACQ.fastq.trim.log
trimmomatic.ENCFF049ACQ.fastq.010ccb14fa7261eac73342cefe094780.mugqic.done
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=24G -N 1 -n 6 | grep "[0-9]" | cut -d\  -f4)
echo "$trimmomatic_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: trimmomatic_8_JOB_ID: trimmomatic.ENCFF000BWK.fastq
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.ENCFF000BWK.fastq
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.ENCFF000BWK.fastq.4370027caf4657c4bb387fd206a273f3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.ENCFF000BWK.fastq.4370027caf4657c4bb387fd206a273f3.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/K562_BRD4_WCE_rep1 && \
`cat > trim/K562_BRD4_WCE_rep1/ENCFF000BWK.fastq.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /scratch/chris11/sb_cofactor/K562/raw/chip-seq/ENCFF000BWK.fastq \
  trim/K562_BRD4_WCE_rep1/ENCFF000BWK.fastq.trim.single.fastq.gz \
  ILLUMINACLIP:trim/K562_BRD4_WCE_rep1/ENCFF000BWK.fastq.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/K562_BRD4_WCE_rep1/ENCFF000BWK.fastq.trim.log
trimmomatic.ENCFF000BWK.fastq.4370027caf4657c4bb387fd206a273f3.mugqic.done
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=24G -N 1 -n 6 | grep "[0-9]" | cut -d\  -f4)
echo "$trimmomatic_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# JOB: trimmomatic_9_JOB_ID: trimmomatic.ENCFF994FIB.fastq
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.ENCFF994FIB.fastq
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.ENCFF994FIB.fastq.1ec8ae04230191271e3568dfe43ddc19.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.ENCFF994FIB.fastq.1ec8ae04230191271e3568dfe43ddc19.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/K562_BRD4_WCE_rep2 && \
`cat > trim/K562_BRD4_WCE_rep2/ENCFF994FIB.fastq.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /scratch/chris11/sb_cofactor/K562/raw/chip-seq/ENCFF994FIB.fastq \
  trim/K562_BRD4_WCE_rep2/ENCFF994FIB.fastq.trim.single.fastq.gz \
  ILLUMINACLIP:trim/K562_BRD4_WCE_rep2/ENCFF994FIB.fastq.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/K562_BRD4_WCE_rep2/ENCFF994FIB.fastq.trim.log
trimmomatic.ENCFF994FIB.fastq.1ec8ae04230191271e3568dfe43ddc19.mugqic.done
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=24G -N 1 -n 6 | grep "[0-9]" | cut -d\  -f4)
echo "$trimmomatic_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

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
JOB_DEPENDENCIES=$trimmomatic_1_JOB_ID:$trimmomatic_2_JOB_ID:$trimmomatic_3_JOB_ID:$trimmomatic_4_JOB_ID:$trimmomatic_5_JOB_ID:$trimmomatic_6_JOB_ID:$trimmomatic_7_JOB_ID:$trimmomatic_8_JOB_ID:$trimmomatic_9_JOB_ID
JOB_DONE=job_output/merge_trimmomatic_stats/merge_trimmomatic_stats.d97c34b5cda21ab77cb7af33083dfffd.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'merge_trimmomatic_stats.d97c34b5cda21ab77cb7af33083dfffd.mugqic.done'
module load mugqic/pandoc/1.15.2 && \
mkdir -p metrics && \
echo 'Sample	Readset	Raw Single Reads #	Surviving Single Reads #	Surviving Single Reads %' > metrics/trimReadsetTable.tsv && \
grep ^Input trim/K562_WCE_rep1/HI.1385.007.Index_6.K562_WCE_MF_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/K562_WCE_rep1	HI.1385.007.Index_6.K562_WCE_MF_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/K562_CDK9_rep1/HI.1385.007.Index_10.K562_CDK9_MF_ChIP2_4_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/K562_CDK9_rep1	HI.1385.007.Index_10.K562_CDK9_MF_ChIP2_4_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/K562_SMC1A_rep1/HI.1385.007.Index_9.K562_SMC1_MF_ChIP2_3_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/K562_SMC1A_rep1	HI.1385.007.Index_9.K562_SMC1_MF_ChIP2_3_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/K562_MED1_rep1/HI.1385.007.Index_8.K562_MED1_MF_ChIP2_2_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/K562_MED1_rep1	HI.1385.007.Index_8.K562_MED1_MF_ChIP2_2_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/K562_NIPBL_rep1/HI.1385.007.Index_7.K562_NIPBL_MF_ChIP2_1_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/K562_NIPBL_rep1	HI.1385.007.Index_7.K562_NIPBL_MF_ChIP2_1_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/K562_BRD4_rep1/ENCFF335PHG.fastq.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/K562_BRD4_rep1	ENCFF335PHG.fastq	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/K562_BRD4_rep2/ENCFF049ACQ.fastq.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/K562_BRD4_rep2	ENCFF049ACQ.fastq	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/K562_BRD4_WCE_rep1/ENCFF000BWK.fastq.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/K562_BRD4_WCE_rep1	ENCFF000BWK.fastq	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/K562_BRD4_WCE_rep2/ENCFF994FIB.fastq.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/K562_BRD4_WCE_rep2	ENCFF994FIB.fastq	\1	\2/' | \
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
merge_trimmomatic_stats.d97c34b5cda21ab77cb7af33083dfffd.mugqic.done
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
# JOB: bwa_mem_picard_sort_sam_1_JOB_ID: bwa_mem_picard_sort_sam.HI.1385.007.Index_6.K562_WCE_MF_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.1385.007.Index_6.K562_WCE_MF_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_1_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.1385.007.Index_6.K562_WCE_MF_rep1.bam.f55013b02e1f2507a4dd224c03ddc8bf.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.1385.007.Index_6.K562_WCE_MF_rep1.bam.f55013b02e1f2507a4dd224c03ddc8bf.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p alignment/K562_WCE_rep1/HI.1385.007.Index_6.K562_WCE_MF_rep1.bam && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:HI.1385.007.Index_6.K562_WCE_MF_rep1.bam	SM:K562_WCE_rep1	LB:K562_WCE_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/K562_WCE_rep1/HI.1385.007.Index_6.K562_WCE_MF_rep1.bam.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/K562_WCE_rep1/HI.1385.007.Index_6.K562_WCE_MF_rep1.bam/HI.1385.007.Index_6.K562_WCE_MF_rep1.bam.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.HI.1385.007.Index_6.K562_WCE_MF_rep1.bam.f55013b02e1f2507a4dd224c03ddc8bf.mugqic.done
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
# JOB: bwa_mem_picard_sort_sam_2_JOB_ID: bwa_mem_picard_sort_sam.HI.1385.007.Index_10.K562_CDK9_MF_ChIP2_4_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.1385.007.Index_10.K562_CDK9_MF_ChIP2_4_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_2_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.1385.007.Index_10.K562_CDK9_MF_ChIP2_4_rep1.bam.20c2ed73411d2290c6d300b23ad65151.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.1385.007.Index_10.K562_CDK9_MF_ChIP2_4_rep1.bam.20c2ed73411d2290c6d300b23ad65151.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p alignment/K562_CDK9_rep1/HI.1385.007.Index_10.K562_CDK9_MF_ChIP2_4_rep1.bam && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:HI.1385.007.Index_10.K562_CDK9_MF_ChIP2_4_rep1.bam	SM:K562_CDK9_rep1	LB:K562_CDK9_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/K562_CDK9_rep1/HI.1385.007.Index_10.K562_CDK9_MF_ChIP2_4_rep1.bam.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/K562_CDK9_rep1/HI.1385.007.Index_10.K562_CDK9_MF_ChIP2_4_rep1.bam/HI.1385.007.Index_10.K562_CDK9_MF_ChIP2_4_rep1.bam.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.HI.1385.007.Index_10.K562_CDK9_MF_ChIP2_4_rep1.bam.20c2ed73411d2290c6d300b23ad65151.mugqic.done
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
# JOB: bwa_mem_picard_sort_sam_3_JOB_ID: bwa_mem_picard_sort_sam.HI.1385.007.Index_9.K562_SMC1_MF_ChIP2_3_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.1385.007.Index_9.K562_SMC1_MF_ChIP2_3_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_3_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.1385.007.Index_9.K562_SMC1_MF_ChIP2_3_rep1.bam.3211f4dac2f6cc9102fac6de7c21973b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.1385.007.Index_9.K562_SMC1_MF_ChIP2_3_rep1.bam.3211f4dac2f6cc9102fac6de7c21973b.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p alignment/K562_SMC1A_rep1/HI.1385.007.Index_9.K562_SMC1_MF_ChIP2_3_rep1.bam && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:HI.1385.007.Index_9.K562_SMC1_MF_ChIP2_3_rep1.bam	SM:K562_SMC1A_rep1	LB:K562_SMC1A_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/K562_SMC1A_rep1/HI.1385.007.Index_9.K562_SMC1_MF_ChIP2_3_rep1.bam.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/K562_SMC1A_rep1/HI.1385.007.Index_9.K562_SMC1_MF_ChIP2_3_rep1.bam/HI.1385.007.Index_9.K562_SMC1_MF_ChIP2_3_rep1.bam.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.HI.1385.007.Index_9.K562_SMC1_MF_ChIP2_3_rep1.bam.3211f4dac2f6cc9102fac6de7c21973b.mugqic.done
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
# JOB: bwa_mem_picard_sort_sam_4_JOB_ID: bwa_mem_picard_sort_sam.HI.1385.007.Index_8.K562_MED1_MF_ChIP2_2_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.1385.007.Index_8.K562_MED1_MF_ChIP2_2_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_4_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.1385.007.Index_8.K562_MED1_MF_ChIP2_2_rep1.bam.3a0d151e71a166466c6214da4ee1edbe.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.1385.007.Index_8.K562_MED1_MF_ChIP2_2_rep1.bam.3a0d151e71a166466c6214da4ee1edbe.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p alignment/K562_MED1_rep1/HI.1385.007.Index_8.K562_MED1_MF_ChIP2_2_rep1.bam && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:HI.1385.007.Index_8.K562_MED1_MF_ChIP2_2_rep1.bam	SM:K562_MED1_rep1	LB:K562_MED1_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/K562_MED1_rep1/HI.1385.007.Index_8.K562_MED1_MF_ChIP2_2_rep1.bam.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/K562_MED1_rep1/HI.1385.007.Index_8.K562_MED1_MF_ChIP2_2_rep1.bam/HI.1385.007.Index_8.K562_MED1_MF_ChIP2_2_rep1.bam.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.HI.1385.007.Index_8.K562_MED1_MF_ChIP2_2_rep1.bam.3a0d151e71a166466c6214da4ee1edbe.mugqic.done
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
# JOB: bwa_mem_picard_sort_sam_5_JOB_ID: bwa_mem_picard_sort_sam.HI.1385.007.Index_7.K562_NIPBL_MF_ChIP2_1_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.1385.007.Index_7.K562_NIPBL_MF_ChIP2_1_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_5_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.1385.007.Index_7.K562_NIPBL_MF_ChIP2_1_rep1.bam.ad13d60c38c9d91cc9602760f25f0e03.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.1385.007.Index_7.K562_NIPBL_MF_ChIP2_1_rep1.bam.ad13d60c38c9d91cc9602760f25f0e03.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p alignment/K562_NIPBL_rep1/HI.1385.007.Index_7.K562_NIPBL_MF_ChIP2_1_rep1.bam && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:HI.1385.007.Index_7.K562_NIPBL_MF_ChIP2_1_rep1.bam	SM:K562_NIPBL_rep1	LB:K562_NIPBL_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/K562_NIPBL_rep1/HI.1385.007.Index_7.K562_NIPBL_MF_ChIP2_1_rep1.bam.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/K562_NIPBL_rep1/HI.1385.007.Index_7.K562_NIPBL_MF_ChIP2_1_rep1.bam/HI.1385.007.Index_7.K562_NIPBL_MF_ChIP2_1_rep1.bam.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.HI.1385.007.Index_7.K562_NIPBL_MF_ChIP2_1_rep1.bam.ad13d60c38c9d91cc9602760f25f0e03.mugqic.done
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
# JOB: bwa_mem_picard_sort_sam_6_JOB_ID: bwa_mem_picard_sort_sam.ENCFF335PHG.fastq
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.ENCFF335PHG.fastq
JOB_DEPENDENCIES=$trimmomatic_6_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.ENCFF335PHG.fastq.a7c4c89e24e82ad861b77e744cdbf52a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.ENCFF335PHG.fastq.a7c4c89e24e82ad861b77e744cdbf52a.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p alignment/K562_BRD4_rep1/ENCFF335PHG.fastq && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:ENCFF335PHG.fastq	SM:K562_BRD4_rep1	LB:K562_BRD4_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/K562_BRD4_rep1/ENCFF335PHG.fastq.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/K562_BRD4_rep1/ENCFF335PHG.fastq/ENCFF335PHG.fastq.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.ENCFF335PHG.fastq.a7c4c89e24e82ad861b77e744cdbf52a.mugqic.done
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
# JOB: bwa_mem_picard_sort_sam_7_JOB_ID: bwa_mem_picard_sort_sam.ENCFF049ACQ.fastq
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.ENCFF049ACQ.fastq
JOB_DEPENDENCIES=$trimmomatic_7_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.ENCFF049ACQ.fastq.61c3e6979b045a9abb2571d37671e82f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.ENCFF049ACQ.fastq.61c3e6979b045a9abb2571d37671e82f.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p alignment/K562_BRD4_rep2/ENCFF049ACQ.fastq && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:ENCFF049ACQ.fastq	SM:K562_BRD4_rep2	LB:K562_BRD4_rep2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/K562_BRD4_rep2/ENCFF049ACQ.fastq.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/K562_BRD4_rep2/ENCFF049ACQ.fastq/ENCFF049ACQ.fastq.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.ENCFF049ACQ.fastq.61c3e6979b045a9abb2571d37671e82f.mugqic.done
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
# JOB: bwa_mem_picard_sort_sam_8_JOB_ID: bwa_mem_picard_sort_sam.ENCFF000BWK.fastq
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.ENCFF000BWK.fastq
JOB_DEPENDENCIES=$trimmomatic_8_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.ENCFF000BWK.fastq.0a36023e9f275e5b30a8c410b8d20438.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.ENCFF000BWK.fastq.0a36023e9f275e5b30a8c410b8d20438.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p alignment/K562_BRD4_WCE_rep1/ENCFF000BWK.fastq && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:ENCFF000BWK.fastq	SM:K562_BRD4_WCE_rep1	LB:K562_BRD4_WCE_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/K562_BRD4_WCE_rep1/ENCFF000BWK.fastq.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/K562_BRD4_WCE_rep1/ENCFF000BWK.fastq/ENCFF000BWK.fastq.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.ENCFF000BWK.fastq.0a36023e9f275e5b30a8c410b8d20438.mugqic.done
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
# JOB: bwa_mem_picard_sort_sam_9_JOB_ID: bwa_mem_picard_sort_sam.ENCFF994FIB.fastq
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.ENCFF994FIB.fastq
JOB_DEPENDENCIES=$trimmomatic_9_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.ENCFF994FIB.fastq.28cab22f251e9c249a1201b6d19cf0fc.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.ENCFF994FIB.fastq.28cab22f251e9c249a1201b6d19cf0fc.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p alignment/K562_BRD4_WCE_rep2/ENCFF994FIB.fastq && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:ENCFF994FIB.fastq	SM:K562_BRD4_WCE_rep2	LB:K562_BRD4_WCE_rep2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/K562_BRD4_WCE_rep2/ENCFF994FIB.fastq.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/K562_BRD4_WCE_rep2/ENCFF994FIB.fastq/ENCFF994FIB.fastq.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.ENCFF994FIB.fastq.28cab22f251e9c249a1201b6d19cf0fc.mugqic.done
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
# JOB: bwa_mem_picard_sort_sam_10_JOB_ID: bwa_mem_picard_sort_sam_report
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam_report
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_1_JOB_ID:$bwa_mem_picard_sort_sam_2_JOB_ID:$bwa_mem_picard_sort_sam_3_JOB_ID:$bwa_mem_picard_sort_sam_4_JOB_ID:$bwa_mem_picard_sort_sam_5_JOB_ID:$bwa_mem_picard_sort_sam_6_JOB_ID:$bwa_mem_picard_sort_sam_7_JOB_ID:$bwa_mem_picard_sort_sam_8_JOB_ID:$bwa_mem_picard_sort_sam_9_JOB_ID
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$bwa_mem_picard_sort_sam_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.2


#-------------------------------------------------------------------------------
# STEP: samtools_view_filter
#-------------------------------------------------------------------------------
STEP=samtools_view_filter
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_1_JOB_ID: samtools_view_filter.HI.1385.007.Index_6.K562_WCE_MF_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.1385.007.Index_6.K562_WCE_MF_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_1_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.1385.007.Index_6.K562_WCE_MF_rep1.bam.b3d7974778634b8368d1eba27772afb8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.1385.007.Index_6.K562_WCE_MF_rep1.bam.b3d7974778634b8368d1eba27772afb8.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/K562_WCE_rep1/HI.1385.007.Index_6.K562_WCE_MF_rep1.bam/HI.1385.007.Index_6.K562_WCE_MF_rep1.bam.sorted.bam \
  > alignment/K562_WCE_rep1/HI.1385.007.Index_6.K562_WCE_MF_rep1.bam/HI.1385.007.Index_6.K562_WCE_MF_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.1385.007.Index_6.K562_WCE_MF_rep1.bam.b3d7974778634b8368d1eba27772afb8.mugqic.done
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
# JOB: samtools_view_filter_2_JOB_ID: samtools_view_filter.HI.1385.007.Index_10.K562_CDK9_MF_ChIP2_4_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.1385.007.Index_10.K562_CDK9_MF_ChIP2_4_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_2_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.1385.007.Index_10.K562_CDK9_MF_ChIP2_4_rep1.bam.af65fbb519e183bbfdbf10e666080166.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.1385.007.Index_10.K562_CDK9_MF_ChIP2_4_rep1.bam.af65fbb519e183bbfdbf10e666080166.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/K562_CDK9_rep1/HI.1385.007.Index_10.K562_CDK9_MF_ChIP2_4_rep1.bam/HI.1385.007.Index_10.K562_CDK9_MF_ChIP2_4_rep1.bam.sorted.bam \
  > alignment/K562_CDK9_rep1/HI.1385.007.Index_10.K562_CDK9_MF_ChIP2_4_rep1.bam/HI.1385.007.Index_10.K562_CDK9_MF_ChIP2_4_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.1385.007.Index_10.K562_CDK9_MF_ChIP2_4_rep1.bam.af65fbb519e183bbfdbf10e666080166.mugqic.done
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
# JOB: samtools_view_filter_3_JOB_ID: samtools_view_filter.HI.1385.007.Index_9.K562_SMC1_MF_ChIP2_3_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.1385.007.Index_9.K562_SMC1_MF_ChIP2_3_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_3_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.1385.007.Index_9.K562_SMC1_MF_ChIP2_3_rep1.bam.f1a38386ff4a6b1c51d25abdfd739506.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.1385.007.Index_9.K562_SMC1_MF_ChIP2_3_rep1.bam.f1a38386ff4a6b1c51d25abdfd739506.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/K562_SMC1A_rep1/HI.1385.007.Index_9.K562_SMC1_MF_ChIP2_3_rep1.bam/HI.1385.007.Index_9.K562_SMC1_MF_ChIP2_3_rep1.bam.sorted.bam \
  > alignment/K562_SMC1A_rep1/HI.1385.007.Index_9.K562_SMC1_MF_ChIP2_3_rep1.bam/HI.1385.007.Index_9.K562_SMC1_MF_ChIP2_3_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.1385.007.Index_9.K562_SMC1_MF_ChIP2_3_rep1.bam.f1a38386ff4a6b1c51d25abdfd739506.mugqic.done
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
# JOB: samtools_view_filter_4_JOB_ID: samtools_view_filter.HI.1385.007.Index_8.K562_MED1_MF_ChIP2_2_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.1385.007.Index_8.K562_MED1_MF_ChIP2_2_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_4_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.1385.007.Index_8.K562_MED1_MF_ChIP2_2_rep1.bam.63561c1c1285d6ad38d79af1141cb7ce.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.1385.007.Index_8.K562_MED1_MF_ChIP2_2_rep1.bam.63561c1c1285d6ad38d79af1141cb7ce.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/K562_MED1_rep1/HI.1385.007.Index_8.K562_MED1_MF_ChIP2_2_rep1.bam/HI.1385.007.Index_8.K562_MED1_MF_ChIP2_2_rep1.bam.sorted.bam \
  > alignment/K562_MED1_rep1/HI.1385.007.Index_8.K562_MED1_MF_ChIP2_2_rep1.bam/HI.1385.007.Index_8.K562_MED1_MF_ChIP2_2_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.1385.007.Index_8.K562_MED1_MF_ChIP2_2_rep1.bam.63561c1c1285d6ad38d79af1141cb7ce.mugqic.done
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
# JOB: samtools_view_filter_5_JOB_ID: samtools_view_filter.HI.1385.007.Index_7.K562_NIPBL_MF_ChIP2_1_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.1385.007.Index_7.K562_NIPBL_MF_ChIP2_1_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_5_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.1385.007.Index_7.K562_NIPBL_MF_ChIP2_1_rep1.bam.097842410b5ee61bb3f1e18fefef1c2b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.1385.007.Index_7.K562_NIPBL_MF_ChIP2_1_rep1.bam.097842410b5ee61bb3f1e18fefef1c2b.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/K562_NIPBL_rep1/HI.1385.007.Index_7.K562_NIPBL_MF_ChIP2_1_rep1.bam/HI.1385.007.Index_7.K562_NIPBL_MF_ChIP2_1_rep1.bam.sorted.bam \
  > alignment/K562_NIPBL_rep1/HI.1385.007.Index_7.K562_NIPBL_MF_ChIP2_1_rep1.bam/HI.1385.007.Index_7.K562_NIPBL_MF_ChIP2_1_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.1385.007.Index_7.K562_NIPBL_MF_ChIP2_1_rep1.bam.097842410b5ee61bb3f1e18fefef1c2b.mugqic.done
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
# JOB: samtools_view_filter_6_JOB_ID: samtools_view_filter.ENCFF335PHG.fastq
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.ENCFF335PHG.fastq
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_6_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.ENCFF335PHG.fastq.926c7bfb42adde8deed26d257d1c91cc.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.ENCFF335PHG.fastq.926c7bfb42adde8deed26d257d1c91cc.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/K562_BRD4_rep1/ENCFF335PHG.fastq/ENCFF335PHG.fastq.sorted.bam \
  > alignment/K562_BRD4_rep1/ENCFF335PHG.fastq/ENCFF335PHG.fastq.sorted.filtered.bam
samtools_view_filter.ENCFF335PHG.fastq.926c7bfb42adde8deed26d257d1c91cc.mugqic.done
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
# JOB: samtools_view_filter_7_JOB_ID: samtools_view_filter.ENCFF049ACQ.fastq
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.ENCFF049ACQ.fastq
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_7_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.ENCFF049ACQ.fastq.c7c1c80e13efe21ede3b8a367cbf337c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.ENCFF049ACQ.fastq.c7c1c80e13efe21ede3b8a367cbf337c.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/K562_BRD4_rep2/ENCFF049ACQ.fastq/ENCFF049ACQ.fastq.sorted.bam \
  > alignment/K562_BRD4_rep2/ENCFF049ACQ.fastq/ENCFF049ACQ.fastq.sorted.filtered.bam
samtools_view_filter.ENCFF049ACQ.fastq.c7c1c80e13efe21ede3b8a367cbf337c.mugqic.done
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
# JOB: samtools_view_filter_8_JOB_ID: samtools_view_filter.ENCFF000BWK.fastq
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.ENCFF000BWK.fastq
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_8_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.ENCFF000BWK.fastq.8199a90ca9dbbb6107866b82383fde08.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.ENCFF000BWK.fastq.8199a90ca9dbbb6107866b82383fde08.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/K562_BRD4_WCE_rep1/ENCFF000BWK.fastq/ENCFF000BWK.fastq.sorted.bam \
  > alignment/K562_BRD4_WCE_rep1/ENCFF000BWK.fastq/ENCFF000BWK.fastq.sorted.filtered.bam
samtools_view_filter.ENCFF000BWK.fastq.8199a90ca9dbbb6107866b82383fde08.mugqic.done
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
# JOB: samtools_view_filter_9_JOB_ID: samtools_view_filter.ENCFF994FIB.fastq
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.ENCFF994FIB.fastq
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_9_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.ENCFF994FIB.fastq.aa4a463cdca6d932ed5c5f09e8e8cb8c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.ENCFF994FIB.fastq.aa4a463cdca6d932ed5c5f09e8e8cb8c.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/K562_BRD4_WCE_rep2/ENCFF994FIB.fastq/ENCFF994FIB.fastq.sorted.bam \
  > alignment/K562_BRD4_WCE_rep2/ENCFF994FIB.fastq/ENCFF994FIB.fastq.sorted.filtered.bam
samtools_view_filter.ENCFF994FIB.fastq.aa4a463cdca6d932ed5c5f09e8e8cb8c.mugqic.done
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
# JOB: samtools_view_filter_10_JOB_ID: samtools_view_filter_report
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter_report
JOB_DEPENDENCIES=$samtools_view_filter_1_JOB_ID:$samtools_view_filter_2_JOB_ID:$samtools_view_filter_3_JOB_ID:$samtools_view_filter_4_JOB_ID:$samtools_view_filter_5_JOB_ID:$samtools_view_filter_6_JOB_ID:$samtools_view_filter_7_JOB_ID:$samtools_view_filter_8_JOB_ID:$samtools_view_filter_9_JOB_ID
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
# STEP: picard_merge_sam_files
#-------------------------------------------------------------------------------
STEP=picard_merge_sam_files
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_1_JOB_ID: symlink_readset_sample_bam.K562_WCE_rep1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.K562_WCE_rep1
JOB_DEPENDENCIES=$samtools_view_filter_1_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.K562_WCE_rep1.bf960fbd6642191c0a60adbb0ac9a464.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.K562_WCE_rep1.bf960fbd6642191c0a60adbb0ac9a464.mugqic.done'
mkdir -p alignment/K562_WCE_rep1 && \
ln -s -f HI.1385.007.Index_6.K562_WCE_MF_rep1.bam/HI.1385.007.Index_6.K562_WCE_MF_rep1.bam.sorted.filtered.bam alignment/K562_WCE_rep1/K562_WCE_rep1.merged.bam
symlink_readset_sample_bam.K562_WCE_rep1.bf960fbd6642191c0a60adbb0ac9a464.mugqic.done
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
# JOB: picard_merge_sam_files_2_JOB_ID: symlink_readset_sample_bam.K562_CDK9_rep1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.K562_CDK9_rep1
JOB_DEPENDENCIES=$samtools_view_filter_2_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.K562_CDK9_rep1.8b451230905602674b8ade37821516ac.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.K562_CDK9_rep1.8b451230905602674b8ade37821516ac.mugqic.done'
mkdir -p alignment/K562_CDK9_rep1 && \
ln -s -f HI.1385.007.Index_10.K562_CDK9_MF_ChIP2_4_rep1.bam/HI.1385.007.Index_10.K562_CDK9_MF_ChIP2_4_rep1.bam.sorted.filtered.bam alignment/K562_CDK9_rep1/K562_CDK9_rep1.merged.bam
symlink_readset_sample_bam.K562_CDK9_rep1.8b451230905602674b8ade37821516ac.mugqic.done
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
# JOB: picard_merge_sam_files_3_JOB_ID: symlink_readset_sample_bam.K562_SMC1A_rep1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.K562_SMC1A_rep1
JOB_DEPENDENCIES=$samtools_view_filter_3_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.K562_SMC1A_rep1.2db377279a12c7ee3aa5c9c7b616d763.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.K562_SMC1A_rep1.2db377279a12c7ee3aa5c9c7b616d763.mugqic.done'
mkdir -p alignment/K562_SMC1A_rep1 && \
ln -s -f HI.1385.007.Index_9.K562_SMC1_MF_ChIP2_3_rep1.bam/HI.1385.007.Index_9.K562_SMC1_MF_ChIP2_3_rep1.bam.sorted.filtered.bam alignment/K562_SMC1A_rep1/K562_SMC1A_rep1.merged.bam
symlink_readset_sample_bam.K562_SMC1A_rep1.2db377279a12c7ee3aa5c9c7b616d763.mugqic.done
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
# JOB: picard_merge_sam_files_4_JOB_ID: symlink_readset_sample_bam.K562_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.K562_MED1_rep1
JOB_DEPENDENCIES=$samtools_view_filter_4_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.K562_MED1_rep1.a03500a8b995324b9c70cbb2f4dafd23.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.K562_MED1_rep1.a03500a8b995324b9c70cbb2f4dafd23.mugqic.done'
mkdir -p alignment/K562_MED1_rep1 && \
ln -s -f HI.1385.007.Index_8.K562_MED1_MF_ChIP2_2_rep1.bam/HI.1385.007.Index_8.K562_MED1_MF_ChIP2_2_rep1.bam.sorted.filtered.bam alignment/K562_MED1_rep1/K562_MED1_rep1.merged.bam
symlink_readset_sample_bam.K562_MED1_rep1.a03500a8b995324b9c70cbb2f4dafd23.mugqic.done
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
# JOB: picard_merge_sam_files_5_JOB_ID: symlink_readset_sample_bam.K562_NIPBL_rep1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.K562_NIPBL_rep1
JOB_DEPENDENCIES=$samtools_view_filter_5_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.K562_NIPBL_rep1.37a1a566a796c874ec13a0d1ec4f9c1f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.K562_NIPBL_rep1.37a1a566a796c874ec13a0d1ec4f9c1f.mugqic.done'
mkdir -p alignment/K562_NIPBL_rep1 && \
ln -s -f HI.1385.007.Index_7.K562_NIPBL_MF_ChIP2_1_rep1.bam/HI.1385.007.Index_7.K562_NIPBL_MF_ChIP2_1_rep1.bam.sorted.filtered.bam alignment/K562_NIPBL_rep1/K562_NIPBL_rep1.merged.bam
symlink_readset_sample_bam.K562_NIPBL_rep1.37a1a566a796c874ec13a0d1ec4f9c1f.mugqic.done
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
# JOB: picard_merge_sam_files_6_JOB_ID: symlink_readset_sample_bam.K562_BRD4_rep1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.K562_BRD4_rep1
JOB_DEPENDENCIES=$samtools_view_filter_6_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.K562_BRD4_rep1.5a192cca609ac9a9c7931fe2a4147e17.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.K562_BRD4_rep1.5a192cca609ac9a9c7931fe2a4147e17.mugqic.done'
mkdir -p alignment/K562_BRD4_rep1 && \
ln -s -f ENCFF335PHG.fastq/ENCFF335PHG.fastq.sorted.filtered.bam alignment/K562_BRD4_rep1/K562_BRD4_rep1.merged.bam
symlink_readset_sample_bam.K562_BRD4_rep1.5a192cca609ac9a9c7931fe2a4147e17.mugqic.done
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
# JOB: picard_merge_sam_files_7_JOB_ID: symlink_readset_sample_bam.K562_BRD4_rep2
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.K562_BRD4_rep2
JOB_DEPENDENCIES=$samtools_view_filter_7_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.K562_BRD4_rep2.5fdca5aad401250a06ecaf584b26fd86.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.K562_BRD4_rep2.5fdca5aad401250a06ecaf584b26fd86.mugqic.done'
mkdir -p alignment/K562_BRD4_rep2 && \
ln -s -f ENCFF049ACQ.fastq/ENCFF049ACQ.fastq.sorted.filtered.bam alignment/K562_BRD4_rep2/K562_BRD4_rep2.merged.bam
symlink_readset_sample_bam.K562_BRD4_rep2.5fdca5aad401250a06ecaf584b26fd86.mugqic.done
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
# JOB: picard_merge_sam_files_8_JOB_ID: symlink_readset_sample_bam.K562_BRD4_WCE_rep1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.K562_BRD4_WCE_rep1
JOB_DEPENDENCIES=$samtools_view_filter_8_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.K562_BRD4_WCE_rep1.bb9347eae08560bdb64ded97d9e40b5a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.K562_BRD4_WCE_rep1.bb9347eae08560bdb64ded97d9e40b5a.mugqic.done'
mkdir -p alignment/K562_BRD4_WCE_rep1 && \
ln -s -f ENCFF000BWK.fastq/ENCFF000BWK.fastq.sorted.filtered.bam alignment/K562_BRD4_WCE_rep1/K562_BRD4_WCE_rep1.merged.bam
symlink_readset_sample_bam.K562_BRD4_WCE_rep1.bb9347eae08560bdb64ded97d9e40b5a.mugqic.done
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
# JOB: picard_merge_sam_files_9_JOB_ID: symlink_readset_sample_bam.K562_BRD4_WCE_rep2
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.K562_BRD4_WCE_rep2
JOB_DEPENDENCIES=$samtools_view_filter_9_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.K562_BRD4_WCE_rep2.66463fd971f1d234a51db31736bba6c5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.K562_BRD4_WCE_rep2.66463fd971f1d234a51db31736bba6c5.mugqic.done'
mkdir -p alignment/K562_BRD4_WCE_rep2 && \
ln -s -f ENCFF994FIB.fastq/ENCFF994FIB.fastq.sorted.filtered.bam alignment/K562_BRD4_WCE_rep2/K562_BRD4_WCE_rep2.merged.bam
symlink_readset_sample_bam.K562_BRD4_WCE_rep2.66463fd971f1d234a51db31736bba6c5.mugqic.done
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
# STEP: picard_mark_duplicates
#-------------------------------------------------------------------------------
STEP=picard_mark_duplicates
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_1_JOB_ID: picard_mark_duplicates.K562_WCE_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.K562_WCE_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_1_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.K562_WCE_rep1.4bbc5469c055fb6167dd605ad8e3c7c4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.K562_WCE_rep1.4bbc5469c055fb6167dd605ad8e3c7c4.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/K562_WCE_rep1/K562_WCE_rep1.merged.bam \
 OUTPUT=alignment/K562_WCE_rep1/K562_WCE_rep1.sorted.dup.bam \
 METRICS_FILE=alignment/K562_WCE_rep1/K562_WCE_rep1.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.K562_WCE_rep1.4bbc5469c055fb6167dd605ad8e3c7c4.mugqic.done
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
# JOB: picard_mark_duplicates_2_JOB_ID: picard_mark_duplicates.K562_CDK9_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.K562_CDK9_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_2_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.K562_CDK9_rep1.34ace401ec7c0b9738656a42a6e18c73.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.K562_CDK9_rep1.34ace401ec7c0b9738656a42a6e18c73.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/K562_CDK9_rep1/K562_CDK9_rep1.merged.bam \
 OUTPUT=alignment/K562_CDK9_rep1/K562_CDK9_rep1.sorted.dup.bam \
 METRICS_FILE=alignment/K562_CDK9_rep1/K562_CDK9_rep1.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.K562_CDK9_rep1.34ace401ec7c0b9738656a42a6e18c73.mugqic.done
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
# JOB: picard_mark_duplicates_3_JOB_ID: picard_mark_duplicates.K562_SMC1A_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.K562_SMC1A_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_3_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.K562_SMC1A_rep1.46b0455ffd2a63a4d8cd4fa04ddcbc0d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.K562_SMC1A_rep1.46b0455ffd2a63a4d8cd4fa04ddcbc0d.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/K562_SMC1A_rep1/K562_SMC1A_rep1.merged.bam \
 OUTPUT=alignment/K562_SMC1A_rep1/K562_SMC1A_rep1.sorted.dup.bam \
 METRICS_FILE=alignment/K562_SMC1A_rep1/K562_SMC1A_rep1.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.K562_SMC1A_rep1.46b0455ffd2a63a4d8cd4fa04ddcbc0d.mugqic.done
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
# JOB: picard_mark_duplicates_4_JOB_ID: picard_mark_duplicates.K562_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.K562_MED1_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_4_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.K562_MED1_rep1.21a65602d26a0382819e1c7fbd3f67e8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.K562_MED1_rep1.21a65602d26a0382819e1c7fbd3f67e8.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/K562_MED1_rep1/K562_MED1_rep1.merged.bam \
 OUTPUT=alignment/K562_MED1_rep1/K562_MED1_rep1.sorted.dup.bam \
 METRICS_FILE=alignment/K562_MED1_rep1/K562_MED1_rep1.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.K562_MED1_rep1.21a65602d26a0382819e1c7fbd3f67e8.mugqic.done
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
# JOB: picard_mark_duplicates_5_JOB_ID: picard_mark_duplicates.K562_NIPBL_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.K562_NIPBL_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_5_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.K562_NIPBL_rep1.ee68bd8cc1bd5e6b6f4e3cfcd5af4fae.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.K562_NIPBL_rep1.ee68bd8cc1bd5e6b6f4e3cfcd5af4fae.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/K562_NIPBL_rep1/K562_NIPBL_rep1.merged.bam \
 OUTPUT=alignment/K562_NIPBL_rep1/K562_NIPBL_rep1.sorted.dup.bam \
 METRICS_FILE=alignment/K562_NIPBL_rep1/K562_NIPBL_rep1.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.K562_NIPBL_rep1.ee68bd8cc1bd5e6b6f4e3cfcd5af4fae.mugqic.done
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
# JOB: picard_mark_duplicates_6_JOB_ID: picard_mark_duplicates.K562_BRD4_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.K562_BRD4_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_6_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.K562_BRD4_rep1.23361ce05a3f3d6633123a60f24c3843.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.K562_BRD4_rep1.23361ce05a3f3d6633123a60f24c3843.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/K562_BRD4_rep1/K562_BRD4_rep1.merged.bam \
 OUTPUT=alignment/K562_BRD4_rep1/K562_BRD4_rep1.sorted.dup.bam \
 METRICS_FILE=alignment/K562_BRD4_rep1/K562_BRD4_rep1.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.K562_BRD4_rep1.23361ce05a3f3d6633123a60f24c3843.mugqic.done
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
# JOB: picard_mark_duplicates_7_JOB_ID: picard_mark_duplicates.K562_BRD4_rep2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.K562_BRD4_rep2
JOB_DEPENDENCIES=$picard_merge_sam_files_7_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.K562_BRD4_rep2.20b8c128290a0c9e71d4d3430264b929.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.K562_BRD4_rep2.20b8c128290a0c9e71d4d3430264b929.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/K562_BRD4_rep2/K562_BRD4_rep2.merged.bam \
 OUTPUT=alignment/K562_BRD4_rep2/K562_BRD4_rep2.sorted.dup.bam \
 METRICS_FILE=alignment/K562_BRD4_rep2/K562_BRD4_rep2.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.K562_BRD4_rep2.20b8c128290a0c9e71d4d3430264b929.mugqic.done
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
# JOB: picard_mark_duplicates_8_JOB_ID: picard_mark_duplicates.K562_BRD4_WCE_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.K562_BRD4_WCE_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_8_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.K562_BRD4_WCE_rep1.5a88bc43459adc2cc757c01a5aa787b9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.K562_BRD4_WCE_rep1.5a88bc43459adc2cc757c01a5aa787b9.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/K562_BRD4_WCE_rep1/K562_BRD4_WCE_rep1.merged.bam \
 OUTPUT=alignment/K562_BRD4_WCE_rep1/K562_BRD4_WCE_rep1.sorted.dup.bam \
 METRICS_FILE=alignment/K562_BRD4_WCE_rep1/K562_BRD4_WCE_rep1.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.K562_BRD4_WCE_rep1.5a88bc43459adc2cc757c01a5aa787b9.mugqic.done
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
# JOB: picard_mark_duplicates_9_JOB_ID: picard_mark_duplicates.K562_BRD4_WCE_rep2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.K562_BRD4_WCE_rep2
JOB_DEPENDENCIES=$picard_merge_sam_files_9_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.K562_BRD4_WCE_rep2.f300ff5ac0aa7c3b3364d7584a6513f7.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.K562_BRD4_WCE_rep2.f300ff5ac0aa7c3b3364d7584a6513f7.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/K562_BRD4_WCE_rep2/K562_BRD4_WCE_rep2.merged.bam \
 OUTPUT=alignment/K562_BRD4_WCE_rep2/K562_BRD4_WCE_rep2.sorted.dup.bam \
 METRICS_FILE=alignment/K562_BRD4_WCE_rep2/K562_BRD4_WCE_rep2.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.K562_BRD4_WCE_rep2.f300ff5ac0aa7c3b3364d7584a6513f7.mugqic.done
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
# JOB: picard_mark_duplicates_10_JOB_ID: picard_mark_duplicates_report
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates_report
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_2_JOB_ID:$picard_mark_duplicates_3_JOB_ID:$picard_mark_duplicates_4_JOB_ID:$picard_mark_duplicates_5_JOB_ID:$picard_mark_duplicates_6_JOB_ID:$picard_mark_duplicates_7_JOB_ID:$picard_mark_duplicates_8_JOB_ID:$picard_mark_duplicates_9_JOB_ID
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$picard_mark_duplicates_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

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
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_2_JOB_ID:$picard_mark_duplicates_3_JOB_ID:$picard_mark_duplicates_4_JOB_ID:$picard_mark_duplicates_5_JOB_ID:$picard_mark_duplicates_6_JOB_ID:$picard_mark_duplicates_7_JOB_ID:$picard_mark_duplicates_8_JOB_ID:$picard_mark_duplicates_9_JOB_ID
JOB_DONE=job_output/metrics/metrics.flagstat.d44a38ea2904575b6e49131fd8a4339d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'metrics.flagstat.d44a38ea2904575b6e49131fd8a4339d.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools flagstat \
  alignment/K562_WCE_rep1/K562_WCE_rep1.sorted.dup.bam \
  > alignment/K562_WCE_rep1/K562_WCE_rep1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/K562_CDK9_rep1/K562_CDK9_rep1.sorted.dup.bam \
  > alignment/K562_CDK9_rep1/K562_CDK9_rep1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/K562_SMC1A_rep1/K562_SMC1A_rep1.sorted.dup.bam \
  > alignment/K562_SMC1A_rep1/K562_SMC1A_rep1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/K562_MED1_rep1/K562_MED1_rep1.sorted.dup.bam \
  > alignment/K562_MED1_rep1/K562_MED1_rep1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/K562_NIPBL_rep1/K562_NIPBL_rep1.sorted.dup.bam \
  > alignment/K562_NIPBL_rep1/K562_NIPBL_rep1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/K562_BRD4_rep1/K562_BRD4_rep1.sorted.dup.bam \
  > alignment/K562_BRD4_rep1/K562_BRD4_rep1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/K562_BRD4_rep2/K562_BRD4_rep2.sorted.dup.bam \
  > alignment/K562_BRD4_rep2/K562_BRD4_rep2.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/K562_BRD4_WCE_rep1/K562_BRD4_WCE_rep1.sorted.dup.bam \
  > alignment/K562_BRD4_WCE_rep1/K562_BRD4_WCE_rep1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/K562_BRD4_WCE_rep2/K562_BRD4_WCE_rep2.sorted.dup.bam \
  > alignment/K562_BRD4_WCE_rep2/K562_BRD4_WCE_rep2.sorted.dup.bam.flagstat
metrics.flagstat.d44a38ea2904575b6e49131fd8a4339d.mugqic.done
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
JOB_DONE=job_output/metrics/metrics_report.31e844e1c6ddb2ad20f3100863992baf.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'metrics_report.31e844e1c6ddb2ad20f3100863992baf.mugqic.done'
module load mugqic/pandoc/1.15.2 && \
for sample in K562_WCE_rep1 K562_CDK9_rep1 K562_SMC1A_rep1 K562_MED1_rep1 K562_NIPBL_rep1 K562_BRD4_rep1 K562_BRD4_rep2 K562_BRD4_WCE_rep1 K562_BRD4_WCE_rep2
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

metrics_report.31e844e1c6ddb2ad20f3100863992baf.mugqic.done
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
# JOB: homer_make_tag_directory_1_JOB_ID: homer_make_tag_directory.K562_WCE_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.K562_WCE_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.K562_WCE_rep1.95a584620faeecae56f265f324cfa1bf.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.K562_WCE_rep1.95a584620faeecae56f265f324cfa1bf.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/K562_WCE_rep1 \
            alignment/K562_WCE_rep1/K562_WCE_rep1.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.K562_WCE_rep1.95a584620faeecae56f265f324cfa1bf.mugqic.done
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
# JOB: homer_make_tag_directory_2_JOB_ID: homer_make_tag_directory.K562_CDK9_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.K562_CDK9_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_2_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.K562_CDK9_rep1.8e6c77838baa9839a842be5049f503e2.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.K562_CDK9_rep1.8e6c77838baa9839a842be5049f503e2.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/K562_CDK9_rep1 \
            alignment/K562_CDK9_rep1/K562_CDK9_rep1.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.K562_CDK9_rep1.8e6c77838baa9839a842be5049f503e2.mugqic.done
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
# JOB: homer_make_tag_directory_3_JOB_ID: homer_make_tag_directory.K562_SMC1A_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.K562_SMC1A_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_3_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.K562_SMC1A_rep1.be48ad5895d2f046a0398f0233e6015d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.K562_SMC1A_rep1.be48ad5895d2f046a0398f0233e6015d.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/K562_SMC1A_rep1 \
            alignment/K562_SMC1A_rep1/K562_SMC1A_rep1.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.K562_SMC1A_rep1.be48ad5895d2f046a0398f0233e6015d.mugqic.done
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
# JOB: homer_make_tag_directory_4_JOB_ID: homer_make_tag_directory.K562_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.K562_MED1_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_4_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.K562_MED1_rep1.2932bde3ab73accfb1121c3c46628de2.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.K562_MED1_rep1.2932bde3ab73accfb1121c3c46628de2.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/K562_MED1_rep1 \
            alignment/K562_MED1_rep1/K562_MED1_rep1.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.K562_MED1_rep1.2932bde3ab73accfb1121c3c46628de2.mugqic.done
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
# JOB: homer_make_tag_directory_5_JOB_ID: homer_make_tag_directory.K562_NIPBL_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.K562_NIPBL_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_5_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.K562_NIPBL_rep1.7b4da0bd375c35c7d79abb5ef8ed3d6e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.K562_NIPBL_rep1.7b4da0bd375c35c7d79abb5ef8ed3d6e.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/K562_NIPBL_rep1 \
            alignment/K562_NIPBL_rep1/K562_NIPBL_rep1.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.K562_NIPBL_rep1.7b4da0bd375c35c7d79abb5ef8ed3d6e.mugqic.done
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
# JOB: homer_make_tag_directory_6_JOB_ID: homer_make_tag_directory.K562_BRD4_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.K562_BRD4_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_6_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.K562_BRD4_rep1.f0f8a47167cf4ce5cff1f15f71a84c12.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.K562_BRD4_rep1.f0f8a47167cf4ce5cff1f15f71a84c12.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/K562_BRD4_rep1 \
            alignment/K562_BRD4_rep1/K562_BRD4_rep1.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.K562_BRD4_rep1.f0f8a47167cf4ce5cff1f15f71a84c12.mugqic.done
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
# JOB: homer_make_tag_directory_7_JOB_ID: homer_make_tag_directory.K562_BRD4_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.K562_BRD4_rep2
JOB_DEPENDENCIES=$picard_mark_duplicates_7_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.K562_BRD4_rep2.b5d94e93297d0d824daa93cc58a58c3a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.K562_BRD4_rep2.b5d94e93297d0d824daa93cc58a58c3a.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/K562_BRD4_rep2 \
            alignment/K562_BRD4_rep2/K562_BRD4_rep2.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.K562_BRD4_rep2.b5d94e93297d0d824daa93cc58a58c3a.mugqic.done
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
# JOB: homer_make_tag_directory_8_JOB_ID: homer_make_tag_directory.K562_BRD4_WCE_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.K562_BRD4_WCE_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_8_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.K562_BRD4_WCE_rep1.6ec6c5a60b2b9f6bd3ec32388521617c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.K562_BRD4_WCE_rep1.6ec6c5a60b2b9f6bd3ec32388521617c.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/K562_BRD4_WCE_rep1 \
            alignment/K562_BRD4_WCE_rep1/K562_BRD4_WCE_rep1.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.K562_BRD4_WCE_rep1.6ec6c5a60b2b9f6bd3ec32388521617c.mugqic.done
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
# JOB: homer_make_tag_directory_9_JOB_ID: homer_make_tag_directory.K562_BRD4_WCE_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.K562_BRD4_WCE_rep2
JOB_DEPENDENCIES=$picard_mark_duplicates_9_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.K562_BRD4_WCE_rep2.b3245a2f3e841f94346169b568982fe3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.K562_BRD4_WCE_rep2.b3245a2f3e841f94346169b568982fe3.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/K562_BRD4_WCE_rep2 \
            alignment/K562_BRD4_WCE_rep2/K562_BRD4_WCE_rep2.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.K562_BRD4_WCE_rep2.b3245a2f3e841f94346169b568982fe3.mugqic.done
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
# STEP: qc_metrics
#-------------------------------------------------------------------------------
STEP=qc_metrics
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: qc_metrics_1_JOB_ID: qc_plots_R
#-------------------------------------------------------------------------------
JOB_NAME=qc_plots_R
JOB_DEPENDENCIES=$homer_make_tag_directory_1_JOB_ID:$homer_make_tag_directory_2_JOB_ID:$homer_make_tag_directory_3_JOB_ID:$homer_make_tag_directory_4_JOB_ID:$homer_make_tag_directory_5_JOB_ID:$homer_make_tag_directory_6_JOB_ID:$homer_make_tag_directory_7_JOB_ID:$homer_make_tag_directory_8_JOB_ID:$homer_make_tag_directory_9_JOB_ID
JOB_DONE=job_output/qc_metrics/qc_plots_R.733b58d35f4829b341fc2efb3a27c4a9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'qc_plots_R.733b58d35f4829b341fc2efb3a27c4a9.mugqic.done'
module load mugqic/mugqic_tools/2.1.9 mugqic/R_Bioconductor/3.5.0_3.7 && \
mkdir -p graphs && \
Rscript $R_TOOLS/chipSeqGenerateQCMetrics.R \
  ../../raw/chip-seq/design.txt \
  /scratch/chris11/sb_cofactor/K562/output/chip-pipeline-GRCh38 && \
cp /cvmfs/soft.mugqic/CentOS6/software/genpipes/genpipes-3.1.2/bfx/report/ChipSeq.qc_metrics.md report/ChipSeq.qc_metrics.md && \
for sample in K562_WCE_rep1 K562_CDK9_rep1 K562_SMC1A_rep1 K562_MED1_rep1 K562_NIPBL_rep1 K562_BRD4_rep1 K562_BRD4_rep2 K562_BRD4_WCE_rep1 K562_BRD4_WCE_rep2
do
  cp --parents graphs/${sample}_QC_Metrics.ps report/
  convert -rotate 90 graphs/${sample}_QC_Metrics.ps report/graphs/${sample}_QC_Metrics.png
  echo -e "----

![QC Metrics for Sample $sample ([download high-res image](graphs/${sample}_QC_Metrics.ps))](graphs/${sample}_QC_Metrics.png)
" \
  >> report/ChipSeq.qc_metrics.md
done
qc_plots_R.733b58d35f4829b341fc2efb3a27c4a9.mugqic.done
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
# JOB: macs2_callpeak_1_JOB_ID: macs2_callpeak.K562_CDK9_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.K562_CDK9_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_2_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.K562_CDK9_rep1.5c3330dd23c98b7185f22d15c30b2268.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.K562_CDK9_rep1.5c3330dd23c98b7185f22d15c30b2268.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/K562_CDK9_rep1 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/K562_CDK9_rep1/K562_CDK9_rep1.sorted.dup.bam \
  --control \
  alignment/K562_WCE_rep1/K562_WCE_rep1.sorted.dup.bam \
  --name peak_call/K562_CDK9_rep1/K562_CDK9_rep1 \
  >& peak_call/K562_CDK9_rep1/K562_CDK9_rep1.diag.macs.out
macs2_callpeak.K562_CDK9_rep1.5c3330dd23c98b7185f22d15c30b2268.mugqic.done
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
# JOB: macs2_callpeak_2_JOB_ID: macs2_callpeak_bigBed.K562_CDK9_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.K562_CDK9_rep1
JOB_DEPENDENCIES=$macs2_callpeak_1_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.K562_CDK9_rep1.037cf7f3ee05a7a38589f6a35378e729.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.K562_CDK9_rep1.037cf7f3ee05a7a38589f6a35378e729.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/K562_CDK9_rep1/K562_CDK9_rep1_peaks.narrowPeak > peak_call/K562_CDK9_rep1/K562_CDK9_rep1_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/K562_CDK9_rep1/K562_CDK9_rep1_peaks.narrowPeak.bed \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/K562_CDK9_rep1/K562_CDK9_rep1_peaks.narrowPeak.bb
macs2_callpeak_bigBed.K562_CDK9_rep1.037cf7f3ee05a7a38589f6a35378e729.mugqic.done
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
# JOB: macs2_callpeak_3_JOB_ID: macs2_callpeak.K562_SMC1A_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.K562_SMC1A_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_3_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.K562_SMC1A_rep1.83b93e26a2cefb782bb9090bfff1da9d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.K562_SMC1A_rep1.83b93e26a2cefb782bb9090bfff1da9d.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/K562_SMC1A_rep1 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/K562_SMC1A_rep1/K562_SMC1A_rep1.sorted.dup.bam \
  --control \
  alignment/K562_WCE_rep1/K562_WCE_rep1.sorted.dup.bam \
  --name peak_call/K562_SMC1A_rep1/K562_SMC1A_rep1 \
  >& peak_call/K562_SMC1A_rep1/K562_SMC1A_rep1.diag.macs.out
macs2_callpeak.K562_SMC1A_rep1.83b93e26a2cefb782bb9090bfff1da9d.mugqic.done
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
# JOB: macs2_callpeak_4_JOB_ID: macs2_callpeak_bigBed.K562_SMC1A_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.K562_SMC1A_rep1
JOB_DEPENDENCIES=$macs2_callpeak_3_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.K562_SMC1A_rep1.c4dd42126fdf5d77f59751998af4ecb7.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.K562_SMC1A_rep1.c4dd42126fdf5d77f59751998af4ecb7.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/K562_SMC1A_rep1/K562_SMC1A_rep1_peaks.narrowPeak > peak_call/K562_SMC1A_rep1/K562_SMC1A_rep1_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/K562_SMC1A_rep1/K562_SMC1A_rep1_peaks.narrowPeak.bed \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/K562_SMC1A_rep1/K562_SMC1A_rep1_peaks.narrowPeak.bb
macs2_callpeak_bigBed.K562_SMC1A_rep1.c4dd42126fdf5d77f59751998af4ecb7.mugqic.done
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
# JOB: macs2_callpeak_5_JOB_ID: macs2_callpeak.K562_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.K562_MED1_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_4_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.K562_MED1_rep1.cedefc5f8f8f718e0bf0d86becb0fe90.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.K562_MED1_rep1.cedefc5f8f8f718e0bf0d86becb0fe90.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/K562_MED1_rep1 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/K562_MED1_rep1/K562_MED1_rep1.sorted.dup.bam \
  --control \
  alignment/K562_WCE_rep1/K562_WCE_rep1.sorted.dup.bam \
  --name peak_call/K562_MED1_rep1/K562_MED1_rep1 \
  >& peak_call/K562_MED1_rep1/K562_MED1_rep1.diag.macs.out
macs2_callpeak.K562_MED1_rep1.cedefc5f8f8f718e0bf0d86becb0fe90.mugqic.done
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
# JOB: macs2_callpeak_6_JOB_ID: macs2_callpeak_bigBed.K562_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.K562_MED1_rep1
JOB_DEPENDENCIES=$macs2_callpeak_5_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.K562_MED1_rep1.52c53f6bc00401b51d0f1222205c8ba7.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.K562_MED1_rep1.52c53f6bc00401b51d0f1222205c8ba7.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/K562_MED1_rep1/K562_MED1_rep1_peaks.narrowPeak > peak_call/K562_MED1_rep1/K562_MED1_rep1_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/K562_MED1_rep1/K562_MED1_rep1_peaks.narrowPeak.bed \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/K562_MED1_rep1/K562_MED1_rep1_peaks.narrowPeak.bb
macs2_callpeak_bigBed.K562_MED1_rep1.52c53f6bc00401b51d0f1222205c8ba7.mugqic.done
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
# JOB: macs2_callpeak_7_JOB_ID: macs2_callpeak.K562_NIPBL_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.K562_NIPBL_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_5_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.K562_NIPBL_rep1.44d30b48b43d06b9bcc9b7e8c4d53ef0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.K562_NIPBL_rep1.44d30b48b43d06b9bcc9b7e8c4d53ef0.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/K562_NIPBL_rep1 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/K562_NIPBL_rep1/K562_NIPBL_rep1.sorted.dup.bam \
  --control \
  alignment/K562_WCE_rep1/K562_WCE_rep1.sorted.dup.bam \
  --name peak_call/K562_NIPBL_rep1/K562_NIPBL_rep1 \
  >& peak_call/K562_NIPBL_rep1/K562_NIPBL_rep1.diag.macs.out
macs2_callpeak.K562_NIPBL_rep1.44d30b48b43d06b9bcc9b7e8c4d53ef0.mugqic.done
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
# JOB: macs2_callpeak_8_JOB_ID: macs2_callpeak_bigBed.K562_NIPBL_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.K562_NIPBL_rep1
JOB_DEPENDENCIES=$macs2_callpeak_7_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.K562_NIPBL_rep1.ff24b6de71de0054dc341eb769766789.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.K562_NIPBL_rep1.ff24b6de71de0054dc341eb769766789.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/K562_NIPBL_rep1/K562_NIPBL_rep1_peaks.narrowPeak > peak_call/K562_NIPBL_rep1/K562_NIPBL_rep1_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/K562_NIPBL_rep1/K562_NIPBL_rep1_peaks.narrowPeak.bed \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/K562_NIPBL_rep1/K562_NIPBL_rep1_peaks.narrowPeak.bb
macs2_callpeak_bigBed.K562_NIPBL_rep1.ff24b6de71de0054dc341eb769766789.mugqic.done
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
# JOB: macs2_callpeak_9_JOB_ID: macs2_callpeak.K562_BRD4_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.K562_BRD4_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_6_JOB_ID:$picard_mark_duplicates_8_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.K562_BRD4_rep1.111a54d42de0bb48d875c63e62580887.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.K562_BRD4_rep1.111a54d42de0bb48d875c63e62580887.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/K562_BRD4_rep1 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/K562_BRD4_rep1/K562_BRD4_rep1.sorted.dup.bam \
  --control \
  alignment/K562_BRD4_WCE_rep1/K562_BRD4_WCE_rep1.sorted.dup.bam \
  --name peak_call/K562_BRD4_rep1/K562_BRD4_rep1 \
  >& peak_call/K562_BRD4_rep1/K562_BRD4_rep1.diag.macs.out
macs2_callpeak.K562_BRD4_rep1.111a54d42de0bb48d875c63e62580887.mugqic.done
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
# JOB: macs2_callpeak_10_JOB_ID: macs2_callpeak_bigBed.K562_BRD4_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.K562_BRD4_rep1
JOB_DEPENDENCIES=$macs2_callpeak_9_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.K562_BRD4_rep1.05d6222eb69cfd7e350ba7b2a502dd24.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.K562_BRD4_rep1.05d6222eb69cfd7e350ba7b2a502dd24.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/K562_BRD4_rep1/K562_BRD4_rep1_peaks.narrowPeak > peak_call/K562_BRD4_rep1/K562_BRD4_rep1_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/K562_BRD4_rep1/K562_BRD4_rep1_peaks.narrowPeak.bed \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/K562_BRD4_rep1/K562_BRD4_rep1_peaks.narrowPeak.bb
macs2_callpeak_bigBed.K562_BRD4_rep1.05d6222eb69cfd7e350ba7b2a502dd24.mugqic.done
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
# JOB: macs2_callpeak_11_JOB_ID: macs2_callpeak.K562_BRD4_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.K562_BRD4_rep2
JOB_DEPENDENCIES=$picard_mark_duplicates_7_JOB_ID:$picard_mark_duplicates_9_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.K562_BRD4_rep2.525cc46a092917dbe807bfe886558515.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.K562_BRD4_rep2.525cc46a092917dbe807bfe886558515.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/K562_BRD4_rep2 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/K562_BRD4_rep2/K562_BRD4_rep2.sorted.dup.bam \
  --control \
  alignment/K562_BRD4_WCE_rep2/K562_BRD4_WCE_rep2.sorted.dup.bam \
  --name peak_call/K562_BRD4_rep2/K562_BRD4_rep2 \
  >& peak_call/K562_BRD4_rep2/K562_BRD4_rep2.diag.macs.out
macs2_callpeak.K562_BRD4_rep2.525cc46a092917dbe807bfe886558515.mugqic.done
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
# JOB: macs2_callpeak_12_JOB_ID: macs2_callpeak_bigBed.K562_BRD4_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.K562_BRD4_rep2
JOB_DEPENDENCIES=$macs2_callpeak_11_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.K562_BRD4_rep2.a0bf75f54a492009b4cdd44ac620525a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.K562_BRD4_rep2.a0bf75f54a492009b4cdd44ac620525a.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/K562_BRD4_rep2/K562_BRD4_rep2_peaks.narrowPeak > peak_call/K562_BRD4_rep2/K562_BRD4_rep2_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/K562_BRD4_rep2/K562_BRD4_rep2_peaks.narrowPeak.bed \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/K562_BRD4_rep2/K562_BRD4_rep2_peaks.narrowPeak.bb
macs2_callpeak_bigBed.K562_BRD4_rep2.a0bf75f54a492009b4cdd44ac620525a.mugqic.done
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
# JOB: macs2_callpeak_13_JOB_ID: macs2_callpeak_report
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_report
JOB_DEPENDENCIES=$macs2_callpeak_1_JOB_ID:$macs2_callpeak_3_JOB_ID:$macs2_callpeak_5_JOB_ID:$macs2_callpeak_7_JOB_ID:$macs2_callpeak_9_JOB_ID:$macs2_callpeak_11_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_report.fe2ea4051828a6a3789677d5fa97f4c8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_report.fe2ea4051828a6a3789677d5fa97f4c8.mugqic.done'
mkdir -p report && \
cp /cvmfs/soft.mugqic/CentOS6/software/genpipes/genpipes-3.1.2/bfx/report/ChipSeq.macs2_callpeak.md report/ && \
for contrast in K562_CDK9_rep1 K562_SMC1A_rep1 K562_MED1_rep1 K562_NIPBL_rep1 K562_BRD4_rep1 K562_BRD4_rep2
do
  cp -a --parents peak_call/$contrast/ report/ && \
  echo -e "* [Peak Calls File for Design $contrast](peak_call/$contrast/${contrast}_peaks.xls)" \
  >> report/ChipSeq.macs2_callpeak.md
done
macs2_callpeak_report.fe2ea4051828a6a3789677d5fa97f4c8.mugqic.done
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
# Call home with pipeline statistics
#-------------------------------------------------------------------------------
LOG_MD5=$(echo $USER-'206.12.124.2-ChipSeq-K562_SMC1A_rep1.HI.1385.007.Index_9.K562_SMC1_MF_ChIP2_3_rep1.bam,K562_NIPBL_rep1.HI.1385.007.Index_7.K562_NIPBL_MF_ChIP2_1_rep1.bam,K562_BRD4_WCE_rep1.ENCFF000BWK.fastq,K562_WCE_rep1.HI.1385.007.Index_6.K562_WCE_MF_rep1.bam,K562_BRD4_WCE_rep2.ENCFF994FIB.fastq,K562_BRD4_rep2.ENCFF049ACQ.fastq,K562_BRD4_rep1.ENCFF335PHG.fastq,K562_MED1_rep1.HI.1385.007.Index_8.K562_MED1_MF_ChIP2_2_rep1.bam,K562_CDK9_rep1.HI.1385.007.Index_10.K562_CDK9_MF_ChIP2_4_rep1.bam' | md5sum | awk '{ print $1 }')
wget "http://mugqic.hpc.mcgill.ca/cgi-bin/pipeline.cgi?hostname=cedar1.cedar.computecanada.ca&ip=206.12.124.2&pipeline=ChipSeq&steps=picard_sam_to_fastq,trimmomatic,merge_trimmomatic_stats,bwa_mem_picard_sort_sam,samtools_view_filter,picard_merge_sam_files,picard_mark_duplicates,metrics,homer_make_tag_directory,qc_metrics,macs2_callpeak&samples=9&md5=$LOG_MD5" --quiet --output-document=/dev/null

