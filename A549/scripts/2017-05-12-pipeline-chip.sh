#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

#-------------------------------------------------------------------------------
# ChipSeq PBSScheduler Job Submission Bash script
# Version: 2.2.0
# Created on: 2017-05-12T14:40:49
# Steps:
#   picard_sam_to_fastq: 16 jobs
#   trimmomatic: 16 jobs
#   merge_trimmomatic_stats: 1 job
#   bwa_mem_picard_sort_sam: 17 jobs
#   samtools_view_filter: 17 jobs
#   picard_merge_sam_files: 11 jobs
#   picard_mark_duplicates: 12 jobs
#   metrics: 2 jobs
#   homer_make_tag_directory: 11 jobs
#   qc_metrics: 1 job
#   homer_make_ucsc_file: 12 jobs
#   macs2_callpeak: 10 jobs
#   homer_annotate_peaks: 10 jobs
#   homer_find_motifs_genome: 10 jobs
#   annotation_graphs: 1 job
#   TOTAL: 147 jobs
#-------------------------------------------------------------------------------

OUTPUT_DIR=/gs/scratch/efournier/CofactorHR/A549/output/chip-pipeline
JOB_OUTPUT_DIR=$OUTPUT_DIR/job_output
TIMESTAMP=`date +%FT%H.%M.%S`
JOB_LIST=$JOB_OUTPUT_DIR/ChipSeq_job_list_$TIMESTAMP
mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR


#-------------------------------------------------------------------------------
# STEP: picard_sam_to_fastq
#-------------------------------------------------------------------------------
STEP=picard_sam_to_fastq
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_1_JOB_ID: picard_sam_to_fastq.HI.1613.003.Index_1.A549_Dex_CDK9_MF_ChIP1_4_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.1613.003.Index_1.A549_Dex_CDK9_MF_ChIP1_4_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.1613.003.Index_1.A549_Dex_CDK9_MF_ChIP1_4_rep1.bam.e3dda43eb2d16190c596aa563f9e298e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.1613.003.Index_1.A549_Dex_CDK9_MF_ChIP1_4_rep1.bam.e3dda43eb2d16190c596aa563f9e298e.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/SamToFastq.jar \
  VALIDATION_STRINGENCY=LENIENT \
  INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.1613.003.Index_1.A549_Dex_CDK9_MF_ChIP1_4_rep1.bam \
  FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.1613.003.Index_1.A549_Dex_CDK9_MF_ChIP1_4_rep1.single.fastq.gz
picard_sam_to_fastq.HI.1613.003.Index_1.A549_Dex_CDK9_MF_ChIP1_4_rep1.bam.e3dda43eb2d16190c596aa563f9e298e.mugqic.done
)
picard_sam_to_fastq_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=3 | grep "[0-9]")
echo "$picard_sam_to_fastq_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_2_JOB_ID: picard_sam_to_fastq.HI.1613.003.Index_12.A549_Dex_MED1_MF_ChIP1_2_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.1613.003.Index_12.A549_Dex_MED1_MF_ChIP1_2_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.1613.003.Index_12.A549_Dex_MED1_MF_ChIP1_2_rep1.bam.aeafad905dccabd75c553773f5044a44.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.1613.003.Index_12.A549_Dex_MED1_MF_ChIP1_2_rep1.bam.aeafad905dccabd75c553773f5044a44.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/SamToFastq.jar \
  VALIDATION_STRINGENCY=LENIENT \
  INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.1613.003.Index_12.A549_Dex_MED1_MF_ChIP1_2_rep1.bam \
  FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.1613.003.Index_12.A549_Dex_MED1_MF_ChIP1_2_rep1.single.fastq.gz
picard_sam_to_fastq.HI.1613.003.Index_12.A549_Dex_MED1_MF_ChIP1_2_rep1.bam.aeafad905dccabd75c553773f5044a44.mugqic.done
)
picard_sam_to_fastq_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=3 | grep "[0-9]")
echo "$picard_sam_to_fastq_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_3_JOB_ID: picard_sam_to_fastq.HI.1613.003.Index_19.A549_Dex_SMC1_MF_ChIP1_3_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.1613.003.Index_19.A549_Dex_SMC1_MF_ChIP1_3_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.1613.003.Index_19.A549_Dex_SMC1_MF_ChIP1_3_rep1.bam.3c948273b033bfdb7fc008e8b6c4371f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.1613.003.Index_19.A549_Dex_SMC1_MF_ChIP1_3_rep1.bam.3c948273b033bfdb7fc008e8b6c4371f.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/SamToFastq.jar \
  VALIDATION_STRINGENCY=LENIENT \
  INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.1613.003.Index_19.A549_Dex_SMC1_MF_ChIP1_3_rep1.bam \
  FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.1613.003.Index_19.A549_Dex_SMC1_MF_ChIP1_3_rep1.single.fastq.gz
picard_sam_to_fastq.HI.1613.003.Index_19.A549_Dex_SMC1_MF_ChIP1_3_rep1.bam.3c948273b033bfdb7fc008e8b6c4371f.mugqic.done
)
picard_sam_to_fastq_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=3 | grep "[0-9]")
echo "$picard_sam_to_fastq_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_4_JOB_ID: picard_sam_to_fastq.HI.1613.003.Index_4.A549_Dex_WCE_MF_ChIP1_8_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.1613.003.Index_4.A549_Dex_WCE_MF_ChIP1_8_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.1613.003.Index_4.A549_Dex_WCE_MF_ChIP1_8_rep1.bam.4974ade2f734eeb1cadfec5d4944a52d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.1613.003.Index_4.A549_Dex_WCE_MF_ChIP1_8_rep1.bam.4974ade2f734eeb1cadfec5d4944a52d.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/SamToFastq.jar \
  VALIDATION_STRINGENCY=LENIENT \
  INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.1613.003.Index_4.A549_Dex_WCE_MF_ChIP1_8_rep1.bam \
  FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.1613.003.Index_4.A549_Dex_WCE_MF_ChIP1_8_rep1.single.fastq.gz
picard_sam_to_fastq.HI.1613.003.Index_4.A549_Dex_WCE_MF_ChIP1_8_rep1.bam.4974ade2f734eeb1cadfec5d4944a52d.mugqic.done
)
picard_sam_to_fastq_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=3 | grep "[0-9]")
echo "$picard_sam_to_fastq_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_5_JOB_ID: picard_sam_to_fastq.HI.1613.003.Index_5.A549_Dex_NIPBL_MF_ChIP1_1_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.1613.003.Index_5.A549_Dex_NIPBL_MF_ChIP1_1_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.1613.003.Index_5.A549_Dex_NIPBL_MF_ChIP1_1_rep1.bam.b1acebc16a106634a250edb8382255b0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.1613.003.Index_5.A549_Dex_NIPBL_MF_ChIP1_1_rep1.bam.b1acebc16a106634a250edb8382255b0.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/SamToFastq.jar \
  VALIDATION_STRINGENCY=LENIENT \
  INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.1613.003.Index_5.A549_Dex_NIPBL_MF_ChIP1_1_rep1.bam \
  FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.1613.003.Index_5.A549_Dex_NIPBL_MF_ChIP1_1_rep1.single.fastq.gz
picard_sam_to_fastq.HI.1613.003.Index_5.A549_Dex_NIPBL_MF_ChIP1_1_rep1.bam.b1acebc16a106634a250edb8382255b0.mugqic.done
)
picard_sam_to_fastq_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=3 | grep "[0-9]")
echo "$picard_sam_to_fastq_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_6_JOB_ID: picard_sam_to_fastq.HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.4363222ad21d004c0591344a4638c689.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.4363222ad21d004c0591344a4638c689.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/SamToFastq.jar \
  VALIDATION_STRINGENCY=LENIENT \
  INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam \
  FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.single.fastq.gz
picard_sam_to_fastq.HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.4363222ad21d004c0591344a4638c689.mugqic.done
)
picard_sam_to_fastq_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=3 | grep "[0-9]")
echo "$picard_sam_to_fastq_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_7_JOB_ID: picard_sam_to_fastq.HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam.857137972b52eeecbcb9225e90e54d2a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam.857137972b52eeecbcb9225e90e54d2a.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/SamToFastq.jar \
  VALIDATION_STRINGENCY=LENIENT \
  INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam \
  FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.single.fastq.gz
picard_sam_to_fastq.HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam.857137972b52eeecbcb9225e90e54d2a.mugqic.done
)
picard_sam_to_fastq_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=3 | grep "[0-9]")
echo "$picard_sam_to_fastq_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_8_JOB_ID: picard_sam_to_fastq.HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.f24faf81e52f4fa8474e29d379f9c625.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.f24faf81e52f4fa8474e29d379f9c625.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/SamToFastq.jar \
  VALIDATION_STRINGENCY=LENIENT \
  INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam \
  FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.single.fastq.gz
picard_sam_to_fastq.HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.f24faf81e52f4fa8474e29d379f9c625.mugqic.done
)
picard_sam_to_fastq_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=3 | grep "[0-9]")
echo "$picard_sam_to_fastq_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_9_JOB_ID: picard_sam_to_fastq.HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.246716814b88747007205c887bae2549.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.246716814b88747007205c887bae2549.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/SamToFastq.jar \
  VALIDATION_STRINGENCY=LENIENT \
  INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam \
  FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.single.fastq.gz
picard_sam_to_fastq.HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.246716814b88747007205c887bae2549.mugqic.done
)
picard_sam_to_fastq_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=3 | grep "[0-9]")
echo "$picard_sam_to_fastq_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_10_JOB_ID: picard_sam_to_fastq.HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam.0bc54e916d4d14d599c4ed39eb0e99fb.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam.0bc54e916d4d14d599c4ed39eb0e99fb.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/SamToFastq.jar \
  VALIDATION_STRINGENCY=LENIENT \
  INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam \
  FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.single.fastq.gz
picard_sam_to_fastq.HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam.0bc54e916d4d14d599c4ed39eb0e99fb.mugqic.done
)
picard_sam_to_fastq_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=3 | grep "[0-9]")
echo "$picard_sam_to_fastq_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_11_JOB_ID: picard_sam_to_fastq.HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam.80413c1eb46475a53aec994419a1c0c4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam.80413c1eb46475a53aec994419a1c0c4.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/SamToFastq.jar \
  VALIDATION_STRINGENCY=LENIENT \
  INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam \
  FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.single.fastq.gz
picard_sam_to_fastq.HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam.80413c1eb46475a53aec994419a1c0c4.mugqic.done
)
picard_sam_to_fastq_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=3 | grep "[0-9]")
echo "$picard_sam_to_fastq_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_12_JOB_ID: picard_sam_to_fastq.HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.1334ebcd543d074a953370402ece5a60.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.1334ebcd543d074a953370402ece5a60.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/SamToFastq.jar \
  VALIDATION_STRINGENCY=LENIENT \
  INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam \
  FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.single.fastq.gz
picard_sam_to_fastq.HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.1334ebcd543d074a953370402ece5a60.mugqic.done
)
picard_sam_to_fastq_12_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=3 | grep "[0-9]")
echo "$picard_sam_to_fastq_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_13_JOB_ID: picard_sam_to_fastq.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.28b7d4f5b6d75e22eafa388bdcc0e220.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.28b7d4f5b6d75e22eafa388bdcc0e220.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/SamToFastq.jar \
  VALIDATION_STRINGENCY=LENIENT \
  INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam \
  FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.single.fastq.gz
picard_sam_to_fastq.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.28b7d4f5b6d75e22eafa388bdcc0e220.mugqic.done
)
picard_sam_to_fastq_13_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=3 | grep "[0-9]")
echo "$picard_sam_to_fastq_13_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_14_JOB_ID: picard_sam_to_fastq.HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam.47cc58be5246a7fa2f5bf220189c56ed.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam.47cc58be5246a7fa2f5bf220189c56ed.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/SamToFastq.jar \
  VALIDATION_STRINGENCY=LENIENT \
  INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam \
  FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.single.fastq.gz
picard_sam_to_fastq.HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam.47cc58be5246a7fa2f5bf220189c56ed.mugqic.done
)
picard_sam_to_fastq_14_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=3 | grep "[0-9]")
echo "$picard_sam_to_fastq_14_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_15_JOB_ID: picard_sam_to_fastq.HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.38827e4b8c26d6083923c2a98479097d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.38827e4b8c26d6083923c2a98479097d.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/SamToFastq.jar \
  VALIDATION_STRINGENCY=LENIENT \
  INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam \
  FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.single.fastq.gz
picard_sam_to_fastq.HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.38827e4b8c26d6083923c2a98479097d.mugqic.done
)
picard_sam_to_fastq_15_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=3 | grep "[0-9]")
echo "$picard_sam_to_fastq_15_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_16_JOB_ID: picard_sam_to_fastq.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.722e253a4817b4964b7e7df403892398.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.722e253a4817b4964b7e7df403892398.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/SamToFastq.jar \
  VALIDATION_STRINGENCY=LENIENT \
  INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam \
  FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.single.fastq.gz
picard_sam_to_fastq.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.722e253a4817b4964b7e7df403892398.mugqic.done
)
picard_sam_to_fastq_16_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=3 | grep "[0-9]")
echo "$picard_sam_to_fastq_16_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: trimmomatic
#-------------------------------------------------------------------------------
STEP=trimmomatic
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: trimmomatic_1_JOB_ID: trimmomatic.HI.1613.003.Index_1.A549_Dex_CDK9_MF_ChIP1_4_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.1613.003.Index_1.A549_Dex_CDK9_MF_ChIP1_4_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_1_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.1613.003.Index_1.A549_Dex_CDK9_MF_ChIP1_4_rep1.bam.411e021bf1d0a642621e09dbb027db46.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.1613.003.Index_1.A549_Dex_CDK9_MF_ChIP1_4_rep1.bam.411e021bf1d0a642621e09dbb027db46.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/A549_DEX_CDK9_rep1 && \
`cat > trim/A549_DEX_CDK9_rep1/HI.1613.003.Index_1.A549_Dex_CDK9_MF_ChIP1_4_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.1613.003.Index_1.A549_Dex_CDK9_MF_ChIP1_4_rep1.single.fastq.gz \
  trim/A549_DEX_CDK9_rep1/HI.1613.003.Index_1.A549_Dex_CDK9_MF_ChIP1_4_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/A549_DEX_CDK9_rep1/HI.1613.003.Index_1.A549_Dex_CDK9_MF_ChIP1_4_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/A549_DEX_CDK9_rep1/HI.1613.003.Index_1.A549_Dex_CDK9_MF_ChIP1_4_rep1.bam.trim.log
trimmomatic.HI.1613.003.Index_1.A549_Dex_CDK9_MF_ChIP1_4_rep1.bam.411e021bf1d0a642621e09dbb027db46.mugqic.done
)
trimmomatic_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_2_JOB_ID: trimmomatic.HI.1613.003.Index_12.A549_Dex_MED1_MF_ChIP1_2_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.1613.003.Index_12.A549_Dex_MED1_MF_ChIP1_2_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_2_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.1613.003.Index_12.A549_Dex_MED1_MF_ChIP1_2_rep1.bam.2e3cfdcfb8ddf1319dc665adb9bb59ee.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.1613.003.Index_12.A549_Dex_MED1_MF_ChIP1_2_rep1.bam.2e3cfdcfb8ddf1319dc665adb9bb59ee.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/A549_DEX_MED1_rep1 && \
`cat > trim/A549_DEX_MED1_rep1/HI.1613.003.Index_12.A549_Dex_MED1_MF_ChIP1_2_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.1613.003.Index_12.A549_Dex_MED1_MF_ChIP1_2_rep1.single.fastq.gz \
  trim/A549_DEX_MED1_rep1/HI.1613.003.Index_12.A549_Dex_MED1_MF_ChIP1_2_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/A549_DEX_MED1_rep1/HI.1613.003.Index_12.A549_Dex_MED1_MF_ChIP1_2_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/A549_DEX_MED1_rep1/HI.1613.003.Index_12.A549_Dex_MED1_MF_ChIP1_2_rep1.bam.trim.log
trimmomatic.HI.1613.003.Index_12.A549_Dex_MED1_MF_ChIP1_2_rep1.bam.2e3cfdcfb8ddf1319dc665adb9bb59ee.mugqic.done
)
trimmomatic_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_3_JOB_ID: trimmomatic.HI.1613.003.Index_19.A549_Dex_SMC1_MF_ChIP1_3_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.1613.003.Index_19.A549_Dex_SMC1_MF_ChIP1_3_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_3_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.1613.003.Index_19.A549_Dex_SMC1_MF_ChIP1_3_rep1.bam.a536e71bbf61e1ded97d8a4999f08393.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.1613.003.Index_19.A549_Dex_SMC1_MF_ChIP1_3_rep1.bam.a536e71bbf61e1ded97d8a4999f08393.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/A549_DEX_BRD4_rep1 && \
`cat > trim/A549_DEX_BRD4_rep1/HI.1613.003.Index_19.A549_Dex_SMC1_MF_ChIP1_3_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.1613.003.Index_19.A549_Dex_SMC1_MF_ChIP1_3_rep1.single.fastq.gz \
  trim/A549_DEX_BRD4_rep1/HI.1613.003.Index_19.A549_Dex_SMC1_MF_ChIP1_3_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/A549_DEX_BRD4_rep1/HI.1613.003.Index_19.A549_Dex_SMC1_MF_ChIP1_3_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/A549_DEX_BRD4_rep1/HI.1613.003.Index_19.A549_Dex_SMC1_MF_ChIP1_3_rep1.bam.trim.log
trimmomatic.HI.1613.003.Index_19.A549_Dex_SMC1_MF_ChIP1_3_rep1.bam.a536e71bbf61e1ded97d8a4999f08393.mugqic.done
)
trimmomatic_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_4_JOB_ID: trimmomatic.HI.1613.003.Index_4.A549_Dex_WCE_MF_ChIP1_8_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.1613.003.Index_4.A549_Dex_WCE_MF_ChIP1_8_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_4_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.1613.003.Index_4.A549_Dex_WCE_MF_ChIP1_8_rep1.bam.928c328fcdeb1523ce20e6c366f3b279.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.1613.003.Index_4.A549_Dex_WCE_MF_ChIP1_8_rep1.bam.928c328fcdeb1523ce20e6c366f3b279.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/A549_DEX_WCE_rep1 && \
`cat > trim/A549_DEX_WCE_rep1/HI.1613.003.Index_4.A549_Dex_WCE_MF_ChIP1_8_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.1613.003.Index_4.A549_Dex_WCE_MF_ChIP1_8_rep1.single.fastq.gz \
  trim/A549_DEX_WCE_rep1/HI.1613.003.Index_4.A549_Dex_WCE_MF_ChIP1_8_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/A549_DEX_WCE_rep1/HI.1613.003.Index_4.A549_Dex_WCE_MF_ChIP1_8_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/A549_DEX_WCE_rep1/HI.1613.003.Index_4.A549_Dex_WCE_MF_ChIP1_8_rep1.bam.trim.log
trimmomatic.HI.1613.003.Index_4.A549_Dex_WCE_MF_ChIP1_8_rep1.bam.928c328fcdeb1523ce20e6c366f3b279.mugqic.done
)
trimmomatic_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_5_JOB_ID: trimmomatic.HI.1613.003.Index_5.A549_Dex_NIPBL_MF_ChIP1_1_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.1613.003.Index_5.A549_Dex_NIPBL_MF_ChIP1_1_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_5_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.1613.003.Index_5.A549_Dex_NIPBL_MF_ChIP1_1_rep1.bam.737ffa89e0bcf2df7cd7edafba9c970c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.1613.003.Index_5.A549_Dex_NIPBL_MF_ChIP1_1_rep1.bam.737ffa89e0bcf2df7cd7edafba9c970c.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/A549_DEX_NIPBL_rep1 && \
`cat > trim/A549_DEX_NIPBL_rep1/HI.1613.003.Index_5.A549_Dex_NIPBL_MF_ChIP1_1_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.1613.003.Index_5.A549_Dex_NIPBL_MF_ChIP1_1_rep1.single.fastq.gz \
  trim/A549_DEX_NIPBL_rep1/HI.1613.003.Index_5.A549_Dex_NIPBL_MF_ChIP1_1_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/A549_DEX_NIPBL_rep1/HI.1613.003.Index_5.A549_Dex_NIPBL_MF_ChIP1_1_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/A549_DEX_NIPBL_rep1/HI.1613.003.Index_5.A549_Dex_NIPBL_MF_ChIP1_1_rep1.bam.trim.log
trimmomatic.HI.1613.003.Index_5.A549_Dex_NIPBL_MF_ChIP1_1_rep1.bam.737ffa89e0bcf2df7cd7edafba9c970c.mugqic.done
)
trimmomatic_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_6_JOB_ID: trimmomatic.HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_6_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.ec5a6e2c096a3718b921830a77f172b7.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.ec5a6e2c096a3718b921830a77f172b7.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/A549_CTRL_SMC1A_rep1 && \
`cat > trim/A549_CTRL_SMC1A_rep1/HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.single.fastq.gz \
  trim/A549_CTRL_SMC1A_rep1/HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/A549_CTRL_SMC1A_rep1/HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/A549_CTRL_SMC1A_rep1/HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.trim.log
trimmomatic.HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.ec5a6e2c096a3718b921830a77f172b7.mugqic.done
)
trimmomatic_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_7_JOB_ID: trimmomatic.HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_7_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam.246b0bdf886cc89686d4d974357cc0f9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam.246b0bdf886cc89686d4d974357cc0f9.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/A549_CTRL_WCE_rep1 && \
`cat > trim/A549_CTRL_WCE_rep1/HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.single.fastq.gz \
  trim/A549_CTRL_WCE_rep1/HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/A549_CTRL_WCE_rep1/HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/A549_CTRL_WCE_rep1/HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam.trim.log
trimmomatic.HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam.246b0bdf886cc89686d4d974357cc0f9.mugqic.done
)
trimmomatic_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_8_JOB_ID: trimmomatic.HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_8_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.83ef5a7ac837e20bb7ddcb37201fcb38.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.83ef5a7ac837e20bb7ddcb37201fcb38.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/A549_CTRL_CDK9_rep1 && \
`cat > trim/A549_CTRL_CDK9_rep1/HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.single.fastq.gz \
  trim/A549_CTRL_CDK9_rep1/HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/A549_CTRL_CDK9_rep1/HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/A549_CTRL_CDK9_rep1/HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.trim.log
trimmomatic.HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.83ef5a7ac837e20bb7ddcb37201fcb38.mugqic.done
)
trimmomatic_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_9_JOB_ID: trimmomatic.HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_9_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.70fb0b3b574b36afe85303338f357c51.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.70fb0b3b574b36afe85303338f357c51.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/A549_CTRL_MED1_rep1 && \
`cat > trim/A549_CTRL_MED1_rep1/HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.single.fastq.gz \
  trim/A549_CTRL_MED1_rep1/HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/A549_CTRL_MED1_rep1/HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/A549_CTRL_MED1_rep1/HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.trim.log
trimmomatic.HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.70fb0b3b574b36afe85303338f357c51.mugqic.done
)
trimmomatic_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_10_JOB_ID: trimmomatic.HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_10_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam.61963918aaa6de7e0b0962962435b799.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam.61963918aaa6de7e0b0962962435b799.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/A549_CTRL_NIPBL_rep1 && \
`cat > trim/A549_CTRL_NIPBL_rep1/HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.single.fastq.gz \
  trim/A549_CTRL_NIPBL_rep1/HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/A549_CTRL_NIPBL_rep1/HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/A549_CTRL_NIPBL_rep1/HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam.trim.log
trimmomatic.HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam.61963918aaa6de7e0b0962962435b799.mugqic.done
)
trimmomatic_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_11_JOB_ID: trimmomatic.HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_11_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam.161786c0ebff96bf27beff62be245d57.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam.161786c0ebff96bf27beff62be245d57.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/A549_DEX_BRD4_rep1 && \
`cat > trim/A549_DEX_BRD4_rep1/HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.single.fastq.gz \
  trim/A549_DEX_BRD4_rep1/HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/A549_DEX_BRD4_rep1/HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/A549_DEX_BRD4_rep1/HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam.trim.log
trimmomatic.HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam.161786c0ebff96bf27beff62be245d57.mugqic.done
)
trimmomatic_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_12_JOB_ID: trimmomatic.HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_12_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.dfb68a8be21b16424cc857fd336ffe23.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.dfb68a8be21b16424cc857fd336ffe23.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/A549_CTRL_BRD4_rep1 && \
`cat > trim/A549_CTRL_BRD4_rep1/HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.single.fastq.gz \
  trim/A549_CTRL_BRD4_rep1/HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/A549_CTRL_BRD4_rep1/HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/A549_CTRL_BRD4_rep1/HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.trim.log
trimmomatic.HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.dfb68a8be21b16424cc857fd336ffe23.mugqic.done
)
trimmomatic_12_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_13_JOB_ID: trimmomatic.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_13_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.4d0d63552a3c1dee95ceb7fb200fb397.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.4d0d63552a3c1dee95ceb7fb200fb397.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/A549_CTRL_SMC1A_rep1 && \
`cat > trim/A549_CTRL_SMC1A_rep1/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.single.fastq.gz \
  trim/A549_CTRL_SMC1A_rep1/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/A549_CTRL_SMC1A_rep1/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/A549_CTRL_SMC1A_rep1/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.trim.log
trimmomatic.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.4d0d63552a3c1dee95ceb7fb200fb397.mugqic.done
)
trimmomatic_13_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_13_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_14_JOB_ID: trimmomatic.HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_14_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam.8cb620ddd07a121116ad740a1e5cc09c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam.8cb620ddd07a121116ad740a1e5cc09c.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/A549_CTRL_WCE_rep1 && \
`cat > trim/A549_CTRL_WCE_rep1/HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.single.fastq.gz \
  trim/A549_CTRL_WCE_rep1/HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/A549_CTRL_WCE_rep1/HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/A549_CTRL_WCE_rep1/HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam.trim.log
trimmomatic.HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam.8cb620ddd07a121116ad740a1e5cc09c.mugqic.done
)
trimmomatic_14_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_14_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_15_JOB_ID: trimmomatic.HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_15_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.3bc52c771a0dbf164adfd163dda41735.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.3bc52c771a0dbf164adfd163dda41735.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/A549_CTRL_NIPBL_rep1 && \
`cat > trim/A549_CTRL_NIPBL_rep1/HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.single.fastq.gz \
  trim/A549_CTRL_NIPBL_rep1/HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/A549_CTRL_NIPBL_rep1/HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/A549_CTRL_NIPBL_rep1/HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.trim.log
trimmomatic.HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.3bc52c771a0dbf164adfd163dda41735.mugqic.done
)
trimmomatic_15_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_15_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_16_JOB_ID: trimmomatic.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_16_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.1ba111621b49ffd234f8a6f1558258fd.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.1ba111621b49ffd234f8a6f1558258fd.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/A549_CTRL_MED1_rep1 && \
`cat > trim/A549_CTRL_MED1_rep1/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.single.fastq.gz \
  trim/A549_CTRL_MED1_rep1/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/A549_CTRL_MED1_rep1/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/A549_CTRL_MED1_rep1/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.trim.log
trimmomatic.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.1ba111621b49ffd234f8a6f1558258fd.mugqic.done
)
trimmomatic_16_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_16_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: merge_trimmomatic_stats
#-------------------------------------------------------------------------------
STEP=merge_trimmomatic_stats
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: merge_trimmomatic_stats_1_JOB_ID: merge_trimmomatic_stats
#-------------------------------------------------------------------------------
JOB_NAME=merge_trimmomatic_stats
JOB_DEPENDENCIES=$trimmomatic_1_JOB_ID:$trimmomatic_2_JOB_ID:$trimmomatic_3_JOB_ID:$trimmomatic_4_JOB_ID:$trimmomatic_5_JOB_ID:$trimmomatic_6_JOB_ID:$trimmomatic_7_JOB_ID:$trimmomatic_8_JOB_ID:$trimmomatic_9_JOB_ID:$trimmomatic_10_JOB_ID:$trimmomatic_11_JOB_ID:$trimmomatic_12_JOB_ID:$trimmomatic_13_JOB_ID:$trimmomatic_14_JOB_ID:$trimmomatic_15_JOB_ID:$trimmomatic_16_JOB_ID
JOB_DONE=job_output/merge_trimmomatic_stats/merge_trimmomatic_stats.c75ced6e6438f5dd06aaa4ec8fd5ca5a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'merge_trimmomatic_stats.c75ced6e6438f5dd06aaa4ec8fd5ca5a.mugqic.done'
module load mugqic/pandoc/1.15.2 && \
mkdir -p metrics && \
echo 'Sample	Readset	Raw Single Reads #	Surviving Single Reads #	Surviving Single Reads %' > metrics/trimReadsetTable.tsv && \
grep ^Input trim/A549_DEX_CDK9_rep1/HI.1613.003.Index_1.A549_Dex_CDK9_MF_ChIP1_4_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/A549_DEX_CDK9_rep1	HI.1613.003.Index_1.A549_Dex_CDK9_MF_ChIP1_4_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/A549_DEX_MED1_rep1/HI.1613.003.Index_12.A549_Dex_MED1_MF_ChIP1_2_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/A549_DEX_MED1_rep1	HI.1613.003.Index_12.A549_Dex_MED1_MF_ChIP1_2_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/A549_DEX_BRD4_rep1/HI.1613.003.Index_19.A549_Dex_SMC1_MF_ChIP1_3_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/A549_DEX_BRD4_rep1	HI.1613.003.Index_19.A549_Dex_SMC1_MF_ChIP1_3_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/A549_DEX_WCE_rep1/HI.1613.003.Index_4.A549_Dex_WCE_MF_ChIP1_8_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/A549_DEX_WCE_rep1	HI.1613.003.Index_4.A549_Dex_WCE_MF_ChIP1_8_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/A549_DEX_NIPBL_rep1/HI.1613.003.Index_5.A549_Dex_NIPBL_MF_ChIP1_1_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/A549_DEX_NIPBL_rep1	HI.1613.003.Index_5.A549_Dex_NIPBL_MF_ChIP1_1_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/A549_CTRL_SMC1A_rep1/HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/A549_CTRL_SMC1A_rep1	HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/A549_CTRL_WCE_rep1/HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/A549_CTRL_WCE_rep1	HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/A549_CTRL_CDK9_rep1/HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/A549_CTRL_CDK9_rep1	HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/A549_CTRL_MED1_rep1/HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/A549_CTRL_MED1_rep1	HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/A549_CTRL_NIPBL_rep1/HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/A549_CTRL_NIPBL_rep1	HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/A549_DEX_BRD4_rep1/HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/A549_DEX_BRD4_rep1	HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/A549_CTRL_BRD4_rep1/HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/A549_CTRL_BRD4_rep1	HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/A549_CTRL_SMC1A_rep1/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/A549_CTRL_SMC1A_rep1	HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/A549_CTRL_WCE_rep1/HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/A549_CTRL_WCE_rep1	HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/A549_CTRL_NIPBL_rep1/HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/A549_CTRL_NIPBL_rep1	HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/A549_CTRL_MED1_rep1/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/A549_CTRL_MED1_rep1	HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
cut -f1,3- metrics/trimReadsetTable.tsv | awk -F"	" '{OFS="	"; if (NR==1) {if ($2=="Raw Paired Reads #") {paired=1};print "Sample", "Raw Reads #", "Surviving Reads #", "Surviving %"} else {if (paired) {$2=$2*2; $3=$3*2}; raw[$1]+=$2; surviving[$1]+=$3}}END{for (sample in raw){print sample, raw[sample], surviving[sample], surviving[sample] / raw[sample] * 100}}' \
  > metrics/trimSampleTable.tsv && \
mkdir -p report && \
cp metrics/trimReadsetTable.tsv metrics/trimSampleTable.tsv report/ && \
trim_readset_table_md=`LC_NUMERIC=en_CA awk -F "	" '{OFS="|"; if (NR == 1) {$1 = $1; print $0; print "-----|-----|-----:|-----:|-----:"} else {print $1, $2, sprintf("%\47d", $3), sprintf("%\47d", $4), sprintf("%.1f", $5)}}' metrics/trimReadsetTable.tsv` && \
pandoc \
  /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/Illumina.merge_trimmomatic_stats.md \
  --template /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/Illumina.merge_trimmomatic_stats.md \
  --variable trailing_min_quality=30 \
  --variable min_length=50 \
  --variable read_type=Single \
  --variable trim_readset_table="$trim_readset_table_md" \
  --to markdown \
  > report/Illumina.merge_trimmomatic_stats.md
merge_trimmomatic_stats.c75ced6e6438f5dd06aaa4ec8fd5ca5a.mugqic.done
)
merge_trimmomatic_stats_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$merge_trimmomatic_stats_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: bwa_mem_picard_sort_sam
#-------------------------------------------------------------------------------
STEP=bwa_mem_picard_sort_sam
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_1_JOB_ID: bwa_mem_picard_sort_sam.HI.1613.003.Index_1.A549_Dex_CDK9_MF_ChIP1_4_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.1613.003.Index_1.A549_Dex_CDK9_MF_ChIP1_4_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_1_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.1613.003.Index_1.A549_Dex_CDK9_MF_ChIP1_4_rep1.bam.15cb9e8b4bb78c199fa4750ba53f5892.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.1613.003.Index_1.A549_Dex_CDK9_MF_ChIP1_4_rep1.bam.15cb9e8b4bb78c199fa4750ba53f5892.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/A549_DEX_CDK9_rep1/HI.1613.003.Index_1.A549_Dex_CDK9_MF_ChIP1_4_rep1.bam && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:HI.1613.003.Index_1.A549_Dex_CDK9_MF_ChIP1_4_rep1.bam	SM:A549_DEX_CDK9_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/bwa_index/Homo_sapiens.hg19.fa \
  trim/A549_DEX_CDK9_rep1/HI.1613.003.Index_1.A549_Dex_CDK9_MF_ChIP1_4_rep1.bam.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/A549_DEX_CDK9_rep1/HI.1613.003.Index_1.A549_Dex_CDK9_MF_ChIP1_4_rep1.bam/HI.1613.003.Index_1.A549_Dex_CDK9_MF_ChIP1_4_rep1.bam.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.HI.1613.003.Index_1.A549_Dex_CDK9_MF_ChIP1_4_rep1.bam.15cb9e8b4bb78c199fa4750ba53f5892.mugqic.done
)
bwa_mem_picard_sort_sam_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_2_JOB_ID: bwa_mem_picard_sort_sam.HI.1613.003.Index_12.A549_Dex_MED1_MF_ChIP1_2_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.1613.003.Index_12.A549_Dex_MED1_MF_ChIP1_2_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_2_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.1613.003.Index_12.A549_Dex_MED1_MF_ChIP1_2_rep1.bam.071f57288888f04ee206a6079f2a5980.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.1613.003.Index_12.A549_Dex_MED1_MF_ChIP1_2_rep1.bam.071f57288888f04ee206a6079f2a5980.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/A549_DEX_MED1_rep1/HI.1613.003.Index_12.A549_Dex_MED1_MF_ChIP1_2_rep1.bam && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:HI.1613.003.Index_12.A549_Dex_MED1_MF_ChIP1_2_rep1.bam	SM:A549_DEX_MED1_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/bwa_index/Homo_sapiens.hg19.fa \
  trim/A549_DEX_MED1_rep1/HI.1613.003.Index_12.A549_Dex_MED1_MF_ChIP1_2_rep1.bam.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/A549_DEX_MED1_rep1/HI.1613.003.Index_12.A549_Dex_MED1_MF_ChIP1_2_rep1.bam/HI.1613.003.Index_12.A549_Dex_MED1_MF_ChIP1_2_rep1.bam.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.HI.1613.003.Index_12.A549_Dex_MED1_MF_ChIP1_2_rep1.bam.071f57288888f04ee206a6079f2a5980.mugqic.done
)
bwa_mem_picard_sort_sam_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_3_JOB_ID: bwa_mem_picard_sort_sam.HI.1613.003.Index_19.A549_Dex_SMC1_MF_ChIP1_3_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.1613.003.Index_19.A549_Dex_SMC1_MF_ChIP1_3_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_3_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.1613.003.Index_19.A549_Dex_SMC1_MF_ChIP1_3_rep1.bam.2e96d505769ef9680accc8f7a9df12de.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.1613.003.Index_19.A549_Dex_SMC1_MF_ChIP1_3_rep1.bam.2e96d505769ef9680accc8f7a9df12de.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/A549_DEX_BRD4_rep1/HI.1613.003.Index_19.A549_Dex_SMC1_MF_ChIP1_3_rep1.bam && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:HI.1613.003.Index_19.A549_Dex_SMC1_MF_ChIP1_3_rep1.bam	SM:A549_DEX_BRD4_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/bwa_index/Homo_sapiens.hg19.fa \
  trim/A549_DEX_BRD4_rep1/HI.1613.003.Index_19.A549_Dex_SMC1_MF_ChIP1_3_rep1.bam.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/A549_DEX_BRD4_rep1/HI.1613.003.Index_19.A549_Dex_SMC1_MF_ChIP1_3_rep1.bam/HI.1613.003.Index_19.A549_Dex_SMC1_MF_ChIP1_3_rep1.bam.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.HI.1613.003.Index_19.A549_Dex_SMC1_MF_ChIP1_3_rep1.bam.2e96d505769ef9680accc8f7a9df12de.mugqic.done
)
bwa_mem_picard_sort_sam_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_4_JOB_ID: bwa_mem_picard_sort_sam.HI.1613.003.Index_4.A549_Dex_WCE_MF_ChIP1_8_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.1613.003.Index_4.A549_Dex_WCE_MF_ChIP1_8_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_4_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.1613.003.Index_4.A549_Dex_WCE_MF_ChIP1_8_rep1.bam.043b703100b0ba72010bcd31a1809abb.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.1613.003.Index_4.A549_Dex_WCE_MF_ChIP1_8_rep1.bam.043b703100b0ba72010bcd31a1809abb.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/A549_DEX_WCE_rep1/HI.1613.003.Index_4.A549_Dex_WCE_MF_ChIP1_8_rep1.bam && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:HI.1613.003.Index_4.A549_Dex_WCE_MF_ChIP1_8_rep1.bam	SM:A549_DEX_WCE_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/bwa_index/Homo_sapiens.hg19.fa \
  trim/A549_DEX_WCE_rep1/HI.1613.003.Index_4.A549_Dex_WCE_MF_ChIP1_8_rep1.bam.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/A549_DEX_WCE_rep1/HI.1613.003.Index_4.A549_Dex_WCE_MF_ChIP1_8_rep1.bam/HI.1613.003.Index_4.A549_Dex_WCE_MF_ChIP1_8_rep1.bam.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.HI.1613.003.Index_4.A549_Dex_WCE_MF_ChIP1_8_rep1.bam.043b703100b0ba72010bcd31a1809abb.mugqic.done
)
bwa_mem_picard_sort_sam_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_5_JOB_ID: bwa_mem_picard_sort_sam.HI.1613.003.Index_5.A549_Dex_NIPBL_MF_ChIP1_1_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.1613.003.Index_5.A549_Dex_NIPBL_MF_ChIP1_1_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_5_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.1613.003.Index_5.A549_Dex_NIPBL_MF_ChIP1_1_rep1.bam.446a4c65d0d7a0116dcd5d9e809e6567.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.1613.003.Index_5.A549_Dex_NIPBL_MF_ChIP1_1_rep1.bam.446a4c65d0d7a0116dcd5d9e809e6567.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/A549_DEX_NIPBL_rep1/HI.1613.003.Index_5.A549_Dex_NIPBL_MF_ChIP1_1_rep1.bam && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:HI.1613.003.Index_5.A549_Dex_NIPBL_MF_ChIP1_1_rep1.bam	SM:A549_DEX_NIPBL_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/bwa_index/Homo_sapiens.hg19.fa \
  trim/A549_DEX_NIPBL_rep1/HI.1613.003.Index_5.A549_Dex_NIPBL_MF_ChIP1_1_rep1.bam.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/A549_DEX_NIPBL_rep1/HI.1613.003.Index_5.A549_Dex_NIPBL_MF_ChIP1_1_rep1.bam/HI.1613.003.Index_5.A549_Dex_NIPBL_MF_ChIP1_1_rep1.bam.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.HI.1613.003.Index_5.A549_Dex_NIPBL_MF_ChIP1_1_rep1.bam.446a4c65d0d7a0116dcd5d9e809e6567.mugqic.done
)
bwa_mem_picard_sort_sam_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_6_JOB_ID: bwa_mem_picard_sort_sam.HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_6_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.251e7f1d145d0694e7e0e9fce5002bec.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.251e7f1d145d0694e7e0e9fce5002bec.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/A549_CTRL_SMC1A_rep1/HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam	SM:A549_CTRL_SMC1A_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/bwa_index/Homo_sapiens.hg19.fa \
  trim/A549_CTRL_SMC1A_rep1/HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/A549_CTRL_SMC1A_rep1/HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam/HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.251e7f1d145d0694e7e0e9fce5002bec.mugqic.done
)
bwa_mem_picard_sort_sam_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_7_JOB_ID: bwa_mem_picard_sort_sam.HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_7_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam.1bb742fb41726b7320d0172e4a59597e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam.1bb742fb41726b7320d0172e4a59597e.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/A549_CTRL_WCE_rep1/HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam	SM:A549_CTRL_WCE_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/bwa_index/Homo_sapiens.hg19.fa \
  trim/A549_CTRL_WCE_rep1/HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/A549_CTRL_WCE_rep1/HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam/HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam.1bb742fb41726b7320d0172e4a59597e.mugqic.done
)
bwa_mem_picard_sort_sam_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_8_JOB_ID: bwa_mem_picard_sort_sam.HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_8_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.9e0fcf05798352d6b93009829d8effc5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.9e0fcf05798352d6b93009829d8effc5.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/A549_CTRL_CDK9_rep1/HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam	SM:A549_CTRL_CDK9_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/bwa_index/Homo_sapiens.hg19.fa \
  trim/A549_CTRL_CDK9_rep1/HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/A549_CTRL_CDK9_rep1/HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam/HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.9e0fcf05798352d6b93009829d8effc5.mugqic.done
)
bwa_mem_picard_sort_sam_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_9_JOB_ID: bwa_mem_picard_sort_sam.HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_9_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.26c9c481f4692a554fe6eec18341f518.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.26c9c481f4692a554fe6eec18341f518.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/A549_CTRL_MED1_rep1/HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam	SM:A549_CTRL_MED1_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/bwa_index/Homo_sapiens.hg19.fa \
  trim/A549_CTRL_MED1_rep1/HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/A549_CTRL_MED1_rep1/HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam/HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.26c9c481f4692a554fe6eec18341f518.mugqic.done
)
bwa_mem_picard_sort_sam_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_10_JOB_ID: bwa_mem_picard_sort_sam.HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_10_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam.9fe46efde6b37a4739c56f8aad6f9bb4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam.9fe46efde6b37a4739c56f8aad6f9bb4.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/A549_CTRL_NIPBL_rep1/HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam	SM:A549_CTRL_NIPBL_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/bwa_index/Homo_sapiens.hg19.fa \
  trim/A549_CTRL_NIPBL_rep1/HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/A549_CTRL_NIPBL_rep1/HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam/HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam.9fe46efde6b37a4739c56f8aad6f9bb4.mugqic.done
)
bwa_mem_picard_sort_sam_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_11_JOB_ID: bwa_mem_picard_sort_sam.HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_11_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam.0fb358d5591020a5d2a3ec125b42df9c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam.0fb358d5591020a5d2a3ec125b42df9c.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/A549_DEX_BRD4_rep1/HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam	SM:A549_DEX_BRD4_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/bwa_index/Homo_sapiens.hg19.fa \
  trim/A549_DEX_BRD4_rep1/HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/A549_DEX_BRD4_rep1/HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam/HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam.0fb358d5591020a5d2a3ec125b42df9c.mugqic.done
)
bwa_mem_picard_sort_sam_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_12_JOB_ID: bwa_mem_picard_sort_sam.HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_12_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.eafbd9d7f74d5e1f024324f3098ca7c9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.eafbd9d7f74d5e1f024324f3098ca7c9.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/A549_CTRL_BRD4_rep1/HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam	SM:A549_CTRL_BRD4_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/bwa_index/Homo_sapiens.hg19.fa \
  trim/A549_CTRL_BRD4_rep1/HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/A549_CTRL_BRD4_rep1/HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam/HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.eafbd9d7f74d5e1f024324f3098ca7c9.mugqic.done
)
bwa_mem_picard_sort_sam_12_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_13_JOB_ID: bwa_mem_picard_sort_sam.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_13_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.1ea8fb39a9e6497d76082025b44a4be8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.1ea8fb39a9e6497d76082025b44a4be8.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/A549_CTRL_SMC1A_rep1/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam	SM:A549_CTRL_SMC1A_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/bwa_index/Homo_sapiens.hg19.fa \
  trim/A549_CTRL_SMC1A_rep1/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/A549_CTRL_SMC1A_rep1/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.1ea8fb39a9e6497d76082025b44a4be8.mugqic.done
)
bwa_mem_picard_sort_sam_13_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_13_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_14_JOB_ID: bwa_mem_picard_sort_sam.HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_14_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam.10fef2642c848159f35494c4b4ca3e36.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam.10fef2642c848159f35494c4b4ca3e36.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/A549_CTRL_WCE_rep1/HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam	SM:A549_CTRL_WCE_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/bwa_index/Homo_sapiens.hg19.fa \
  trim/A549_CTRL_WCE_rep1/HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/A549_CTRL_WCE_rep1/HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam/HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam.10fef2642c848159f35494c4b4ca3e36.mugqic.done
)
bwa_mem_picard_sort_sam_14_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_14_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_15_JOB_ID: bwa_mem_picard_sort_sam.HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_15_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.d588716ad5c0f58d54bc41b8ec6e96f5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.d588716ad5c0f58d54bc41b8ec6e96f5.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/A549_CTRL_NIPBL_rep1/HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam	SM:A549_CTRL_NIPBL_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/bwa_index/Homo_sapiens.hg19.fa \
  trim/A549_CTRL_NIPBL_rep1/HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/A549_CTRL_NIPBL_rep1/HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam/HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.d588716ad5c0f58d54bc41b8ec6e96f5.mugqic.done
)
bwa_mem_picard_sort_sam_15_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_15_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_16_JOB_ID: bwa_mem_picard_sort_sam.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_16_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.49d89140d128b5c3eb972148b427ee6f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.49d89140d128b5c3eb972148b427ee6f.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/A549_CTRL_MED1_rep1/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam	SM:A549_CTRL_MED1_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/bwa_index/Homo_sapiens.hg19.fa \
  trim/A549_CTRL_MED1_rep1/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/A549_CTRL_MED1_rep1/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.49d89140d128b5c3eb972148b427ee6f.mugqic.done
)
bwa_mem_picard_sort_sam_16_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_16_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_17_JOB_ID: bwa_mem_picard_sort_sam_report
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam_report
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_1_JOB_ID:$bwa_mem_picard_sort_sam_2_JOB_ID:$bwa_mem_picard_sort_sam_3_JOB_ID:$bwa_mem_picard_sort_sam_4_JOB_ID:$bwa_mem_picard_sort_sam_5_JOB_ID:$bwa_mem_picard_sort_sam_6_JOB_ID:$bwa_mem_picard_sort_sam_7_JOB_ID:$bwa_mem_picard_sort_sam_8_JOB_ID:$bwa_mem_picard_sort_sam_9_JOB_ID:$bwa_mem_picard_sort_sam_10_JOB_ID:$bwa_mem_picard_sort_sam_11_JOB_ID:$bwa_mem_picard_sort_sam_12_JOB_ID:$bwa_mem_picard_sort_sam_13_JOB_ID:$bwa_mem_picard_sort_sam_14_JOB_ID:$bwa_mem_picard_sort_sam_15_JOB_ID:$bwa_mem_picard_sort_sam_16_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam_report.ccf108459067256486a8f50c1a551880.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam_report.ccf108459067256486a8f50c1a551880.mugqic.done'
module load mugqic/pandoc/1.15.2 && \
mkdir -p report && \
pandoc --to=markdown \
  --template /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/DnaSeq.bwa_mem_picard_sort_sam.md \
  --variable scientific_name="Homo_sapiens" \
  --variable assembly="hg19" \
  /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/DnaSeq.bwa_mem_picard_sort_sam.md \
  > report/DnaSeq.bwa_mem_picard_sort_sam.md
bwa_mem_picard_sort_sam_report.ccf108459067256486a8f50c1a551880.mugqic.done
)
bwa_mem_picard_sort_sam_17_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_17_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: samtools_view_filter
#-------------------------------------------------------------------------------
STEP=samtools_view_filter
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_1_JOB_ID: samtools_view_filter.HI.1613.003.Index_1.A549_Dex_CDK9_MF_ChIP1_4_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.1613.003.Index_1.A549_Dex_CDK9_MF_ChIP1_4_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_1_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.1613.003.Index_1.A549_Dex_CDK9_MF_ChIP1_4_rep1.bam.3e1d707c5123fa47709068b4e818c3b2.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.1613.003.Index_1.A549_Dex_CDK9_MF_ChIP1_4_rep1.bam.3e1d707c5123fa47709068b4e818c3b2.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/A549_DEX_CDK9_rep1/HI.1613.003.Index_1.A549_Dex_CDK9_MF_ChIP1_4_rep1.bam/HI.1613.003.Index_1.A549_Dex_CDK9_MF_ChIP1_4_rep1.bam.sorted.bam \
  > alignment/A549_DEX_CDK9_rep1/HI.1613.003.Index_1.A549_Dex_CDK9_MF_ChIP1_4_rep1.bam/HI.1613.003.Index_1.A549_Dex_CDK9_MF_ChIP1_4_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.1613.003.Index_1.A549_Dex_CDK9_MF_ChIP1_4_rep1.bam.3e1d707c5123fa47709068b4e818c3b2.mugqic.done
)
samtools_view_filter_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_2_JOB_ID: samtools_view_filter.HI.1613.003.Index_12.A549_Dex_MED1_MF_ChIP1_2_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.1613.003.Index_12.A549_Dex_MED1_MF_ChIP1_2_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_2_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.1613.003.Index_12.A549_Dex_MED1_MF_ChIP1_2_rep1.bam.ba6a8605bf68b44c494a4b6ce3aa80d5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.1613.003.Index_12.A549_Dex_MED1_MF_ChIP1_2_rep1.bam.ba6a8605bf68b44c494a4b6ce3aa80d5.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/A549_DEX_MED1_rep1/HI.1613.003.Index_12.A549_Dex_MED1_MF_ChIP1_2_rep1.bam/HI.1613.003.Index_12.A549_Dex_MED1_MF_ChIP1_2_rep1.bam.sorted.bam \
  > alignment/A549_DEX_MED1_rep1/HI.1613.003.Index_12.A549_Dex_MED1_MF_ChIP1_2_rep1.bam/HI.1613.003.Index_12.A549_Dex_MED1_MF_ChIP1_2_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.1613.003.Index_12.A549_Dex_MED1_MF_ChIP1_2_rep1.bam.ba6a8605bf68b44c494a4b6ce3aa80d5.mugqic.done
)
samtools_view_filter_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_3_JOB_ID: samtools_view_filter.HI.1613.003.Index_19.A549_Dex_SMC1_MF_ChIP1_3_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.1613.003.Index_19.A549_Dex_SMC1_MF_ChIP1_3_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_3_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.1613.003.Index_19.A549_Dex_SMC1_MF_ChIP1_3_rep1.bam.30ecccf09ab25b604c3cb8404efd7ff4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.1613.003.Index_19.A549_Dex_SMC1_MF_ChIP1_3_rep1.bam.30ecccf09ab25b604c3cb8404efd7ff4.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/A549_DEX_BRD4_rep1/HI.1613.003.Index_19.A549_Dex_SMC1_MF_ChIP1_3_rep1.bam/HI.1613.003.Index_19.A549_Dex_SMC1_MF_ChIP1_3_rep1.bam.sorted.bam \
  > alignment/A549_DEX_BRD4_rep1/HI.1613.003.Index_19.A549_Dex_SMC1_MF_ChIP1_3_rep1.bam/HI.1613.003.Index_19.A549_Dex_SMC1_MF_ChIP1_3_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.1613.003.Index_19.A549_Dex_SMC1_MF_ChIP1_3_rep1.bam.30ecccf09ab25b604c3cb8404efd7ff4.mugqic.done
)
samtools_view_filter_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_4_JOB_ID: samtools_view_filter.HI.1613.003.Index_4.A549_Dex_WCE_MF_ChIP1_8_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.1613.003.Index_4.A549_Dex_WCE_MF_ChIP1_8_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_4_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.1613.003.Index_4.A549_Dex_WCE_MF_ChIP1_8_rep1.bam.85ad48df33c57f9f72c054acc18867a8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.1613.003.Index_4.A549_Dex_WCE_MF_ChIP1_8_rep1.bam.85ad48df33c57f9f72c054acc18867a8.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/A549_DEX_WCE_rep1/HI.1613.003.Index_4.A549_Dex_WCE_MF_ChIP1_8_rep1.bam/HI.1613.003.Index_4.A549_Dex_WCE_MF_ChIP1_8_rep1.bam.sorted.bam \
  > alignment/A549_DEX_WCE_rep1/HI.1613.003.Index_4.A549_Dex_WCE_MF_ChIP1_8_rep1.bam/HI.1613.003.Index_4.A549_Dex_WCE_MF_ChIP1_8_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.1613.003.Index_4.A549_Dex_WCE_MF_ChIP1_8_rep1.bam.85ad48df33c57f9f72c054acc18867a8.mugqic.done
)
samtools_view_filter_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_5_JOB_ID: samtools_view_filter.HI.1613.003.Index_5.A549_Dex_NIPBL_MF_ChIP1_1_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.1613.003.Index_5.A549_Dex_NIPBL_MF_ChIP1_1_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_5_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.1613.003.Index_5.A549_Dex_NIPBL_MF_ChIP1_1_rep1.bam.8b03b45b867b14b801a504cb5b7559de.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.1613.003.Index_5.A549_Dex_NIPBL_MF_ChIP1_1_rep1.bam.8b03b45b867b14b801a504cb5b7559de.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/A549_DEX_NIPBL_rep1/HI.1613.003.Index_5.A549_Dex_NIPBL_MF_ChIP1_1_rep1.bam/HI.1613.003.Index_5.A549_Dex_NIPBL_MF_ChIP1_1_rep1.bam.sorted.bam \
  > alignment/A549_DEX_NIPBL_rep1/HI.1613.003.Index_5.A549_Dex_NIPBL_MF_ChIP1_1_rep1.bam/HI.1613.003.Index_5.A549_Dex_NIPBL_MF_ChIP1_1_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.1613.003.Index_5.A549_Dex_NIPBL_MF_ChIP1_1_rep1.bam.8b03b45b867b14b801a504cb5b7559de.mugqic.done
)
samtools_view_filter_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_6_JOB_ID: samtools_view_filter.HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_6_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.b0e841ee34f2eb16e40809cc42ad36d9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.b0e841ee34f2eb16e40809cc42ad36d9.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/A549_CTRL_SMC1A_rep1/HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam/HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.sorted.bam \
  > alignment/A549_CTRL_SMC1A_rep1/HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam/HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.b0e841ee34f2eb16e40809cc42ad36d9.mugqic.done
)
samtools_view_filter_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_7_JOB_ID: samtools_view_filter.HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_7_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam.2f5933ae79487b5e7676f7b8cd5bfcc5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam.2f5933ae79487b5e7676f7b8cd5bfcc5.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/A549_CTRL_WCE_rep1/HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam/HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam.sorted.bam \
  > alignment/A549_CTRL_WCE_rep1/HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam/HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam.2f5933ae79487b5e7676f7b8cd5bfcc5.mugqic.done
)
samtools_view_filter_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_8_JOB_ID: samtools_view_filter.HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_8_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.96934114352ce56c8bed7b9c9a5295f4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.96934114352ce56c8bed7b9c9a5295f4.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/A549_CTRL_CDK9_rep1/HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam/HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.sorted.bam \
  > alignment/A549_CTRL_CDK9_rep1/HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam/HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.96934114352ce56c8bed7b9c9a5295f4.mugqic.done
)
samtools_view_filter_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_9_JOB_ID: samtools_view_filter.HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_9_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.2625a6d17ef75e85c27bcfdf10162108.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.2625a6d17ef75e85c27bcfdf10162108.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/A549_CTRL_MED1_rep1/HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam/HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.sorted.bam \
  > alignment/A549_CTRL_MED1_rep1/HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam/HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.2625a6d17ef75e85c27bcfdf10162108.mugqic.done
)
samtools_view_filter_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_10_JOB_ID: samtools_view_filter.HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_10_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam.135dad8fa59620817f06d155c5723d92.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam.135dad8fa59620817f06d155c5723d92.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/A549_CTRL_NIPBL_rep1/HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam/HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam.sorted.bam \
  > alignment/A549_CTRL_NIPBL_rep1/HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam/HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam.135dad8fa59620817f06d155c5723d92.mugqic.done
)
samtools_view_filter_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_11_JOB_ID: samtools_view_filter.HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_11_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam.7d044d15571ce51808eb20c4e680edfc.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam.7d044d15571ce51808eb20c4e680edfc.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/A549_DEX_BRD4_rep1/HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam/HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam.sorted.bam \
  > alignment/A549_DEX_BRD4_rep1/HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam/HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam.7d044d15571ce51808eb20c4e680edfc.mugqic.done
)
samtools_view_filter_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_12_JOB_ID: samtools_view_filter.HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_12_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.2ee4205352822ad3d64cf3702f32c3e6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.2ee4205352822ad3d64cf3702f32c3e6.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/A549_CTRL_BRD4_rep1/HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam/HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.sorted.bam \
  > alignment/A549_CTRL_BRD4_rep1/HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam/HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.2ee4205352822ad3d64cf3702f32c3e6.mugqic.done
)
samtools_view_filter_12_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_13_JOB_ID: samtools_view_filter.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_13_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.38bd9e8611b7aaa10e8b55d369b37497.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.38bd9e8611b7aaa10e8b55d369b37497.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/A549_CTRL_SMC1A_rep1/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.sorted.bam \
  > alignment/A549_CTRL_SMC1A_rep1/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.38bd9e8611b7aaa10e8b55d369b37497.mugqic.done
)
samtools_view_filter_13_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_13_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_14_JOB_ID: samtools_view_filter.HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_14_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam.07388b169f0c55e103cc2620c0a3154f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam.07388b169f0c55e103cc2620c0a3154f.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/A549_CTRL_WCE_rep1/HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam/HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam.sorted.bam \
  > alignment/A549_CTRL_WCE_rep1/HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam/HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam.07388b169f0c55e103cc2620c0a3154f.mugqic.done
)
samtools_view_filter_14_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_14_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_15_JOB_ID: samtools_view_filter.HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_15_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.f5f1450baddebc5e22ad2892a71bd0e9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.f5f1450baddebc5e22ad2892a71bd0e9.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/A549_CTRL_NIPBL_rep1/HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam/HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.sorted.bam \
  > alignment/A549_CTRL_NIPBL_rep1/HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam/HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.f5f1450baddebc5e22ad2892a71bd0e9.mugqic.done
)
samtools_view_filter_15_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_15_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_16_JOB_ID: samtools_view_filter.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_16_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.f8b99fe53ccc0508c21526388cd451f0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.f8b99fe53ccc0508c21526388cd451f0.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/A549_CTRL_MED1_rep1/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.sorted.bam \
  > alignment/A549_CTRL_MED1_rep1/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.f8b99fe53ccc0508c21526388cd451f0.mugqic.done
)
samtools_view_filter_16_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_16_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_17_JOB_ID: samtools_view_filter_report
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter_report
JOB_DEPENDENCIES=$samtools_view_filter_1_JOB_ID:$samtools_view_filter_2_JOB_ID:$samtools_view_filter_3_JOB_ID:$samtools_view_filter_4_JOB_ID:$samtools_view_filter_5_JOB_ID:$samtools_view_filter_6_JOB_ID:$samtools_view_filter_7_JOB_ID:$samtools_view_filter_8_JOB_ID:$samtools_view_filter_9_JOB_ID:$samtools_view_filter_10_JOB_ID:$samtools_view_filter_11_JOB_ID:$samtools_view_filter_12_JOB_ID:$samtools_view_filter_13_JOB_ID:$samtools_view_filter_14_JOB_ID:$samtools_view_filter_15_JOB_ID:$samtools_view_filter_16_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter_report.b5cac64c2a4ab579900138ddfbed49a6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter_report.b5cac64c2a4ab579900138ddfbed49a6.mugqic.done'
module load mugqic/pandoc/1.15.2 && \
mkdir -p report && \
pandoc --to=markdown \
  --template /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/ChipSeq.samtools_view_filter.md \
  --variable min_mapq="20" \
  /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/ChipSeq.samtools_view_filter.md \
  > report/ChipSeq.samtools_view_filter.md
samtools_view_filter_report.b5cac64c2a4ab579900138ddfbed49a6.mugqic.done
)
samtools_view_filter_17_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_17_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: picard_merge_sam_files
#-------------------------------------------------------------------------------
STEP=picard_merge_sam_files
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_1_JOB_ID: symlink_readset_sample_bam.A549_DEX_CDK9_rep1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.A549_DEX_CDK9_rep1
JOB_DEPENDENCIES=$samtools_view_filter_1_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.A549_DEX_CDK9_rep1.3447c9fbb2a9ecc02d966ba20ebcdc20.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.A549_DEX_CDK9_rep1.3447c9fbb2a9ecc02d966ba20ebcdc20.mugqic.done'
mkdir -p alignment/A549_DEX_CDK9_rep1 && \
ln -s -f HI.1613.003.Index_1.A549_Dex_CDK9_MF_ChIP1_4_rep1.bam/HI.1613.003.Index_1.A549_Dex_CDK9_MF_ChIP1_4_rep1.bam.sorted.filtered.bam alignment/A549_DEX_CDK9_rep1/A549_DEX_CDK9_rep1.merged.bam
symlink_readset_sample_bam.A549_DEX_CDK9_rep1.3447c9fbb2a9ecc02d966ba20ebcdc20.mugqic.done
)
picard_merge_sam_files_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_2_JOB_ID: symlink_readset_sample_bam.A549_DEX_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.A549_DEX_MED1_rep1
JOB_DEPENDENCIES=$samtools_view_filter_2_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.A549_DEX_MED1_rep1.7d4a74f099219c155aac16e5ae6009e3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.A549_DEX_MED1_rep1.7d4a74f099219c155aac16e5ae6009e3.mugqic.done'
mkdir -p alignment/A549_DEX_MED1_rep1 && \
ln -s -f HI.1613.003.Index_12.A549_Dex_MED1_MF_ChIP1_2_rep1.bam/HI.1613.003.Index_12.A549_Dex_MED1_MF_ChIP1_2_rep1.bam.sorted.filtered.bam alignment/A549_DEX_MED1_rep1/A549_DEX_MED1_rep1.merged.bam
symlink_readset_sample_bam.A549_DEX_MED1_rep1.7d4a74f099219c155aac16e5ae6009e3.mugqic.done
)
picard_merge_sam_files_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_3_JOB_ID: picard_merge_sam_files.A549_DEX_BRD4_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_merge_sam_files.A549_DEX_BRD4_rep1
JOB_DEPENDENCIES=$samtools_view_filter_3_JOB_ID:$samtools_view_filter_11_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/picard_merge_sam_files.A549_DEX_BRD4_rep1.aa936666da78b8079526b016ab15e21e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_merge_sam_files.A549_DEX_BRD4_rep1.aa936666da78b8079526b016ab15e21e.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/A549_DEX_BRD4_rep1 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx1700M -jar $PICARD_HOME/MergeSamFiles.jar \
  VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/A549_DEX_BRD4_rep1/HI.1613.003.Index_19.A549_Dex_SMC1_MF_ChIP1_3_rep1.bam/HI.1613.003.Index_19.A549_Dex_SMC1_MF_ChIP1_3_rep1.bam.sorted.filtered.bam \
  INPUT=alignment/A549_DEX_BRD4_rep1/HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam/HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam.sorted.filtered.bam \
  OUTPUT=alignment/A549_DEX_BRD4_rep1/A549_DEX_BRD4_rep1.merged.bam \
  MAX_RECORDS_IN_RAM=250000
picard_merge_sam_files.A549_DEX_BRD4_rep1.aa936666da78b8079526b016ab15e21e.mugqic.done
)
picard_merge_sam_files_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=35:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_4_JOB_ID: symlink_readset_sample_bam.A549_DEX_WCE_rep1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.A549_DEX_WCE_rep1
JOB_DEPENDENCIES=$samtools_view_filter_4_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.A549_DEX_WCE_rep1.2dbb1807fbf5b3dc1f7ef8ad42259788.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.A549_DEX_WCE_rep1.2dbb1807fbf5b3dc1f7ef8ad42259788.mugqic.done'
mkdir -p alignment/A549_DEX_WCE_rep1 && \
ln -s -f HI.1613.003.Index_4.A549_Dex_WCE_MF_ChIP1_8_rep1.bam/HI.1613.003.Index_4.A549_Dex_WCE_MF_ChIP1_8_rep1.bam.sorted.filtered.bam alignment/A549_DEX_WCE_rep1/A549_DEX_WCE_rep1.merged.bam
symlink_readset_sample_bam.A549_DEX_WCE_rep1.2dbb1807fbf5b3dc1f7ef8ad42259788.mugqic.done
)
picard_merge_sam_files_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_5_JOB_ID: symlink_readset_sample_bam.A549_DEX_NIPBL_rep1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.A549_DEX_NIPBL_rep1
JOB_DEPENDENCIES=$samtools_view_filter_5_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.A549_DEX_NIPBL_rep1.7facbcf391e5d1f2b87fde11c5113108.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.A549_DEX_NIPBL_rep1.7facbcf391e5d1f2b87fde11c5113108.mugqic.done'
mkdir -p alignment/A549_DEX_NIPBL_rep1 && \
ln -s -f HI.1613.003.Index_5.A549_Dex_NIPBL_MF_ChIP1_1_rep1.bam/HI.1613.003.Index_5.A549_Dex_NIPBL_MF_ChIP1_1_rep1.bam.sorted.filtered.bam alignment/A549_DEX_NIPBL_rep1/A549_DEX_NIPBL_rep1.merged.bam
symlink_readset_sample_bam.A549_DEX_NIPBL_rep1.7facbcf391e5d1f2b87fde11c5113108.mugqic.done
)
picard_merge_sam_files_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_6_JOB_ID: picard_merge_sam_files.A549_CTRL_SMC1A_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_merge_sam_files.A549_CTRL_SMC1A_rep1
JOB_DEPENDENCIES=$samtools_view_filter_6_JOB_ID:$samtools_view_filter_13_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/picard_merge_sam_files.A549_CTRL_SMC1A_rep1.03abe1e5b1422d6fd986837aa75a4baf.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_merge_sam_files.A549_CTRL_SMC1A_rep1.03abe1e5b1422d6fd986837aa75a4baf.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/A549_CTRL_SMC1A_rep1 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx1700M -jar $PICARD_HOME/MergeSamFiles.jar \
  VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/A549_CTRL_SMC1A_rep1/HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam/HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.sorted.filtered.bam \
  INPUT=alignment/A549_CTRL_SMC1A_rep1/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.sorted.filtered.bam \
  OUTPUT=alignment/A549_CTRL_SMC1A_rep1/A549_CTRL_SMC1A_rep1.merged.bam \
  MAX_RECORDS_IN_RAM=250000
picard_merge_sam_files.A549_CTRL_SMC1A_rep1.03abe1e5b1422d6fd986837aa75a4baf.mugqic.done
)
picard_merge_sam_files_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=35:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_7_JOB_ID: picard_merge_sam_files.A549_CTRL_WCE_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_merge_sam_files.A549_CTRL_WCE_rep1
JOB_DEPENDENCIES=$samtools_view_filter_7_JOB_ID:$samtools_view_filter_14_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/picard_merge_sam_files.A549_CTRL_WCE_rep1.1f7700b1962f6cf43ddc6ec69c4a5ee2.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_merge_sam_files.A549_CTRL_WCE_rep1.1f7700b1962f6cf43ddc6ec69c4a5ee2.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/A549_CTRL_WCE_rep1 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx1700M -jar $PICARD_HOME/MergeSamFiles.jar \
  VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/A549_CTRL_WCE_rep1/HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam/HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam.sorted.filtered.bam \
  INPUT=alignment/A549_CTRL_WCE_rep1/HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam/HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam.sorted.filtered.bam \
  OUTPUT=alignment/A549_CTRL_WCE_rep1/A549_CTRL_WCE_rep1.merged.bam \
  MAX_RECORDS_IN_RAM=250000
picard_merge_sam_files.A549_CTRL_WCE_rep1.1f7700b1962f6cf43ddc6ec69c4a5ee2.mugqic.done
)
picard_merge_sam_files_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=35:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_8_JOB_ID: symlink_readset_sample_bam.A549_CTRL_CDK9_rep1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.A549_CTRL_CDK9_rep1
JOB_DEPENDENCIES=$samtools_view_filter_8_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.A549_CTRL_CDK9_rep1.e37c74f70c72e58da1e581952f7395b1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.A549_CTRL_CDK9_rep1.e37c74f70c72e58da1e581952f7395b1.mugqic.done'
mkdir -p alignment/A549_CTRL_CDK9_rep1 && \
ln -s -f HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam/HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.sorted.filtered.bam alignment/A549_CTRL_CDK9_rep1/A549_CTRL_CDK9_rep1.merged.bam
symlink_readset_sample_bam.A549_CTRL_CDK9_rep1.e37c74f70c72e58da1e581952f7395b1.mugqic.done
)
picard_merge_sam_files_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_9_JOB_ID: picard_merge_sam_files.A549_CTRL_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_merge_sam_files.A549_CTRL_MED1_rep1
JOB_DEPENDENCIES=$samtools_view_filter_9_JOB_ID:$samtools_view_filter_16_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/picard_merge_sam_files.A549_CTRL_MED1_rep1.e260f599c1419f9d46d2fcbaa6284eb7.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_merge_sam_files.A549_CTRL_MED1_rep1.e260f599c1419f9d46d2fcbaa6284eb7.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/A549_CTRL_MED1_rep1 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx1700M -jar $PICARD_HOME/MergeSamFiles.jar \
  VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/A549_CTRL_MED1_rep1/HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam/HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.sorted.filtered.bam \
  INPUT=alignment/A549_CTRL_MED1_rep1/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.sorted.filtered.bam \
  OUTPUT=alignment/A549_CTRL_MED1_rep1/A549_CTRL_MED1_rep1.merged.bam \
  MAX_RECORDS_IN_RAM=250000
picard_merge_sam_files.A549_CTRL_MED1_rep1.e260f599c1419f9d46d2fcbaa6284eb7.mugqic.done
)
picard_merge_sam_files_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=35:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_10_JOB_ID: picard_merge_sam_files.A549_CTRL_NIPBL_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_merge_sam_files.A549_CTRL_NIPBL_rep1
JOB_DEPENDENCIES=$samtools_view_filter_10_JOB_ID:$samtools_view_filter_15_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/picard_merge_sam_files.A549_CTRL_NIPBL_rep1.4fd1b7098d1ba73b7da49b5f41db9c26.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_merge_sam_files.A549_CTRL_NIPBL_rep1.4fd1b7098d1ba73b7da49b5f41db9c26.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/A549_CTRL_NIPBL_rep1 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx1700M -jar $PICARD_HOME/MergeSamFiles.jar \
  VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/A549_CTRL_NIPBL_rep1/HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam/HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam.sorted.filtered.bam \
  INPUT=alignment/A549_CTRL_NIPBL_rep1/HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam/HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.sorted.filtered.bam \
  OUTPUT=alignment/A549_CTRL_NIPBL_rep1/A549_CTRL_NIPBL_rep1.merged.bam \
  MAX_RECORDS_IN_RAM=250000
picard_merge_sam_files.A549_CTRL_NIPBL_rep1.4fd1b7098d1ba73b7da49b5f41db9c26.mugqic.done
)
picard_merge_sam_files_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=35:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_11_JOB_ID: symlink_readset_sample_bam.A549_CTRL_BRD4_rep1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.A549_CTRL_BRD4_rep1
JOB_DEPENDENCIES=$samtools_view_filter_12_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.A549_CTRL_BRD4_rep1.122782b49de4447e23add94e5b7d5d28.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.A549_CTRL_BRD4_rep1.122782b49de4447e23add94e5b7d5d28.mugqic.done'
mkdir -p alignment/A549_CTRL_BRD4_rep1 && \
ln -s -f HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam/HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.sorted.filtered.bam alignment/A549_CTRL_BRD4_rep1/A549_CTRL_BRD4_rep1.merged.bam
symlink_readset_sample_bam.A549_CTRL_BRD4_rep1.122782b49de4447e23add94e5b7d5d28.mugqic.done
)
picard_merge_sam_files_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: picard_mark_duplicates
#-------------------------------------------------------------------------------
STEP=picard_mark_duplicates
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_1_JOB_ID: picard_mark_duplicates.A549_DEX_CDK9_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.A549_DEX_CDK9_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_1_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.A549_DEX_CDK9_rep1.f5714b2f0ffb03277e2922a5a1f0ad7e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.A549_DEX_CDK9_rep1.f5714b2f0ffb03277e2922a5a1f0ad7e.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/A549_DEX_CDK9_rep1/A549_DEX_CDK9_rep1.merged.bam \
  OUTPUT=alignment/A549_DEX_CDK9_rep1/A549_DEX_CDK9_rep1.sorted.dup.bam \
  METRICS_FILE=alignment/A549_DEX_CDK9_rep1/A549_DEX_CDK9_rep1.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.A549_DEX_CDK9_rep1.f5714b2f0ffb03277e2922a5a1f0ad7e.mugqic.done
)
picard_mark_duplicates_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_2_JOB_ID: picard_mark_duplicates.A549_DEX_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.A549_DEX_MED1_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_2_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.A549_DEX_MED1_rep1.d73c0ee78f800e438f17a431d740a592.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.A549_DEX_MED1_rep1.d73c0ee78f800e438f17a431d740a592.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/A549_DEX_MED1_rep1/A549_DEX_MED1_rep1.merged.bam \
  OUTPUT=alignment/A549_DEX_MED1_rep1/A549_DEX_MED1_rep1.sorted.dup.bam \
  METRICS_FILE=alignment/A549_DEX_MED1_rep1/A549_DEX_MED1_rep1.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.A549_DEX_MED1_rep1.d73c0ee78f800e438f17a431d740a592.mugqic.done
)
picard_mark_duplicates_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_3_JOB_ID: picard_mark_duplicates.A549_DEX_BRD4_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.A549_DEX_BRD4_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_3_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.A549_DEX_BRD4_rep1.e6bdc5d06caf3788c703e49ebb08697a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.A549_DEX_BRD4_rep1.e6bdc5d06caf3788c703e49ebb08697a.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/A549_DEX_BRD4_rep1/A549_DEX_BRD4_rep1.merged.bam \
  OUTPUT=alignment/A549_DEX_BRD4_rep1/A549_DEX_BRD4_rep1.sorted.dup.bam \
  METRICS_FILE=alignment/A549_DEX_BRD4_rep1/A549_DEX_BRD4_rep1.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.A549_DEX_BRD4_rep1.e6bdc5d06caf3788c703e49ebb08697a.mugqic.done
)
picard_mark_duplicates_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_4_JOB_ID: picard_mark_duplicates.A549_DEX_WCE_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.A549_DEX_WCE_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_4_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.A549_DEX_WCE_rep1.8119f6e610d9a14eb827885664d13265.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.A549_DEX_WCE_rep1.8119f6e610d9a14eb827885664d13265.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/A549_DEX_WCE_rep1/A549_DEX_WCE_rep1.merged.bam \
  OUTPUT=alignment/A549_DEX_WCE_rep1/A549_DEX_WCE_rep1.sorted.dup.bam \
  METRICS_FILE=alignment/A549_DEX_WCE_rep1/A549_DEX_WCE_rep1.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.A549_DEX_WCE_rep1.8119f6e610d9a14eb827885664d13265.mugqic.done
)
picard_mark_duplicates_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_5_JOB_ID: picard_mark_duplicates.A549_DEX_NIPBL_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.A549_DEX_NIPBL_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_5_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.A549_DEX_NIPBL_rep1.8478805c67a35486dd1df36cddba2f5b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.A549_DEX_NIPBL_rep1.8478805c67a35486dd1df36cddba2f5b.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/A549_DEX_NIPBL_rep1/A549_DEX_NIPBL_rep1.merged.bam \
  OUTPUT=alignment/A549_DEX_NIPBL_rep1/A549_DEX_NIPBL_rep1.sorted.dup.bam \
  METRICS_FILE=alignment/A549_DEX_NIPBL_rep1/A549_DEX_NIPBL_rep1.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.A549_DEX_NIPBL_rep1.8478805c67a35486dd1df36cddba2f5b.mugqic.done
)
picard_mark_duplicates_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_6_JOB_ID: picard_mark_duplicates.A549_CTRL_SMC1A_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.A549_CTRL_SMC1A_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_6_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.A549_CTRL_SMC1A_rep1.6993e854018e8adf50539dde283e2e22.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.A549_CTRL_SMC1A_rep1.6993e854018e8adf50539dde283e2e22.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/A549_CTRL_SMC1A_rep1/A549_CTRL_SMC1A_rep1.merged.bam \
  OUTPUT=alignment/A549_CTRL_SMC1A_rep1/A549_CTRL_SMC1A_rep1.sorted.dup.bam \
  METRICS_FILE=alignment/A549_CTRL_SMC1A_rep1/A549_CTRL_SMC1A_rep1.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.A549_CTRL_SMC1A_rep1.6993e854018e8adf50539dde283e2e22.mugqic.done
)
picard_mark_duplicates_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_7_JOB_ID: picard_mark_duplicates.A549_CTRL_WCE_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.A549_CTRL_WCE_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_7_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.A549_CTRL_WCE_rep1.cf95d58d7f2eb0af41907dae0769b507.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.A549_CTRL_WCE_rep1.cf95d58d7f2eb0af41907dae0769b507.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/A549_CTRL_WCE_rep1/A549_CTRL_WCE_rep1.merged.bam \
  OUTPUT=alignment/A549_CTRL_WCE_rep1/A549_CTRL_WCE_rep1.sorted.dup.bam \
  METRICS_FILE=alignment/A549_CTRL_WCE_rep1/A549_CTRL_WCE_rep1.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.A549_CTRL_WCE_rep1.cf95d58d7f2eb0af41907dae0769b507.mugqic.done
)
picard_mark_duplicates_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_8_JOB_ID: picard_mark_duplicates.A549_CTRL_CDK9_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.A549_CTRL_CDK9_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_8_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.A549_CTRL_CDK9_rep1.0116686c4d68835a723d7cfa67b81910.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.A549_CTRL_CDK9_rep1.0116686c4d68835a723d7cfa67b81910.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/A549_CTRL_CDK9_rep1/A549_CTRL_CDK9_rep1.merged.bam \
  OUTPUT=alignment/A549_CTRL_CDK9_rep1/A549_CTRL_CDK9_rep1.sorted.dup.bam \
  METRICS_FILE=alignment/A549_CTRL_CDK9_rep1/A549_CTRL_CDK9_rep1.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.A549_CTRL_CDK9_rep1.0116686c4d68835a723d7cfa67b81910.mugqic.done
)
picard_mark_duplicates_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_9_JOB_ID: picard_mark_duplicates.A549_CTRL_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.A549_CTRL_MED1_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_9_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.A549_CTRL_MED1_rep1.84b1db20b25f0858b7e39dbadaded85d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.A549_CTRL_MED1_rep1.84b1db20b25f0858b7e39dbadaded85d.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/A549_CTRL_MED1_rep1/A549_CTRL_MED1_rep1.merged.bam \
  OUTPUT=alignment/A549_CTRL_MED1_rep1/A549_CTRL_MED1_rep1.sorted.dup.bam \
  METRICS_FILE=alignment/A549_CTRL_MED1_rep1/A549_CTRL_MED1_rep1.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.A549_CTRL_MED1_rep1.84b1db20b25f0858b7e39dbadaded85d.mugqic.done
)
picard_mark_duplicates_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_10_JOB_ID: picard_mark_duplicates.A549_CTRL_NIPBL_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.A549_CTRL_NIPBL_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_10_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.A549_CTRL_NIPBL_rep1.ef0493545db1cb6bae52a33c6b22888e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.A549_CTRL_NIPBL_rep1.ef0493545db1cb6bae52a33c6b22888e.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/A549_CTRL_NIPBL_rep1/A549_CTRL_NIPBL_rep1.merged.bam \
  OUTPUT=alignment/A549_CTRL_NIPBL_rep1/A549_CTRL_NIPBL_rep1.sorted.dup.bam \
  METRICS_FILE=alignment/A549_CTRL_NIPBL_rep1/A549_CTRL_NIPBL_rep1.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.A549_CTRL_NIPBL_rep1.ef0493545db1cb6bae52a33c6b22888e.mugqic.done
)
picard_mark_duplicates_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_11_JOB_ID: picard_mark_duplicates.A549_CTRL_BRD4_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.A549_CTRL_BRD4_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_11_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.A549_CTRL_BRD4_rep1.b91fc383835a2112da1bf75596efd0cc.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.A549_CTRL_BRD4_rep1.b91fc383835a2112da1bf75596efd0cc.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/A549_CTRL_BRD4_rep1/A549_CTRL_BRD4_rep1.merged.bam \
  OUTPUT=alignment/A549_CTRL_BRD4_rep1/A549_CTRL_BRD4_rep1.sorted.dup.bam \
  METRICS_FILE=alignment/A549_CTRL_BRD4_rep1/A549_CTRL_BRD4_rep1.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.A549_CTRL_BRD4_rep1.b91fc383835a2112da1bf75596efd0cc.mugqic.done
)
picard_mark_duplicates_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_12_JOB_ID: picard_mark_duplicates_report
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates_report
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_2_JOB_ID:$picard_mark_duplicates_3_JOB_ID:$picard_mark_duplicates_4_JOB_ID:$picard_mark_duplicates_5_JOB_ID:$picard_mark_duplicates_6_JOB_ID:$picard_mark_duplicates_7_JOB_ID:$picard_mark_duplicates_8_JOB_ID:$picard_mark_duplicates_9_JOB_ID:$picard_mark_duplicates_10_JOB_ID:$picard_mark_duplicates_11_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates_report.c80ab57aaab7999422d809ae0583b4b5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates_report.c80ab57aaab7999422d809ae0583b4b5.mugqic.done'
mkdir -p report && \
cp \
  /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/ChipSeq.picard_mark_duplicates.md \
  report/ChipSeq.picard_mark_duplicates.md
picard_mark_duplicates_report.c80ab57aaab7999422d809ae0583b4b5.mugqic.done
)
picard_mark_duplicates_12_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: metrics
#-------------------------------------------------------------------------------
STEP=metrics
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: metrics_1_JOB_ID: metrics.flagstat
#-------------------------------------------------------------------------------
JOB_NAME=metrics.flagstat
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_2_JOB_ID:$picard_mark_duplicates_3_JOB_ID:$picard_mark_duplicates_4_JOB_ID:$picard_mark_duplicates_5_JOB_ID:$picard_mark_duplicates_6_JOB_ID:$picard_mark_duplicates_7_JOB_ID:$picard_mark_duplicates_8_JOB_ID:$picard_mark_duplicates_9_JOB_ID:$picard_mark_duplicates_10_JOB_ID:$picard_mark_duplicates_11_JOB_ID
JOB_DONE=job_output/metrics/metrics.flagstat.8b93f83bb70095fa0717cfe942629368.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'metrics.flagstat.8b93f83bb70095fa0717cfe942629368.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools flagstat \
  alignment/A549_DEX_CDK9_rep1/A549_DEX_CDK9_rep1.sorted.dup.bam \
  > alignment/A549_DEX_CDK9_rep1/A549_DEX_CDK9_rep1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/A549_DEX_MED1_rep1/A549_DEX_MED1_rep1.sorted.dup.bam \
  > alignment/A549_DEX_MED1_rep1/A549_DEX_MED1_rep1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/A549_DEX_BRD4_rep1/A549_DEX_BRD4_rep1.sorted.dup.bam \
  > alignment/A549_DEX_BRD4_rep1/A549_DEX_BRD4_rep1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/A549_DEX_WCE_rep1/A549_DEX_WCE_rep1.sorted.dup.bam \
  > alignment/A549_DEX_WCE_rep1/A549_DEX_WCE_rep1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/A549_DEX_NIPBL_rep1/A549_DEX_NIPBL_rep1.sorted.dup.bam \
  > alignment/A549_DEX_NIPBL_rep1/A549_DEX_NIPBL_rep1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/A549_CTRL_SMC1A_rep1/A549_CTRL_SMC1A_rep1.sorted.dup.bam \
  > alignment/A549_CTRL_SMC1A_rep1/A549_CTRL_SMC1A_rep1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/A549_CTRL_WCE_rep1/A549_CTRL_WCE_rep1.sorted.dup.bam \
  > alignment/A549_CTRL_WCE_rep1/A549_CTRL_WCE_rep1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/A549_CTRL_CDK9_rep1/A549_CTRL_CDK9_rep1.sorted.dup.bam \
  > alignment/A549_CTRL_CDK9_rep1/A549_CTRL_CDK9_rep1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/A549_CTRL_MED1_rep1/A549_CTRL_MED1_rep1.sorted.dup.bam \
  > alignment/A549_CTRL_MED1_rep1/A549_CTRL_MED1_rep1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/A549_CTRL_NIPBL_rep1/A549_CTRL_NIPBL_rep1.sorted.dup.bam \
  > alignment/A549_CTRL_NIPBL_rep1/A549_CTRL_NIPBL_rep1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/A549_CTRL_BRD4_rep1/A549_CTRL_BRD4_rep1.sorted.dup.bam \
  > alignment/A549_CTRL_BRD4_rep1/A549_CTRL_BRD4_rep1.sorted.dup.bam.flagstat
metrics.flagstat.8b93f83bb70095fa0717cfe942629368.mugqic.done
)
metrics_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$metrics_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: metrics_2_JOB_ID: metrics_report
#-------------------------------------------------------------------------------
JOB_NAME=metrics_report
JOB_DEPENDENCIES=$metrics_1_JOB_ID
JOB_DONE=job_output/metrics/metrics_report.be5c58bfc5691d63fd93a41c63843afc.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'metrics_report.be5c58bfc5691d63fd93a41c63843afc.mugqic.done'
module load mugqic/pandoc/1.15.2 && \
for sample in A549_DEX_CDK9_rep1 A549_DEX_MED1_rep1 A549_DEX_BRD4_rep1 A549_DEX_WCE_rep1 A549_DEX_NIPBL_rep1 A549_CTRL_SMC1A_rep1 A549_CTRL_WCE_rep1 A549_CTRL_CDK9_rep1 A549_CTRL_MED1_rep1 A549_CTRL_NIPBL_rep1 A549_CTRL_BRD4_rep1
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
  --template /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/ChipSeq.metrics.md \
  --variable trim_mem_sample_table="$trim_mem_sample_table" \
  /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/ChipSeq.metrics.md \
  > report/ChipSeq.metrics.md

metrics_report.be5c58bfc5691d63fd93a41c63843afc.mugqic.done
)
metrics_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$metrics_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: homer_make_tag_directory
#-------------------------------------------------------------------------------
STEP=homer_make_tag_directory
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_1_JOB_ID: homer_make_tag_directory.A549_DEX_CDK9_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.A549_DEX_CDK9_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.A549_DEX_CDK9_rep1.6535012dc035c6036fe1f4f0654a50f8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.A549_DEX_CDK9_rep1.6535012dc035c6036fe1f4f0654a50f8.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/A549_DEX_CDK9_rep1 \
  alignment/A549_DEX_CDK9_rep1/A549_DEX_CDK9_rep1.sorted.dup.bam \
  -checkGC -genome hg19
homer_make_tag_directory.A549_DEX_CDK9_rep1.6535012dc035c6036fe1f4f0654a50f8.mugqic.done
)
homer_make_tag_directory_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_2_JOB_ID: homer_make_tag_directory.A549_DEX_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.A549_DEX_MED1_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_2_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.A549_DEX_MED1_rep1.99ce93388b73a1a7dfd48620b2364742.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.A549_DEX_MED1_rep1.99ce93388b73a1a7dfd48620b2364742.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/A549_DEX_MED1_rep1 \
  alignment/A549_DEX_MED1_rep1/A549_DEX_MED1_rep1.sorted.dup.bam \
  -checkGC -genome hg19
homer_make_tag_directory.A549_DEX_MED1_rep1.99ce93388b73a1a7dfd48620b2364742.mugqic.done
)
homer_make_tag_directory_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_3_JOB_ID: homer_make_tag_directory.A549_DEX_BRD4_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.A549_DEX_BRD4_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_3_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.A549_DEX_BRD4_rep1.eaabf9dcaae2438f577493e4ba9a7f04.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.A549_DEX_BRD4_rep1.eaabf9dcaae2438f577493e4ba9a7f04.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/A549_DEX_BRD4_rep1 \
  alignment/A549_DEX_BRD4_rep1/A549_DEX_BRD4_rep1.sorted.dup.bam \
  -checkGC -genome hg19
homer_make_tag_directory.A549_DEX_BRD4_rep1.eaabf9dcaae2438f577493e4ba9a7f04.mugqic.done
)
homer_make_tag_directory_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_4_JOB_ID: homer_make_tag_directory.A549_DEX_WCE_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.A549_DEX_WCE_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_4_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.A549_DEX_WCE_rep1.6c47a6c5535695d17254d9e1120ccb26.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.A549_DEX_WCE_rep1.6c47a6c5535695d17254d9e1120ccb26.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/A549_DEX_WCE_rep1 \
  alignment/A549_DEX_WCE_rep1/A549_DEX_WCE_rep1.sorted.dup.bam \
  -checkGC -genome hg19
homer_make_tag_directory.A549_DEX_WCE_rep1.6c47a6c5535695d17254d9e1120ccb26.mugqic.done
)
homer_make_tag_directory_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_5_JOB_ID: homer_make_tag_directory.A549_DEX_NIPBL_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.A549_DEX_NIPBL_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_5_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.A549_DEX_NIPBL_rep1.954ef2b4b01905a0b80ef683d8cbf832.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.A549_DEX_NIPBL_rep1.954ef2b4b01905a0b80ef683d8cbf832.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/A549_DEX_NIPBL_rep1 \
  alignment/A549_DEX_NIPBL_rep1/A549_DEX_NIPBL_rep1.sorted.dup.bam \
  -checkGC -genome hg19
homer_make_tag_directory.A549_DEX_NIPBL_rep1.954ef2b4b01905a0b80ef683d8cbf832.mugqic.done
)
homer_make_tag_directory_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_6_JOB_ID: homer_make_tag_directory.A549_CTRL_SMC1A_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.A549_CTRL_SMC1A_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_6_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.A549_CTRL_SMC1A_rep1.2cbe0f483235ee177217c4daf606e189.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.A549_CTRL_SMC1A_rep1.2cbe0f483235ee177217c4daf606e189.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/A549_CTRL_SMC1A_rep1 \
  alignment/A549_CTRL_SMC1A_rep1/A549_CTRL_SMC1A_rep1.sorted.dup.bam \
  -checkGC -genome hg19
homer_make_tag_directory.A549_CTRL_SMC1A_rep1.2cbe0f483235ee177217c4daf606e189.mugqic.done
)
homer_make_tag_directory_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_7_JOB_ID: homer_make_tag_directory.A549_CTRL_WCE_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.A549_CTRL_WCE_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_7_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.A549_CTRL_WCE_rep1.f6c767605058fb49efb0726a622d1759.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.A549_CTRL_WCE_rep1.f6c767605058fb49efb0726a622d1759.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/A549_CTRL_WCE_rep1 \
  alignment/A549_CTRL_WCE_rep1/A549_CTRL_WCE_rep1.sorted.dup.bam \
  -checkGC -genome hg19
homer_make_tag_directory.A549_CTRL_WCE_rep1.f6c767605058fb49efb0726a622d1759.mugqic.done
)
homer_make_tag_directory_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_8_JOB_ID: homer_make_tag_directory.A549_CTRL_CDK9_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.A549_CTRL_CDK9_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_8_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.A549_CTRL_CDK9_rep1.0219d0b5dc6f20479b3d34b3b13c2aba.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.A549_CTRL_CDK9_rep1.0219d0b5dc6f20479b3d34b3b13c2aba.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/A549_CTRL_CDK9_rep1 \
  alignment/A549_CTRL_CDK9_rep1/A549_CTRL_CDK9_rep1.sorted.dup.bam \
  -checkGC -genome hg19
homer_make_tag_directory.A549_CTRL_CDK9_rep1.0219d0b5dc6f20479b3d34b3b13c2aba.mugqic.done
)
homer_make_tag_directory_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_9_JOB_ID: homer_make_tag_directory.A549_CTRL_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.A549_CTRL_MED1_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_9_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.A549_CTRL_MED1_rep1.114af6b33d97f0f99bb65115921f6fd3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.A549_CTRL_MED1_rep1.114af6b33d97f0f99bb65115921f6fd3.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/A549_CTRL_MED1_rep1 \
  alignment/A549_CTRL_MED1_rep1/A549_CTRL_MED1_rep1.sorted.dup.bam \
  -checkGC -genome hg19
homer_make_tag_directory.A549_CTRL_MED1_rep1.114af6b33d97f0f99bb65115921f6fd3.mugqic.done
)
homer_make_tag_directory_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_10_JOB_ID: homer_make_tag_directory.A549_CTRL_NIPBL_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.A549_CTRL_NIPBL_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_10_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.A549_CTRL_NIPBL_rep1.9c68ce060fb1de4e81a288844da7e160.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.A549_CTRL_NIPBL_rep1.9c68ce060fb1de4e81a288844da7e160.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/A549_CTRL_NIPBL_rep1 \
  alignment/A549_CTRL_NIPBL_rep1/A549_CTRL_NIPBL_rep1.sorted.dup.bam \
  -checkGC -genome hg19
homer_make_tag_directory.A549_CTRL_NIPBL_rep1.9c68ce060fb1de4e81a288844da7e160.mugqic.done
)
homer_make_tag_directory_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_11_JOB_ID: homer_make_tag_directory.A549_CTRL_BRD4_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.A549_CTRL_BRD4_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_11_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.A549_CTRL_BRD4_rep1.2f0eb96fabf6b853fc186b2d09b5b88d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.A549_CTRL_BRD4_rep1.2f0eb96fabf6b853fc186b2d09b5b88d.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/A549_CTRL_BRD4_rep1 \
  alignment/A549_CTRL_BRD4_rep1/A549_CTRL_BRD4_rep1.sorted.dup.bam \
  -checkGC -genome hg19
homer_make_tag_directory.A549_CTRL_BRD4_rep1.2f0eb96fabf6b853fc186b2d09b5b88d.mugqic.done
)
homer_make_tag_directory_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: qc_metrics
#-------------------------------------------------------------------------------
STEP=qc_metrics
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: qc_metrics_1_JOB_ID: qc_plots_R
#-------------------------------------------------------------------------------
JOB_NAME=qc_plots_R
JOB_DEPENDENCIES=$homer_make_tag_directory_1_JOB_ID:$homer_make_tag_directory_2_JOB_ID:$homer_make_tag_directory_3_JOB_ID:$homer_make_tag_directory_4_JOB_ID:$homer_make_tag_directory_5_JOB_ID:$homer_make_tag_directory_6_JOB_ID:$homer_make_tag_directory_7_JOB_ID:$homer_make_tag_directory_8_JOB_ID:$homer_make_tag_directory_9_JOB_ID:$homer_make_tag_directory_10_JOB_ID:$homer_make_tag_directory_11_JOB_ID
JOB_DONE=job_output/qc_metrics/qc_plots_R.947ca4b9fd1de01b4fbbe04f27d2d113.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'qc_plots_R.947ca4b9fd1de01b4fbbe04f27d2d113.mugqic.done'
module load mugqic/mugqic_tools/2.1.5 mugqic/R_Bioconductor/3.2.3_3.2 && \
mkdir -p graphs && \
Rscript $R_TOOLS/chipSeqGenerateQCMetrics.R \
  ../../raw/chip-seq/design.txt \
  /gs/scratch/efournier/CofactorHR/A549/output/chip-pipeline && \
cp /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/ChipSeq.qc_metrics.md report/ChipSeq.qc_metrics.md && \
for sample in A549_DEX_CDK9_rep1 A549_DEX_MED1_rep1 A549_DEX_BRD4_rep1 A549_DEX_WCE_rep1 A549_DEX_NIPBL_rep1 A549_CTRL_SMC1A_rep1 A549_CTRL_WCE_rep1 A549_CTRL_CDK9_rep1 A549_CTRL_MED1_rep1 A549_CTRL_NIPBL_rep1 A549_CTRL_BRD4_rep1
do
  cp --parents graphs/${sample}_QC_Metrics.ps report/
  convert -rotate 90 graphs/${sample}_QC_Metrics.ps report/graphs/${sample}_QC_Metrics.png
  echo -e "----

![QC Metrics for Sample $sample ([download high-res image](graphs/${sample}_QC_Metrics.ps))](graphs/${sample}_QC_Metrics.png)
" \
  >> report/ChipSeq.qc_metrics.md
done
qc_plots_R.947ca4b9fd1de01b4fbbe04f27d2d113.mugqic.done
)
qc_metrics_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$qc_metrics_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: homer_make_ucsc_file
#-------------------------------------------------------------------------------
STEP=homer_make_ucsc_file
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_1_JOB_ID: homer_make_ucsc_file.A549_DEX_CDK9_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.A549_DEX_CDK9_rep1
JOB_DEPENDENCIES=$homer_make_tag_directory_1_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.A549_DEX_CDK9_rep1.9e28a01d2eb36b3f6fea71fd6c950401.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.A549_DEX_CDK9_rep1.9e28a01d2eb36b3f6fea71fd6c950401.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/A549_DEX_CDK9_rep1 && \
makeUCSCfile \
  tags/A549_DEX_CDK9_rep1 | \
gzip -1 -c > tracks/A549_DEX_CDK9_rep1/A549_DEX_CDK9_rep1.ucsc.bedGraph.gz
homer_make_ucsc_file.A549_DEX_CDK9_rep1.9e28a01d2eb36b3f6fea71fd6c950401.mugqic.done
)
homer_make_ucsc_file_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_2_JOB_ID: homer_make_ucsc_file.A549_DEX_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.A549_DEX_MED1_rep1
JOB_DEPENDENCIES=$homer_make_tag_directory_2_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.A549_DEX_MED1_rep1.be3464040605eddc317f4c2131f4a9b7.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.A549_DEX_MED1_rep1.be3464040605eddc317f4c2131f4a9b7.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/A549_DEX_MED1_rep1 && \
makeUCSCfile \
  tags/A549_DEX_MED1_rep1 | \
gzip -1 -c > tracks/A549_DEX_MED1_rep1/A549_DEX_MED1_rep1.ucsc.bedGraph.gz
homer_make_ucsc_file.A549_DEX_MED1_rep1.be3464040605eddc317f4c2131f4a9b7.mugqic.done
)
homer_make_ucsc_file_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_3_JOB_ID: homer_make_ucsc_file.A549_DEX_BRD4_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.A549_DEX_BRD4_rep1
JOB_DEPENDENCIES=$homer_make_tag_directory_3_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.A549_DEX_BRD4_rep1.ff579ffbbd8962d7ef6841cf6b7da54d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.A549_DEX_BRD4_rep1.ff579ffbbd8962d7ef6841cf6b7da54d.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/A549_DEX_BRD4_rep1 && \
makeUCSCfile \
  tags/A549_DEX_BRD4_rep1 | \
gzip -1 -c > tracks/A549_DEX_BRD4_rep1/A549_DEX_BRD4_rep1.ucsc.bedGraph.gz
homer_make_ucsc_file.A549_DEX_BRD4_rep1.ff579ffbbd8962d7ef6841cf6b7da54d.mugqic.done
)
homer_make_ucsc_file_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_4_JOB_ID: homer_make_ucsc_file.A549_DEX_WCE_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.A549_DEX_WCE_rep1
JOB_DEPENDENCIES=$homer_make_tag_directory_4_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.A549_DEX_WCE_rep1.db19b39709de948f522547e0899a0cf9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.A549_DEX_WCE_rep1.db19b39709de948f522547e0899a0cf9.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/A549_DEX_WCE_rep1 && \
makeUCSCfile \
  tags/A549_DEX_WCE_rep1 | \
gzip -1 -c > tracks/A549_DEX_WCE_rep1/A549_DEX_WCE_rep1.ucsc.bedGraph.gz
homer_make_ucsc_file.A549_DEX_WCE_rep1.db19b39709de948f522547e0899a0cf9.mugqic.done
)
homer_make_ucsc_file_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_5_JOB_ID: homer_make_ucsc_file.A549_DEX_NIPBL_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.A549_DEX_NIPBL_rep1
JOB_DEPENDENCIES=$homer_make_tag_directory_5_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.A549_DEX_NIPBL_rep1.c53fef7fd643e7035253540005cf417c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.A549_DEX_NIPBL_rep1.c53fef7fd643e7035253540005cf417c.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/A549_DEX_NIPBL_rep1 && \
makeUCSCfile \
  tags/A549_DEX_NIPBL_rep1 | \
gzip -1 -c > tracks/A549_DEX_NIPBL_rep1/A549_DEX_NIPBL_rep1.ucsc.bedGraph.gz
homer_make_ucsc_file.A549_DEX_NIPBL_rep1.c53fef7fd643e7035253540005cf417c.mugqic.done
)
homer_make_ucsc_file_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_6_JOB_ID: homer_make_ucsc_file.A549_CTRL_SMC1A_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.A549_CTRL_SMC1A_rep1
JOB_DEPENDENCIES=$homer_make_tag_directory_6_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.A549_CTRL_SMC1A_rep1.8151a087d4dbd42f0cb8b2641709cbf1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.A549_CTRL_SMC1A_rep1.8151a087d4dbd42f0cb8b2641709cbf1.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/A549_CTRL_SMC1A_rep1 && \
makeUCSCfile \
  tags/A549_CTRL_SMC1A_rep1 | \
gzip -1 -c > tracks/A549_CTRL_SMC1A_rep1/A549_CTRL_SMC1A_rep1.ucsc.bedGraph.gz
homer_make_ucsc_file.A549_CTRL_SMC1A_rep1.8151a087d4dbd42f0cb8b2641709cbf1.mugqic.done
)
homer_make_ucsc_file_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_7_JOB_ID: homer_make_ucsc_file.A549_CTRL_WCE_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.A549_CTRL_WCE_rep1
JOB_DEPENDENCIES=$homer_make_tag_directory_7_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.A549_CTRL_WCE_rep1.892ada9cd3c63d8f21f1d7de934cc2c9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.A549_CTRL_WCE_rep1.892ada9cd3c63d8f21f1d7de934cc2c9.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/A549_CTRL_WCE_rep1 && \
makeUCSCfile \
  tags/A549_CTRL_WCE_rep1 | \
gzip -1 -c > tracks/A549_CTRL_WCE_rep1/A549_CTRL_WCE_rep1.ucsc.bedGraph.gz
homer_make_ucsc_file.A549_CTRL_WCE_rep1.892ada9cd3c63d8f21f1d7de934cc2c9.mugqic.done
)
homer_make_ucsc_file_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_8_JOB_ID: homer_make_ucsc_file.A549_CTRL_CDK9_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.A549_CTRL_CDK9_rep1
JOB_DEPENDENCIES=$homer_make_tag_directory_8_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.A549_CTRL_CDK9_rep1.8b92e915c752dd8076578cc0f7a009ac.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.A549_CTRL_CDK9_rep1.8b92e915c752dd8076578cc0f7a009ac.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/A549_CTRL_CDK9_rep1 && \
makeUCSCfile \
  tags/A549_CTRL_CDK9_rep1 | \
gzip -1 -c > tracks/A549_CTRL_CDK9_rep1/A549_CTRL_CDK9_rep1.ucsc.bedGraph.gz
homer_make_ucsc_file.A549_CTRL_CDK9_rep1.8b92e915c752dd8076578cc0f7a009ac.mugqic.done
)
homer_make_ucsc_file_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_9_JOB_ID: homer_make_ucsc_file.A549_CTRL_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.A549_CTRL_MED1_rep1
JOB_DEPENDENCIES=$homer_make_tag_directory_9_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.A549_CTRL_MED1_rep1.e26c72e9a6cd08c337dc5efb8c8f0ba1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.A549_CTRL_MED1_rep1.e26c72e9a6cd08c337dc5efb8c8f0ba1.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/A549_CTRL_MED1_rep1 && \
makeUCSCfile \
  tags/A549_CTRL_MED1_rep1 | \
gzip -1 -c > tracks/A549_CTRL_MED1_rep1/A549_CTRL_MED1_rep1.ucsc.bedGraph.gz
homer_make_ucsc_file.A549_CTRL_MED1_rep1.e26c72e9a6cd08c337dc5efb8c8f0ba1.mugqic.done
)
homer_make_ucsc_file_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_10_JOB_ID: homer_make_ucsc_file.A549_CTRL_NIPBL_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.A549_CTRL_NIPBL_rep1
JOB_DEPENDENCIES=$homer_make_tag_directory_10_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.A549_CTRL_NIPBL_rep1.3afbc9614dba78040cbd37864e45c63b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.A549_CTRL_NIPBL_rep1.3afbc9614dba78040cbd37864e45c63b.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/A549_CTRL_NIPBL_rep1 && \
makeUCSCfile \
  tags/A549_CTRL_NIPBL_rep1 | \
gzip -1 -c > tracks/A549_CTRL_NIPBL_rep1/A549_CTRL_NIPBL_rep1.ucsc.bedGraph.gz
homer_make_ucsc_file.A549_CTRL_NIPBL_rep1.3afbc9614dba78040cbd37864e45c63b.mugqic.done
)
homer_make_ucsc_file_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_11_JOB_ID: homer_make_ucsc_file.A549_CTRL_BRD4_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.A549_CTRL_BRD4_rep1
JOB_DEPENDENCIES=$homer_make_tag_directory_11_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.A549_CTRL_BRD4_rep1.907640651b97c9e62556ef346b21599f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.A549_CTRL_BRD4_rep1.907640651b97c9e62556ef346b21599f.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/A549_CTRL_BRD4_rep1 && \
makeUCSCfile \
  tags/A549_CTRL_BRD4_rep1 | \
gzip -1 -c > tracks/A549_CTRL_BRD4_rep1/A549_CTRL_BRD4_rep1.ucsc.bedGraph.gz
homer_make_ucsc_file.A549_CTRL_BRD4_rep1.907640651b97c9e62556ef346b21599f.mugqic.done
)
homer_make_ucsc_file_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_12_JOB_ID: homer_make_ucsc_file_report
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file_report
JOB_DEPENDENCIES=$homer_make_ucsc_file_11_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file_report.c9cac04b3bc71f34cbcde042cf5b419d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file_report.c9cac04b3bc71f34cbcde042cf5b419d.mugqic.done'
mkdir -p report && \
zip -r report/tracks.zip tracks/*/*.ucsc.bedGraph.gz && \
cp /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/ChipSeq.homer_make_ucsc_file.md report/
homer_make_ucsc_file_report.c9cac04b3bc71f34cbcde042cf5b419d.mugqic.done
)
homer_make_ucsc_file_12_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: macs2_callpeak
#-------------------------------------------------------------------------------
STEP=macs2_callpeak
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_1_JOB_ID: macs2_callpeak.A549_CTRL_BRD4_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_CTRL_BRD4_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_7_JOB_ID:$picard_mark_duplicates_11_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_CTRL_BRD4_rep1.a126c32befb8514bcb999136f0755168.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_CTRL_BRD4_rep1.a126c32befb8514bcb999136f0755168.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/A549_CTRL_BRD4_rep1 && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2509729011.2 \
  --treatment \
  alignment/A549_CTRL_BRD4_rep1/A549_CTRL_BRD4_rep1.sorted.dup.bam \
  --control \
  alignment/A549_CTRL_WCE_rep1/A549_CTRL_WCE_rep1.sorted.dup.bam \
  --name peak_call/A549_CTRL_BRD4_rep1/A549_CTRL_BRD4_rep1 \
  >& peak_call/A549_CTRL_BRD4_rep1/A549_CTRL_BRD4_rep1.diag.macs.out
macs2_callpeak.A549_CTRL_BRD4_rep1.a126c32befb8514bcb999136f0755168.mugqic.done
)
macs2_callpeak_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_2_JOB_ID: macs2_callpeak.A549_CTRL_CDK9_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_CTRL_CDK9_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_7_JOB_ID:$picard_mark_duplicates_8_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_CTRL_CDK9_rep1.b64fb6bfc8d4b25da549fc370a44bb41.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_CTRL_CDK9_rep1.b64fb6bfc8d4b25da549fc370a44bb41.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/A549_CTRL_CDK9_rep1 && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2509729011.2 \
  --treatment \
  alignment/A549_CTRL_CDK9_rep1/A549_CTRL_CDK9_rep1.sorted.dup.bam \
  --control \
  alignment/A549_CTRL_WCE_rep1/A549_CTRL_WCE_rep1.sorted.dup.bam \
  --name peak_call/A549_CTRL_CDK9_rep1/A549_CTRL_CDK9_rep1 \
  >& peak_call/A549_CTRL_CDK9_rep1/A549_CTRL_CDK9_rep1.diag.macs.out
macs2_callpeak.A549_CTRL_CDK9_rep1.b64fb6bfc8d4b25da549fc370a44bb41.mugqic.done
)
macs2_callpeak_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_3_JOB_ID: macs2_callpeak.A549_CTRL_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_CTRL_MED1_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_7_JOB_ID:$picard_mark_duplicates_9_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_CTRL_MED1_rep1.cb1adc1433375c57a1613e1b783534ee.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_CTRL_MED1_rep1.cb1adc1433375c57a1613e1b783534ee.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/A549_CTRL_MED1_rep1 && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2509729011.2 \
  --treatment \
  alignment/A549_CTRL_MED1_rep1/A549_CTRL_MED1_rep1.sorted.dup.bam \
  --control \
  alignment/A549_CTRL_WCE_rep1/A549_CTRL_WCE_rep1.sorted.dup.bam \
  --name peak_call/A549_CTRL_MED1_rep1/A549_CTRL_MED1_rep1 \
  >& peak_call/A549_CTRL_MED1_rep1/A549_CTRL_MED1_rep1.diag.macs.out
macs2_callpeak.A549_CTRL_MED1_rep1.cb1adc1433375c57a1613e1b783534ee.mugqic.done
)
macs2_callpeak_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_4_JOB_ID: macs2_callpeak.A549_CTRL_NIPBL_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_CTRL_NIPBL_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_7_JOB_ID:$picard_mark_duplicates_10_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_CTRL_NIPBL_rep1.de69b783536660cbcfcbf21e23255f3f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_CTRL_NIPBL_rep1.de69b783536660cbcfcbf21e23255f3f.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/A549_CTRL_NIPBL_rep1 && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2509729011.2 \
  --treatment \
  alignment/A549_CTRL_NIPBL_rep1/A549_CTRL_NIPBL_rep1.sorted.dup.bam \
  --control \
  alignment/A549_CTRL_WCE_rep1/A549_CTRL_WCE_rep1.sorted.dup.bam \
  --name peak_call/A549_CTRL_NIPBL_rep1/A549_CTRL_NIPBL_rep1 \
  >& peak_call/A549_CTRL_NIPBL_rep1/A549_CTRL_NIPBL_rep1.diag.macs.out
macs2_callpeak.A549_CTRL_NIPBL_rep1.de69b783536660cbcfcbf21e23255f3f.mugqic.done
)
macs2_callpeak_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_5_JOB_ID: macs2_callpeak.A549_CTRL_SMC1A_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_CTRL_SMC1A_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_6_JOB_ID:$picard_mark_duplicates_7_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_CTRL_SMC1A_rep1.8cb5b263fc86cc07927aa047913df479.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_CTRL_SMC1A_rep1.8cb5b263fc86cc07927aa047913df479.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/A549_CTRL_SMC1A_rep1 && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2509729011.2 \
  --treatment \
  alignment/A549_CTRL_SMC1A_rep1/A549_CTRL_SMC1A_rep1.sorted.dup.bam \
  --control \
  alignment/A549_CTRL_WCE_rep1/A549_CTRL_WCE_rep1.sorted.dup.bam \
  --name peak_call/A549_CTRL_SMC1A_rep1/A549_CTRL_SMC1A_rep1 \
  >& peak_call/A549_CTRL_SMC1A_rep1/A549_CTRL_SMC1A_rep1.diag.macs.out
macs2_callpeak.A549_CTRL_SMC1A_rep1.8cb5b263fc86cc07927aa047913df479.mugqic.done
)
macs2_callpeak_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_6_JOB_ID: macs2_callpeak.A549_DEX_BRD4_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_DEX_BRD4_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_3_JOB_ID:$picard_mark_duplicates_4_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_DEX_BRD4_rep1.e23eff6c938f8fe7630bbf112ce6777a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_DEX_BRD4_rep1.e23eff6c938f8fe7630bbf112ce6777a.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/A549_DEX_BRD4_rep1 && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2509729011.2 \
  --treatment \
  alignment/A549_DEX_BRD4_rep1/A549_DEX_BRD4_rep1.sorted.dup.bam \
  --control \
  alignment/A549_DEX_WCE_rep1/A549_DEX_WCE_rep1.sorted.dup.bam \
  --name peak_call/A549_DEX_BRD4_rep1/A549_DEX_BRD4_rep1 \
  >& peak_call/A549_DEX_BRD4_rep1/A549_DEX_BRD4_rep1.diag.macs.out
macs2_callpeak.A549_DEX_BRD4_rep1.e23eff6c938f8fe7630bbf112ce6777a.mugqic.done
)
macs2_callpeak_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_7_JOB_ID: macs2_callpeak.A549_DEX_CDK9_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_DEX_CDK9_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_4_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_DEX_CDK9_rep1.e4d827c05e1de219aa0c6c0cac2b7673.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_DEX_CDK9_rep1.e4d827c05e1de219aa0c6c0cac2b7673.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/A549_DEX_CDK9_rep1 && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2509729011.2 \
  --treatment \
  alignment/A549_DEX_CDK9_rep1/A549_DEX_CDK9_rep1.sorted.dup.bam \
  --control \
  alignment/A549_DEX_WCE_rep1/A549_DEX_WCE_rep1.sorted.dup.bam \
  --name peak_call/A549_DEX_CDK9_rep1/A549_DEX_CDK9_rep1 \
  >& peak_call/A549_DEX_CDK9_rep1/A549_DEX_CDK9_rep1.diag.macs.out
macs2_callpeak.A549_DEX_CDK9_rep1.e4d827c05e1de219aa0c6c0cac2b7673.mugqic.done
)
macs2_callpeak_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_8_JOB_ID: macs2_callpeak.A549_DEX_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_DEX_MED1_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_2_JOB_ID:$picard_mark_duplicates_4_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_DEX_MED1_rep1.c6eab73ccdaf34bf6b27aa3578afbd3e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_DEX_MED1_rep1.c6eab73ccdaf34bf6b27aa3578afbd3e.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/A549_DEX_MED1_rep1 && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2509729011.2 \
  --treatment \
  alignment/A549_DEX_MED1_rep1/A549_DEX_MED1_rep1.sorted.dup.bam \
  --control \
  alignment/A549_DEX_WCE_rep1/A549_DEX_WCE_rep1.sorted.dup.bam \
  --name peak_call/A549_DEX_MED1_rep1/A549_DEX_MED1_rep1 \
  >& peak_call/A549_DEX_MED1_rep1/A549_DEX_MED1_rep1.diag.macs.out
macs2_callpeak.A549_DEX_MED1_rep1.c6eab73ccdaf34bf6b27aa3578afbd3e.mugqic.done
)
macs2_callpeak_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_9_JOB_ID: macs2_callpeak.A549_DEX_NIPBL_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_DEX_NIPBL_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_4_JOB_ID:$picard_mark_duplicates_5_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_DEX_NIPBL_rep1.c04e3b3507f31fe9aca332d743aee4b1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_DEX_NIPBL_rep1.c04e3b3507f31fe9aca332d743aee4b1.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/A549_DEX_NIPBL_rep1 && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2509729011.2 \
  --treatment \
  alignment/A549_DEX_NIPBL_rep1/A549_DEX_NIPBL_rep1.sorted.dup.bam \
  --control \
  alignment/A549_DEX_WCE_rep1/A549_DEX_WCE_rep1.sorted.dup.bam \
  --name peak_call/A549_DEX_NIPBL_rep1/A549_DEX_NIPBL_rep1 \
  >& peak_call/A549_DEX_NIPBL_rep1/A549_DEX_NIPBL_rep1.diag.macs.out
macs2_callpeak.A549_DEX_NIPBL_rep1.c04e3b3507f31fe9aca332d743aee4b1.mugqic.done
)
macs2_callpeak_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_10_JOB_ID: macs2_callpeak_report
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_report
JOB_DEPENDENCIES=$macs2_callpeak_1_JOB_ID:$macs2_callpeak_2_JOB_ID:$macs2_callpeak_3_JOB_ID:$macs2_callpeak_4_JOB_ID:$macs2_callpeak_5_JOB_ID:$macs2_callpeak_6_JOB_ID:$macs2_callpeak_7_JOB_ID:$macs2_callpeak_8_JOB_ID:$macs2_callpeak_9_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_report.18d6627390a9896da322d4cc4fbfb157.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_report.18d6627390a9896da322d4cc4fbfb157.mugqic.done'
mkdir -p report && \
cp /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/ChipSeq.macs2_callpeak.md report/ && \
for contrast in A549_CTRL_BRD4_rep1 A549_CTRL_CDK9_rep1 A549_CTRL_MED1_rep1 A549_CTRL_NIPBL_rep1 A549_CTRL_SMC1A_rep1 A549_DEX_BRD4_rep1 A549_DEX_CDK9_rep1 A549_DEX_MED1_rep1 A549_DEX_NIPBL_rep1
do
  cp -a --parents peak_call/$contrast/ report/ && \
  echo -e "* [Peak Calls File for Design $contrast](peak_call/$contrast/${contrast}_peaks.xls)" \
  >> report/ChipSeq.macs2_callpeak.md
done
macs2_callpeak_report.18d6627390a9896da322d4cc4fbfb157.mugqic.done
)
macs2_callpeak_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: homer_annotate_peaks
#-------------------------------------------------------------------------------
STEP=homer_annotate_peaks
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_1_JOB_ID: homer_annotate_peaks.A549_CTRL_BRD4_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.A549_CTRL_BRD4_rep1
JOB_DEPENDENCIES=$macs2_callpeak_1_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.A549_CTRL_BRD4_rep1.d53dfdcff6cbe442957e357fe3595760.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.A549_CTRL_BRD4_rep1.d53dfdcff6cbe442957e357fe3595760.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/A549_CTRL_BRD4_rep1/A549_CTRL_BRD4_rep1 && \
annotatePeaks.pl \
  peak_call/A549_CTRL_BRD4_rep1/A549_CTRL_BRD4_rep1_peaks.narrowPeak \
  hg19 \
  -gsize hg19 \
  -cons -CpG \
  -go annotation/A549_CTRL_BRD4_rep1/A549_CTRL_BRD4_rep1 \
  -genomeOntology annotation/A549_CTRL_BRD4_rep1/A549_CTRL_BRD4_rep1 \
  > annotation/A549_CTRL_BRD4_rep1/A549_CTRL_BRD4_rep1.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/A549_CTRL_BRD4_rep1/A549_CTRL_BRD4_rep1.annotated.csv",
  "annotation/A549_CTRL_BRD4_rep1/A549_CTRL_BRD4_rep1",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.A549_CTRL_BRD4_rep1.d53dfdcff6cbe442957e357fe3595760.mugqic.done
)
homer_annotate_peaks_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_2_JOB_ID: homer_annotate_peaks.A549_CTRL_CDK9_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.A549_CTRL_CDK9_rep1
JOB_DEPENDENCIES=$macs2_callpeak_2_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.A549_CTRL_CDK9_rep1.66b1daa6efe1174f841b2a8cb3315943.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.A549_CTRL_CDK9_rep1.66b1daa6efe1174f841b2a8cb3315943.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/A549_CTRL_CDK9_rep1/A549_CTRL_CDK9_rep1 && \
annotatePeaks.pl \
  peak_call/A549_CTRL_CDK9_rep1/A549_CTRL_CDK9_rep1_peaks.narrowPeak \
  hg19 \
  -gsize hg19 \
  -cons -CpG \
  -go annotation/A549_CTRL_CDK9_rep1/A549_CTRL_CDK9_rep1 \
  -genomeOntology annotation/A549_CTRL_CDK9_rep1/A549_CTRL_CDK9_rep1 \
  > annotation/A549_CTRL_CDK9_rep1/A549_CTRL_CDK9_rep1.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/A549_CTRL_CDK9_rep1/A549_CTRL_CDK9_rep1.annotated.csv",
  "annotation/A549_CTRL_CDK9_rep1/A549_CTRL_CDK9_rep1",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.A549_CTRL_CDK9_rep1.66b1daa6efe1174f841b2a8cb3315943.mugqic.done
)
homer_annotate_peaks_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_3_JOB_ID: homer_annotate_peaks.A549_CTRL_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.A549_CTRL_MED1_rep1
JOB_DEPENDENCIES=$macs2_callpeak_3_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.A549_CTRL_MED1_rep1.8d835b28b8f0e89a09e3be0e6bc7c913.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.A549_CTRL_MED1_rep1.8d835b28b8f0e89a09e3be0e6bc7c913.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/A549_CTRL_MED1_rep1/A549_CTRL_MED1_rep1 && \
annotatePeaks.pl \
  peak_call/A549_CTRL_MED1_rep1/A549_CTRL_MED1_rep1_peaks.narrowPeak \
  hg19 \
  -gsize hg19 \
  -cons -CpG \
  -go annotation/A549_CTRL_MED1_rep1/A549_CTRL_MED1_rep1 \
  -genomeOntology annotation/A549_CTRL_MED1_rep1/A549_CTRL_MED1_rep1 \
  > annotation/A549_CTRL_MED1_rep1/A549_CTRL_MED1_rep1.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/A549_CTRL_MED1_rep1/A549_CTRL_MED1_rep1.annotated.csv",
  "annotation/A549_CTRL_MED1_rep1/A549_CTRL_MED1_rep1",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.A549_CTRL_MED1_rep1.8d835b28b8f0e89a09e3be0e6bc7c913.mugqic.done
)
homer_annotate_peaks_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_4_JOB_ID: homer_annotate_peaks.A549_CTRL_NIPBL_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.A549_CTRL_NIPBL_rep1
JOB_DEPENDENCIES=$macs2_callpeak_4_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.A549_CTRL_NIPBL_rep1.90e50ce5172f2df2609fd94b4feee43a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.A549_CTRL_NIPBL_rep1.90e50ce5172f2df2609fd94b4feee43a.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/A549_CTRL_NIPBL_rep1/A549_CTRL_NIPBL_rep1 && \
annotatePeaks.pl \
  peak_call/A549_CTRL_NIPBL_rep1/A549_CTRL_NIPBL_rep1_peaks.narrowPeak \
  hg19 \
  -gsize hg19 \
  -cons -CpG \
  -go annotation/A549_CTRL_NIPBL_rep1/A549_CTRL_NIPBL_rep1 \
  -genomeOntology annotation/A549_CTRL_NIPBL_rep1/A549_CTRL_NIPBL_rep1 \
  > annotation/A549_CTRL_NIPBL_rep1/A549_CTRL_NIPBL_rep1.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/A549_CTRL_NIPBL_rep1/A549_CTRL_NIPBL_rep1.annotated.csv",
  "annotation/A549_CTRL_NIPBL_rep1/A549_CTRL_NIPBL_rep1",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.A549_CTRL_NIPBL_rep1.90e50ce5172f2df2609fd94b4feee43a.mugqic.done
)
homer_annotate_peaks_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_5_JOB_ID: homer_annotate_peaks.A549_CTRL_SMC1A_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.A549_CTRL_SMC1A_rep1
JOB_DEPENDENCIES=$macs2_callpeak_5_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.A549_CTRL_SMC1A_rep1.feb66823aeb67b1b52cfd0714924de26.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.A549_CTRL_SMC1A_rep1.feb66823aeb67b1b52cfd0714924de26.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/A549_CTRL_SMC1A_rep1/A549_CTRL_SMC1A_rep1 && \
annotatePeaks.pl \
  peak_call/A549_CTRL_SMC1A_rep1/A549_CTRL_SMC1A_rep1_peaks.narrowPeak \
  hg19 \
  -gsize hg19 \
  -cons -CpG \
  -go annotation/A549_CTRL_SMC1A_rep1/A549_CTRL_SMC1A_rep1 \
  -genomeOntology annotation/A549_CTRL_SMC1A_rep1/A549_CTRL_SMC1A_rep1 \
  > annotation/A549_CTRL_SMC1A_rep1/A549_CTRL_SMC1A_rep1.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/A549_CTRL_SMC1A_rep1/A549_CTRL_SMC1A_rep1.annotated.csv",
  "annotation/A549_CTRL_SMC1A_rep1/A549_CTRL_SMC1A_rep1",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.A549_CTRL_SMC1A_rep1.feb66823aeb67b1b52cfd0714924de26.mugqic.done
)
homer_annotate_peaks_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_6_JOB_ID: homer_annotate_peaks.A549_DEX_BRD4_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.A549_DEX_BRD4_rep1
JOB_DEPENDENCIES=$macs2_callpeak_6_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.A549_DEX_BRD4_rep1.fd58af76f0c10f0abd4de2caa85d0872.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.A549_DEX_BRD4_rep1.fd58af76f0c10f0abd4de2caa85d0872.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/A549_DEX_BRD4_rep1/A549_DEX_BRD4_rep1 && \
annotatePeaks.pl \
  peak_call/A549_DEX_BRD4_rep1/A549_DEX_BRD4_rep1_peaks.narrowPeak \
  hg19 \
  -gsize hg19 \
  -cons -CpG \
  -go annotation/A549_DEX_BRD4_rep1/A549_DEX_BRD4_rep1 \
  -genomeOntology annotation/A549_DEX_BRD4_rep1/A549_DEX_BRD4_rep1 \
  > annotation/A549_DEX_BRD4_rep1/A549_DEX_BRD4_rep1.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/A549_DEX_BRD4_rep1/A549_DEX_BRD4_rep1.annotated.csv",
  "annotation/A549_DEX_BRD4_rep1/A549_DEX_BRD4_rep1",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.A549_DEX_BRD4_rep1.fd58af76f0c10f0abd4de2caa85d0872.mugqic.done
)
homer_annotate_peaks_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_7_JOB_ID: homer_annotate_peaks.A549_DEX_CDK9_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.A549_DEX_CDK9_rep1
JOB_DEPENDENCIES=$macs2_callpeak_7_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.A549_DEX_CDK9_rep1.bfe0edf9a0e08696f870c6cd73e24540.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.A549_DEX_CDK9_rep1.bfe0edf9a0e08696f870c6cd73e24540.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/A549_DEX_CDK9_rep1/A549_DEX_CDK9_rep1 && \
annotatePeaks.pl \
  peak_call/A549_DEX_CDK9_rep1/A549_DEX_CDK9_rep1_peaks.narrowPeak \
  hg19 \
  -gsize hg19 \
  -cons -CpG \
  -go annotation/A549_DEX_CDK9_rep1/A549_DEX_CDK9_rep1 \
  -genomeOntology annotation/A549_DEX_CDK9_rep1/A549_DEX_CDK9_rep1 \
  > annotation/A549_DEX_CDK9_rep1/A549_DEX_CDK9_rep1.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/A549_DEX_CDK9_rep1/A549_DEX_CDK9_rep1.annotated.csv",
  "annotation/A549_DEX_CDK9_rep1/A549_DEX_CDK9_rep1",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.A549_DEX_CDK9_rep1.bfe0edf9a0e08696f870c6cd73e24540.mugqic.done
)
homer_annotate_peaks_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_8_JOB_ID: homer_annotate_peaks.A549_DEX_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.A549_DEX_MED1_rep1
JOB_DEPENDENCIES=$macs2_callpeak_8_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.A549_DEX_MED1_rep1.c667a43c5b8d87dcdf0f86a539d9037b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.A549_DEX_MED1_rep1.c667a43c5b8d87dcdf0f86a539d9037b.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/A549_DEX_MED1_rep1/A549_DEX_MED1_rep1 && \
annotatePeaks.pl \
  peak_call/A549_DEX_MED1_rep1/A549_DEX_MED1_rep1_peaks.narrowPeak \
  hg19 \
  -gsize hg19 \
  -cons -CpG \
  -go annotation/A549_DEX_MED1_rep1/A549_DEX_MED1_rep1 \
  -genomeOntology annotation/A549_DEX_MED1_rep1/A549_DEX_MED1_rep1 \
  > annotation/A549_DEX_MED1_rep1/A549_DEX_MED1_rep1.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/A549_DEX_MED1_rep1/A549_DEX_MED1_rep1.annotated.csv",
  "annotation/A549_DEX_MED1_rep1/A549_DEX_MED1_rep1",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.A549_DEX_MED1_rep1.c667a43c5b8d87dcdf0f86a539d9037b.mugqic.done
)
homer_annotate_peaks_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_9_JOB_ID: homer_annotate_peaks.A549_DEX_NIPBL_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.A549_DEX_NIPBL_rep1
JOB_DEPENDENCIES=$macs2_callpeak_9_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.A549_DEX_NIPBL_rep1.cd89116013d7ee5b098ba8cc478fcf92.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.A549_DEX_NIPBL_rep1.cd89116013d7ee5b098ba8cc478fcf92.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/A549_DEX_NIPBL_rep1/A549_DEX_NIPBL_rep1 && \
annotatePeaks.pl \
  peak_call/A549_DEX_NIPBL_rep1/A549_DEX_NIPBL_rep1_peaks.narrowPeak \
  hg19 \
  -gsize hg19 \
  -cons -CpG \
  -go annotation/A549_DEX_NIPBL_rep1/A549_DEX_NIPBL_rep1 \
  -genomeOntology annotation/A549_DEX_NIPBL_rep1/A549_DEX_NIPBL_rep1 \
  > annotation/A549_DEX_NIPBL_rep1/A549_DEX_NIPBL_rep1.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/A549_DEX_NIPBL_rep1/A549_DEX_NIPBL_rep1.annotated.csv",
  "annotation/A549_DEX_NIPBL_rep1/A549_DEX_NIPBL_rep1",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.A549_DEX_NIPBL_rep1.cd89116013d7ee5b098ba8cc478fcf92.mugqic.done
)
homer_annotate_peaks_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_10_JOB_ID: homer_annotate_peaks_report
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks_report
JOB_DEPENDENCIES=$homer_annotate_peaks_1_JOB_ID:$homer_annotate_peaks_2_JOB_ID:$homer_annotate_peaks_3_JOB_ID:$homer_annotate_peaks_4_JOB_ID:$homer_annotate_peaks_5_JOB_ID:$homer_annotate_peaks_6_JOB_ID:$homer_annotate_peaks_7_JOB_ID:$homer_annotate_peaks_8_JOB_ID:$homer_annotate_peaks_9_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks_report.94b242b078b0b2c79758ad717be2ac32.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks_report.94b242b078b0b2c79758ad717be2ac32.mugqic.done'
mkdir -p report/annotation/ && \
cp /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/ChipSeq.homer_annotate_peaks.md report/ && \
for contrast in A549_CTRL_BRD4_rep1 A549_CTRL_CDK9_rep1 A549_CTRL_MED1_rep1 A549_CTRL_NIPBL_rep1 A549_CTRL_SMC1A_rep1 A549_DEX_BRD4_rep1 A549_DEX_CDK9_rep1 A549_DEX_MED1_rep1 A549_DEX_NIPBL_rep1
do
  rsync -avP annotation/$contrast report/annotation/ && \
  echo -e "* [Gene Annotations for Design $contrast](annotation/$contrast/${contrast}.annotated.csv)
* [HOMER Gene Ontology Annotations for Design $contrast](annotation/$contrast/$contrast/geneOntology.html)
* [HOMER Genome Ontology Annotations for Design $contrast](annotation/$contrast/$contrast/GenomeOntology.html)" \
  >> report/ChipSeq.homer_annotate_peaks.md
done
homer_annotate_peaks_report.94b242b078b0b2c79758ad717be2ac32.mugqic.done
)
homer_annotate_peaks_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: homer_find_motifs_genome
#-------------------------------------------------------------------------------
STEP=homer_find_motifs_genome
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_1_JOB_ID: homer_find_motifs_genome.A549_CTRL_BRD4_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.A549_CTRL_BRD4_rep1
JOB_DEPENDENCIES=$macs2_callpeak_1_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.A549_CTRL_BRD4_rep1.e4c1c91e8faa4484caf40387a4bb2371.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.A549_CTRL_BRD4_rep1.e4c1c91e8faa4484caf40387a4bb2371.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/A549_CTRL_BRD4_rep1/A549_CTRL_BRD4_rep1 && \
findMotifsGenome.pl \
  peak_call/A549_CTRL_BRD4_rep1/A549_CTRL_BRD4_rep1_peaks.narrowPeak \
  hg19 \
  annotation/A549_CTRL_BRD4_rep1/A549_CTRL_BRD4_rep1 \
  -preparsedDir annotation/A549_CTRL_BRD4_rep1/A549_CTRL_BRD4_rep1/preparsed \
  -p 4
homer_find_motifs_genome.A549_CTRL_BRD4_rep1.e4c1c91e8faa4484caf40387a4bb2371.mugqic.done
)
homer_find_motifs_genome_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_2_JOB_ID: homer_find_motifs_genome.A549_CTRL_CDK9_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.A549_CTRL_CDK9_rep1
JOB_DEPENDENCIES=$macs2_callpeak_2_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.A549_CTRL_CDK9_rep1.01146e983b8fe20acd6f02359378c2a6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.A549_CTRL_CDK9_rep1.01146e983b8fe20acd6f02359378c2a6.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/A549_CTRL_CDK9_rep1/A549_CTRL_CDK9_rep1 && \
findMotifsGenome.pl \
  peak_call/A549_CTRL_CDK9_rep1/A549_CTRL_CDK9_rep1_peaks.narrowPeak \
  hg19 \
  annotation/A549_CTRL_CDK9_rep1/A549_CTRL_CDK9_rep1 \
  -preparsedDir annotation/A549_CTRL_CDK9_rep1/A549_CTRL_CDK9_rep1/preparsed \
  -p 4
homer_find_motifs_genome.A549_CTRL_CDK9_rep1.01146e983b8fe20acd6f02359378c2a6.mugqic.done
)
homer_find_motifs_genome_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_3_JOB_ID: homer_find_motifs_genome.A549_CTRL_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.A549_CTRL_MED1_rep1
JOB_DEPENDENCIES=$macs2_callpeak_3_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.A549_CTRL_MED1_rep1.b49327b684f6d0124c54b14580f03df3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.A549_CTRL_MED1_rep1.b49327b684f6d0124c54b14580f03df3.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/A549_CTRL_MED1_rep1/A549_CTRL_MED1_rep1 && \
findMotifsGenome.pl \
  peak_call/A549_CTRL_MED1_rep1/A549_CTRL_MED1_rep1_peaks.narrowPeak \
  hg19 \
  annotation/A549_CTRL_MED1_rep1/A549_CTRL_MED1_rep1 \
  -preparsedDir annotation/A549_CTRL_MED1_rep1/A549_CTRL_MED1_rep1/preparsed \
  -p 4
homer_find_motifs_genome.A549_CTRL_MED1_rep1.b49327b684f6d0124c54b14580f03df3.mugqic.done
)
homer_find_motifs_genome_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_4_JOB_ID: homer_find_motifs_genome.A549_CTRL_NIPBL_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.A549_CTRL_NIPBL_rep1
JOB_DEPENDENCIES=$macs2_callpeak_4_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.A549_CTRL_NIPBL_rep1.5da3e2ba37a4b27accbe1ae10d3dac37.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.A549_CTRL_NIPBL_rep1.5da3e2ba37a4b27accbe1ae10d3dac37.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/A549_CTRL_NIPBL_rep1/A549_CTRL_NIPBL_rep1 && \
findMotifsGenome.pl \
  peak_call/A549_CTRL_NIPBL_rep1/A549_CTRL_NIPBL_rep1_peaks.narrowPeak \
  hg19 \
  annotation/A549_CTRL_NIPBL_rep1/A549_CTRL_NIPBL_rep1 \
  -preparsedDir annotation/A549_CTRL_NIPBL_rep1/A549_CTRL_NIPBL_rep1/preparsed \
  -p 4
homer_find_motifs_genome.A549_CTRL_NIPBL_rep1.5da3e2ba37a4b27accbe1ae10d3dac37.mugqic.done
)
homer_find_motifs_genome_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_5_JOB_ID: homer_find_motifs_genome.A549_CTRL_SMC1A_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.A549_CTRL_SMC1A_rep1
JOB_DEPENDENCIES=$macs2_callpeak_5_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.A549_CTRL_SMC1A_rep1.a86095e677d1cc87fe24e672dfded8d7.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.A549_CTRL_SMC1A_rep1.a86095e677d1cc87fe24e672dfded8d7.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/A549_CTRL_SMC1A_rep1/A549_CTRL_SMC1A_rep1 && \
findMotifsGenome.pl \
  peak_call/A549_CTRL_SMC1A_rep1/A549_CTRL_SMC1A_rep1_peaks.narrowPeak \
  hg19 \
  annotation/A549_CTRL_SMC1A_rep1/A549_CTRL_SMC1A_rep1 \
  -preparsedDir annotation/A549_CTRL_SMC1A_rep1/A549_CTRL_SMC1A_rep1/preparsed \
  -p 4
homer_find_motifs_genome.A549_CTRL_SMC1A_rep1.a86095e677d1cc87fe24e672dfded8d7.mugqic.done
)
homer_find_motifs_genome_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_6_JOB_ID: homer_find_motifs_genome.A549_DEX_BRD4_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.A549_DEX_BRD4_rep1
JOB_DEPENDENCIES=$macs2_callpeak_6_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.A549_DEX_BRD4_rep1.f567b9ea298ba6e7c86b880fd8547aaf.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.A549_DEX_BRD4_rep1.f567b9ea298ba6e7c86b880fd8547aaf.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/A549_DEX_BRD4_rep1/A549_DEX_BRD4_rep1 && \
findMotifsGenome.pl \
  peak_call/A549_DEX_BRD4_rep1/A549_DEX_BRD4_rep1_peaks.narrowPeak \
  hg19 \
  annotation/A549_DEX_BRD4_rep1/A549_DEX_BRD4_rep1 \
  -preparsedDir annotation/A549_DEX_BRD4_rep1/A549_DEX_BRD4_rep1/preparsed \
  -p 4
homer_find_motifs_genome.A549_DEX_BRD4_rep1.f567b9ea298ba6e7c86b880fd8547aaf.mugqic.done
)
homer_find_motifs_genome_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_7_JOB_ID: homer_find_motifs_genome.A549_DEX_CDK9_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.A549_DEX_CDK9_rep1
JOB_DEPENDENCIES=$macs2_callpeak_7_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.A549_DEX_CDK9_rep1.dfa7e7cc72ac65c0bde93d8fe1a4a8d4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.A549_DEX_CDK9_rep1.dfa7e7cc72ac65c0bde93d8fe1a4a8d4.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/A549_DEX_CDK9_rep1/A549_DEX_CDK9_rep1 && \
findMotifsGenome.pl \
  peak_call/A549_DEX_CDK9_rep1/A549_DEX_CDK9_rep1_peaks.narrowPeak \
  hg19 \
  annotation/A549_DEX_CDK9_rep1/A549_DEX_CDK9_rep1 \
  -preparsedDir annotation/A549_DEX_CDK9_rep1/A549_DEX_CDK9_rep1/preparsed \
  -p 4
homer_find_motifs_genome.A549_DEX_CDK9_rep1.dfa7e7cc72ac65c0bde93d8fe1a4a8d4.mugqic.done
)
homer_find_motifs_genome_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_8_JOB_ID: homer_find_motifs_genome.A549_DEX_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.A549_DEX_MED1_rep1
JOB_DEPENDENCIES=$macs2_callpeak_8_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.A549_DEX_MED1_rep1.605ee47205fb9573e717ff8017500fe5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.A549_DEX_MED1_rep1.605ee47205fb9573e717ff8017500fe5.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/A549_DEX_MED1_rep1/A549_DEX_MED1_rep1 && \
findMotifsGenome.pl \
  peak_call/A549_DEX_MED1_rep1/A549_DEX_MED1_rep1_peaks.narrowPeak \
  hg19 \
  annotation/A549_DEX_MED1_rep1/A549_DEX_MED1_rep1 \
  -preparsedDir annotation/A549_DEX_MED1_rep1/A549_DEX_MED1_rep1/preparsed \
  -p 4
homer_find_motifs_genome.A549_DEX_MED1_rep1.605ee47205fb9573e717ff8017500fe5.mugqic.done
)
homer_find_motifs_genome_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_9_JOB_ID: homer_find_motifs_genome.A549_DEX_NIPBL_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.A549_DEX_NIPBL_rep1
JOB_DEPENDENCIES=$macs2_callpeak_9_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.A549_DEX_NIPBL_rep1.4097e2558a66b5a5369f4423e44149b6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.A549_DEX_NIPBL_rep1.4097e2558a66b5a5369f4423e44149b6.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/A549_DEX_NIPBL_rep1/A549_DEX_NIPBL_rep1 && \
findMotifsGenome.pl \
  peak_call/A549_DEX_NIPBL_rep1/A549_DEX_NIPBL_rep1_peaks.narrowPeak \
  hg19 \
  annotation/A549_DEX_NIPBL_rep1/A549_DEX_NIPBL_rep1 \
  -preparsedDir annotation/A549_DEX_NIPBL_rep1/A549_DEX_NIPBL_rep1/preparsed \
  -p 4
homer_find_motifs_genome.A549_DEX_NIPBL_rep1.4097e2558a66b5a5369f4423e44149b6.mugqic.done
)
homer_find_motifs_genome_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_10_JOB_ID: homer_find_motifs_genome_report
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome_report
JOB_DEPENDENCIES=$homer_find_motifs_genome_1_JOB_ID:$homer_find_motifs_genome_2_JOB_ID:$homer_find_motifs_genome_3_JOB_ID:$homer_find_motifs_genome_4_JOB_ID:$homer_find_motifs_genome_5_JOB_ID:$homer_find_motifs_genome_6_JOB_ID:$homer_find_motifs_genome_7_JOB_ID:$homer_find_motifs_genome_8_JOB_ID:$homer_find_motifs_genome_9_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome_report.4e9be54f739d3cffdce4a3946b2baa36.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome_report.4e9be54f739d3cffdce4a3946b2baa36.mugqic.done'
mkdir -p report/annotation/ && \
cp /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/ChipSeq.homer_find_motifs_genome.md report/ && \
for contrast in A549_CTRL_BRD4_rep1 A549_CTRL_CDK9_rep1 A549_CTRL_MED1_rep1 A549_CTRL_NIPBL_rep1 A549_CTRL_SMC1A_rep1 A549_DEX_BRD4_rep1 A549_DEX_CDK9_rep1 A549_DEX_MED1_rep1 A549_DEX_NIPBL_rep1
do
  rsync -avP annotation/$contrast report/annotation/ && \
  echo -e "* [HOMER _De Novo_ Motif Results for Design $contrast](annotation/$contrast/$contrast/homerResults.html)
* [HOMER Known Motif Results for Design $contrast](annotation/$contrast/$contrast/knownResults.html)" \
  >> report/ChipSeq.homer_find_motifs_genome.md
done
homer_find_motifs_genome_report.4e9be54f739d3cffdce4a3946b2baa36.mugqic.done
)
homer_find_motifs_genome_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: annotation_graphs
#-------------------------------------------------------------------------------
STEP=annotation_graphs
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: annotation_graphs_1_JOB_ID: annotation_graphs
#-------------------------------------------------------------------------------
JOB_NAME=annotation_graphs
JOB_DEPENDENCIES=$homer_annotate_peaks_1_JOB_ID:$homer_annotate_peaks_2_JOB_ID:$homer_annotate_peaks_3_JOB_ID:$homer_annotate_peaks_4_JOB_ID:$homer_annotate_peaks_5_JOB_ID:$homer_annotate_peaks_6_JOB_ID:$homer_annotate_peaks_7_JOB_ID:$homer_annotate_peaks_8_JOB_ID:$homer_annotate_peaks_9_JOB_ID
JOB_DONE=job_output/annotation_graphs/annotation_graphs.4cb816b15a519602c91a0915f3c4926b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'annotation_graphs.4cb816b15a519602c91a0915f3c4926b.mugqic.done'
module load mugqic/mugqic_tools/2.1.5 mugqic/R_Bioconductor/3.2.3_3.2 mugqic/pandoc/1.15.2 && \
mkdir -p graphs && \
Rscript $R_TOOLS/chipSeqgenerateAnnotationGraphs.R \
  ../../raw/chip-seq/design.txt \
  /gs/scratch/efournier/CofactorHR/A549/output/chip-pipeline && \
mkdir -p report/annotation/ && \
if [[ -f annotation/peak_stats.csv ]]
then
  cp annotation/peak_stats.csv report/annotation/
peak_stats_table=`LC_NUMERIC=en_CA awk -F "," '{OFS="|"; if (NR == 1) {$1 = $1; print $0; print "-----|-----|-----:|-----:|-----:|-----:|-----:|-----:"} else {print $1, $2,  sprintf("%\47d", $3), $4, sprintf("%\47.1f", $5), sprintf("%\47.1f", $6), sprintf("%\47.1f", $7), sprintf("%\47.1f", $8)}}' annotation/peak_stats.csv`
else
  peak_stats_table=""
fi
pandoc --to=markdown \
  --template /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/ChipSeq.annotation_graphs.md \
  --variable peak_stats_table="$peak_stats_table" \
  --variable proximal_distance="2" \
  --variable distal_distance="10" \
  --variable distance5d_lower="10" \
  --variable distance5d_upper="100" \
  --variable gene_desert_size="100" \
  /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/ChipSeq.annotation_graphs.md \
  > report/ChipSeq.annotation_graphs.md && \
for contrast in A549_CTRL_BRD4_rep1 A549_CTRL_CDK9_rep1 A549_CTRL_MED1_rep1 A549_CTRL_NIPBL_rep1 A549_CTRL_SMC1A_rep1 A549_DEX_BRD4_rep1 A549_DEX_CDK9_rep1 A549_DEX_MED1_rep1 A549_DEX_NIPBL_rep1
do
  cp --parents graphs/${contrast}_Misc_Graphs.ps report/
  convert -rotate 90 graphs/${contrast}_Misc_Graphs.ps report/graphs/${contrast}_Misc_Graphs.png
  echo -e "----

![Annotation Statistics for Design $contrast ([download high-res image](graphs/${contrast}_Misc_Graphs.ps))](graphs/${contrast}_Misc_Graphs.png)
" \
  >> report/ChipSeq.annotation_graphs.md
done
annotation_graphs.4cb816b15a519602c91a0915f3c4926b.mugqic.done
)
annotation_graphs_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$annotation_graphs_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# Call home with pipeline statistics
#-------------------------------------------------------------------------------
wget "http://mugqic.hpc.mcgill.ca/cgi-bin/pipeline.cgi?hostname=lg-1r17-n03&ip=10.241.129.13&pipeline=ChipSeq&steps=picard_sam_to_fastq,trimmomatic,merge_trimmomatic_stats,bwa_mem_picard_sort_sam,samtools_view_filter,picard_merge_sam_files,picard_mark_duplicates,metrics,homer_make_tag_directory,qc_metrics,homer_make_ucsc_file,macs2_callpeak,homer_annotate_peaks,homer_find_motifs_genome,annotation_graphs&samples=11" --quiet --output-document=/dev/null

