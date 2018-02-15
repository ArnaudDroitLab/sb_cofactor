#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

#-------------------------------------------------------------------------------
# ChipSeq PBSScheduler Job Submission Bash script
# Version: 2.2.0
# Created on: 2017-08-28T20:04:25
# Steps:
#   picard_sam_to_fastq: 8 jobs
#   trimmomatic: 8 jobs
#   merge_trimmomatic_stats: 1 job
#   bwa_mem_picard_sort_sam: 9 jobs
#   samtools_view_filter: 9 jobs
#   picard_merge_sam_files: 8 jobs
#   picard_mark_duplicates: 9 jobs
#   metrics: 2 jobs
#   homer_make_tag_directory: 8 jobs
#   qc_metrics: 1 job
#   homer_make_ucsc_file: 9 jobs
#   macs2_callpeak: 9 jobs
#   homer_annotate_peaks: 9 jobs
#   homer_find_motifs_genome: 9 jobs
#   annotation_graphs: 1 job
#   TOTAL: 100 jobs
#-------------------------------------------------------------------------------

OUTPUT_DIR=/gs/scratch/efournier/sb_cofactor_hr/A549/output/chip-pipeline-GRCh38
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
# JOB: picard_sam_to_fastq_1_JOB_ID: picard_sam_to_fastq.HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.f7ede0ee6a5175d27ae574adfe8145e5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.f7ede0ee6a5175d27ae574adfe8145e5.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/SamToFastq.jar \
  VALIDATION_STRINGENCY=LENIENT \
  INPUT=/gs/scratch/efournier/sb_cofactor_hr/A549/raw/chip-seq/HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam \
  FASTQ=/gs/scratch/efournier/sb_cofactor_hr/A549/raw/chip-seq/HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.single.fastq.gz
picard_sam_to_fastq.HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.f7ede0ee6a5175d27ae574adfe8145e5.mugqic.done
)
picard_sam_to_fastq_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=3 | grep "[0-9]")
echo "$picard_sam_to_fastq_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_2_JOB_ID: picard_sam_to_fastq.HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.472c4b9d873ba986450a4af1323fafb3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.472c4b9d873ba986450a4af1323fafb3.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/SamToFastq.jar \
  VALIDATION_STRINGENCY=LENIENT \
  INPUT=/gs/scratch/efournier/sb_cofactor_hr/A549/raw/chip-seq/HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam \
  FASTQ=/gs/scratch/efournier/sb_cofactor_hr/A549/raw/chip-seq/HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.single.fastq.gz
picard_sam_to_fastq.HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.472c4b9d873ba986450a4af1323fafb3.mugqic.done
)
picard_sam_to_fastq_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=3 | grep "[0-9]")
echo "$picard_sam_to_fastq_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_3_JOB_ID: picard_sam_to_fastq.HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.db7b005b9b2cc485e28cba004c9c1b47.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.db7b005b9b2cc485e28cba004c9c1b47.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/SamToFastq.jar \
  VALIDATION_STRINGENCY=LENIENT \
  INPUT=/gs/scratch/efournier/sb_cofactor_hr/A549/raw/chip-seq/HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam \
  FASTQ=/gs/scratch/efournier/sb_cofactor_hr/A549/raw/chip-seq/HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.single.fastq.gz
picard_sam_to_fastq.HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.db7b005b9b2cc485e28cba004c9c1b47.mugqic.done
)
picard_sam_to_fastq_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=3 | grep "[0-9]")
echo "$picard_sam_to_fastq_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_4_JOB_ID: picard_sam_to_fastq.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.7840c4b632a1fe188e786c306d8d65c8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.7840c4b632a1fe188e786c306d8d65c8.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/SamToFastq.jar \
  VALIDATION_STRINGENCY=LENIENT \
  INPUT=/gs/scratch/efournier/sb_cofactor_hr/A549/raw/chip-seq/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam \
  FASTQ=/gs/scratch/efournier/sb_cofactor_hr/A549/raw/chip-seq/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.single.fastq.gz
picard_sam_to_fastq.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.7840c4b632a1fe188e786c306d8d65c8.mugqic.done
)
picard_sam_to_fastq_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=3 | grep "[0-9]")
echo "$picard_sam_to_fastq_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_5_JOB_ID: picard_sam_to_fastq.HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.4874473ba3b139dee7242ae2f3b4af1c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.4874473ba3b139dee7242ae2f3b4af1c.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/SamToFastq.jar \
  VALIDATION_STRINGENCY=LENIENT \
  INPUT=/gs/scratch/efournier/sb_cofactor_hr/A549/raw/chip-seq/HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam \
  FASTQ=/gs/scratch/efournier/sb_cofactor_hr/A549/raw/chip-seq/HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.single.fastq.gz
picard_sam_to_fastq.HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.4874473ba3b139dee7242ae2f3b4af1c.mugqic.done
)
picard_sam_to_fastq_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=3 | grep "[0-9]")
echo "$picard_sam_to_fastq_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_6_JOB_ID: picard_sam_to_fastq.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.7f36c3719b4f2804418ba8ef84bd9c10.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.7f36c3719b4f2804418ba8ef84bd9c10.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/SamToFastq.jar \
  VALIDATION_STRINGENCY=LENIENT \
  INPUT=/gs/scratch/efournier/sb_cofactor_hr/A549/raw/chip-seq/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam \
  FASTQ=/gs/scratch/efournier/sb_cofactor_hr/A549/raw/chip-seq/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.single.fastq.gz
picard_sam_to_fastq.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.7f36c3719b4f2804418ba8ef84bd9c10.mugqic.done
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
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam.cf3a93760b5cd9a94839fe24d214c999.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam.cf3a93760b5cd9a94839fe24d214c999.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/SamToFastq.jar \
  VALIDATION_STRINGENCY=LENIENT \
  INPUT=/gs/scratch/efournier/sb_cofactor_hr/A549/raw/chip-seq/HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam \
  FASTQ=/gs/scratch/efournier/sb_cofactor_hr/A549/raw/chip-seq/HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.single.fastq.gz
picard_sam_to_fastq.HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam.cf3a93760b5cd9a94839fe24d214c999.mugqic.done
)
picard_sam_to_fastq_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=3 | grep "[0-9]")
echo "$picard_sam_to_fastq_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_8_JOB_ID: picard_sam_to_fastq.HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam.3adfed36bd9af721dd8057d10b46c7b4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam.3adfed36bd9af721dd8057d10b46c7b4.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/SamToFastq.jar \
  VALIDATION_STRINGENCY=LENIENT \
  INPUT=/gs/scratch/efournier/sb_cofactor_hr/A549/raw/chip-seq/HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam \
  FASTQ=/gs/scratch/efournier/sb_cofactor_hr/A549/raw/chip-seq/HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.single.fastq.gz
picard_sam_to_fastq.HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam.3adfed36bd9af721dd8057d10b46c7b4.mugqic.done
)
picard_sam_to_fastq_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=3 | grep "[0-9]")
echo "$picard_sam_to_fastq_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


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
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.bccd4034ca68da99f08b95ca918618d5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.bccd4034ca68da99f08b95ca918618d5.mugqic.done'
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
  /gs/scratch/efournier/sb_cofactor_hr/A549/raw/chip-seq/HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.single.fastq.gz \
  trim/A549_CTRL_BRD4_rep1/HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/A549_CTRL_BRD4_rep1/HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/A549_CTRL_BRD4_rep1/HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.trim.log
trimmomatic.HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.bccd4034ca68da99f08b95ca918618d5.mugqic.done
)
trimmomatic_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_2_JOB_ID: trimmomatic.HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_2_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.92e08fb7f5899de183452841605ebf98.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.92e08fb7f5899de183452841605ebf98.mugqic.done'
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
  /gs/scratch/efournier/sb_cofactor_hr/A549/raw/chip-seq/HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.single.fastq.gz \
  trim/A549_CTRL_CDK9_rep1/HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/A549_CTRL_CDK9_rep1/HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/A549_CTRL_CDK9_rep1/HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.trim.log
trimmomatic.HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.92e08fb7f5899de183452841605ebf98.mugqic.done
)
trimmomatic_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_3_JOB_ID: trimmomatic.HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_3_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.91d89521af31294c3755ad4f80ee952a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.91d89521af31294c3755ad4f80ee952a.mugqic.done'
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
  /gs/scratch/efournier/sb_cofactor_hr/A549/raw/chip-seq/HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.single.fastq.gz \
  trim/A549_CTRL_MED1_rep1/HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/A549_CTRL_MED1_rep1/HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/A549_CTRL_MED1_rep1/HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.trim.log
trimmomatic.HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.91d89521af31294c3755ad4f80ee952a.mugqic.done
)
trimmomatic_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_4_JOB_ID: trimmomatic.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_4_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.ae9569d021d42a8e045c3243cab91c27.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.ae9569d021d42a8e045c3243cab91c27.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/A549_CTRL_MED1_rep2 && \
`cat > trim/A549_CTRL_MED1_rep2/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/sb_cofactor_hr/A549/raw/chip-seq/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.single.fastq.gz \
  trim/A549_CTRL_MED1_rep2/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/A549_CTRL_MED1_rep2/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/A549_CTRL_MED1_rep2/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.trim.log
trimmomatic.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.ae9569d021d42a8e045c3243cab91c27.mugqic.done
)
trimmomatic_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_5_JOB_ID: trimmomatic.HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_5_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.d4351975a49dbd5495d3c480d6757fdf.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.d4351975a49dbd5495d3c480d6757fdf.mugqic.done'
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
  /gs/scratch/efournier/sb_cofactor_hr/A549/raw/chip-seq/HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.single.fastq.gz \
  trim/A549_CTRL_SMC1A_rep1/HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/A549_CTRL_SMC1A_rep1/HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/A549_CTRL_SMC1A_rep1/HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.trim.log
trimmomatic.HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.d4351975a49dbd5495d3c480d6757fdf.mugqic.done
)
trimmomatic_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_6_JOB_ID: trimmomatic.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_6_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.174c08b796f5f4f098931d5ef384c3e0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.174c08b796f5f4f098931d5ef384c3e0.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/A549_CTRL_SMC1A_rep2 && \
`cat > trim/A549_CTRL_SMC1A_rep2/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/sb_cofactor_hr/A549/raw/chip-seq/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.single.fastq.gz \
  trim/A549_CTRL_SMC1A_rep2/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/A549_CTRL_SMC1A_rep2/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/A549_CTRL_SMC1A_rep2/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.trim.log
trimmomatic.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.174c08b796f5f4f098931d5ef384c3e0.mugqic.done
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
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam.4c92bdc610214b629727fafaea059445.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam.4c92bdc610214b629727fafaea059445.mugqic.done'
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
  /gs/scratch/efournier/sb_cofactor_hr/A549/raw/chip-seq/HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.single.fastq.gz \
  trim/A549_CTRL_WCE_rep1/HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/A549_CTRL_WCE_rep1/HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/A549_CTRL_WCE_rep1/HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam.trim.log
trimmomatic.HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam.4c92bdc610214b629727fafaea059445.mugqic.done
)
trimmomatic_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_8_JOB_ID: trimmomatic.HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_8_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam.6c02d671e74462627f0a445c2584f0e6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam.6c02d671e74462627f0a445c2584f0e6.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/A549_CTRL_WCE_rep2 && \
`cat > trim/A549_CTRL_WCE_rep2/HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/sb_cofactor_hr/A549/raw/chip-seq/HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.single.fastq.gz \
  trim/A549_CTRL_WCE_rep2/HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/A549_CTRL_WCE_rep2/HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/A549_CTRL_WCE_rep2/HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam.trim.log
trimmomatic.HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam.6c02d671e74462627f0a445c2584f0e6.mugqic.done
)
trimmomatic_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: merge_trimmomatic_stats
#-------------------------------------------------------------------------------
STEP=merge_trimmomatic_stats
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: merge_trimmomatic_stats_1_JOB_ID: merge_trimmomatic_stats
#-------------------------------------------------------------------------------
JOB_NAME=merge_trimmomatic_stats
JOB_DEPENDENCIES=$trimmomatic_1_JOB_ID:$trimmomatic_2_JOB_ID:$trimmomatic_3_JOB_ID:$trimmomatic_4_JOB_ID:$trimmomatic_5_JOB_ID:$trimmomatic_6_JOB_ID:$trimmomatic_7_JOB_ID:$trimmomatic_8_JOB_ID
JOB_DONE=job_output/merge_trimmomatic_stats/merge_trimmomatic_stats.14536a4b5acc24eea6ee277440416bbc.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'merge_trimmomatic_stats.14536a4b5acc24eea6ee277440416bbc.mugqic.done'
module load mugqic/pandoc/1.15.2 && \
mkdir -p metrics && \
echo 'Sample	Readset	Raw Single Reads #	Surviving Single Reads #	Surviving Single Reads %' > metrics/trimReadsetTable.tsv && \
grep ^Input trim/A549_CTRL_BRD4_rep1/HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/A549_CTRL_BRD4_rep1	HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam	\1	\2/' | \
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
grep ^Input trim/A549_CTRL_MED1_rep2/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/A549_CTRL_MED1_rep2	HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/A549_CTRL_SMC1A_rep1/HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/A549_CTRL_SMC1A_rep1	HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/A549_CTRL_SMC1A_rep2/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/A549_CTRL_SMC1A_rep2	HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/A549_CTRL_WCE_rep1/HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/A549_CTRL_WCE_rep1	HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/A549_CTRL_WCE_rep2/HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/A549_CTRL_WCE_rep2	HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam	\1	\2/' | \
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
merge_trimmomatic_stats.14536a4b5acc24eea6ee277440416bbc.mugqic.done
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
# JOB: bwa_mem_picard_sort_sam_1_JOB_ID: bwa_mem_picard_sort_sam.HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_1_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.fe08c6d21b5b753763e5df242594fb7f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.fe08c6d21b5b753763e5df242594fb7f.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/A549_CTRL_BRD4_rep1/HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam	SM:A549_CTRL_BRD4_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/A549_CTRL_BRD4_rep1/HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/A549_CTRL_BRD4_rep1/HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam/HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.fe08c6d21b5b753763e5df242594fb7f.mugqic.done
)
bwa_mem_picard_sort_sam_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_2_JOB_ID: bwa_mem_picard_sort_sam.HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_2_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.1fb594a27602e282fa8ac81a1631bf7b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.1fb594a27602e282fa8ac81a1631bf7b.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/A549_CTRL_CDK9_rep1/HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam	SM:A549_CTRL_CDK9_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/A549_CTRL_CDK9_rep1/HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/A549_CTRL_CDK9_rep1/HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam/HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.1fb594a27602e282fa8ac81a1631bf7b.mugqic.done
)
bwa_mem_picard_sort_sam_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_3_JOB_ID: bwa_mem_picard_sort_sam.HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_3_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.0771535ca4b38281dc7fe76506b9c90c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.0771535ca4b38281dc7fe76506b9c90c.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/A549_CTRL_MED1_rep1/HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam	SM:A549_CTRL_MED1_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/A549_CTRL_MED1_rep1/HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/A549_CTRL_MED1_rep1/HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam/HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.0771535ca4b38281dc7fe76506b9c90c.mugqic.done
)
bwa_mem_picard_sort_sam_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_4_JOB_ID: bwa_mem_picard_sort_sam.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_4_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.e43b16d6ee37c9f02bd502ba302fab37.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.e43b16d6ee37c9f02bd502ba302fab37.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/A549_CTRL_MED1_rep2/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam	SM:A549_CTRL_MED1_rep2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/A549_CTRL_MED1_rep2/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/A549_CTRL_MED1_rep2/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.e43b16d6ee37c9f02bd502ba302fab37.mugqic.done
)
bwa_mem_picard_sort_sam_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_5_JOB_ID: bwa_mem_picard_sort_sam.HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_5_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.0fb97cd2e1e85cd9844d34f1816efe86.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.0fb97cd2e1e85cd9844d34f1816efe86.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/A549_CTRL_SMC1A_rep1/HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam	SM:A549_CTRL_SMC1A_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/A549_CTRL_SMC1A_rep1/HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/A549_CTRL_SMC1A_rep1/HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam/HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.0fb97cd2e1e85cd9844d34f1816efe86.mugqic.done
)
bwa_mem_picard_sort_sam_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_6_JOB_ID: bwa_mem_picard_sort_sam.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_6_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.e687a77a8b331d608bccaafbf891f857.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.e687a77a8b331d608bccaafbf891f857.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/A549_CTRL_SMC1A_rep2/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam	SM:A549_CTRL_SMC1A_rep2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/A549_CTRL_SMC1A_rep2/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/A549_CTRL_SMC1A_rep2/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.e687a77a8b331d608bccaafbf891f857.mugqic.done
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
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam.9e88a087f1e50cef607a6494c536cad8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam.9e88a087f1e50cef607a6494c536cad8.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/A549_CTRL_WCE_rep1/HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam	SM:A549_CTRL_WCE_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/A549_CTRL_WCE_rep1/HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/A549_CTRL_WCE_rep1/HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam/HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam.9e88a087f1e50cef607a6494c536cad8.mugqic.done
)
bwa_mem_picard_sort_sam_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_8_JOB_ID: bwa_mem_picard_sort_sam.HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_8_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam.8c6fb2684ce810b322d6682b9b769093.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam.8c6fb2684ce810b322d6682b9b769093.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/A549_CTRL_WCE_rep2/HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam	SM:A549_CTRL_WCE_rep2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/A549_CTRL_WCE_rep2/HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/A549_CTRL_WCE_rep2/HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam/HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam.8c6fb2684ce810b322d6682b9b769093.mugqic.done
)
bwa_mem_picard_sort_sam_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_9_JOB_ID: bwa_mem_picard_sort_sam_report
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam_report
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_1_JOB_ID:$bwa_mem_picard_sort_sam_2_JOB_ID:$bwa_mem_picard_sort_sam_3_JOB_ID:$bwa_mem_picard_sort_sam_4_JOB_ID:$bwa_mem_picard_sort_sam_5_JOB_ID:$bwa_mem_picard_sort_sam_6_JOB_ID:$bwa_mem_picard_sort_sam_7_JOB_ID:$bwa_mem_picard_sort_sam_8_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam_report.b10d0f7eacf0b8b59354726b32e5aa9c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam_report.b10d0f7eacf0b8b59354726b32e5aa9c.mugqic.done'
module load mugqic/pandoc/1.15.2 && \
mkdir -p report && \
pandoc --to=markdown \
  --template /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/DnaSeq.bwa_mem_picard_sort_sam.md \
  --variable scientific_name="Homo_sapiens" \
  --variable assembly="GRCh38" \
  /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/DnaSeq.bwa_mem_picard_sort_sam.md \
  > report/DnaSeq.bwa_mem_picard_sort_sam.md
bwa_mem_picard_sort_sam_report.b10d0f7eacf0b8b59354726b32e5aa9c.mugqic.done
)
bwa_mem_picard_sort_sam_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


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
samtools_view_filter_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_2_JOB_ID: samtools_view_filter.HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_2_JOB_ID
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
samtools_view_filter_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_3_JOB_ID: samtools_view_filter.HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_3_JOB_ID
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
samtools_view_filter_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_4_JOB_ID: samtools_view_filter.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_4_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.fc1d0f403db60848b20373f3c938aacc.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.fc1d0f403db60848b20373f3c938aacc.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/A549_CTRL_MED1_rep2/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.sorted.bam \
  > alignment/A549_CTRL_MED1_rep2/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.fc1d0f403db60848b20373f3c938aacc.mugqic.done
)
samtools_view_filter_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_5_JOB_ID: samtools_view_filter.HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_5_JOB_ID
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
samtools_view_filter_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_6_JOB_ID: samtools_view_filter.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_6_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.98753d32e43009e83dd33caf52999ea6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.98753d32e43009e83dd33caf52999ea6.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/A549_CTRL_SMC1A_rep2/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.sorted.bam \
  > alignment/A549_CTRL_SMC1A_rep2/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.98753d32e43009e83dd33caf52999ea6.mugqic.done
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
# JOB: samtools_view_filter_8_JOB_ID: samtools_view_filter.HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_8_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam.aebd3a9429add22730c8af7809028c3d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam.aebd3a9429add22730c8af7809028c3d.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/A549_CTRL_WCE_rep2/HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam/HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam.sorted.bam \
  > alignment/A549_CTRL_WCE_rep2/HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam/HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam.aebd3a9429add22730c8af7809028c3d.mugqic.done
)
samtools_view_filter_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_9_JOB_ID: samtools_view_filter_report
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter_report
JOB_DEPENDENCIES=$samtools_view_filter_1_JOB_ID:$samtools_view_filter_2_JOB_ID:$samtools_view_filter_3_JOB_ID:$samtools_view_filter_4_JOB_ID:$samtools_view_filter_5_JOB_ID:$samtools_view_filter_6_JOB_ID:$samtools_view_filter_7_JOB_ID:$samtools_view_filter_8_JOB_ID
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
samtools_view_filter_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: picard_merge_sam_files
#-------------------------------------------------------------------------------
STEP=picard_merge_sam_files
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_1_JOB_ID: symlink_readset_sample_bam.A549_CTRL_BRD4_rep1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.A549_CTRL_BRD4_rep1
JOB_DEPENDENCIES=$samtools_view_filter_1_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.A549_CTRL_BRD4_rep1.122782b49de4447e23add94e5b7d5d28.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.A549_CTRL_BRD4_rep1.122782b49de4447e23add94e5b7d5d28.mugqic.done'
mkdir -p alignment/A549_CTRL_BRD4_rep1 && \
ln -s -f HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam/HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.sorted.filtered.bam alignment/A549_CTRL_BRD4_rep1/A549_CTRL_BRD4_rep1.merged.bam
symlink_readset_sample_bam.A549_CTRL_BRD4_rep1.122782b49de4447e23add94e5b7d5d28.mugqic.done
)
picard_merge_sam_files_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_2_JOB_ID: symlink_readset_sample_bam.A549_CTRL_CDK9_rep1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.A549_CTRL_CDK9_rep1
JOB_DEPENDENCIES=$samtools_view_filter_2_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.A549_CTRL_CDK9_rep1.e37c74f70c72e58da1e581952f7395b1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.A549_CTRL_CDK9_rep1.e37c74f70c72e58da1e581952f7395b1.mugqic.done'
mkdir -p alignment/A549_CTRL_CDK9_rep1 && \
ln -s -f HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam/HI.1613.004.Index_16.A549_ctrl_CDK9_MF_ChIP1_4_rep1.bam.sorted.filtered.bam alignment/A549_CTRL_CDK9_rep1/A549_CTRL_CDK9_rep1.merged.bam
symlink_readset_sample_bam.A549_CTRL_CDK9_rep1.e37c74f70c72e58da1e581952f7395b1.mugqic.done
)
picard_merge_sam_files_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_3_JOB_ID: symlink_readset_sample_bam.A549_CTRL_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.A549_CTRL_MED1_rep1
JOB_DEPENDENCIES=$samtools_view_filter_3_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.A549_CTRL_MED1_rep1.93b78a91a30a2b05a8ec84948942bd09.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.A549_CTRL_MED1_rep1.93b78a91a30a2b05a8ec84948942bd09.mugqic.done'
mkdir -p alignment/A549_CTRL_MED1_rep1 && \
ln -s -f HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam/HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.sorted.filtered.bam alignment/A549_CTRL_MED1_rep1/A549_CTRL_MED1_rep1.merged.bam
symlink_readset_sample_bam.A549_CTRL_MED1_rep1.93b78a91a30a2b05a8ec84948942bd09.mugqic.done
)
picard_merge_sam_files_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_4_JOB_ID: symlink_readset_sample_bam.A549_CTRL_MED1_rep2
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.A549_CTRL_MED1_rep2
JOB_DEPENDENCIES=$samtools_view_filter_4_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.A549_CTRL_MED1_rep2.0dd9318287ce9f4696134ffa8bca3173.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.A549_CTRL_MED1_rep2.0dd9318287ce9f4696134ffa8bca3173.mugqic.done'
mkdir -p alignment/A549_CTRL_MED1_rep2 && \
ln -s -f HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.sorted.filtered.bam alignment/A549_CTRL_MED1_rep2/A549_CTRL_MED1_rep2.merged.bam
symlink_readset_sample_bam.A549_CTRL_MED1_rep2.0dd9318287ce9f4696134ffa8bca3173.mugqic.done
)
picard_merge_sam_files_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_5_JOB_ID: symlink_readset_sample_bam.A549_CTRL_SMC1A_rep1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.A549_CTRL_SMC1A_rep1
JOB_DEPENDENCIES=$samtools_view_filter_5_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.A549_CTRL_SMC1A_rep1.301863dc2fbaabe3fa9463fd1a0ab344.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.A549_CTRL_SMC1A_rep1.301863dc2fbaabe3fa9463fd1a0ab344.mugqic.done'
mkdir -p alignment/A549_CTRL_SMC1A_rep1 && \
ln -s -f HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam/HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.sorted.filtered.bam alignment/A549_CTRL_SMC1A_rep1/A549_CTRL_SMC1A_rep1.merged.bam
symlink_readset_sample_bam.A549_CTRL_SMC1A_rep1.301863dc2fbaabe3fa9463fd1a0ab344.mugqic.done
)
picard_merge_sam_files_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_6_JOB_ID: symlink_readset_sample_bam.A549_CTRL_SMC1A_rep2
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.A549_CTRL_SMC1A_rep2
JOB_DEPENDENCIES=$samtools_view_filter_6_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.A549_CTRL_SMC1A_rep2.e14f5c2700d1c3b7065db2d93c834ed5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.A549_CTRL_SMC1A_rep2.e14f5c2700d1c3b7065db2d93c834ed5.mugqic.done'
mkdir -p alignment/A549_CTRL_SMC1A_rep2 && \
ln -s -f HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.sorted.filtered.bam alignment/A549_CTRL_SMC1A_rep2/A549_CTRL_SMC1A_rep2.merged.bam
symlink_readset_sample_bam.A549_CTRL_SMC1A_rep2.e14f5c2700d1c3b7065db2d93c834ed5.mugqic.done
)
picard_merge_sam_files_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_7_JOB_ID: symlink_readset_sample_bam.A549_CTRL_WCE_rep1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.A549_CTRL_WCE_rep1
JOB_DEPENDENCIES=$samtools_view_filter_7_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.A549_CTRL_WCE_rep1.55133cc6a5e4bf8d3079533ed1955e28.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.A549_CTRL_WCE_rep1.55133cc6a5e4bf8d3079533ed1955e28.mugqic.done'
mkdir -p alignment/A549_CTRL_WCE_rep1 && \
ln -s -f HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam/HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam.sorted.filtered.bam alignment/A549_CTRL_WCE_rep1/A549_CTRL_WCE_rep1.merged.bam
symlink_readset_sample_bam.A549_CTRL_WCE_rep1.55133cc6a5e4bf8d3079533ed1955e28.mugqic.done
)
picard_merge_sam_files_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_8_JOB_ID: symlink_readset_sample_bam.A549_CTRL_WCE_rep2
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.A549_CTRL_WCE_rep2
JOB_DEPENDENCIES=$samtools_view_filter_8_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.A549_CTRL_WCE_rep2.7a3c1b2ed4bf50855ad5df5f445f6478.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.A549_CTRL_WCE_rep2.7a3c1b2ed4bf50855ad5df5f445f6478.mugqic.done'
mkdir -p alignment/A549_CTRL_WCE_rep2 && \
ln -s -f HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam/HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam.sorted.filtered.bam alignment/A549_CTRL_WCE_rep2/A549_CTRL_WCE_rep2.merged.bam
symlink_readset_sample_bam.A549_CTRL_WCE_rep2.7a3c1b2ed4bf50855ad5df5f445f6478.mugqic.done
)
picard_merge_sam_files_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: picard_mark_duplicates
#-------------------------------------------------------------------------------
STEP=picard_mark_duplicates
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_1_JOB_ID: picard_mark_duplicates.A549_CTRL_BRD4_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.A549_CTRL_BRD4_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_1_JOB_ID
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
picard_mark_duplicates_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_2_JOB_ID: picard_mark_duplicates.A549_CTRL_CDK9_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.A549_CTRL_CDK9_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_2_JOB_ID
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
picard_mark_duplicates_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_3_JOB_ID: picard_mark_duplicates.A549_CTRL_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.A549_CTRL_MED1_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_3_JOB_ID
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
picard_mark_duplicates_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_4_JOB_ID: picard_mark_duplicates.A549_CTRL_MED1_rep2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.A549_CTRL_MED1_rep2
JOB_DEPENDENCIES=$picard_merge_sam_files_4_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.A549_CTRL_MED1_rep2.6f3905b4ed731df98d72958be933143e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.A549_CTRL_MED1_rep2.6f3905b4ed731df98d72958be933143e.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/A549_CTRL_MED1_rep2/A549_CTRL_MED1_rep2.merged.bam \
  OUTPUT=alignment/A549_CTRL_MED1_rep2/A549_CTRL_MED1_rep2.sorted.dup.bam \
  METRICS_FILE=alignment/A549_CTRL_MED1_rep2/A549_CTRL_MED1_rep2.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.A549_CTRL_MED1_rep2.6f3905b4ed731df98d72958be933143e.mugqic.done
)
picard_mark_duplicates_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_5_JOB_ID: picard_mark_duplicates.A549_CTRL_SMC1A_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.A549_CTRL_SMC1A_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_5_JOB_ID
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
picard_mark_duplicates_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_6_JOB_ID: picard_mark_duplicates.A549_CTRL_SMC1A_rep2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.A549_CTRL_SMC1A_rep2
JOB_DEPENDENCIES=$picard_merge_sam_files_6_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.A549_CTRL_SMC1A_rep2.f9546d32bb0406861c6d37db880c35ea.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.A549_CTRL_SMC1A_rep2.f9546d32bb0406861c6d37db880c35ea.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/A549_CTRL_SMC1A_rep2/A549_CTRL_SMC1A_rep2.merged.bam \
  OUTPUT=alignment/A549_CTRL_SMC1A_rep2/A549_CTRL_SMC1A_rep2.sorted.dup.bam \
  METRICS_FILE=alignment/A549_CTRL_SMC1A_rep2/A549_CTRL_SMC1A_rep2.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.A549_CTRL_SMC1A_rep2.f9546d32bb0406861c6d37db880c35ea.mugqic.done
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
# JOB: picard_mark_duplicates_8_JOB_ID: picard_mark_duplicates.A549_CTRL_WCE_rep2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.A549_CTRL_WCE_rep2
JOB_DEPENDENCIES=$picard_merge_sam_files_8_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.A549_CTRL_WCE_rep2.3aa5c0ad4c853af4f6b6c8ee3b0e9a46.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.A549_CTRL_WCE_rep2.3aa5c0ad4c853af4f6b6c8ee3b0e9a46.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/A549_CTRL_WCE_rep2/A549_CTRL_WCE_rep2.merged.bam \
  OUTPUT=alignment/A549_CTRL_WCE_rep2/A549_CTRL_WCE_rep2.sorted.dup.bam \
  METRICS_FILE=alignment/A549_CTRL_WCE_rep2/A549_CTRL_WCE_rep2.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.A549_CTRL_WCE_rep2.3aa5c0ad4c853af4f6b6c8ee3b0e9a46.mugqic.done
)
picard_mark_duplicates_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_9_JOB_ID: picard_mark_duplicates_report
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates_report
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_2_JOB_ID:$picard_mark_duplicates_3_JOB_ID:$picard_mark_duplicates_4_JOB_ID:$picard_mark_duplicates_5_JOB_ID:$picard_mark_duplicates_6_JOB_ID:$picard_mark_duplicates_7_JOB_ID:$picard_mark_duplicates_8_JOB_ID
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
picard_mark_duplicates_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: metrics
#-------------------------------------------------------------------------------
STEP=metrics
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: metrics_1_JOB_ID: metrics.flagstat
#-------------------------------------------------------------------------------
JOB_NAME=metrics.flagstat
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_2_JOB_ID:$picard_mark_duplicates_3_JOB_ID:$picard_mark_duplicates_4_JOB_ID:$picard_mark_duplicates_5_JOB_ID:$picard_mark_duplicates_6_JOB_ID:$picard_mark_duplicates_7_JOB_ID:$picard_mark_duplicates_8_JOB_ID
JOB_DONE=job_output/metrics/metrics.flagstat.4b66a58c70c4e7131b439ccc81e3dfc9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'metrics.flagstat.4b66a58c70c4e7131b439ccc81e3dfc9.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools flagstat \
  alignment/A549_CTRL_BRD4_rep1/A549_CTRL_BRD4_rep1.sorted.dup.bam \
  > alignment/A549_CTRL_BRD4_rep1/A549_CTRL_BRD4_rep1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/A549_CTRL_CDK9_rep1/A549_CTRL_CDK9_rep1.sorted.dup.bam \
  > alignment/A549_CTRL_CDK9_rep1/A549_CTRL_CDK9_rep1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/A549_CTRL_MED1_rep1/A549_CTRL_MED1_rep1.sorted.dup.bam \
  > alignment/A549_CTRL_MED1_rep1/A549_CTRL_MED1_rep1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/A549_CTRL_MED1_rep2/A549_CTRL_MED1_rep2.sorted.dup.bam \
  > alignment/A549_CTRL_MED1_rep2/A549_CTRL_MED1_rep2.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/A549_CTRL_SMC1A_rep1/A549_CTRL_SMC1A_rep1.sorted.dup.bam \
  > alignment/A549_CTRL_SMC1A_rep1/A549_CTRL_SMC1A_rep1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/A549_CTRL_SMC1A_rep2/A549_CTRL_SMC1A_rep2.sorted.dup.bam \
  > alignment/A549_CTRL_SMC1A_rep2/A549_CTRL_SMC1A_rep2.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/A549_CTRL_WCE_rep1/A549_CTRL_WCE_rep1.sorted.dup.bam \
  > alignment/A549_CTRL_WCE_rep1/A549_CTRL_WCE_rep1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/A549_CTRL_WCE_rep2/A549_CTRL_WCE_rep2.sorted.dup.bam \
  > alignment/A549_CTRL_WCE_rep2/A549_CTRL_WCE_rep2.sorted.dup.bam.flagstat
metrics.flagstat.4b66a58c70c4e7131b439ccc81e3dfc9.mugqic.done
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
JOB_DONE=job_output/metrics/metrics_report.5c4f199f340c312dbc00abce711dc56e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'metrics_report.5c4f199f340c312dbc00abce711dc56e.mugqic.done'
module load mugqic/pandoc/1.15.2 && \
for sample in A549_CTRL_BRD4_rep1 A549_CTRL_CDK9_rep1 A549_CTRL_MED1_rep1 A549_CTRL_MED1_rep2 A549_CTRL_SMC1A_rep1 A549_CTRL_SMC1A_rep2 A549_CTRL_WCE_rep1 A549_CTRL_WCE_rep2
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

metrics_report.5c4f199f340c312dbc00abce711dc56e.mugqic.done
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
# JOB: homer_make_tag_directory_1_JOB_ID: homer_make_tag_directory.A549_CTRL_BRD4_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.A549_CTRL_BRD4_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.A549_CTRL_BRD4_rep1.802187dfbf60113ffbd3f42bc58cd8db.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.A549_CTRL_BRD4_rep1.802187dfbf60113ffbd3f42bc58cd8db.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/A549_CTRL_BRD4_rep1 \
  alignment/A549_CTRL_BRD4_rep1/A549_CTRL_BRD4_rep1.sorted.dup.bam \
  -checkGC -genome GRCh38
homer_make_tag_directory.A549_CTRL_BRD4_rep1.802187dfbf60113ffbd3f42bc58cd8db.mugqic.done
)
homer_make_tag_directory_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_2_JOB_ID: homer_make_tag_directory.A549_CTRL_CDK9_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.A549_CTRL_CDK9_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_2_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.A549_CTRL_CDK9_rep1.15aef0cac5224951f47fae753cc060f7.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.A549_CTRL_CDK9_rep1.15aef0cac5224951f47fae753cc060f7.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/A549_CTRL_CDK9_rep1 \
  alignment/A549_CTRL_CDK9_rep1/A549_CTRL_CDK9_rep1.sorted.dup.bam \
  -checkGC -genome GRCh38
homer_make_tag_directory.A549_CTRL_CDK9_rep1.15aef0cac5224951f47fae753cc060f7.mugqic.done
)
homer_make_tag_directory_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_3_JOB_ID: homer_make_tag_directory.A549_CTRL_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.A549_CTRL_MED1_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_3_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.A549_CTRL_MED1_rep1.8ca61a68a9a240f2a7adb896a7011671.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.A549_CTRL_MED1_rep1.8ca61a68a9a240f2a7adb896a7011671.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/A549_CTRL_MED1_rep1 \
  alignment/A549_CTRL_MED1_rep1/A549_CTRL_MED1_rep1.sorted.dup.bam \
  -checkGC -genome GRCh38
homer_make_tag_directory.A549_CTRL_MED1_rep1.8ca61a68a9a240f2a7adb896a7011671.mugqic.done
)
homer_make_tag_directory_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_4_JOB_ID: homer_make_tag_directory.A549_CTRL_MED1_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.A549_CTRL_MED1_rep2
JOB_DEPENDENCIES=$picard_mark_duplicates_4_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.A549_CTRL_MED1_rep2.f8bd857eb1cb06ca5fcde04785ecc3ac.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.A549_CTRL_MED1_rep2.f8bd857eb1cb06ca5fcde04785ecc3ac.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/A549_CTRL_MED1_rep2 \
  alignment/A549_CTRL_MED1_rep2/A549_CTRL_MED1_rep2.sorted.dup.bam \
  -checkGC -genome GRCh38
homer_make_tag_directory.A549_CTRL_MED1_rep2.f8bd857eb1cb06ca5fcde04785ecc3ac.mugqic.done
)
homer_make_tag_directory_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_5_JOB_ID: homer_make_tag_directory.A549_CTRL_SMC1A_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.A549_CTRL_SMC1A_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_5_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.A549_CTRL_SMC1A_rep1.156ca44eddec7c1a7e653d82944162d3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.A549_CTRL_SMC1A_rep1.156ca44eddec7c1a7e653d82944162d3.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/A549_CTRL_SMC1A_rep1 \
  alignment/A549_CTRL_SMC1A_rep1/A549_CTRL_SMC1A_rep1.sorted.dup.bam \
  -checkGC -genome GRCh38
homer_make_tag_directory.A549_CTRL_SMC1A_rep1.156ca44eddec7c1a7e653d82944162d3.mugqic.done
)
homer_make_tag_directory_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_6_JOB_ID: homer_make_tag_directory.A549_CTRL_SMC1A_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.A549_CTRL_SMC1A_rep2
JOB_DEPENDENCIES=$picard_mark_duplicates_6_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.A549_CTRL_SMC1A_rep2.76e5e1e1de1e5f98cc4fb092c9ca1512.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.A549_CTRL_SMC1A_rep2.76e5e1e1de1e5f98cc4fb092c9ca1512.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/A549_CTRL_SMC1A_rep2 \
  alignment/A549_CTRL_SMC1A_rep2/A549_CTRL_SMC1A_rep2.sorted.dup.bam \
  -checkGC -genome GRCh38
homer_make_tag_directory.A549_CTRL_SMC1A_rep2.76e5e1e1de1e5f98cc4fb092c9ca1512.mugqic.done
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
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.A549_CTRL_WCE_rep1.d5dcede8e2bc38efd89b15d795b7c4a4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.A549_CTRL_WCE_rep1.d5dcede8e2bc38efd89b15d795b7c4a4.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/A549_CTRL_WCE_rep1 \
  alignment/A549_CTRL_WCE_rep1/A549_CTRL_WCE_rep1.sorted.dup.bam \
  -checkGC -genome GRCh38
homer_make_tag_directory.A549_CTRL_WCE_rep1.d5dcede8e2bc38efd89b15d795b7c4a4.mugqic.done
)
homer_make_tag_directory_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_8_JOB_ID: homer_make_tag_directory.A549_CTRL_WCE_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.A549_CTRL_WCE_rep2
JOB_DEPENDENCIES=$picard_mark_duplicates_8_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.A549_CTRL_WCE_rep2.6264c2b85603257df48234695f74eb40.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.A549_CTRL_WCE_rep2.6264c2b85603257df48234695f74eb40.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/A549_CTRL_WCE_rep2 \
  alignment/A549_CTRL_WCE_rep2/A549_CTRL_WCE_rep2.sorted.dup.bam \
  -checkGC -genome GRCh38
homer_make_tag_directory.A549_CTRL_WCE_rep2.6264c2b85603257df48234695f74eb40.mugqic.done
)
homer_make_tag_directory_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: qc_metrics
#-------------------------------------------------------------------------------
STEP=qc_metrics
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: qc_metrics_1_JOB_ID: qc_plots_R
#-------------------------------------------------------------------------------
JOB_NAME=qc_plots_R
JOB_DEPENDENCIES=$homer_make_tag_directory_1_JOB_ID:$homer_make_tag_directory_2_JOB_ID:$homer_make_tag_directory_3_JOB_ID:$homer_make_tag_directory_4_JOB_ID:$homer_make_tag_directory_5_JOB_ID:$homer_make_tag_directory_6_JOB_ID:$homer_make_tag_directory_7_JOB_ID:$homer_make_tag_directory_8_JOB_ID
JOB_DONE=job_output/qc_metrics/qc_plots_R.cc9b2ae86f6165897050416a90b26f68.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'qc_plots_R.cc9b2ae86f6165897050416a90b26f68.mugqic.done'
module load mugqic/mugqic_tools/2.1.5 mugqic/R_Bioconductor/3.2.3_3.2 && \
mkdir -p graphs && \
Rscript $R_TOOLS/chipSeqGenerateQCMetrics.R \
  ../../raw/chip-seq/design_partial.txt \
  /gs/scratch/efournier/sb_cofactor_hr/A549/output/chip-pipeline-GRCh38 && \
cp /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/ChipSeq.qc_metrics.md report/ChipSeq.qc_metrics.md && \
for sample in A549_CTRL_BRD4_rep1 A549_CTRL_CDK9_rep1 A549_CTRL_MED1_rep1 A549_CTRL_MED1_rep2 A549_CTRL_SMC1A_rep1 A549_CTRL_SMC1A_rep2 A549_CTRL_WCE_rep1 A549_CTRL_WCE_rep2
do
  cp --parents graphs/${sample}_QC_Metrics.ps report/
  convert -rotate 90 graphs/${sample}_QC_Metrics.ps report/graphs/${sample}_QC_Metrics.png
  echo -e "----

![QC Metrics for Sample $sample ([download high-res image](graphs/${sample}_QC_Metrics.ps))](graphs/${sample}_QC_Metrics.png)
" \
  >> report/ChipSeq.qc_metrics.md
done
qc_plots_R.cc9b2ae86f6165897050416a90b26f68.mugqic.done
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
# JOB: homer_make_ucsc_file_1_JOB_ID: homer_make_ucsc_file.A549_CTRL_BRD4_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.A549_CTRL_BRD4_rep1
JOB_DEPENDENCIES=$homer_make_tag_directory_1_JOB_ID
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
homer_make_ucsc_file_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_2_JOB_ID: homer_make_ucsc_file.A549_CTRL_CDK9_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.A549_CTRL_CDK9_rep1
JOB_DEPENDENCIES=$homer_make_tag_directory_2_JOB_ID
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
homer_make_ucsc_file_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_3_JOB_ID: homer_make_ucsc_file.A549_CTRL_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.A549_CTRL_MED1_rep1
JOB_DEPENDENCIES=$homer_make_tag_directory_3_JOB_ID
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
homer_make_ucsc_file_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_4_JOB_ID: homer_make_ucsc_file.A549_CTRL_MED1_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.A549_CTRL_MED1_rep2
JOB_DEPENDENCIES=$homer_make_tag_directory_4_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.A549_CTRL_MED1_rep2.aff034630438187f951ef9e8ac40be28.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.A549_CTRL_MED1_rep2.aff034630438187f951ef9e8ac40be28.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/A549_CTRL_MED1_rep2 && \
makeUCSCfile \
  tags/A549_CTRL_MED1_rep2 | \
gzip -1 -c > tracks/A549_CTRL_MED1_rep2/A549_CTRL_MED1_rep2.ucsc.bedGraph.gz
homer_make_ucsc_file.A549_CTRL_MED1_rep2.aff034630438187f951ef9e8ac40be28.mugqic.done
)
homer_make_ucsc_file_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_5_JOB_ID: homer_make_ucsc_file.A549_CTRL_SMC1A_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.A549_CTRL_SMC1A_rep1
JOB_DEPENDENCIES=$homer_make_tag_directory_5_JOB_ID
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
homer_make_ucsc_file_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_6_JOB_ID: homer_make_ucsc_file.A549_CTRL_SMC1A_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.A549_CTRL_SMC1A_rep2
JOB_DEPENDENCIES=$homer_make_tag_directory_6_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.A549_CTRL_SMC1A_rep2.e0e3a9cb9f8131fd0590b5bceff272da.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.A549_CTRL_SMC1A_rep2.e0e3a9cb9f8131fd0590b5bceff272da.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/A549_CTRL_SMC1A_rep2 && \
makeUCSCfile \
  tags/A549_CTRL_SMC1A_rep2 | \
gzip -1 -c > tracks/A549_CTRL_SMC1A_rep2/A549_CTRL_SMC1A_rep2.ucsc.bedGraph.gz
homer_make_ucsc_file.A549_CTRL_SMC1A_rep2.e0e3a9cb9f8131fd0590b5bceff272da.mugqic.done
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
# JOB: homer_make_ucsc_file_8_JOB_ID: homer_make_ucsc_file.A549_CTRL_WCE_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.A549_CTRL_WCE_rep2
JOB_DEPENDENCIES=$homer_make_tag_directory_8_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.A549_CTRL_WCE_rep2.1d135660d7ac619d28b5337014333765.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.A549_CTRL_WCE_rep2.1d135660d7ac619d28b5337014333765.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/A549_CTRL_WCE_rep2 && \
makeUCSCfile \
  tags/A549_CTRL_WCE_rep2 | \
gzip -1 -c > tracks/A549_CTRL_WCE_rep2/A549_CTRL_WCE_rep2.ucsc.bedGraph.gz
homer_make_ucsc_file.A549_CTRL_WCE_rep2.1d135660d7ac619d28b5337014333765.mugqic.done
)
homer_make_ucsc_file_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_9_JOB_ID: homer_make_ucsc_file_report
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file_report
JOB_DEPENDENCIES=$homer_make_ucsc_file_8_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file_report.c9cac04b3bc71f34cbcde042cf5b419d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file_report.c9cac04b3bc71f34cbcde042cf5b419d.mugqic.done'
mkdir -p report && \
zip -r report/tracks.zip tracks/*/*.ucsc.bedGraph.gz && \
cp /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/ChipSeq.homer_make_ucsc_file.md report/
homer_make_ucsc_file_report.c9cac04b3bc71f34cbcde042cf5b419d.mugqic.done
)
homer_make_ucsc_file_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: macs2_callpeak
#-------------------------------------------------------------------------------
STEP=macs2_callpeak
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_1_JOB_ID: macs2_callpeak.A549_CTRL_BRD4
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_CTRL_BRD4
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_7_JOB_ID:$picard_mark_duplicates_8_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_CTRL_BRD4.d734f63d867cf2c4ed1403d1dffe13ac.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_CTRL_BRD4.d734f63d867cf2c4ed1403d1dffe13ac.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/A549_CTRL_BRD4 && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2479938032.8 \
  --treatment \
  alignment/A549_CTRL_BRD4_rep1/A549_CTRL_BRD4_rep1.sorted.dup.bam \
  --control \
  alignment/A549_CTRL_WCE_rep1/A549_CTRL_WCE_rep1.sorted.dup.bam \
  alignment/A549_CTRL_WCE_rep2/A549_CTRL_WCE_rep2.sorted.dup.bam \
  --name peak_call/A549_CTRL_BRD4/A549_CTRL_BRD4 \
  >& peak_call/A549_CTRL_BRD4/A549_CTRL_BRD4.diag.macs.out
macs2_callpeak.A549_CTRL_BRD4.d734f63d867cf2c4ed1403d1dffe13ac.mugqic.done
)
macs2_callpeak_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_2_JOB_ID: macs2_callpeak.A549_CTRL_CDK9
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_CTRL_CDK9
JOB_DEPENDENCIES=$picard_mark_duplicates_2_JOB_ID:$picard_mark_duplicates_7_JOB_ID:$picard_mark_duplicates_8_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_CTRL_CDK9.1755e2909855cc55c4eb2b93ec1d760a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_CTRL_CDK9.1755e2909855cc55c4eb2b93ec1d760a.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/A549_CTRL_CDK9 && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2479938032.8 \
  --treatment \
  alignment/A549_CTRL_CDK9_rep1/A549_CTRL_CDK9_rep1.sorted.dup.bam \
  --control \
  alignment/A549_CTRL_WCE_rep1/A549_CTRL_WCE_rep1.sorted.dup.bam \
  alignment/A549_CTRL_WCE_rep2/A549_CTRL_WCE_rep2.sorted.dup.bam \
  --name peak_call/A549_CTRL_CDK9/A549_CTRL_CDK9 \
  >& peak_call/A549_CTRL_CDK9/A549_CTRL_CDK9.diag.macs.out
macs2_callpeak.A549_CTRL_CDK9.1755e2909855cc55c4eb2b93ec1d760a.mugqic.done
)
macs2_callpeak_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_3_JOB_ID: macs2_callpeak.A549_CTRL_MED1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_CTRL_MED1
JOB_DEPENDENCIES=$picard_mark_duplicates_3_JOB_ID:$picard_mark_duplicates_4_JOB_ID:$picard_mark_duplicates_7_JOB_ID:$picard_mark_duplicates_8_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_CTRL_MED1.0d3bd9f3e882228e39185d195aa32aff.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_CTRL_MED1.0d3bd9f3e882228e39185d195aa32aff.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/A549_CTRL_MED1 && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2479938032.8 \
  --treatment \
  alignment/A549_CTRL_MED1_rep1/A549_CTRL_MED1_rep1.sorted.dup.bam \
  alignment/A549_CTRL_MED1_rep2/A549_CTRL_MED1_rep2.sorted.dup.bam \
  --control \
  alignment/A549_CTRL_WCE_rep1/A549_CTRL_WCE_rep1.sorted.dup.bam \
  alignment/A549_CTRL_WCE_rep2/A549_CTRL_WCE_rep2.sorted.dup.bam \
  --name peak_call/A549_CTRL_MED1/A549_CTRL_MED1 \
  >& peak_call/A549_CTRL_MED1/A549_CTRL_MED1.diag.macs.out
macs2_callpeak.A549_CTRL_MED1.0d3bd9f3e882228e39185d195aa32aff.mugqic.done
)
macs2_callpeak_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_4_JOB_ID: macs2_callpeak.A549_CTRL_SMC1A
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_CTRL_SMC1A
JOB_DEPENDENCIES=$picard_mark_duplicates_5_JOB_ID:$picard_mark_duplicates_6_JOB_ID:$picard_mark_duplicates_7_JOB_ID:$picard_mark_duplicates_8_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_CTRL_SMC1A.b6feef6614901f18e40dd89a9bcd3dc7.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_CTRL_SMC1A.b6feef6614901f18e40dd89a9bcd3dc7.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/A549_CTRL_SMC1A && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2479938032.8 \
  --treatment \
  alignment/A549_CTRL_SMC1A_rep1/A549_CTRL_SMC1A_rep1.sorted.dup.bam \
  alignment/A549_CTRL_SMC1A_rep2/A549_CTRL_SMC1A_rep2.sorted.dup.bam \
  --control \
  alignment/A549_CTRL_WCE_rep1/A549_CTRL_WCE_rep1.sorted.dup.bam \
  alignment/A549_CTRL_WCE_rep2/A549_CTRL_WCE_rep2.sorted.dup.bam \
  --name peak_call/A549_CTRL_SMC1A/A549_CTRL_SMC1A \
  >& peak_call/A549_CTRL_SMC1A/A549_CTRL_SMC1A.diag.macs.out
macs2_callpeak.A549_CTRL_SMC1A.b6feef6614901f18e40dd89a9bcd3dc7.mugqic.done
)
macs2_callpeak_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_5_JOB_ID: macs2_callpeak.A549_CTRL_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_CTRL_MED1_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_3_JOB_ID:$picard_mark_duplicates_7_JOB_ID:$picard_mark_duplicates_8_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_CTRL_MED1_rep1.94f4aa0182548a088fb322b7a021b8ed.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_CTRL_MED1_rep1.94f4aa0182548a088fb322b7a021b8ed.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/A549_CTRL_MED1_rep1 && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2479938032.8 \
  --treatment \
  alignment/A549_CTRL_MED1_rep1/A549_CTRL_MED1_rep1.sorted.dup.bam \
  --control \
  alignment/A549_CTRL_WCE_rep1/A549_CTRL_WCE_rep1.sorted.dup.bam \
  alignment/A549_CTRL_WCE_rep2/A549_CTRL_WCE_rep2.sorted.dup.bam \
  --name peak_call/A549_CTRL_MED1_rep1/A549_CTRL_MED1_rep1 \
  >& peak_call/A549_CTRL_MED1_rep1/A549_CTRL_MED1_rep1.diag.macs.out
macs2_callpeak.A549_CTRL_MED1_rep1.94f4aa0182548a088fb322b7a021b8ed.mugqic.done
)
macs2_callpeak_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_6_JOB_ID: macs2_callpeak.A549_CTRL_SMC1A_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_CTRL_SMC1A_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_5_JOB_ID:$picard_mark_duplicates_7_JOB_ID:$picard_mark_duplicates_8_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_CTRL_SMC1A_rep1.5405b6d488e76cfc51aecda63958297d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_CTRL_SMC1A_rep1.5405b6d488e76cfc51aecda63958297d.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/A549_CTRL_SMC1A_rep1 && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2479938032.8 \
  --treatment \
  alignment/A549_CTRL_SMC1A_rep1/A549_CTRL_SMC1A_rep1.sorted.dup.bam \
  --control \
  alignment/A549_CTRL_WCE_rep1/A549_CTRL_WCE_rep1.sorted.dup.bam \
  alignment/A549_CTRL_WCE_rep2/A549_CTRL_WCE_rep2.sorted.dup.bam \
  --name peak_call/A549_CTRL_SMC1A_rep1/A549_CTRL_SMC1A_rep1 \
  >& peak_call/A549_CTRL_SMC1A_rep1/A549_CTRL_SMC1A_rep1.diag.macs.out
macs2_callpeak.A549_CTRL_SMC1A_rep1.5405b6d488e76cfc51aecda63958297d.mugqic.done
)
macs2_callpeak_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_7_JOB_ID: macs2_callpeak.A549_CTRL_MED1_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_CTRL_MED1_rep2
JOB_DEPENDENCIES=$picard_mark_duplicates_4_JOB_ID:$picard_mark_duplicates_7_JOB_ID:$picard_mark_duplicates_8_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_CTRL_MED1_rep2.49f9ec70193388d80bcabb811ccaae6f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_CTRL_MED1_rep2.49f9ec70193388d80bcabb811ccaae6f.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/A549_CTRL_MED1_rep2 && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2479938032.8 \
  --treatment \
  alignment/A549_CTRL_MED1_rep2/A549_CTRL_MED1_rep2.sorted.dup.bam \
  --control \
  alignment/A549_CTRL_WCE_rep1/A549_CTRL_WCE_rep1.sorted.dup.bam \
  alignment/A549_CTRL_WCE_rep2/A549_CTRL_WCE_rep2.sorted.dup.bam \
  --name peak_call/A549_CTRL_MED1_rep2/A549_CTRL_MED1_rep2 \
  >& peak_call/A549_CTRL_MED1_rep2/A549_CTRL_MED1_rep2.diag.macs.out
macs2_callpeak.A549_CTRL_MED1_rep2.49f9ec70193388d80bcabb811ccaae6f.mugqic.done
)
macs2_callpeak_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_8_JOB_ID: macs2_callpeak.A549_CTRL_SMC1A_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_CTRL_SMC1A_rep2
JOB_DEPENDENCIES=$picard_mark_duplicates_6_JOB_ID:$picard_mark_duplicates_7_JOB_ID:$picard_mark_duplicates_8_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_CTRL_SMC1A_rep2.47493f30b72a0f94f8b5b9a4c5068e56.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_CTRL_SMC1A_rep2.47493f30b72a0f94f8b5b9a4c5068e56.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/A549_CTRL_SMC1A_rep2 && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2479938032.8 \
  --treatment \
  alignment/A549_CTRL_SMC1A_rep2/A549_CTRL_SMC1A_rep2.sorted.dup.bam \
  --control \
  alignment/A549_CTRL_WCE_rep1/A549_CTRL_WCE_rep1.sorted.dup.bam \
  alignment/A549_CTRL_WCE_rep2/A549_CTRL_WCE_rep2.sorted.dup.bam \
  --name peak_call/A549_CTRL_SMC1A_rep2/A549_CTRL_SMC1A_rep2 \
  >& peak_call/A549_CTRL_SMC1A_rep2/A549_CTRL_SMC1A_rep2.diag.macs.out
macs2_callpeak.A549_CTRL_SMC1A_rep2.47493f30b72a0f94f8b5b9a4c5068e56.mugqic.done
)
macs2_callpeak_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_9_JOB_ID: macs2_callpeak_report
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_report
JOB_DEPENDENCIES=$macs2_callpeak_1_JOB_ID:$macs2_callpeak_2_JOB_ID:$macs2_callpeak_3_JOB_ID:$macs2_callpeak_4_JOB_ID:$macs2_callpeak_5_JOB_ID:$macs2_callpeak_6_JOB_ID:$macs2_callpeak_7_JOB_ID:$macs2_callpeak_8_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_report.620eb83c7693ca9689779238b126cfd3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_report.620eb83c7693ca9689779238b126cfd3.mugqic.done'
mkdir -p report && \
cp /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/ChipSeq.macs2_callpeak.md report/ && \
for contrast in A549_CTRL_BRD4 A549_CTRL_CDK9 A549_CTRL_MED1 A549_CTRL_SMC1A A549_CTRL_MED1_rep1 A549_CTRL_SMC1A_rep1 A549_CTRL_MED1_rep2 A549_CTRL_SMC1A_rep2
do
  cp -a --parents peak_call/$contrast/ report/ && \
  echo -e "* [Peak Calls File for Design $contrast](peak_call/$contrast/${contrast}_peaks.xls)" \
  >> report/ChipSeq.macs2_callpeak.md
done
macs2_callpeak_report.620eb83c7693ca9689779238b126cfd3.mugqic.done
)
macs2_callpeak_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: homer_annotate_peaks
#-------------------------------------------------------------------------------
STEP=homer_annotate_peaks
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_1_JOB_ID: homer_annotate_peaks.A549_CTRL_BRD4
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.A549_CTRL_BRD4
JOB_DEPENDENCIES=$macs2_callpeak_1_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.A549_CTRL_BRD4.75f688521827f0101ff99802a4359299.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.A549_CTRL_BRD4.75f688521827f0101ff99802a4359299.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/A549_CTRL_BRD4/A549_CTRL_BRD4 && \
annotatePeaks.pl \
  peak_call/A549_CTRL_BRD4/A549_CTRL_BRD4_peaks.narrowPeak \
  GRCh38 \
  -gsize GRCh38 \
  -cons -CpG \
  -go annotation/A549_CTRL_BRD4/A549_CTRL_BRD4 \
  -genomeOntology annotation/A549_CTRL_BRD4/A549_CTRL_BRD4 \
  > annotation/A549_CTRL_BRD4/A549_CTRL_BRD4.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/A549_CTRL_BRD4/A549_CTRL_BRD4.annotated.csv",
  "annotation/A549_CTRL_BRD4/A549_CTRL_BRD4",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.A549_CTRL_BRD4.75f688521827f0101ff99802a4359299.mugqic.done
)
homer_annotate_peaks_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_2_JOB_ID: homer_annotate_peaks.A549_CTRL_CDK9
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.A549_CTRL_CDK9
JOB_DEPENDENCIES=$macs2_callpeak_2_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.A549_CTRL_CDK9.e1d1c317355ac2a3e03dc93f52472334.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.A549_CTRL_CDK9.e1d1c317355ac2a3e03dc93f52472334.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/A549_CTRL_CDK9/A549_CTRL_CDK9 && \
annotatePeaks.pl \
  peak_call/A549_CTRL_CDK9/A549_CTRL_CDK9_peaks.narrowPeak \
  GRCh38 \
  -gsize GRCh38 \
  -cons -CpG \
  -go annotation/A549_CTRL_CDK9/A549_CTRL_CDK9 \
  -genomeOntology annotation/A549_CTRL_CDK9/A549_CTRL_CDK9 \
  > annotation/A549_CTRL_CDK9/A549_CTRL_CDK9.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/A549_CTRL_CDK9/A549_CTRL_CDK9.annotated.csv",
  "annotation/A549_CTRL_CDK9/A549_CTRL_CDK9",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.A549_CTRL_CDK9.e1d1c317355ac2a3e03dc93f52472334.mugqic.done
)
homer_annotate_peaks_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_3_JOB_ID: homer_annotate_peaks.A549_CTRL_MED1
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.A549_CTRL_MED1
JOB_DEPENDENCIES=$macs2_callpeak_3_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.A549_CTRL_MED1.115719e9dbaebfb8c918f2a81907c679.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.A549_CTRL_MED1.115719e9dbaebfb8c918f2a81907c679.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/A549_CTRL_MED1/A549_CTRL_MED1 && \
annotatePeaks.pl \
  peak_call/A549_CTRL_MED1/A549_CTRL_MED1_peaks.narrowPeak \
  GRCh38 \
  -gsize GRCh38 \
  -cons -CpG \
  -go annotation/A549_CTRL_MED1/A549_CTRL_MED1 \
  -genomeOntology annotation/A549_CTRL_MED1/A549_CTRL_MED1 \
  > annotation/A549_CTRL_MED1/A549_CTRL_MED1.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/A549_CTRL_MED1/A549_CTRL_MED1.annotated.csv",
  "annotation/A549_CTRL_MED1/A549_CTRL_MED1",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.A549_CTRL_MED1.115719e9dbaebfb8c918f2a81907c679.mugqic.done
)
homer_annotate_peaks_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_4_JOB_ID: homer_annotate_peaks.A549_CTRL_SMC1A
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.A549_CTRL_SMC1A
JOB_DEPENDENCIES=$macs2_callpeak_4_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.A549_CTRL_SMC1A.215a110f4bf8e140510bd655d9a0535f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.A549_CTRL_SMC1A.215a110f4bf8e140510bd655d9a0535f.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/A549_CTRL_SMC1A/A549_CTRL_SMC1A && \
annotatePeaks.pl \
  peak_call/A549_CTRL_SMC1A/A549_CTRL_SMC1A_peaks.narrowPeak \
  GRCh38 \
  -gsize GRCh38 \
  -cons -CpG \
  -go annotation/A549_CTRL_SMC1A/A549_CTRL_SMC1A \
  -genomeOntology annotation/A549_CTRL_SMC1A/A549_CTRL_SMC1A \
  > annotation/A549_CTRL_SMC1A/A549_CTRL_SMC1A.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/A549_CTRL_SMC1A/A549_CTRL_SMC1A.annotated.csv",
  "annotation/A549_CTRL_SMC1A/A549_CTRL_SMC1A",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.A549_CTRL_SMC1A.215a110f4bf8e140510bd655d9a0535f.mugqic.done
)
homer_annotate_peaks_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_5_JOB_ID: homer_annotate_peaks.A549_CTRL_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.A549_CTRL_MED1_rep1
JOB_DEPENDENCIES=$macs2_callpeak_5_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.A549_CTRL_MED1_rep1.71a5c3d923ac0d644415fb06cd751228.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.A549_CTRL_MED1_rep1.71a5c3d923ac0d644415fb06cd751228.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/A549_CTRL_MED1_rep1/A549_CTRL_MED1_rep1 && \
annotatePeaks.pl \
  peak_call/A549_CTRL_MED1_rep1/A549_CTRL_MED1_rep1_peaks.narrowPeak \
  GRCh38 \
  -gsize GRCh38 \
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
homer_annotate_peaks.A549_CTRL_MED1_rep1.71a5c3d923ac0d644415fb06cd751228.mugqic.done
)
homer_annotate_peaks_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_6_JOB_ID: homer_annotate_peaks.A549_CTRL_SMC1A_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.A549_CTRL_SMC1A_rep1
JOB_DEPENDENCIES=$macs2_callpeak_6_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.A549_CTRL_SMC1A_rep1.87dfee4af8ab0db8023d97ba81e13fb8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.A549_CTRL_SMC1A_rep1.87dfee4af8ab0db8023d97ba81e13fb8.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/A549_CTRL_SMC1A_rep1/A549_CTRL_SMC1A_rep1 && \
annotatePeaks.pl \
  peak_call/A549_CTRL_SMC1A_rep1/A549_CTRL_SMC1A_rep1_peaks.narrowPeak \
  GRCh38 \
  -gsize GRCh38 \
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
homer_annotate_peaks.A549_CTRL_SMC1A_rep1.87dfee4af8ab0db8023d97ba81e13fb8.mugqic.done
)
homer_annotate_peaks_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_7_JOB_ID: homer_annotate_peaks.A549_CTRL_MED1_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.A549_CTRL_MED1_rep2
JOB_DEPENDENCIES=$macs2_callpeak_7_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.A549_CTRL_MED1_rep2.d844582fd097725761e73b3913cb37e3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.A549_CTRL_MED1_rep2.d844582fd097725761e73b3913cb37e3.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/A549_CTRL_MED1_rep2/A549_CTRL_MED1_rep2 && \
annotatePeaks.pl \
  peak_call/A549_CTRL_MED1_rep2/A549_CTRL_MED1_rep2_peaks.narrowPeak \
  GRCh38 \
  -gsize GRCh38 \
  -cons -CpG \
  -go annotation/A549_CTRL_MED1_rep2/A549_CTRL_MED1_rep2 \
  -genomeOntology annotation/A549_CTRL_MED1_rep2/A549_CTRL_MED1_rep2 \
  > annotation/A549_CTRL_MED1_rep2/A549_CTRL_MED1_rep2.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/A549_CTRL_MED1_rep2/A549_CTRL_MED1_rep2.annotated.csv",
  "annotation/A549_CTRL_MED1_rep2/A549_CTRL_MED1_rep2",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.A549_CTRL_MED1_rep2.d844582fd097725761e73b3913cb37e3.mugqic.done
)
homer_annotate_peaks_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_8_JOB_ID: homer_annotate_peaks.A549_CTRL_SMC1A_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.A549_CTRL_SMC1A_rep2
JOB_DEPENDENCIES=$macs2_callpeak_8_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.A549_CTRL_SMC1A_rep2.649ed9fd6f222a45cadc705b45a41e52.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.A549_CTRL_SMC1A_rep2.649ed9fd6f222a45cadc705b45a41e52.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/A549_CTRL_SMC1A_rep2/A549_CTRL_SMC1A_rep2 && \
annotatePeaks.pl \
  peak_call/A549_CTRL_SMC1A_rep2/A549_CTRL_SMC1A_rep2_peaks.narrowPeak \
  GRCh38 \
  -gsize GRCh38 \
  -cons -CpG \
  -go annotation/A549_CTRL_SMC1A_rep2/A549_CTRL_SMC1A_rep2 \
  -genomeOntology annotation/A549_CTRL_SMC1A_rep2/A549_CTRL_SMC1A_rep2 \
  > annotation/A549_CTRL_SMC1A_rep2/A549_CTRL_SMC1A_rep2.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/A549_CTRL_SMC1A_rep2/A549_CTRL_SMC1A_rep2.annotated.csv",
  "annotation/A549_CTRL_SMC1A_rep2/A549_CTRL_SMC1A_rep2",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.A549_CTRL_SMC1A_rep2.649ed9fd6f222a45cadc705b45a41e52.mugqic.done
)
homer_annotate_peaks_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_9_JOB_ID: homer_annotate_peaks_report
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks_report
JOB_DEPENDENCIES=$homer_annotate_peaks_1_JOB_ID:$homer_annotate_peaks_2_JOB_ID:$homer_annotate_peaks_3_JOB_ID:$homer_annotate_peaks_4_JOB_ID:$homer_annotate_peaks_5_JOB_ID:$homer_annotate_peaks_6_JOB_ID:$homer_annotate_peaks_7_JOB_ID:$homer_annotate_peaks_8_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks_report.d14cab46c0db8ab092063448b51e57c2.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks_report.d14cab46c0db8ab092063448b51e57c2.mugqic.done'
mkdir -p report/annotation/ && \
cp /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/ChipSeq.homer_annotate_peaks.md report/ && \
for contrast in A549_CTRL_BRD4 A549_CTRL_CDK9 A549_CTRL_MED1 A549_CTRL_SMC1A A549_CTRL_MED1_rep1 A549_CTRL_SMC1A_rep1 A549_CTRL_MED1_rep2 A549_CTRL_SMC1A_rep2
do
  rsync -avP annotation/$contrast report/annotation/ && \
  echo -e "* [Gene Annotations for Design $contrast](annotation/$contrast/${contrast}.annotated.csv)
* [HOMER Gene Ontology Annotations for Design $contrast](annotation/$contrast/$contrast/geneOntology.html)
* [HOMER Genome Ontology Annotations for Design $contrast](annotation/$contrast/$contrast/GenomeOntology.html)" \
  >> report/ChipSeq.homer_annotate_peaks.md
done
homer_annotate_peaks_report.d14cab46c0db8ab092063448b51e57c2.mugqic.done
)
homer_annotate_peaks_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: homer_find_motifs_genome
#-------------------------------------------------------------------------------
STEP=homer_find_motifs_genome
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_1_JOB_ID: homer_find_motifs_genome.A549_CTRL_BRD4
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.A549_CTRL_BRD4
JOB_DEPENDENCIES=$macs2_callpeak_1_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.A549_CTRL_BRD4.54e758e6b8815cdfe7bcb9bbfa5f0469.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.A549_CTRL_BRD4.54e758e6b8815cdfe7bcb9bbfa5f0469.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/A549_CTRL_BRD4/A549_CTRL_BRD4 && \
findMotifsGenome.pl \
  peak_call/A549_CTRL_BRD4/A549_CTRL_BRD4_peaks.narrowPeak \
  GRCh38 \
  annotation/A549_CTRL_BRD4/A549_CTRL_BRD4 \
  -preparsedDir annotation/A549_CTRL_BRD4/A549_CTRL_BRD4/preparsed \
  -p 4
homer_find_motifs_genome.A549_CTRL_BRD4.54e758e6b8815cdfe7bcb9bbfa5f0469.mugqic.done
)
homer_find_motifs_genome_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_2_JOB_ID: homer_find_motifs_genome.A549_CTRL_CDK9
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.A549_CTRL_CDK9
JOB_DEPENDENCIES=$macs2_callpeak_2_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.A549_CTRL_CDK9.53b55dc53a644a6849e0e5eb76bbc000.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.A549_CTRL_CDK9.53b55dc53a644a6849e0e5eb76bbc000.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/A549_CTRL_CDK9/A549_CTRL_CDK9 && \
findMotifsGenome.pl \
  peak_call/A549_CTRL_CDK9/A549_CTRL_CDK9_peaks.narrowPeak \
  GRCh38 \
  annotation/A549_CTRL_CDK9/A549_CTRL_CDK9 \
  -preparsedDir annotation/A549_CTRL_CDK9/A549_CTRL_CDK9/preparsed \
  -p 4
homer_find_motifs_genome.A549_CTRL_CDK9.53b55dc53a644a6849e0e5eb76bbc000.mugqic.done
)
homer_find_motifs_genome_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_3_JOB_ID: homer_find_motifs_genome.A549_CTRL_MED1
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.A549_CTRL_MED1
JOB_DEPENDENCIES=$macs2_callpeak_3_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.A549_CTRL_MED1.54b0b50e2b5a6912f8c61166fc442409.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.A549_CTRL_MED1.54b0b50e2b5a6912f8c61166fc442409.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/A549_CTRL_MED1/A549_CTRL_MED1 && \
findMotifsGenome.pl \
  peak_call/A549_CTRL_MED1/A549_CTRL_MED1_peaks.narrowPeak \
  GRCh38 \
  annotation/A549_CTRL_MED1/A549_CTRL_MED1 \
  -preparsedDir annotation/A549_CTRL_MED1/A549_CTRL_MED1/preparsed \
  -p 4
homer_find_motifs_genome.A549_CTRL_MED1.54b0b50e2b5a6912f8c61166fc442409.mugqic.done
)
homer_find_motifs_genome_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_4_JOB_ID: homer_find_motifs_genome.A549_CTRL_SMC1A
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.A549_CTRL_SMC1A
JOB_DEPENDENCIES=$macs2_callpeak_4_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.A549_CTRL_SMC1A.c7455fae160965a155ee64107e6ca56a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.A549_CTRL_SMC1A.c7455fae160965a155ee64107e6ca56a.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/A549_CTRL_SMC1A/A549_CTRL_SMC1A && \
findMotifsGenome.pl \
  peak_call/A549_CTRL_SMC1A/A549_CTRL_SMC1A_peaks.narrowPeak \
  GRCh38 \
  annotation/A549_CTRL_SMC1A/A549_CTRL_SMC1A \
  -preparsedDir annotation/A549_CTRL_SMC1A/A549_CTRL_SMC1A/preparsed \
  -p 4
homer_find_motifs_genome.A549_CTRL_SMC1A.c7455fae160965a155ee64107e6ca56a.mugqic.done
)
homer_find_motifs_genome_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_5_JOB_ID: homer_find_motifs_genome.A549_CTRL_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.A549_CTRL_MED1_rep1
JOB_DEPENDENCIES=$macs2_callpeak_5_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.A549_CTRL_MED1_rep1.1a467ed565a0061a8ce66040494e38e7.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.A549_CTRL_MED1_rep1.1a467ed565a0061a8ce66040494e38e7.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/A549_CTRL_MED1_rep1/A549_CTRL_MED1_rep1 && \
findMotifsGenome.pl \
  peak_call/A549_CTRL_MED1_rep1/A549_CTRL_MED1_rep1_peaks.narrowPeak \
  GRCh38 \
  annotation/A549_CTRL_MED1_rep1/A549_CTRL_MED1_rep1 \
  -preparsedDir annotation/A549_CTRL_MED1_rep1/A549_CTRL_MED1_rep1/preparsed \
  -p 4
homer_find_motifs_genome.A549_CTRL_MED1_rep1.1a467ed565a0061a8ce66040494e38e7.mugqic.done
)
homer_find_motifs_genome_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_6_JOB_ID: homer_find_motifs_genome.A549_CTRL_SMC1A_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.A549_CTRL_SMC1A_rep1
JOB_DEPENDENCIES=$macs2_callpeak_6_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.A549_CTRL_SMC1A_rep1.0ca9a7dd922a9f0f230dce07a8c90aff.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.A549_CTRL_SMC1A_rep1.0ca9a7dd922a9f0f230dce07a8c90aff.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/A549_CTRL_SMC1A_rep1/A549_CTRL_SMC1A_rep1 && \
findMotifsGenome.pl \
  peak_call/A549_CTRL_SMC1A_rep1/A549_CTRL_SMC1A_rep1_peaks.narrowPeak \
  GRCh38 \
  annotation/A549_CTRL_SMC1A_rep1/A549_CTRL_SMC1A_rep1 \
  -preparsedDir annotation/A549_CTRL_SMC1A_rep1/A549_CTRL_SMC1A_rep1/preparsed \
  -p 4
homer_find_motifs_genome.A549_CTRL_SMC1A_rep1.0ca9a7dd922a9f0f230dce07a8c90aff.mugqic.done
)
homer_find_motifs_genome_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_7_JOB_ID: homer_find_motifs_genome.A549_CTRL_MED1_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.A549_CTRL_MED1_rep2
JOB_DEPENDENCIES=$macs2_callpeak_7_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.A549_CTRL_MED1_rep2.f61803ee3640c5e8a20c85c4c216daef.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.A549_CTRL_MED1_rep2.f61803ee3640c5e8a20c85c4c216daef.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/A549_CTRL_MED1_rep2/A549_CTRL_MED1_rep2 && \
findMotifsGenome.pl \
  peak_call/A549_CTRL_MED1_rep2/A549_CTRL_MED1_rep2_peaks.narrowPeak \
  GRCh38 \
  annotation/A549_CTRL_MED1_rep2/A549_CTRL_MED1_rep2 \
  -preparsedDir annotation/A549_CTRL_MED1_rep2/A549_CTRL_MED1_rep2/preparsed \
  -p 4
homer_find_motifs_genome.A549_CTRL_MED1_rep2.f61803ee3640c5e8a20c85c4c216daef.mugqic.done
)
homer_find_motifs_genome_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_8_JOB_ID: homer_find_motifs_genome.A549_CTRL_SMC1A_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.A549_CTRL_SMC1A_rep2
JOB_DEPENDENCIES=$macs2_callpeak_8_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.A549_CTRL_SMC1A_rep2.275897bcbb598113d8a46ca57a6c9150.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.A549_CTRL_SMC1A_rep2.275897bcbb598113d8a46ca57a6c9150.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/A549_CTRL_SMC1A_rep2/A549_CTRL_SMC1A_rep2 && \
findMotifsGenome.pl \
  peak_call/A549_CTRL_SMC1A_rep2/A549_CTRL_SMC1A_rep2_peaks.narrowPeak \
  GRCh38 \
  annotation/A549_CTRL_SMC1A_rep2/A549_CTRL_SMC1A_rep2 \
  -preparsedDir annotation/A549_CTRL_SMC1A_rep2/A549_CTRL_SMC1A_rep2/preparsed \
  -p 4
homer_find_motifs_genome.A549_CTRL_SMC1A_rep2.275897bcbb598113d8a46ca57a6c9150.mugqic.done
)
homer_find_motifs_genome_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_9_JOB_ID: homer_find_motifs_genome_report
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome_report
JOB_DEPENDENCIES=$homer_find_motifs_genome_1_JOB_ID:$homer_find_motifs_genome_2_JOB_ID:$homer_find_motifs_genome_3_JOB_ID:$homer_find_motifs_genome_4_JOB_ID:$homer_find_motifs_genome_5_JOB_ID:$homer_find_motifs_genome_6_JOB_ID:$homer_find_motifs_genome_7_JOB_ID:$homer_find_motifs_genome_8_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome_report.e77fdfbaf52f832e95b64b95184e1e3a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome_report.e77fdfbaf52f832e95b64b95184e1e3a.mugqic.done'
mkdir -p report/annotation/ && \
cp /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/ChipSeq.homer_find_motifs_genome.md report/ && \
for contrast in A549_CTRL_BRD4 A549_CTRL_CDK9 A549_CTRL_MED1 A549_CTRL_SMC1A A549_CTRL_MED1_rep1 A549_CTRL_SMC1A_rep1 A549_CTRL_MED1_rep2 A549_CTRL_SMC1A_rep2
do
  rsync -avP annotation/$contrast report/annotation/ && \
  echo -e "* [HOMER _De Novo_ Motif Results for Design $contrast](annotation/$contrast/$contrast/homerResults.html)
* [HOMER Known Motif Results for Design $contrast](annotation/$contrast/$contrast/knownResults.html)" \
  >> report/ChipSeq.homer_find_motifs_genome.md
done
homer_find_motifs_genome_report.e77fdfbaf52f832e95b64b95184e1e3a.mugqic.done
)
homer_find_motifs_genome_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: annotation_graphs
#-------------------------------------------------------------------------------
STEP=annotation_graphs
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: annotation_graphs_1_JOB_ID: annotation_graphs
#-------------------------------------------------------------------------------
JOB_NAME=annotation_graphs
JOB_DEPENDENCIES=$homer_annotate_peaks_1_JOB_ID:$homer_annotate_peaks_2_JOB_ID:$homer_annotate_peaks_3_JOB_ID:$homer_annotate_peaks_4_JOB_ID:$homer_annotate_peaks_5_JOB_ID:$homer_annotate_peaks_6_JOB_ID:$homer_annotate_peaks_7_JOB_ID:$homer_annotate_peaks_8_JOB_ID
JOB_DONE=job_output/annotation_graphs/annotation_graphs.7c90459ac0753def350af11f1472e1a0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'annotation_graphs.7c90459ac0753def350af11f1472e1a0.mugqic.done'
module load mugqic/mugqic_tools/2.1.5 mugqic/R_Bioconductor/3.2.3_3.2 mugqic/pandoc/1.15.2 && \
mkdir -p graphs && \
Rscript $R_TOOLS/chipSeqgenerateAnnotationGraphs.R \
  ../../raw/chip-seq/design_partial.txt \
  /gs/scratch/efournier/sb_cofactor_hr/A549/output/chip-pipeline-GRCh38 && \
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
for contrast in A549_CTRL_BRD4 A549_CTRL_CDK9 A549_CTRL_MED1 A549_CTRL_SMC1A A549_CTRL_MED1_rep1 A549_CTRL_SMC1A_rep1 A549_CTRL_MED1_rep2 A549_CTRL_SMC1A_rep2
do
  cp --parents graphs/${contrast}_Misc_Graphs.ps report/
  convert -rotate 90 graphs/${contrast}_Misc_Graphs.ps report/graphs/${contrast}_Misc_Graphs.png
  echo -e "----

![Annotation Statistics for Design $contrast ([download high-res image](graphs/${contrast}_Misc_Graphs.ps))](graphs/${contrast}_Misc_Graphs.png)
" \
  >> report/ChipSeq.annotation_graphs.md
done
annotation_graphs.7c90459ac0753def350af11f1472e1a0.mugqic.done
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
wget "http://mugqic.hpc.mcgill.ca/cgi-bin/pipeline.cgi?hostname=lg-1r17-n02&ip=10.241.129.12&pipeline=ChipSeq&steps=picard_sam_to_fastq,trimmomatic,merge_trimmomatic_stats,bwa_mem_picard_sort_sam,samtools_view_filter,picard_merge_sam_files,picard_mark_duplicates,metrics,homer_make_tag_directory,qc_metrics,homer_make_ucsc_file,macs2_callpeak,homer_annotate_peaks,homer_find_motifs_genome,annotation_graphs&samples=8" --quiet --output-document=/dev/null

