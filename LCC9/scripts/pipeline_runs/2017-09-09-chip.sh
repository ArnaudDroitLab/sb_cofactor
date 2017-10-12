#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

#-------------------------------------------------------------------------------
# ChipSeq PBSScheduler Job Submission Bash script
# Version: 2.2.0
# Created on: 2017-09-10T02:40:56
# Steps:
#   picard_sam_to_fastq: 12 jobs
#   trimmomatic: 12 jobs
#   merge_trimmomatic_stats: 1 job
#   bwa_mem_picard_sort_sam: 13 jobs
#   samtools_view_filter: 13 jobs
#   picard_merge_sam_files: 12 jobs
#   picard_mark_duplicates: 13 jobs
#   metrics: 2 jobs
#   homer_make_tag_directory: 12 jobs
#   qc_metrics: 1 job
#   homer_make_ucsc_file: 13 jobs
#   macs2_callpeak: 11 jobs
#   homer_annotate_peaks: 11 jobs
#   homer_find_motifs_genome: 9 jobs
#   annotation_graphs: 1 job
#   TOTAL: 136 jobs
#-------------------------------------------------------------------------------

OUTPUT_DIR=/gs/scratch/efournier/sb_cofactor_hr/LCC9/output/chip-pipeline-GRCh38
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
# JOB: picard_sam_to_fastq_1_JOB_ID: picard_sam_to_fastq.HI.1943.008.Index_14.LCC9_CTRL_MED1_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.1943.008.Index_14.LCC9_CTRL_MED1_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.1943.008.Index_14.LCC9_CTRL_MED1_ChIP8_GB_rep1.bam.13cd0117348203faccdc0c786a99683a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.1943.008.Index_14.LCC9_CTRL_MED1_ChIP8_GB_rep1.bam.13cd0117348203faccdc0c786a99683a.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/SamToFastq.jar \
  VALIDATION_STRINGENCY=LENIENT \
  INPUT=/gs/scratch/efournier/sb_cofactor_hr/LCC9/raw/chip-seq/HI.1943.008.Index_14.LCC9_CTRL_MED1_ChIP8_GB_rep1.bam \
  FASTQ=/gs/scratch/efournier/sb_cofactor_hr/LCC9/raw/chip-seq/HI.1943.008.Index_14.LCC9_CTRL_MED1_ChIP8_GB_rep1.single.fastq.gz
picard_sam_to_fastq.HI.1943.008.Index_14.LCC9_CTRL_MED1_ChIP8_GB_rep1.bam.13cd0117348203faccdc0c786a99683a.mugqic.done
)
picard_sam_to_fastq_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=3 | grep "[0-9]")
echo "$picard_sam_to_fastq_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_2_JOB_ID: picard_sam_to_fastq.HI.1943.008.Index_15.LCC9_E2_WCE_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.1943.008.Index_15.LCC9_E2_WCE_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.1943.008.Index_15.LCC9_E2_WCE_ChIP8_GB_rep1.bam.51ef5e58bb33e2b4d8c571c87314dfc1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.1943.008.Index_15.LCC9_E2_WCE_ChIP8_GB_rep1.bam.51ef5e58bb33e2b4d8c571c87314dfc1.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/SamToFastq.jar \
  VALIDATION_STRINGENCY=LENIENT \
  INPUT=/gs/scratch/efournier/sb_cofactor_hr/LCC9/raw/chip-seq/HI.1943.008.Index_15.LCC9_E2_WCE_ChIP8_GB_rep1.bam \
  FASTQ=/gs/scratch/efournier/sb_cofactor_hr/LCC9/raw/chip-seq/HI.1943.008.Index_15.LCC9_E2_WCE_ChIP8_GB_rep1.single.fastq.gz
picard_sam_to_fastq.HI.1943.008.Index_15.LCC9_E2_WCE_ChIP8_GB_rep1.bam.51ef5e58bb33e2b4d8c571c87314dfc1.mugqic.done
)
picard_sam_to_fastq_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=3 | grep "[0-9]")
echo "$picard_sam_to_fastq_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_3_JOB_ID: picard_sam_to_fastq.HI.1943.008.Index_16.LCC9_E2_MED1_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.1943.008.Index_16.LCC9_E2_MED1_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.1943.008.Index_16.LCC9_E2_MED1_ChIP8_GB_rep1.bam.f389aa1a667385ffd794691718c354f1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.1943.008.Index_16.LCC9_E2_MED1_ChIP8_GB_rep1.bam.f389aa1a667385ffd794691718c354f1.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/SamToFastq.jar \
  VALIDATION_STRINGENCY=LENIENT \
  INPUT=/gs/scratch/efournier/sb_cofactor_hr/LCC9/raw/chip-seq/HI.1943.008.Index_16.LCC9_E2_MED1_ChIP8_GB_rep1.bam \
  FASTQ=/gs/scratch/efournier/sb_cofactor_hr/LCC9/raw/chip-seq/HI.1943.008.Index_16.LCC9_E2_MED1_ChIP8_GB_rep1.single.fastq.gz
picard_sam_to_fastq.HI.1943.008.Index_16.LCC9_E2_MED1_ChIP8_GB_rep1.bam.f389aa1a667385ffd794691718c354f1.mugqic.done
)
picard_sam_to_fastq_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=3 | grep "[0-9]")
echo "$picard_sam_to_fastq_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_4_JOB_ID: picard_sam_to_fastq.HI.1943.008.Index_6.LCC9_CTRL_WCE_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.1943.008.Index_6.LCC9_CTRL_WCE_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.1943.008.Index_6.LCC9_CTRL_WCE_ChIP8_GB_rep1.bam.14de190662c4724c1d61c3165bf3d28c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.1943.008.Index_6.LCC9_CTRL_WCE_ChIP8_GB_rep1.bam.14de190662c4724c1d61c3165bf3d28c.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/SamToFastq.jar \
  VALIDATION_STRINGENCY=LENIENT \
  INPUT=/gs/scratch/efournier/sb_cofactor_hr/LCC9/raw/chip-seq/HI.1943.008.Index_6.LCC9_CTRL_WCE_ChIP8_GB_rep1.bam \
  FASTQ=/gs/scratch/efournier/sb_cofactor_hr/LCC9/raw/chip-seq/HI.1943.008.Index_6.LCC9_CTRL_WCE_ChIP8_GB_rep1.single.fastq.gz
picard_sam_to_fastq.HI.1943.008.Index_6.LCC9_CTRL_WCE_ChIP8_GB_rep1.bam.14de190662c4724c1d61c3165bf3d28c.mugqic.done
)
picard_sam_to_fastq_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=3 | grep "[0-9]")
echo "$picard_sam_to_fastq_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_5_JOB_ID: picard_sam_to_fastq.HI.1997.004.Index_12.LCC9_CTRL_SMC1_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.1997.004.Index_12.LCC9_CTRL_SMC1_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.1997.004.Index_12.LCC9_CTRL_SMC1_ChIP8_GB_rep1.bam.588745bd50e5d3279a8d159e21220236.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.1997.004.Index_12.LCC9_CTRL_SMC1_ChIP8_GB_rep1.bam.588745bd50e5d3279a8d159e21220236.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/SamToFastq.jar \
  VALIDATION_STRINGENCY=LENIENT \
  INPUT=/gs/scratch/efournier/sb_cofactor_hr/LCC9/raw/chip-seq/HI.1997.004.Index_12.LCC9_CTRL_SMC1_ChIP8_GB_rep1.bam \
  FASTQ=/gs/scratch/efournier/sb_cofactor_hr/LCC9/raw/chip-seq/HI.1997.004.Index_12.LCC9_CTRL_SMC1_ChIP8_GB_rep1.single.fastq.gz
picard_sam_to_fastq.HI.1997.004.Index_12.LCC9_CTRL_SMC1_ChIP8_GB_rep1.bam.588745bd50e5d3279a8d159e21220236.mugqic.done
)
picard_sam_to_fastq_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=3 | grep "[0-9]")
echo "$picard_sam_to_fastq_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_6_JOB_ID: picard_sam_to_fastq.HI.1997.004.Index_19.LCC9_E2_SMC1_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.1997.004.Index_19.LCC9_E2_SMC1_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.1997.004.Index_19.LCC9_E2_SMC1_ChIP8_GB_rep1.bam.425f149fb23430c84030592d002103e8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.1997.004.Index_19.LCC9_E2_SMC1_ChIP8_GB_rep1.bam.425f149fb23430c84030592d002103e8.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/SamToFastq.jar \
  VALIDATION_STRINGENCY=LENIENT \
  INPUT=/gs/scratch/efournier/sb_cofactor_hr/LCC9/raw/chip-seq/HI.1997.004.Index_19.LCC9_E2_SMC1_ChIP8_GB_rep1.bam \
  FASTQ=/gs/scratch/efournier/sb_cofactor_hr/LCC9/raw/chip-seq/HI.1997.004.Index_19.LCC9_E2_SMC1_ChIP8_GB_rep1.single.fastq.gz
picard_sam_to_fastq.HI.1997.004.Index_19.LCC9_E2_SMC1_ChIP8_GB_rep1.bam.425f149fb23430c84030592d002103e8.mugqic.done
)
picard_sam_to_fastq_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=3 | grep "[0-9]")
echo "$picard_sam_to_fastq_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_7_JOB_ID: picard_sam_to_fastq.HI.1997.004.Index_8.LCC9_E2_NIPBL_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.1997.004.Index_8.LCC9_E2_NIPBL_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.1997.004.Index_8.LCC9_E2_NIPBL_ChIP8_GB_rep1.bam.51d554655b10e54249763679d9d8037d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.1997.004.Index_8.LCC9_E2_NIPBL_ChIP8_GB_rep1.bam.51d554655b10e54249763679d9d8037d.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/SamToFastq.jar \
  VALIDATION_STRINGENCY=LENIENT \
  INPUT=/gs/scratch/efournier/sb_cofactor_hr/LCC9/raw/chip-seq/HI.1997.004.Index_8.LCC9_E2_NIPBL_ChIP8_GB_rep1.bam \
  FASTQ=/gs/scratch/efournier/sb_cofactor_hr/LCC9/raw/chip-seq/HI.1997.004.Index_8.LCC9_E2_NIPBL_ChIP8_GB_rep1.single.fastq.gz
picard_sam_to_fastq.HI.1997.004.Index_8.LCC9_E2_NIPBL_ChIP8_GB_rep1.bam.51d554655b10e54249763679d9d8037d.mugqic.done
)
picard_sam_to_fastq_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=3 | grep "[0-9]")
echo "$picard_sam_to_fastq_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_8_JOB_ID: picard_sam_to_fastq.HI.1997.004.Index_9.LCC9_CTRL_NIPBL_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.1997.004.Index_9.LCC9_CTRL_NIPBL_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.1997.004.Index_9.LCC9_CTRL_NIPBL_ChIP8_GB_rep1.bam.3ce1c3f17b922e8cfdb2c0b06db26aef.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.1997.004.Index_9.LCC9_CTRL_NIPBL_ChIP8_GB_rep1.bam.3ce1c3f17b922e8cfdb2c0b06db26aef.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/SamToFastq.jar \
  VALIDATION_STRINGENCY=LENIENT \
  INPUT=/gs/scratch/efournier/sb_cofactor_hr/LCC9/raw/chip-seq/HI.1997.004.Index_9.LCC9_CTRL_NIPBL_ChIP8_GB_rep1.bam \
  FASTQ=/gs/scratch/efournier/sb_cofactor_hr/LCC9/raw/chip-seq/HI.1997.004.Index_9.LCC9_CTRL_NIPBL_ChIP8_GB_rep1.single.fastq.gz
picard_sam_to_fastq.HI.1997.004.Index_9.LCC9_CTRL_NIPBL_ChIP8_GB_rep1.bam.3ce1c3f17b922e8cfdb2c0b06db26aef.mugqic.done
)
picard_sam_to_fastq_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=3 | grep "[0-9]")
echo "$picard_sam_to_fastq_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_9_JOB_ID: picard_sam_to_fastq.HI.1997.005.Index_20.LCC9_CTRL_POL2_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.1997.005.Index_20.LCC9_CTRL_POL2_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.1997.005.Index_20.LCC9_CTRL_POL2_ChIP8_GB_rep1.bam.616635861d3653b18b17b7b89a5232af.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.1997.005.Index_20.LCC9_CTRL_POL2_ChIP8_GB_rep1.bam.616635861d3653b18b17b7b89a5232af.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/SamToFastq.jar \
  VALIDATION_STRINGENCY=LENIENT \
  INPUT=/gs/scratch/efournier/sb_cofactor_hr/LCC9/raw/chip-seq/HI.1997.005.Index_20.LCC9_CTRL_POL2_ChIP8_GB_rep1.bam \
  FASTQ=/gs/scratch/efournier/sb_cofactor_hr/LCC9/raw/chip-seq/HI.1997.005.Index_20.LCC9_CTRL_POL2_ChIP8_GB_rep1.single.fastq.gz
picard_sam_to_fastq.HI.1997.005.Index_20.LCC9_CTRL_POL2_ChIP8_GB_rep1.bam.616635861d3653b18b17b7b89a5232af.mugqic.done
)
picard_sam_to_fastq_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=3 | grep "[0-9]")
echo "$picard_sam_to_fastq_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_10_JOB_ID: picard_sam_to_fastq.HI.1997.005.Index_22.LCC9_E2_POL2_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.1997.005.Index_22.LCC9_E2_POL2_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.1997.005.Index_22.LCC9_E2_POL2_ChIP8_GB_rep1.bam.033b20d7a3a9409144ce77334fc745d0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.1997.005.Index_22.LCC9_E2_POL2_ChIP8_GB_rep1.bam.033b20d7a3a9409144ce77334fc745d0.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/SamToFastq.jar \
  VALIDATION_STRINGENCY=LENIENT \
  INPUT=/gs/scratch/efournier/sb_cofactor_hr/LCC9/raw/chip-seq/HI.1997.005.Index_22.LCC9_E2_POL2_ChIP8_GB_rep1.bam \
  FASTQ=/gs/scratch/efournier/sb_cofactor_hr/LCC9/raw/chip-seq/HI.1997.005.Index_22.LCC9_E2_POL2_ChIP8_GB_rep1.single.fastq.gz
picard_sam_to_fastq.HI.1997.005.Index_22.LCC9_E2_POL2_ChIP8_GB_rep1.bam.033b20d7a3a9409144ce77334fc745d0.mugqic.done
)
picard_sam_to_fastq_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=3 | grep "[0-9]")
echo "$picard_sam_to_fastq_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_11_JOB_ID: picard_sam_to_fastq.HI.1997.005.Index_23.LCC9_CTRL_ERA_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.1997.005.Index_23.LCC9_CTRL_ERA_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.1997.005.Index_23.LCC9_CTRL_ERA_ChIP8_GB_rep1.bam.1743d69e11f478b36360b07934e7bc67.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.1997.005.Index_23.LCC9_CTRL_ERA_ChIP8_GB_rep1.bam.1743d69e11f478b36360b07934e7bc67.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/SamToFastq.jar \
  VALIDATION_STRINGENCY=LENIENT \
  INPUT=/gs/scratch/efournier/sb_cofactor_hr/LCC9/raw/chip-seq/HI.1997.005.Index_23.LCC9_CTRL_ERA_ChIP8_GB_rep1.bam \
  FASTQ=/gs/scratch/efournier/sb_cofactor_hr/LCC9/raw/chip-seq/HI.1997.005.Index_23.LCC9_CTRL_ERA_ChIP8_GB_rep1.single.fastq.gz
picard_sam_to_fastq.HI.1997.005.Index_23.LCC9_CTRL_ERA_ChIP8_GB_rep1.bam.1743d69e11f478b36360b07934e7bc67.mugqic.done
)
picard_sam_to_fastq_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=3 | grep "[0-9]")
echo "$picard_sam_to_fastq_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_12_JOB_ID: picard_sam_to_fastq.HI.1997.005.Index_27.LCC9_E2_ERA_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.HI.1997.005.Index_27.LCC9_E2_ERA_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.HI.1997.005.Index_27.LCC9_E2_ERA_ChIP8_GB_rep1.bam.6d452076fd188c57568170ba30e474e4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.HI.1997.005.Index_27.LCC9_E2_ERA_ChIP8_GB_rep1.bam.6d452076fd188c57568170ba30e474e4.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/SamToFastq.jar \
  VALIDATION_STRINGENCY=LENIENT \
  INPUT=/gs/scratch/efournier/sb_cofactor_hr/LCC9/raw/chip-seq/HI.1997.005.Index_27.LCC9_E2_ERA_ChIP8_GB_rep1.bam \
  FASTQ=/gs/scratch/efournier/sb_cofactor_hr/LCC9/raw/chip-seq/HI.1997.005.Index_27.LCC9_E2_ERA_ChIP8_GB_rep1.single.fastq.gz
picard_sam_to_fastq.HI.1997.005.Index_27.LCC9_E2_ERA_ChIP8_GB_rep1.bam.6d452076fd188c57568170ba30e474e4.mugqic.done
)
picard_sam_to_fastq_12_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=3 | grep "[0-9]")
echo "$picard_sam_to_fastq_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: trimmomatic
#-------------------------------------------------------------------------------
STEP=trimmomatic
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: trimmomatic_1_JOB_ID: trimmomatic.HI.1943.008.Index_14.LCC9_CTRL_MED1_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.1943.008.Index_14.LCC9_CTRL_MED1_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_1_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.1943.008.Index_14.LCC9_CTRL_MED1_ChIP8_GB_rep1.bam.6a3c29023a1d54a02bd9d461e9dc225c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.1943.008.Index_14.LCC9_CTRL_MED1_ChIP8_GB_rep1.bam.6a3c29023a1d54a02bd9d461e9dc225c.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/LCC9_CTRL_MED1 && \
`cat > trim/LCC9_CTRL_MED1/HI.1943.008.Index_14.LCC9_CTRL_MED1_ChIP8_GB_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/sb_cofactor_hr/LCC9/raw/chip-seq/HI.1943.008.Index_14.LCC9_CTRL_MED1_ChIP8_GB_rep1.single.fastq.gz \
  trim/LCC9_CTRL_MED1/HI.1943.008.Index_14.LCC9_CTRL_MED1_ChIP8_GB_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/LCC9_CTRL_MED1/HI.1943.008.Index_14.LCC9_CTRL_MED1_ChIP8_GB_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/LCC9_CTRL_MED1/HI.1943.008.Index_14.LCC9_CTRL_MED1_ChIP8_GB_rep1.bam.trim.log
trimmomatic.HI.1943.008.Index_14.LCC9_CTRL_MED1_ChIP8_GB_rep1.bam.6a3c29023a1d54a02bd9d461e9dc225c.mugqic.done
)
trimmomatic_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_2_JOB_ID: trimmomatic.HI.1943.008.Index_15.LCC9_E2_WCE_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.1943.008.Index_15.LCC9_E2_WCE_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_2_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.1943.008.Index_15.LCC9_E2_WCE_ChIP8_GB_rep1.bam.0da3a9bda17f5ad59eb6d0940c99c651.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.1943.008.Index_15.LCC9_E2_WCE_ChIP8_GB_rep1.bam.0da3a9bda17f5ad59eb6d0940c99c651.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/LCC9_E2_WCE && \
`cat > trim/LCC9_E2_WCE/HI.1943.008.Index_15.LCC9_E2_WCE_ChIP8_GB_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/sb_cofactor_hr/LCC9/raw/chip-seq/HI.1943.008.Index_15.LCC9_E2_WCE_ChIP8_GB_rep1.single.fastq.gz \
  trim/LCC9_E2_WCE/HI.1943.008.Index_15.LCC9_E2_WCE_ChIP8_GB_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/LCC9_E2_WCE/HI.1943.008.Index_15.LCC9_E2_WCE_ChIP8_GB_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/LCC9_E2_WCE/HI.1943.008.Index_15.LCC9_E2_WCE_ChIP8_GB_rep1.bam.trim.log
trimmomatic.HI.1943.008.Index_15.LCC9_E2_WCE_ChIP8_GB_rep1.bam.0da3a9bda17f5ad59eb6d0940c99c651.mugqic.done
)
trimmomatic_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_3_JOB_ID: trimmomatic.HI.1943.008.Index_16.LCC9_E2_MED1_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.1943.008.Index_16.LCC9_E2_MED1_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_3_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.1943.008.Index_16.LCC9_E2_MED1_ChIP8_GB_rep1.bam.2f46cb784efee4e003786782d5504b7f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.1943.008.Index_16.LCC9_E2_MED1_ChIP8_GB_rep1.bam.2f46cb784efee4e003786782d5504b7f.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/LCC9_E2_MED1 && \
`cat > trim/LCC9_E2_MED1/HI.1943.008.Index_16.LCC9_E2_MED1_ChIP8_GB_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/sb_cofactor_hr/LCC9/raw/chip-seq/HI.1943.008.Index_16.LCC9_E2_MED1_ChIP8_GB_rep1.single.fastq.gz \
  trim/LCC9_E2_MED1/HI.1943.008.Index_16.LCC9_E2_MED1_ChIP8_GB_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/LCC9_E2_MED1/HI.1943.008.Index_16.LCC9_E2_MED1_ChIP8_GB_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/LCC9_E2_MED1/HI.1943.008.Index_16.LCC9_E2_MED1_ChIP8_GB_rep1.bam.trim.log
trimmomatic.HI.1943.008.Index_16.LCC9_E2_MED1_ChIP8_GB_rep1.bam.2f46cb784efee4e003786782d5504b7f.mugqic.done
)
trimmomatic_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_4_JOB_ID: trimmomatic.HI.1943.008.Index_6.LCC9_CTRL_WCE_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.1943.008.Index_6.LCC9_CTRL_WCE_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_4_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.1943.008.Index_6.LCC9_CTRL_WCE_ChIP8_GB_rep1.bam.f4ff259b4b9c26d2b04e8881e1691b15.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.1943.008.Index_6.LCC9_CTRL_WCE_ChIP8_GB_rep1.bam.f4ff259b4b9c26d2b04e8881e1691b15.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/LCC9_CTRL_WCE && \
`cat > trim/LCC9_CTRL_WCE/HI.1943.008.Index_6.LCC9_CTRL_WCE_ChIP8_GB_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/sb_cofactor_hr/LCC9/raw/chip-seq/HI.1943.008.Index_6.LCC9_CTRL_WCE_ChIP8_GB_rep1.single.fastq.gz \
  trim/LCC9_CTRL_WCE/HI.1943.008.Index_6.LCC9_CTRL_WCE_ChIP8_GB_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/LCC9_CTRL_WCE/HI.1943.008.Index_6.LCC9_CTRL_WCE_ChIP8_GB_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/LCC9_CTRL_WCE/HI.1943.008.Index_6.LCC9_CTRL_WCE_ChIP8_GB_rep1.bam.trim.log
trimmomatic.HI.1943.008.Index_6.LCC9_CTRL_WCE_ChIP8_GB_rep1.bam.f4ff259b4b9c26d2b04e8881e1691b15.mugqic.done
)
trimmomatic_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_5_JOB_ID: trimmomatic.HI.1997.004.Index_12.LCC9_CTRL_SMC1_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.1997.004.Index_12.LCC9_CTRL_SMC1_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_5_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.1997.004.Index_12.LCC9_CTRL_SMC1_ChIP8_GB_rep1.bam.b6dafbf3c74dd0ea45c1d81043bb9fe0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.1997.004.Index_12.LCC9_CTRL_SMC1_ChIP8_GB_rep1.bam.b6dafbf3c74dd0ea45c1d81043bb9fe0.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/LCC9_CTRL_SMC1 && \
`cat > trim/LCC9_CTRL_SMC1/HI.1997.004.Index_12.LCC9_CTRL_SMC1_ChIP8_GB_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/sb_cofactor_hr/LCC9/raw/chip-seq/HI.1997.004.Index_12.LCC9_CTRL_SMC1_ChIP8_GB_rep1.single.fastq.gz \
  trim/LCC9_CTRL_SMC1/HI.1997.004.Index_12.LCC9_CTRL_SMC1_ChIP8_GB_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/LCC9_CTRL_SMC1/HI.1997.004.Index_12.LCC9_CTRL_SMC1_ChIP8_GB_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/LCC9_CTRL_SMC1/HI.1997.004.Index_12.LCC9_CTRL_SMC1_ChIP8_GB_rep1.bam.trim.log
trimmomatic.HI.1997.004.Index_12.LCC9_CTRL_SMC1_ChIP8_GB_rep1.bam.b6dafbf3c74dd0ea45c1d81043bb9fe0.mugqic.done
)
trimmomatic_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_6_JOB_ID: trimmomatic.HI.1997.004.Index_19.LCC9_E2_SMC1_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.1997.004.Index_19.LCC9_E2_SMC1_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_6_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.1997.004.Index_19.LCC9_E2_SMC1_ChIP8_GB_rep1.bam.f66120e3eaec1f3825acf9b694301458.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.1997.004.Index_19.LCC9_E2_SMC1_ChIP8_GB_rep1.bam.f66120e3eaec1f3825acf9b694301458.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/LCC9_E2_SMC1 && \
`cat > trim/LCC9_E2_SMC1/HI.1997.004.Index_19.LCC9_E2_SMC1_ChIP8_GB_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/sb_cofactor_hr/LCC9/raw/chip-seq/HI.1997.004.Index_19.LCC9_E2_SMC1_ChIP8_GB_rep1.single.fastq.gz \
  trim/LCC9_E2_SMC1/HI.1997.004.Index_19.LCC9_E2_SMC1_ChIP8_GB_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/LCC9_E2_SMC1/HI.1997.004.Index_19.LCC9_E2_SMC1_ChIP8_GB_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/LCC9_E2_SMC1/HI.1997.004.Index_19.LCC9_E2_SMC1_ChIP8_GB_rep1.bam.trim.log
trimmomatic.HI.1997.004.Index_19.LCC9_E2_SMC1_ChIP8_GB_rep1.bam.f66120e3eaec1f3825acf9b694301458.mugqic.done
)
trimmomatic_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_7_JOB_ID: trimmomatic.HI.1997.004.Index_8.LCC9_E2_NIPBL_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.1997.004.Index_8.LCC9_E2_NIPBL_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_7_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.1997.004.Index_8.LCC9_E2_NIPBL_ChIP8_GB_rep1.bam.a926b8da92c518d8cddeb9017ace1cf5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.1997.004.Index_8.LCC9_E2_NIPBL_ChIP8_GB_rep1.bam.a926b8da92c518d8cddeb9017ace1cf5.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/LCC9_E2_NIPBL && \
`cat > trim/LCC9_E2_NIPBL/HI.1997.004.Index_8.LCC9_E2_NIPBL_ChIP8_GB_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/sb_cofactor_hr/LCC9/raw/chip-seq/HI.1997.004.Index_8.LCC9_E2_NIPBL_ChIP8_GB_rep1.single.fastq.gz \
  trim/LCC9_E2_NIPBL/HI.1997.004.Index_8.LCC9_E2_NIPBL_ChIP8_GB_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/LCC9_E2_NIPBL/HI.1997.004.Index_8.LCC9_E2_NIPBL_ChIP8_GB_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/LCC9_E2_NIPBL/HI.1997.004.Index_8.LCC9_E2_NIPBL_ChIP8_GB_rep1.bam.trim.log
trimmomatic.HI.1997.004.Index_8.LCC9_E2_NIPBL_ChIP8_GB_rep1.bam.a926b8da92c518d8cddeb9017ace1cf5.mugqic.done
)
trimmomatic_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_8_JOB_ID: trimmomatic.HI.1997.004.Index_9.LCC9_CTRL_NIPBL_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.1997.004.Index_9.LCC9_CTRL_NIPBL_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_8_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.1997.004.Index_9.LCC9_CTRL_NIPBL_ChIP8_GB_rep1.bam.e5fe85f1921afddcb5e57b232a201f15.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.1997.004.Index_9.LCC9_CTRL_NIPBL_ChIP8_GB_rep1.bam.e5fe85f1921afddcb5e57b232a201f15.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/LCC9_CTRL_NIPBL && \
`cat > trim/LCC9_CTRL_NIPBL/HI.1997.004.Index_9.LCC9_CTRL_NIPBL_ChIP8_GB_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/sb_cofactor_hr/LCC9/raw/chip-seq/HI.1997.004.Index_9.LCC9_CTRL_NIPBL_ChIP8_GB_rep1.single.fastq.gz \
  trim/LCC9_CTRL_NIPBL/HI.1997.004.Index_9.LCC9_CTRL_NIPBL_ChIP8_GB_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/LCC9_CTRL_NIPBL/HI.1997.004.Index_9.LCC9_CTRL_NIPBL_ChIP8_GB_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/LCC9_CTRL_NIPBL/HI.1997.004.Index_9.LCC9_CTRL_NIPBL_ChIP8_GB_rep1.bam.trim.log
trimmomatic.HI.1997.004.Index_9.LCC9_CTRL_NIPBL_ChIP8_GB_rep1.bam.e5fe85f1921afddcb5e57b232a201f15.mugqic.done
)
trimmomatic_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_9_JOB_ID: trimmomatic.HI.1997.005.Index_20.LCC9_CTRL_POL2_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.1997.005.Index_20.LCC9_CTRL_POL2_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_9_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.1997.005.Index_20.LCC9_CTRL_POL2_ChIP8_GB_rep1.bam.bd95910b563b879b1596f3c0c602f024.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.1997.005.Index_20.LCC9_CTRL_POL2_ChIP8_GB_rep1.bam.bd95910b563b879b1596f3c0c602f024.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/LCC9_CTRL_POL2 && \
`cat > trim/LCC9_CTRL_POL2/HI.1997.005.Index_20.LCC9_CTRL_POL2_ChIP8_GB_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/sb_cofactor_hr/LCC9/raw/chip-seq/HI.1997.005.Index_20.LCC9_CTRL_POL2_ChIP8_GB_rep1.single.fastq.gz \
  trim/LCC9_CTRL_POL2/HI.1997.005.Index_20.LCC9_CTRL_POL2_ChIP8_GB_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/LCC9_CTRL_POL2/HI.1997.005.Index_20.LCC9_CTRL_POL2_ChIP8_GB_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/LCC9_CTRL_POL2/HI.1997.005.Index_20.LCC9_CTRL_POL2_ChIP8_GB_rep1.bam.trim.log
trimmomatic.HI.1997.005.Index_20.LCC9_CTRL_POL2_ChIP8_GB_rep1.bam.bd95910b563b879b1596f3c0c602f024.mugqic.done
)
trimmomatic_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_10_JOB_ID: trimmomatic.HI.1997.005.Index_22.LCC9_E2_POL2_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.1997.005.Index_22.LCC9_E2_POL2_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_10_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.1997.005.Index_22.LCC9_E2_POL2_ChIP8_GB_rep1.bam.f364434bd675c4acb417f7f7ba5ff705.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.1997.005.Index_22.LCC9_E2_POL2_ChIP8_GB_rep1.bam.f364434bd675c4acb417f7f7ba5ff705.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/LCC9_E2_POL2 && \
`cat > trim/LCC9_E2_POL2/HI.1997.005.Index_22.LCC9_E2_POL2_ChIP8_GB_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/sb_cofactor_hr/LCC9/raw/chip-seq/HI.1997.005.Index_22.LCC9_E2_POL2_ChIP8_GB_rep1.single.fastq.gz \
  trim/LCC9_E2_POL2/HI.1997.005.Index_22.LCC9_E2_POL2_ChIP8_GB_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/LCC9_E2_POL2/HI.1997.005.Index_22.LCC9_E2_POL2_ChIP8_GB_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/LCC9_E2_POL2/HI.1997.005.Index_22.LCC9_E2_POL2_ChIP8_GB_rep1.bam.trim.log
trimmomatic.HI.1997.005.Index_22.LCC9_E2_POL2_ChIP8_GB_rep1.bam.f364434bd675c4acb417f7f7ba5ff705.mugqic.done
)
trimmomatic_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_11_JOB_ID: trimmomatic.HI.1997.005.Index_23.LCC9_CTRL_ERA_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.1997.005.Index_23.LCC9_CTRL_ERA_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_11_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.1997.005.Index_23.LCC9_CTRL_ERA_ChIP8_GB_rep1.bam.ee63b1f5c0f6f34fc97ca34926d344ad.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.1997.005.Index_23.LCC9_CTRL_ERA_ChIP8_GB_rep1.bam.ee63b1f5c0f6f34fc97ca34926d344ad.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/LCC9_CTRL_ERA && \
`cat > trim/LCC9_CTRL_ERA/HI.1997.005.Index_23.LCC9_CTRL_ERA_ChIP8_GB_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/sb_cofactor_hr/LCC9/raw/chip-seq/HI.1997.005.Index_23.LCC9_CTRL_ERA_ChIP8_GB_rep1.single.fastq.gz \
  trim/LCC9_CTRL_ERA/HI.1997.005.Index_23.LCC9_CTRL_ERA_ChIP8_GB_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/LCC9_CTRL_ERA/HI.1997.005.Index_23.LCC9_CTRL_ERA_ChIP8_GB_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/LCC9_CTRL_ERA/HI.1997.005.Index_23.LCC9_CTRL_ERA_ChIP8_GB_rep1.bam.trim.log
trimmomatic.HI.1997.005.Index_23.LCC9_CTRL_ERA_ChIP8_GB_rep1.bam.ee63b1f5c0f6f34fc97ca34926d344ad.mugqic.done
)
trimmomatic_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_12_JOB_ID: trimmomatic.HI.1997.005.Index_27.LCC9_E2_ERA_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.1997.005.Index_27.LCC9_E2_ERA_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$picard_sam_to_fastq_12_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.1997.005.Index_27.LCC9_E2_ERA_ChIP8_GB_rep1.bam.21ef91ee372cd971f97f306f6e713553.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.1997.005.Index_27.LCC9_E2_ERA_ChIP8_GB_rep1.bam.21ef91ee372cd971f97f306f6e713553.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/LCC9_E2_ERA && \
`cat > trim/LCC9_E2_ERA/HI.1997.005.Index_27.LCC9_E2_ERA_ChIP8_GB_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/sb_cofactor_hr/LCC9/raw/chip-seq/HI.1997.005.Index_27.LCC9_E2_ERA_ChIP8_GB_rep1.single.fastq.gz \
  trim/LCC9_E2_ERA/HI.1997.005.Index_27.LCC9_E2_ERA_ChIP8_GB_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/LCC9_E2_ERA/HI.1997.005.Index_27.LCC9_E2_ERA_ChIP8_GB_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/LCC9_E2_ERA/HI.1997.005.Index_27.LCC9_E2_ERA_ChIP8_GB_rep1.bam.trim.log
trimmomatic.HI.1997.005.Index_27.LCC9_E2_ERA_ChIP8_GB_rep1.bam.21ef91ee372cd971f97f306f6e713553.mugqic.done
)
trimmomatic_12_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: merge_trimmomatic_stats
#-------------------------------------------------------------------------------
STEP=merge_trimmomatic_stats
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: merge_trimmomatic_stats_1_JOB_ID: merge_trimmomatic_stats
#-------------------------------------------------------------------------------
JOB_NAME=merge_trimmomatic_stats
JOB_DEPENDENCIES=$trimmomatic_1_JOB_ID:$trimmomatic_2_JOB_ID:$trimmomatic_3_JOB_ID:$trimmomatic_4_JOB_ID:$trimmomatic_5_JOB_ID:$trimmomatic_6_JOB_ID:$trimmomatic_7_JOB_ID:$trimmomatic_8_JOB_ID:$trimmomatic_9_JOB_ID:$trimmomatic_10_JOB_ID:$trimmomatic_11_JOB_ID:$trimmomatic_12_JOB_ID
JOB_DONE=job_output/merge_trimmomatic_stats/merge_trimmomatic_stats.3693ad33276f6e2c997c732de3a9789a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'merge_trimmomatic_stats.3693ad33276f6e2c997c732de3a9789a.mugqic.done'
module load mugqic/pandoc/1.15.2 && \
mkdir -p metrics && \
echo 'Sample	Readset	Raw Single Reads #	Surviving Single Reads #	Surviving Single Reads %' > metrics/trimReadsetTable.tsv && \
grep ^Input trim/LCC9_CTRL_MED1/HI.1943.008.Index_14.LCC9_CTRL_MED1_ChIP8_GB_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/LCC9_CTRL_MED1	HI.1943.008.Index_14.LCC9_CTRL_MED1_ChIP8_GB_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/LCC9_E2_WCE/HI.1943.008.Index_15.LCC9_E2_WCE_ChIP8_GB_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/LCC9_E2_WCE	HI.1943.008.Index_15.LCC9_E2_WCE_ChIP8_GB_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/LCC9_E2_MED1/HI.1943.008.Index_16.LCC9_E2_MED1_ChIP8_GB_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/LCC9_E2_MED1	HI.1943.008.Index_16.LCC9_E2_MED1_ChIP8_GB_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/LCC9_CTRL_WCE/HI.1943.008.Index_6.LCC9_CTRL_WCE_ChIP8_GB_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/LCC9_CTRL_WCE	HI.1943.008.Index_6.LCC9_CTRL_WCE_ChIP8_GB_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/LCC9_CTRL_SMC1/HI.1997.004.Index_12.LCC9_CTRL_SMC1_ChIP8_GB_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/LCC9_CTRL_SMC1	HI.1997.004.Index_12.LCC9_CTRL_SMC1_ChIP8_GB_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/LCC9_E2_SMC1/HI.1997.004.Index_19.LCC9_E2_SMC1_ChIP8_GB_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/LCC9_E2_SMC1	HI.1997.004.Index_19.LCC9_E2_SMC1_ChIP8_GB_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/LCC9_E2_NIPBL/HI.1997.004.Index_8.LCC9_E2_NIPBL_ChIP8_GB_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/LCC9_E2_NIPBL	HI.1997.004.Index_8.LCC9_E2_NIPBL_ChIP8_GB_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/LCC9_CTRL_NIPBL/HI.1997.004.Index_9.LCC9_CTRL_NIPBL_ChIP8_GB_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/LCC9_CTRL_NIPBL	HI.1997.004.Index_9.LCC9_CTRL_NIPBL_ChIP8_GB_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/LCC9_CTRL_POL2/HI.1997.005.Index_20.LCC9_CTRL_POL2_ChIP8_GB_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/LCC9_CTRL_POL2	HI.1997.005.Index_20.LCC9_CTRL_POL2_ChIP8_GB_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/LCC9_E2_POL2/HI.1997.005.Index_22.LCC9_E2_POL2_ChIP8_GB_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/LCC9_E2_POL2	HI.1997.005.Index_22.LCC9_E2_POL2_ChIP8_GB_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/LCC9_CTRL_ERA/HI.1997.005.Index_23.LCC9_CTRL_ERA_ChIP8_GB_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/LCC9_CTRL_ERA	HI.1997.005.Index_23.LCC9_CTRL_ERA_ChIP8_GB_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/LCC9_E2_ERA/HI.1997.005.Index_27.LCC9_E2_ERA_ChIP8_GB_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/LCC9_E2_ERA	HI.1997.005.Index_27.LCC9_E2_ERA_ChIP8_GB_rep1.bam	\1	\2/' | \
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
merge_trimmomatic_stats.3693ad33276f6e2c997c732de3a9789a.mugqic.done
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
# JOB: bwa_mem_picard_sort_sam_1_JOB_ID: bwa_mem_picard_sort_sam.HI.1943.008.Index_14.LCC9_CTRL_MED1_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.1943.008.Index_14.LCC9_CTRL_MED1_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_1_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.1943.008.Index_14.LCC9_CTRL_MED1_ChIP8_GB_rep1.bam.9c9d2c9a1b6179233a32a407398ff692.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.1943.008.Index_14.LCC9_CTRL_MED1_ChIP8_GB_rep1.bam.9c9d2c9a1b6179233a32a407398ff692.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/LCC9_CTRL_MED1/HI.1943.008.Index_14.LCC9_CTRL_MED1_ChIP8_GB_rep1.bam && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:HI.1943.008.Index_14.LCC9_CTRL_MED1_ChIP8_GB_rep1.bam	SM:LCC9_CTRL_MED1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/LCC9_CTRL_MED1/HI.1943.008.Index_14.LCC9_CTRL_MED1_ChIP8_GB_rep1.bam.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/LCC9_CTRL_MED1/HI.1943.008.Index_14.LCC9_CTRL_MED1_ChIP8_GB_rep1.bam/HI.1943.008.Index_14.LCC9_CTRL_MED1_ChIP8_GB_rep1.bam.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.HI.1943.008.Index_14.LCC9_CTRL_MED1_ChIP8_GB_rep1.bam.9c9d2c9a1b6179233a32a407398ff692.mugqic.done
)
bwa_mem_picard_sort_sam_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_2_JOB_ID: bwa_mem_picard_sort_sam.HI.1943.008.Index_15.LCC9_E2_WCE_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.1943.008.Index_15.LCC9_E2_WCE_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_2_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.1943.008.Index_15.LCC9_E2_WCE_ChIP8_GB_rep1.bam.8121573482cb6876cb1e6a302e645629.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.1943.008.Index_15.LCC9_E2_WCE_ChIP8_GB_rep1.bam.8121573482cb6876cb1e6a302e645629.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/LCC9_E2_WCE/HI.1943.008.Index_15.LCC9_E2_WCE_ChIP8_GB_rep1.bam && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:HI.1943.008.Index_15.LCC9_E2_WCE_ChIP8_GB_rep1.bam	SM:LCC9_E2_WCE	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/LCC9_E2_WCE/HI.1943.008.Index_15.LCC9_E2_WCE_ChIP8_GB_rep1.bam.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/LCC9_E2_WCE/HI.1943.008.Index_15.LCC9_E2_WCE_ChIP8_GB_rep1.bam/HI.1943.008.Index_15.LCC9_E2_WCE_ChIP8_GB_rep1.bam.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.HI.1943.008.Index_15.LCC9_E2_WCE_ChIP8_GB_rep1.bam.8121573482cb6876cb1e6a302e645629.mugqic.done
)
bwa_mem_picard_sort_sam_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_3_JOB_ID: bwa_mem_picard_sort_sam.HI.1943.008.Index_16.LCC9_E2_MED1_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.1943.008.Index_16.LCC9_E2_MED1_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_3_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.1943.008.Index_16.LCC9_E2_MED1_ChIP8_GB_rep1.bam.503c3072d0c049ae36c8f81086a84d0c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.1943.008.Index_16.LCC9_E2_MED1_ChIP8_GB_rep1.bam.503c3072d0c049ae36c8f81086a84d0c.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/LCC9_E2_MED1/HI.1943.008.Index_16.LCC9_E2_MED1_ChIP8_GB_rep1.bam && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:HI.1943.008.Index_16.LCC9_E2_MED1_ChIP8_GB_rep1.bam	SM:LCC9_E2_MED1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/LCC9_E2_MED1/HI.1943.008.Index_16.LCC9_E2_MED1_ChIP8_GB_rep1.bam.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/LCC9_E2_MED1/HI.1943.008.Index_16.LCC9_E2_MED1_ChIP8_GB_rep1.bam/HI.1943.008.Index_16.LCC9_E2_MED1_ChIP8_GB_rep1.bam.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.HI.1943.008.Index_16.LCC9_E2_MED1_ChIP8_GB_rep1.bam.503c3072d0c049ae36c8f81086a84d0c.mugqic.done
)
bwa_mem_picard_sort_sam_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_4_JOB_ID: bwa_mem_picard_sort_sam.HI.1943.008.Index_6.LCC9_CTRL_WCE_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.1943.008.Index_6.LCC9_CTRL_WCE_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_4_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.1943.008.Index_6.LCC9_CTRL_WCE_ChIP8_GB_rep1.bam.cc0a2dedbd6fd893879b71b2daa04744.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.1943.008.Index_6.LCC9_CTRL_WCE_ChIP8_GB_rep1.bam.cc0a2dedbd6fd893879b71b2daa04744.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/LCC9_CTRL_WCE/HI.1943.008.Index_6.LCC9_CTRL_WCE_ChIP8_GB_rep1.bam && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:HI.1943.008.Index_6.LCC9_CTRL_WCE_ChIP8_GB_rep1.bam	SM:LCC9_CTRL_WCE	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/LCC9_CTRL_WCE/HI.1943.008.Index_6.LCC9_CTRL_WCE_ChIP8_GB_rep1.bam.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/LCC9_CTRL_WCE/HI.1943.008.Index_6.LCC9_CTRL_WCE_ChIP8_GB_rep1.bam/HI.1943.008.Index_6.LCC9_CTRL_WCE_ChIP8_GB_rep1.bam.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.HI.1943.008.Index_6.LCC9_CTRL_WCE_ChIP8_GB_rep1.bam.cc0a2dedbd6fd893879b71b2daa04744.mugqic.done
)
bwa_mem_picard_sort_sam_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_5_JOB_ID: bwa_mem_picard_sort_sam.HI.1997.004.Index_12.LCC9_CTRL_SMC1_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.1997.004.Index_12.LCC9_CTRL_SMC1_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_5_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.1997.004.Index_12.LCC9_CTRL_SMC1_ChIP8_GB_rep1.bam.87cdd45636cdd37e561ceb5ec6c99af3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.1997.004.Index_12.LCC9_CTRL_SMC1_ChIP8_GB_rep1.bam.87cdd45636cdd37e561ceb5ec6c99af3.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/LCC9_CTRL_SMC1/HI.1997.004.Index_12.LCC9_CTRL_SMC1_ChIP8_GB_rep1.bam && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:HI.1997.004.Index_12.LCC9_CTRL_SMC1_ChIP8_GB_rep1.bam	SM:LCC9_CTRL_SMC1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/LCC9_CTRL_SMC1/HI.1997.004.Index_12.LCC9_CTRL_SMC1_ChIP8_GB_rep1.bam.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/LCC9_CTRL_SMC1/HI.1997.004.Index_12.LCC9_CTRL_SMC1_ChIP8_GB_rep1.bam/HI.1997.004.Index_12.LCC9_CTRL_SMC1_ChIP8_GB_rep1.bam.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.HI.1997.004.Index_12.LCC9_CTRL_SMC1_ChIP8_GB_rep1.bam.87cdd45636cdd37e561ceb5ec6c99af3.mugqic.done
)
bwa_mem_picard_sort_sam_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_6_JOB_ID: bwa_mem_picard_sort_sam.HI.1997.004.Index_19.LCC9_E2_SMC1_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.1997.004.Index_19.LCC9_E2_SMC1_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_6_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.1997.004.Index_19.LCC9_E2_SMC1_ChIP8_GB_rep1.bam.9b2ee855731359dabd447eabba432391.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.1997.004.Index_19.LCC9_E2_SMC1_ChIP8_GB_rep1.bam.9b2ee855731359dabd447eabba432391.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/LCC9_E2_SMC1/HI.1997.004.Index_19.LCC9_E2_SMC1_ChIP8_GB_rep1.bam && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:HI.1997.004.Index_19.LCC9_E2_SMC1_ChIP8_GB_rep1.bam	SM:LCC9_E2_SMC1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/LCC9_E2_SMC1/HI.1997.004.Index_19.LCC9_E2_SMC1_ChIP8_GB_rep1.bam.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/LCC9_E2_SMC1/HI.1997.004.Index_19.LCC9_E2_SMC1_ChIP8_GB_rep1.bam/HI.1997.004.Index_19.LCC9_E2_SMC1_ChIP8_GB_rep1.bam.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.HI.1997.004.Index_19.LCC9_E2_SMC1_ChIP8_GB_rep1.bam.9b2ee855731359dabd447eabba432391.mugqic.done
)
bwa_mem_picard_sort_sam_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_7_JOB_ID: bwa_mem_picard_sort_sam.HI.1997.004.Index_8.LCC9_E2_NIPBL_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.1997.004.Index_8.LCC9_E2_NIPBL_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_7_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.1997.004.Index_8.LCC9_E2_NIPBL_ChIP8_GB_rep1.bam.888db2e539bb745f81233e4bbeef33ef.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.1997.004.Index_8.LCC9_E2_NIPBL_ChIP8_GB_rep1.bam.888db2e539bb745f81233e4bbeef33ef.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/LCC9_E2_NIPBL/HI.1997.004.Index_8.LCC9_E2_NIPBL_ChIP8_GB_rep1.bam && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:HI.1997.004.Index_8.LCC9_E2_NIPBL_ChIP8_GB_rep1.bam	SM:LCC9_E2_NIPBL	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/LCC9_E2_NIPBL/HI.1997.004.Index_8.LCC9_E2_NIPBL_ChIP8_GB_rep1.bam.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/LCC9_E2_NIPBL/HI.1997.004.Index_8.LCC9_E2_NIPBL_ChIP8_GB_rep1.bam/HI.1997.004.Index_8.LCC9_E2_NIPBL_ChIP8_GB_rep1.bam.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.HI.1997.004.Index_8.LCC9_E2_NIPBL_ChIP8_GB_rep1.bam.888db2e539bb745f81233e4bbeef33ef.mugqic.done
)
bwa_mem_picard_sort_sam_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_8_JOB_ID: bwa_mem_picard_sort_sam.HI.1997.004.Index_9.LCC9_CTRL_NIPBL_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.1997.004.Index_9.LCC9_CTRL_NIPBL_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_8_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.1997.004.Index_9.LCC9_CTRL_NIPBL_ChIP8_GB_rep1.bam.1d7895afddafdb18f6b1ae821db94ced.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.1997.004.Index_9.LCC9_CTRL_NIPBL_ChIP8_GB_rep1.bam.1d7895afddafdb18f6b1ae821db94ced.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/LCC9_CTRL_NIPBL/HI.1997.004.Index_9.LCC9_CTRL_NIPBL_ChIP8_GB_rep1.bam && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:HI.1997.004.Index_9.LCC9_CTRL_NIPBL_ChIP8_GB_rep1.bam	SM:LCC9_CTRL_NIPBL	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/LCC9_CTRL_NIPBL/HI.1997.004.Index_9.LCC9_CTRL_NIPBL_ChIP8_GB_rep1.bam.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/LCC9_CTRL_NIPBL/HI.1997.004.Index_9.LCC9_CTRL_NIPBL_ChIP8_GB_rep1.bam/HI.1997.004.Index_9.LCC9_CTRL_NIPBL_ChIP8_GB_rep1.bam.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.HI.1997.004.Index_9.LCC9_CTRL_NIPBL_ChIP8_GB_rep1.bam.1d7895afddafdb18f6b1ae821db94ced.mugqic.done
)
bwa_mem_picard_sort_sam_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_9_JOB_ID: bwa_mem_picard_sort_sam.HI.1997.005.Index_20.LCC9_CTRL_POL2_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.1997.005.Index_20.LCC9_CTRL_POL2_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_9_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.1997.005.Index_20.LCC9_CTRL_POL2_ChIP8_GB_rep1.bam.890be09ea5f37262c1edc3b592b79c1d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.1997.005.Index_20.LCC9_CTRL_POL2_ChIP8_GB_rep1.bam.890be09ea5f37262c1edc3b592b79c1d.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/LCC9_CTRL_POL2/HI.1997.005.Index_20.LCC9_CTRL_POL2_ChIP8_GB_rep1.bam && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:HI.1997.005.Index_20.LCC9_CTRL_POL2_ChIP8_GB_rep1.bam	SM:LCC9_CTRL_POL2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/LCC9_CTRL_POL2/HI.1997.005.Index_20.LCC9_CTRL_POL2_ChIP8_GB_rep1.bam.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/LCC9_CTRL_POL2/HI.1997.005.Index_20.LCC9_CTRL_POL2_ChIP8_GB_rep1.bam/HI.1997.005.Index_20.LCC9_CTRL_POL2_ChIP8_GB_rep1.bam.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.HI.1997.005.Index_20.LCC9_CTRL_POL2_ChIP8_GB_rep1.bam.890be09ea5f37262c1edc3b592b79c1d.mugqic.done
)
bwa_mem_picard_sort_sam_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_10_JOB_ID: bwa_mem_picard_sort_sam.HI.1997.005.Index_22.LCC9_E2_POL2_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.1997.005.Index_22.LCC9_E2_POL2_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_10_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.1997.005.Index_22.LCC9_E2_POL2_ChIP8_GB_rep1.bam.0a524bfb99f4733bb0490f84e4fa3b10.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.1997.005.Index_22.LCC9_E2_POL2_ChIP8_GB_rep1.bam.0a524bfb99f4733bb0490f84e4fa3b10.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/LCC9_E2_POL2/HI.1997.005.Index_22.LCC9_E2_POL2_ChIP8_GB_rep1.bam && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:HI.1997.005.Index_22.LCC9_E2_POL2_ChIP8_GB_rep1.bam	SM:LCC9_E2_POL2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/LCC9_E2_POL2/HI.1997.005.Index_22.LCC9_E2_POL2_ChIP8_GB_rep1.bam.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/LCC9_E2_POL2/HI.1997.005.Index_22.LCC9_E2_POL2_ChIP8_GB_rep1.bam/HI.1997.005.Index_22.LCC9_E2_POL2_ChIP8_GB_rep1.bam.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.HI.1997.005.Index_22.LCC9_E2_POL2_ChIP8_GB_rep1.bam.0a524bfb99f4733bb0490f84e4fa3b10.mugqic.done
)
bwa_mem_picard_sort_sam_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_11_JOB_ID: bwa_mem_picard_sort_sam.HI.1997.005.Index_23.LCC9_CTRL_ERA_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.1997.005.Index_23.LCC9_CTRL_ERA_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_11_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.1997.005.Index_23.LCC9_CTRL_ERA_ChIP8_GB_rep1.bam.da8f68303ec4099c629092efad6416e2.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.1997.005.Index_23.LCC9_CTRL_ERA_ChIP8_GB_rep1.bam.da8f68303ec4099c629092efad6416e2.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/LCC9_CTRL_ERA/HI.1997.005.Index_23.LCC9_CTRL_ERA_ChIP8_GB_rep1.bam && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:HI.1997.005.Index_23.LCC9_CTRL_ERA_ChIP8_GB_rep1.bam	SM:LCC9_CTRL_ERA	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/LCC9_CTRL_ERA/HI.1997.005.Index_23.LCC9_CTRL_ERA_ChIP8_GB_rep1.bam.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/LCC9_CTRL_ERA/HI.1997.005.Index_23.LCC9_CTRL_ERA_ChIP8_GB_rep1.bam/HI.1997.005.Index_23.LCC9_CTRL_ERA_ChIP8_GB_rep1.bam.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.HI.1997.005.Index_23.LCC9_CTRL_ERA_ChIP8_GB_rep1.bam.da8f68303ec4099c629092efad6416e2.mugqic.done
)
bwa_mem_picard_sort_sam_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_12_JOB_ID: bwa_mem_picard_sort_sam.HI.1997.005.Index_27.LCC9_E2_ERA_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.1997.005.Index_27.LCC9_E2_ERA_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_12_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.1997.005.Index_27.LCC9_E2_ERA_ChIP8_GB_rep1.bam.94a3ca82e541fc09a8c9b6c05563bd83.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.1997.005.Index_27.LCC9_E2_ERA_ChIP8_GB_rep1.bam.94a3ca82e541fc09a8c9b6c05563bd83.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/LCC9_E2_ERA/HI.1997.005.Index_27.LCC9_E2_ERA_ChIP8_GB_rep1.bam && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:HI.1997.005.Index_27.LCC9_E2_ERA_ChIP8_GB_rep1.bam	SM:LCC9_E2_ERA	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/LCC9_E2_ERA/HI.1997.005.Index_27.LCC9_E2_ERA_ChIP8_GB_rep1.bam.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/LCC9_E2_ERA/HI.1997.005.Index_27.LCC9_E2_ERA_ChIP8_GB_rep1.bam/HI.1997.005.Index_27.LCC9_E2_ERA_ChIP8_GB_rep1.bam.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.HI.1997.005.Index_27.LCC9_E2_ERA_ChIP8_GB_rep1.bam.94a3ca82e541fc09a8c9b6c05563bd83.mugqic.done
)
bwa_mem_picard_sort_sam_12_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_13_JOB_ID: bwa_mem_picard_sort_sam_report
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam_report
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_1_JOB_ID:$bwa_mem_picard_sort_sam_2_JOB_ID:$bwa_mem_picard_sort_sam_3_JOB_ID:$bwa_mem_picard_sort_sam_4_JOB_ID:$bwa_mem_picard_sort_sam_5_JOB_ID:$bwa_mem_picard_sort_sam_6_JOB_ID:$bwa_mem_picard_sort_sam_7_JOB_ID:$bwa_mem_picard_sort_sam_8_JOB_ID:$bwa_mem_picard_sort_sam_9_JOB_ID:$bwa_mem_picard_sort_sam_10_JOB_ID:$bwa_mem_picard_sort_sam_11_JOB_ID:$bwa_mem_picard_sort_sam_12_JOB_ID
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
bwa_mem_picard_sort_sam_13_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_13_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: samtools_view_filter
#-------------------------------------------------------------------------------
STEP=samtools_view_filter
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_1_JOB_ID: samtools_view_filter.HI.1943.008.Index_14.LCC9_CTRL_MED1_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.1943.008.Index_14.LCC9_CTRL_MED1_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_1_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.1943.008.Index_14.LCC9_CTRL_MED1_ChIP8_GB_rep1.bam.3cb20451b210b47fbdf7ecbbe1a475a5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.1943.008.Index_14.LCC9_CTRL_MED1_ChIP8_GB_rep1.bam.3cb20451b210b47fbdf7ecbbe1a475a5.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/LCC9_CTRL_MED1/HI.1943.008.Index_14.LCC9_CTRL_MED1_ChIP8_GB_rep1.bam/HI.1943.008.Index_14.LCC9_CTRL_MED1_ChIP8_GB_rep1.bam.sorted.bam \
  > alignment/LCC9_CTRL_MED1/HI.1943.008.Index_14.LCC9_CTRL_MED1_ChIP8_GB_rep1.bam/HI.1943.008.Index_14.LCC9_CTRL_MED1_ChIP8_GB_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.1943.008.Index_14.LCC9_CTRL_MED1_ChIP8_GB_rep1.bam.3cb20451b210b47fbdf7ecbbe1a475a5.mugqic.done
)
samtools_view_filter_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_2_JOB_ID: samtools_view_filter.HI.1943.008.Index_15.LCC9_E2_WCE_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.1943.008.Index_15.LCC9_E2_WCE_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_2_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.1943.008.Index_15.LCC9_E2_WCE_ChIP8_GB_rep1.bam.b21494ad8272d2b66ac22316806aa94a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.1943.008.Index_15.LCC9_E2_WCE_ChIP8_GB_rep1.bam.b21494ad8272d2b66ac22316806aa94a.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/LCC9_E2_WCE/HI.1943.008.Index_15.LCC9_E2_WCE_ChIP8_GB_rep1.bam/HI.1943.008.Index_15.LCC9_E2_WCE_ChIP8_GB_rep1.bam.sorted.bam \
  > alignment/LCC9_E2_WCE/HI.1943.008.Index_15.LCC9_E2_WCE_ChIP8_GB_rep1.bam/HI.1943.008.Index_15.LCC9_E2_WCE_ChIP8_GB_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.1943.008.Index_15.LCC9_E2_WCE_ChIP8_GB_rep1.bam.b21494ad8272d2b66ac22316806aa94a.mugqic.done
)
samtools_view_filter_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_3_JOB_ID: samtools_view_filter.HI.1943.008.Index_16.LCC9_E2_MED1_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.1943.008.Index_16.LCC9_E2_MED1_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_3_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.1943.008.Index_16.LCC9_E2_MED1_ChIP8_GB_rep1.bam.ce0419cad836bbd010898266a67cf799.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.1943.008.Index_16.LCC9_E2_MED1_ChIP8_GB_rep1.bam.ce0419cad836bbd010898266a67cf799.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/LCC9_E2_MED1/HI.1943.008.Index_16.LCC9_E2_MED1_ChIP8_GB_rep1.bam/HI.1943.008.Index_16.LCC9_E2_MED1_ChIP8_GB_rep1.bam.sorted.bam \
  > alignment/LCC9_E2_MED1/HI.1943.008.Index_16.LCC9_E2_MED1_ChIP8_GB_rep1.bam/HI.1943.008.Index_16.LCC9_E2_MED1_ChIP8_GB_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.1943.008.Index_16.LCC9_E2_MED1_ChIP8_GB_rep1.bam.ce0419cad836bbd010898266a67cf799.mugqic.done
)
samtools_view_filter_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_4_JOB_ID: samtools_view_filter.HI.1943.008.Index_6.LCC9_CTRL_WCE_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.1943.008.Index_6.LCC9_CTRL_WCE_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_4_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.1943.008.Index_6.LCC9_CTRL_WCE_ChIP8_GB_rep1.bam.7cfb3e514a73e3aaf2eb307885318811.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.1943.008.Index_6.LCC9_CTRL_WCE_ChIP8_GB_rep1.bam.7cfb3e514a73e3aaf2eb307885318811.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/LCC9_CTRL_WCE/HI.1943.008.Index_6.LCC9_CTRL_WCE_ChIP8_GB_rep1.bam/HI.1943.008.Index_6.LCC9_CTRL_WCE_ChIP8_GB_rep1.bam.sorted.bam \
  > alignment/LCC9_CTRL_WCE/HI.1943.008.Index_6.LCC9_CTRL_WCE_ChIP8_GB_rep1.bam/HI.1943.008.Index_6.LCC9_CTRL_WCE_ChIP8_GB_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.1943.008.Index_6.LCC9_CTRL_WCE_ChIP8_GB_rep1.bam.7cfb3e514a73e3aaf2eb307885318811.mugqic.done
)
samtools_view_filter_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_5_JOB_ID: samtools_view_filter.HI.1997.004.Index_12.LCC9_CTRL_SMC1_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.1997.004.Index_12.LCC9_CTRL_SMC1_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_5_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.1997.004.Index_12.LCC9_CTRL_SMC1_ChIP8_GB_rep1.bam.e0a58edcc0c8042dfda23022bbd4abd5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.1997.004.Index_12.LCC9_CTRL_SMC1_ChIP8_GB_rep1.bam.e0a58edcc0c8042dfda23022bbd4abd5.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/LCC9_CTRL_SMC1/HI.1997.004.Index_12.LCC9_CTRL_SMC1_ChIP8_GB_rep1.bam/HI.1997.004.Index_12.LCC9_CTRL_SMC1_ChIP8_GB_rep1.bam.sorted.bam \
  > alignment/LCC9_CTRL_SMC1/HI.1997.004.Index_12.LCC9_CTRL_SMC1_ChIP8_GB_rep1.bam/HI.1997.004.Index_12.LCC9_CTRL_SMC1_ChIP8_GB_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.1997.004.Index_12.LCC9_CTRL_SMC1_ChIP8_GB_rep1.bam.e0a58edcc0c8042dfda23022bbd4abd5.mugqic.done
)
samtools_view_filter_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_6_JOB_ID: samtools_view_filter.HI.1997.004.Index_19.LCC9_E2_SMC1_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.1997.004.Index_19.LCC9_E2_SMC1_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_6_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.1997.004.Index_19.LCC9_E2_SMC1_ChIP8_GB_rep1.bam.0170628414c512cf38cfbb2b8a5316ea.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.1997.004.Index_19.LCC9_E2_SMC1_ChIP8_GB_rep1.bam.0170628414c512cf38cfbb2b8a5316ea.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/LCC9_E2_SMC1/HI.1997.004.Index_19.LCC9_E2_SMC1_ChIP8_GB_rep1.bam/HI.1997.004.Index_19.LCC9_E2_SMC1_ChIP8_GB_rep1.bam.sorted.bam \
  > alignment/LCC9_E2_SMC1/HI.1997.004.Index_19.LCC9_E2_SMC1_ChIP8_GB_rep1.bam/HI.1997.004.Index_19.LCC9_E2_SMC1_ChIP8_GB_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.1997.004.Index_19.LCC9_E2_SMC1_ChIP8_GB_rep1.bam.0170628414c512cf38cfbb2b8a5316ea.mugqic.done
)
samtools_view_filter_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_7_JOB_ID: samtools_view_filter.HI.1997.004.Index_8.LCC9_E2_NIPBL_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.1997.004.Index_8.LCC9_E2_NIPBL_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_7_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.1997.004.Index_8.LCC9_E2_NIPBL_ChIP8_GB_rep1.bam.fa46a1f5b9f17a220ca0c3154636fa7d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.1997.004.Index_8.LCC9_E2_NIPBL_ChIP8_GB_rep1.bam.fa46a1f5b9f17a220ca0c3154636fa7d.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/LCC9_E2_NIPBL/HI.1997.004.Index_8.LCC9_E2_NIPBL_ChIP8_GB_rep1.bam/HI.1997.004.Index_8.LCC9_E2_NIPBL_ChIP8_GB_rep1.bam.sorted.bam \
  > alignment/LCC9_E2_NIPBL/HI.1997.004.Index_8.LCC9_E2_NIPBL_ChIP8_GB_rep1.bam/HI.1997.004.Index_8.LCC9_E2_NIPBL_ChIP8_GB_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.1997.004.Index_8.LCC9_E2_NIPBL_ChIP8_GB_rep1.bam.fa46a1f5b9f17a220ca0c3154636fa7d.mugqic.done
)
samtools_view_filter_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_8_JOB_ID: samtools_view_filter.HI.1997.004.Index_9.LCC9_CTRL_NIPBL_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.1997.004.Index_9.LCC9_CTRL_NIPBL_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_8_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.1997.004.Index_9.LCC9_CTRL_NIPBL_ChIP8_GB_rep1.bam.35a8014d1a4682d98eb3225a734f5d71.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.1997.004.Index_9.LCC9_CTRL_NIPBL_ChIP8_GB_rep1.bam.35a8014d1a4682d98eb3225a734f5d71.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/LCC9_CTRL_NIPBL/HI.1997.004.Index_9.LCC9_CTRL_NIPBL_ChIP8_GB_rep1.bam/HI.1997.004.Index_9.LCC9_CTRL_NIPBL_ChIP8_GB_rep1.bam.sorted.bam \
  > alignment/LCC9_CTRL_NIPBL/HI.1997.004.Index_9.LCC9_CTRL_NIPBL_ChIP8_GB_rep1.bam/HI.1997.004.Index_9.LCC9_CTRL_NIPBL_ChIP8_GB_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.1997.004.Index_9.LCC9_CTRL_NIPBL_ChIP8_GB_rep1.bam.35a8014d1a4682d98eb3225a734f5d71.mugqic.done
)
samtools_view_filter_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_9_JOB_ID: samtools_view_filter.HI.1997.005.Index_20.LCC9_CTRL_POL2_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.1997.005.Index_20.LCC9_CTRL_POL2_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_9_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.1997.005.Index_20.LCC9_CTRL_POL2_ChIP8_GB_rep1.bam.f2bf79d58d28028386124301b86336b5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.1997.005.Index_20.LCC9_CTRL_POL2_ChIP8_GB_rep1.bam.f2bf79d58d28028386124301b86336b5.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/LCC9_CTRL_POL2/HI.1997.005.Index_20.LCC9_CTRL_POL2_ChIP8_GB_rep1.bam/HI.1997.005.Index_20.LCC9_CTRL_POL2_ChIP8_GB_rep1.bam.sorted.bam \
  > alignment/LCC9_CTRL_POL2/HI.1997.005.Index_20.LCC9_CTRL_POL2_ChIP8_GB_rep1.bam/HI.1997.005.Index_20.LCC9_CTRL_POL2_ChIP8_GB_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.1997.005.Index_20.LCC9_CTRL_POL2_ChIP8_GB_rep1.bam.f2bf79d58d28028386124301b86336b5.mugqic.done
)
samtools_view_filter_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_10_JOB_ID: samtools_view_filter.HI.1997.005.Index_22.LCC9_E2_POL2_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.1997.005.Index_22.LCC9_E2_POL2_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_10_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.1997.005.Index_22.LCC9_E2_POL2_ChIP8_GB_rep1.bam.c5ed0ae485c35276e1927a354fcf1a60.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.1997.005.Index_22.LCC9_E2_POL2_ChIP8_GB_rep1.bam.c5ed0ae485c35276e1927a354fcf1a60.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/LCC9_E2_POL2/HI.1997.005.Index_22.LCC9_E2_POL2_ChIP8_GB_rep1.bam/HI.1997.005.Index_22.LCC9_E2_POL2_ChIP8_GB_rep1.bam.sorted.bam \
  > alignment/LCC9_E2_POL2/HI.1997.005.Index_22.LCC9_E2_POL2_ChIP8_GB_rep1.bam/HI.1997.005.Index_22.LCC9_E2_POL2_ChIP8_GB_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.1997.005.Index_22.LCC9_E2_POL2_ChIP8_GB_rep1.bam.c5ed0ae485c35276e1927a354fcf1a60.mugqic.done
)
samtools_view_filter_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_11_JOB_ID: samtools_view_filter.HI.1997.005.Index_23.LCC9_CTRL_ERA_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.1997.005.Index_23.LCC9_CTRL_ERA_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_11_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.1997.005.Index_23.LCC9_CTRL_ERA_ChIP8_GB_rep1.bam.da139f186e9e0495e881b0d139fba2b4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.1997.005.Index_23.LCC9_CTRL_ERA_ChIP8_GB_rep1.bam.da139f186e9e0495e881b0d139fba2b4.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/LCC9_CTRL_ERA/HI.1997.005.Index_23.LCC9_CTRL_ERA_ChIP8_GB_rep1.bam/HI.1997.005.Index_23.LCC9_CTRL_ERA_ChIP8_GB_rep1.bam.sorted.bam \
  > alignment/LCC9_CTRL_ERA/HI.1997.005.Index_23.LCC9_CTRL_ERA_ChIP8_GB_rep1.bam/HI.1997.005.Index_23.LCC9_CTRL_ERA_ChIP8_GB_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.1997.005.Index_23.LCC9_CTRL_ERA_ChIP8_GB_rep1.bam.da139f186e9e0495e881b0d139fba2b4.mugqic.done
)
samtools_view_filter_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_12_JOB_ID: samtools_view_filter.HI.1997.005.Index_27.LCC9_E2_ERA_ChIP8_GB_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.1997.005.Index_27.LCC9_E2_ERA_ChIP8_GB_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_12_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.1997.005.Index_27.LCC9_E2_ERA_ChIP8_GB_rep1.bam.db26b746577864cf1408ef306b5e8bec.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.1997.005.Index_27.LCC9_E2_ERA_ChIP8_GB_rep1.bam.db26b746577864cf1408ef306b5e8bec.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/LCC9_E2_ERA/HI.1997.005.Index_27.LCC9_E2_ERA_ChIP8_GB_rep1.bam/HI.1997.005.Index_27.LCC9_E2_ERA_ChIP8_GB_rep1.bam.sorted.bam \
  > alignment/LCC9_E2_ERA/HI.1997.005.Index_27.LCC9_E2_ERA_ChIP8_GB_rep1.bam/HI.1997.005.Index_27.LCC9_E2_ERA_ChIP8_GB_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.1997.005.Index_27.LCC9_E2_ERA_ChIP8_GB_rep1.bam.db26b746577864cf1408ef306b5e8bec.mugqic.done
)
samtools_view_filter_12_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_13_JOB_ID: samtools_view_filter_report
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter_report
JOB_DEPENDENCIES=$samtools_view_filter_1_JOB_ID:$samtools_view_filter_2_JOB_ID:$samtools_view_filter_3_JOB_ID:$samtools_view_filter_4_JOB_ID:$samtools_view_filter_5_JOB_ID:$samtools_view_filter_6_JOB_ID:$samtools_view_filter_7_JOB_ID:$samtools_view_filter_8_JOB_ID:$samtools_view_filter_9_JOB_ID:$samtools_view_filter_10_JOB_ID:$samtools_view_filter_11_JOB_ID:$samtools_view_filter_12_JOB_ID
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
samtools_view_filter_13_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_13_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: picard_merge_sam_files
#-------------------------------------------------------------------------------
STEP=picard_merge_sam_files
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_1_JOB_ID: symlink_readset_sample_bam.LCC9_CTRL_MED1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.LCC9_CTRL_MED1
JOB_DEPENDENCIES=$samtools_view_filter_1_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.LCC9_CTRL_MED1.c7fef4fcfe76eb03b3aa9e5523676f06.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.LCC9_CTRL_MED1.c7fef4fcfe76eb03b3aa9e5523676f06.mugqic.done'
mkdir -p alignment/LCC9_CTRL_MED1 && \
ln -s -f HI.1943.008.Index_14.LCC9_CTRL_MED1_ChIP8_GB_rep1.bam/HI.1943.008.Index_14.LCC9_CTRL_MED1_ChIP8_GB_rep1.bam.sorted.filtered.bam alignment/LCC9_CTRL_MED1/LCC9_CTRL_MED1.merged.bam
symlink_readset_sample_bam.LCC9_CTRL_MED1.c7fef4fcfe76eb03b3aa9e5523676f06.mugqic.done
)
picard_merge_sam_files_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_2_JOB_ID: symlink_readset_sample_bam.LCC9_E2_WCE
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.LCC9_E2_WCE
JOB_DEPENDENCIES=$samtools_view_filter_2_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.LCC9_E2_WCE.6525d08201c34c4b4e782a0c69d14e84.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.LCC9_E2_WCE.6525d08201c34c4b4e782a0c69d14e84.mugqic.done'
mkdir -p alignment/LCC9_E2_WCE && \
ln -s -f HI.1943.008.Index_15.LCC9_E2_WCE_ChIP8_GB_rep1.bam/HI.1943.008.Index_15.LCC9_E2_WCE_ChIP8_GB_rep1.bam.sorted.filtered.bam alignment/LCC9_E2_WCE/LCC9_E2_WCE.merged.bam
symlink_readset_sample_bam.LCC9_E2_WCE.6525d08201c34c4b4e782a0c69d14e84.mugqic.done
)
picard_merge_sam_files_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_3_JOB_ID: symlink_readset_sample_bam.LCC9_E2_MED1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.LCC9_E2_MED1
JOB_DEPENDENCIES=$samtools_view_filter_3_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.LCC9_E2_MED1.bec47288d4941ec5ac31b6da8350f171.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.LCC9_E2_MED1.bec47288d4941ec5ac31b6da8350f171.mugqic.done'
mkdir -p alignment/LCC9_E2_MED1 && \
ln -s -f HI.1943.008.Index_16.LCC9_E2_MED1_ChIP8_GB_rep1.bam/HI.1943.008.Index_16.LCC9_E2_MED1_ChIP8_GB_rep1.bam.sorted.filtered.bam alignment/LCC9_E2_MED1/LCC9_E2_MED1.merged.bam
symlink_readset_sample_bam.LCC9_E2_MED1.bec47288d4941ec5ac31b6da8350f171.mugqic.done
)
picard_merge_sam_files_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_4_JOB_ID: symlink_readset_sample_bam.LCC9_CTRL_WCE
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.LCC9_CTRL_WCE
JOB_DEPENDENCIES=$samtools_view_filter_4_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.LCC9_CTRL_WCE.8a217add0381eee389f32e444e3363ee.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.LCC9_CTRL_WCE.8a217add0381eee389f32e444e3363ee.mugqic.done'
mkdir -p alignment/LCC9_CTRL_WCE && \
ln -s -f HI.1943.008.Index_6.LCC9_CTRL_WCE_ChIP8_GB_rep1.bam/HI.1943.008.Index_6.LCC9_CTRL_WCE_ChIP8_GB_rep1.bam.sorted.filtered.bam alignment/LCC9_CTRL_WCE/LCC9_CTRL_WCE.merged.bam
symlink_readset_sample_bam.LCC9_CTRL_WCE.8a217add0381eee389f32e444e3363ee.mugqic.done
)
picard_merge_sam_files_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_5_JOB_ID: symlink_readset_sample_bam.LCC9_CTRL_SMC1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.LCC9_CTRL_SMC1
JOB_DEPENDENCIES=$samtools_view_filter_5_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.LCC9_CTRL_SMC1.f88b797eeafbf28ee058f86588363740.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.LCC9_CTRL_SMC1.f88b797eeafbf28ee058f86588363740.mugqic.done'
mkdir -p alignment/LCC9_CTRL_SMC1 && \
ln -s -f HI.1997.004.Index_12.LCC9_CTRL_SMC1_ChIP8_GB_rep1.bam/HI.1997.004.Index_12.LCC9_CTRL_SMC1_ChIP8_GB_rep1.bam.sorted.filtered.bam alignment/LCC9_CTRL_SMC1/LCC9_CTRL_SMC1.merged.bam
symlink_readset_sample_bam.LCC9_CTRL_SMC1.f88b797eeafbf28ee058f86588363740.mugqic.done
)
picard_merge_sam_files_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_6_JOB_ID: symlink_readset_sample_bam.LCC9_E2_SMC1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.LCC9_E2_SMC1
JOB_DEPENDENCIES=$samtools_view_filter_6_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.LCC9_E2_SMC1.7b229605b87c5076f29dbcae823e79ba.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.LCC9_E2_SMC1.7b229605b87c5076f29dbcae823e79ba.mugqic.done'
mkdir -p alignment/LCC9_E2_SMC1 && \
ln -s -f HI.1997.004.Index_19.LCC9_E2_SMC1_ChIP8_GB_rep1.bam/HI.1997.004.Index_19.LCC9_E2_SMC1_ChIP8_GB_rep1.bam.sorted.filtered.bam alignment/LCC9_E2_SMC1/LCC9_E2_SMC1.merged.bam
symlink_readset_sample_bam.LCC9_E2_SMC1.7b229605b87c5076f29dbcae823e79ba.mugqic.done
)
picard_merge_sam_files_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_7_JOB_ID: symlink_readset_sample_bam.LCC9_E2_NIPBL
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.LCC9_E2_NIPBL
JOB_DEPENDENCIES=$samtools_view_filter_7_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.LCC9_E2_NIPBL.3c63bc0c10ff17112f7b4b15d34ecead.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.LCC9_E2_NIPBL.3c63bc0c10ff17112f7b4b15d34ecead.mugqic.done'
mkdir -p alignment/LCC9_E2_NIPBL && \
ln -s -f HI.1997.004.Index_8.LCC9_E2_NIPBL_ChIP8_GB_rep1.bam/HI.1997.004.Index_8.LCC9_E2_NIPBL_ChIP8_GB_rep1.bam.sorted.filtered.bam alignment/LCC9_E2_NIPBL/LCC9_E2_NIPBL.merged.bam
symlink_readset_sample_bam.LCC9_E2_NIPBL.3c63bc0c10ff17112f7b4b15d34ecead.mugqic.done
)
picard_merge_sam_files_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_8_JOB_ID: symlink_readset_sample_bam.LCC9_CTRL_NIPBL
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.LCC9_CTRL_NIPBL
JOB_DEPENDENCIES=$samtools_view_filter_8_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.LCC9_CTRL_NIPBL.1ab1060995df79508cb6e310056bfbfe.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.LCC9_CTRL_NIPBL.1ab1060995df79508cb6e310056bfbfe.mugqic.done'
mkdir -p alignment/LCC9_CTRL_NIPBL && \
ln -s -f HI.1997.004.Index_9.LCC9_CTRL_NIPBL_ChIP8_GB_rep1.bam/HI.1997.004.Index_9.LCC9_CTRL_NIPBL_ChIP8_GB_rep1.bam.sorted.filtered.bam alignment/LCC9_CTRL_NIPBL/LCC9_CTRL_NIPBL.merged.bam
symlink_readset_sample_bam.LCC9_CTRL_NIPBL.1ab1060995df79508cb6e310056bfbfe.mugqic.done
)
picard_merge_sam_files_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_9_JOB_ID: symlink_readset_sample_bam.LCC9_CTRL_POL2
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.LCC9_CTRL_POL2
JOB_DEPENDENCIES=$samtools_view_filter_9_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.LCC9_CTRL_POL2.03a2dbc51acc67b335233728bd306acf.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.LCC9_CTRL_POL2.03a2dbc51acc67b335233728bd306acf.mugqic.done'
mkdir -p alignment/LCC9_CTRL_POL2 && \
ln -s -f HI.1997.005.Index_20.LCC9_CTRL_POL2_ChIP8_GB_rep1.bam/HI.1997.005.Index_20.LCC9_CTRL_POL2_ChIP8_GB_rep1.bam.sorted.filtered.bam alignment/LCC9_CTRL_POL2/LCC9_CTRL_POL2.merged.bam
symlink_readset_sample_bam.LCC9_CTRL_POL2.03a2dbc51acc67b335233728bd306acf.mugqic.done
)
picard_merge_sam_files_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_10_JOB_ID: symlink_readset_sample_bam.LCC9_E2_POL2
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.LCC9_E2_POL2
JOB_DEPENDENCIES=$samtools_view_filter_10_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.LCC9_E2_POL2.e0b5fbb6733549fad997c2a2ffc61b48.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.LCC9_E2_POL2.e0b5fbb6733549fad997c2a2ffc61b48.mugqic.done'
mkdir -p alignment/LCC9_E2_POL2 && \
ln -s -f HI.1997.005.Index_22.LCC9_E2_POL2_ChIP8_GB_rep1.bam/HI.1997.005.Index_22.LCC9_E2_POL2_ChIP8_GB_rep1.bam.sorted.filtered.bam alignment/LCC9_E2_POL2/LCC9_E2_POL2.merged.bam
symlink_readset_sample_bam.LCC9_E2_POL2.e0b5fbb6733549fad997c2a2ffc61b48.mugqic.done
)
picard_merge_sam_files_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_11_JOB_ID: symlink_readset_sample_bam.LCC9_CTRL_ERA
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.LCC9_CTRL_ERA
JOB_DEPENDENCIES=$samtools_view_filter_11_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.LCC9_CTRL_ERA.360fcf289cfdee96167e7bbc7bd3ca65.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.LCC9_CTRL_ERA.360fcf289cfdee96167e7bbc7bd3ca65.mugqic.done'
mkdir -p alignment/LCC9_CTRL_ERA && \
ln -s -f HI.1997.005.Index_23.LCC9_CTRL_ERA_ChIP8_GB_rep1.bam/HI.1997.005.Index_23.LCC9_CTRL_ERA_ChIP8_GB_rep1.bam.sorted.filtered.bam alignment/LCC9_CTRL_ERA/LCC9_CTRL_ERA.merged.bam
symlink_readset_sample_bam.LCC9_CTRL_ERA.360fcf289cfdee96167e7bbc7bd3ca65.mugqic.done
)
picard_merge_sam_files_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_12_JOB_ID: symlink_readset_sample_bam.LCC9_E2_ERA
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.LCC9_E2_ERA
JOB_DEPENDENCIES=$samtools_view_filter_12_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.LCC9_E2_ERA.64cf905eefdd570858db086e164b92ee.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.LCC9_E2_ERA.64cf905eefdd570858db086e164b92ee.mugqic.done'
mkdir -p alignment/LCC9_E2_ERA && \
ln -s -f HI.1997.005.Index_27.LCC9_E2_ERA_ChIP8_GB_rep1.bam/HI.1997.005.Index_27.LCC9_E2_ERA_ChIP8_GB_rep1.bam.sorted.filtered.bam alignment/LCC9_E2_ERA/LCC9_E2_ERA.merged.bam
symlink_readset_sample_bam.LCC9_E2_ERA.64cf905eefdd570858db086e164b92ee.mugqic.done
)
picard_merge_sam_files_12_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: picard_mark_duplicates
#-------------------------------------------------------------------------------
STEP=picard_mark_duplicates
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_1_JOB_ID: picard_mark_duplicates.LCC9_CTRL_MED1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.LCC9_CTRL_MED1
JOB_DEPENDENCIES=$picard_merge_sam_files_1_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.LCC9_CTRL_MED1.b3ebe92532fd49cc6276f2b22c1c3161.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.LCC9_CTRL_MED1.b3ebe92532fd49cc6276f2b22c1c3161.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/LCC9_CTRL_MED1/LCC9_CTRL_MED1.merged.bam \
  OUTPUT=alignment/LCC9_CTRL_MED1/LCC9_CTRL_MED1.sorted.dup.bam \
  METRICS_FILE=alignment/LCC9_CTRL_MED1/LCC9_CTRL_MED1.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.LCC9_CTRL_MED1.b3ebe92532fd49cc6276f2b22c1c3161.mugqic.done
)
picard_mark_duplicates_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_2_JOB_ID: picard_mark_duplicates.LCC9_E2_WCE
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.LCC9_E2_WCE
JOB_DEPENDENCIES=$picard_merge_sam_files_2_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.LCC9_E2_WCE.91a8fc98acd7607c9f42990a93df8e17.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.LCC9_E2_WCE.91a8fc98acd7607c9f42990a93df8e17.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/LCC9_E2_WCE/LCC9_E2_WCE.merged.bam \
  OUTPUT=alignment/LCC9_E2_WCE/LCC9_E2_WCE.sorted.dup.bam \
  METRICS_FILE=alignment/LCC9_E2_WCE/LCC9_E2_WCE.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.LCC9_E2_WCE.91a8fc98acd7607c9f42990a93df8e17.mugqic.done
)
picard_mark_duplicates_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_3_JOB_ID: picard_mark_duplicates.LCC9_E2_MED1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.LCC9_E2_MED1
JOB_DEPENDENCIES=$picard_merge_sam_files_3_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.LCC9_E2_MED1.a88ce9d7bf4a435ddf9aa462ca5f8aa2.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.LCC9_E2_MED1.a88ce9d7bf4a435ddf9aa462ca5f8aa2.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/LCC9_E2_MED1/LCC9_E2_MED1.merged.bam \
  OUTPUT=alignment/LCC9_E2_MED1/LCC9_E2_MED1.sorted.dup.bam \
  METRICS_FILE=alignment/LCC9_E2_MED1/LCC9_E2_MED1.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.LCC9_E2_MED1.a88ce9d7bf4a435ddf9aa462ca5f8aa2.mugqic.done
)
picard_mark_duplicates_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_4_JOB_ID: picard_mark_duplicates.LCC9_CTRL_WCE
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.LCC9_CTRL_WCE
JOB_DEPENDENCIES=$picard_merge_sam_files_4_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.LCC9_CTRL_WCE.8e480d0f217a3feee0d5957bc55f7c85.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.LCC9_CTRL_WCE.8e480d0f217a3feee0d5957bc55f7c85.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/LCC9_CTRL_WCE/LCC9_CTRL_WCE.merged.bam \
  OUTPUT=alignment/LCC9_CTRL_WCE/LCC9_CTRL_WCE.sorted.dup.bam \
  METRICS_FILE=alignment/LCC9_CTRL_WCE/LCC9_CTRL_WCE.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.LCC9_CTRL_WCE.8e480d0f217a3feee0d5957bc55f7c85.mugqic.done
)
picard_mark_duplicates_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_5_JOB_ID: picard_mark_duplicates.LCC9_CTRL_SMC1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.LCC9_CTRL_SMC1
JOB_DEPENDENCIES=$picard_merge_sam_files_5_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.LCC9_CTRL_SMC1.e0b265c507689f18cac0b382dbd80d8f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.LCC9_CTRL_SMC1.e0b265c507689f18cac0b382dbd80d8f.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/LCC9_CTRL_SMC1/LCC9_CTRL_SMC1.merged.bam \
  OUTPUT=alignment/LCC9_CTRL_SMC1/LCC9_CTRL_SMC1.sorted.dup.bam \
  METRICS_FILE=alignment/LCC9_CTRL_SMC1/LCC9_CTRL_SMC1.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.LCC9_CTRL_SMC1.e0b265c507689f18cac0b382dbd80d8f.mugqic.done
)
picard_mark_duplicates_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_6_JOB_ID: picard_mark_duplicates.LCC9_E2_SMC1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.LCC9_E2_SMC1
JOB_DEPENDENCIES=$picard_merge_sam_files_6_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.LCC9_E2_SMC1.5079854d22e3487009288123010003cc.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.LCC9_E2_SMC1.5079854d22e3487009288123010003cc.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/LCC9_E2_SMC1/LCC9_E2_SMC1.merged.bam \
  OUTPUT=alignment/LCC9_E2_SMC1/LCC9_E2_SMC1.sorted.dup.bam \
  METRICS_FILE=alignment/LCC9_E2_SMC1/LCC9_E2_SMC1.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.LCC9_E2_SMC1.5079854d22e3487009288123010003cc.mugqic.done
)
picard_mark_duplicates_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_7_JOB_ID: picard_mark_duplicates.LCC9_E2_NIPBL
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.LCC9_E2_NIPBL
JOB_DEPENDENCIES=$picard_merge_sam_files_7_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.LCC9_E2_NIPBL.059d0e9711b3f25cbf3fbd871c1004ed.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.LCC9_E2_NIPBL.059d0e9711b3f25cbf3fbd871c1004ed.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/LCC9_E2_NIPBL/LCC9_E2_NIPBL.merged.bam \
  OUTPUT=alignment/LCC9_E2_NIPBL/LCC9_E2_NIPBL.sorted.dup.bam \
  METRICS_FILE=alignment/LCC9_E2_NIPBL/LCC9_E2_NIPBL.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.LCC9_E2_NIPBL.059d0e9711b3f25cbf3fbd871c1004ed.mugqic.done
)
picard_mark_duplicates_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_8_JOB_ID: picard_mark_duplicates.LCC9_CTRL_NIPBL
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.LCC9_CTRL_NIPBL
JOB_DEPENDENCIES=$picard_merge_sam_files_8_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.LCC9_CTRL_NIPBL.3fb3a711b5c31725d0996dfa50eddf42.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.LCC9_CTRL_NIPBL.3fb3a711b5c31725d0996dfa50eddf42.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/LCC9_CTRL_NIPBL/LCC9_CTRL_NIPBL.merged.bam \
  OUTPUT=alignment/LCC9_CTRL_NIPBL/LCC9_CTRL_NIPBL.sorted.dup.bam \
  METRICS_FILE=alignment/LCC9_CTRL_NIPBL/LCC9_CTRL_NIPBL.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.LCC9_CTRL_NIPBL.3fb3a711b5c31725d0996dfa50eddf42.mugqic.done
)
picard_mark_duplicates_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_9_JOB_ID: picard_mark_duplicates.LCC9_CTRL_POL2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.LCC9_CTRL_POL2
JOB_DEPENDENCIES=$picard_merge_sam_files_9_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.LCC9_CTRL_POL2.b1ca12fcb54398055be9b7a07a009b2a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.LCC9_CTRL_POL2.b1ca12fcb54398055be9b7a07a009b2a.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/LCC9_CTRL_POL2/LCC9_CTRL_POL2.merged.bam \
  OUTPUT=alignment/LCC9_CTRL_POL2/LCC9_CTRL_POL2.sorted.dup.bam \
  METRICS_FILE=alignment/LCC9_CTRL_POL2/LCC9_CTRL_POL2.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.LCC9_CTRL_POL2.b1ca12fcb54398055be9b7a07a009b2a.mugqic.done
)
picard_mark_duplicates_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_10_JOB_ID: picard_mark_duplicates.LCC9_E2_POL2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.LCC9_E2_POL2
JOB_DEPENDENCIES=$picard_merge_sam_files_10_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.LCC9_E2_POL2.ce57ae4ec180ec9d97ab02f52343fc3b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.LCC9_E2_POL2.ce57ae4ec180ec9d97ab02f52343fc3b.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/LCC9_E2_POL2/LCC9_E2_POL2.merged.bam \
  OUTPUT=alignment/LCC9_E2_POL2/LCC9_E2_POL2.sorted.dup.bam \
  METRICS_FILE=alignment/LCC9_E2_POL2/LCC9_E2_POL2.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.LCC9_E2_POL2.ce57ae4ec180ec9d97ab02f52343fc3b.mugqic.done
)
picard_mark_duplicates_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_11_JOB_ID: picard_mark_duplicates.LCC9_CTRL_ERA
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.LCC9_CTRL_ERA
JOB_DEPENDENCIES=$picard_merge_sam_files_11_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.LCC9_CTRL_ERA.d0eb9c84d624b54d1d0d8b2533077e21.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.LCC9_CTRL_ERA.d0eb9c84d624b54d1d0d8b2533077e21.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/LCC9_CTRL_ERA/LCC9_CTRL_ERA.merged.bam \
  OUTPUT=alignment/LCC9_CTRL_ERA/LCC9_CTRL_ERA.sorted.dup.bam \
  METRICS_FILE=alignment/LCC9_CTRL_ERA/LCC9_CTRL_ERA.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.LCC9_CTRL_ERA.d0eb9c84d624b54d1d0d8b2533077e21.mugqic.done
)
picard_mark_duplicates_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_12_JOB_ID: picard_mark_duplicates.LCC9_E2_ERA
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.LCC9_E2_ERA
JOB_DEPENDENCIES=$picard_merge_sam_files_12_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.LCC9_E2_ERA.adb9cdd3eddd278ac46a74b49b0e0a42.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.LCC9_E2_ERA.adb9cdd3eddd278ac46a74b49b0e0a42.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/LCC9_E2_ERA/LCC9_E2_ERA.merged.bam \
  OUTPUT=alignment/LCC9_E2_ERA/LCC9_E2_ERA.sorted.dup.bam \
  METRICS_FILE=alignment/LCC9_E2_ERA/LCC9_E2_ERA.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.LCC9_E2_ERA.adb9cdd3eddd278ac46a74b49b0e0a42.mugqic.done
)
picard_mark_duplicates_12_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_13_JOB_ID: picard_mark_duplicates_report
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates_report
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_2_JOB_ID:$picard_mark_duplicates_3_JOB_ID:$picard_mark_duplicates_4_JOB_ID:$picard_mark_duplicates_5_JOB_ID:$picard_mark_duplicates_6_JOB_ID:$picard_mark_duplicates_7_JOB_ID:$picard_mark_duplicates_8_JOB_ID:$picard_mark_duplicates_9_JOB_ID:$picard_mark_duplicates_10_JOB_ID:$picard_mark_duplicates_11_JOB_ID:$picard_mark_duplicates_12_JOB_ID
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
picard_mark_duplicates_13_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_13_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: metrics
#-------------------------------------------------------------------------------
STEP=metrics
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: metrics_1_JOB_ID: metrics.flagstat
#-------------------------------------------------------------------------------
JOB_NAME=metrics.flagstat
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_2_JOB_ID:$picard_mark_duplicates_3_JOB_ID:$picard_mark_duplicates_4_JOB_ID:$picard_mark_duplicates_5_JOB_ID:$picard_mark_duplicates_6_JOB_ID:$picard_mark_duplicates_7_JOB_ID:$picard_mark_duplicates_8_JOB_ID:$picard_mark_duplicates_9_JOB_ID:$picard_mark_duplicates_10_JOB_ID:$picard_mark_duplicates_11_JOB_ID:$picard_mark_duplicates_12_JOB_ID
JOB_DONE=job_output/metrics/metrics.flagstat.1d563895d90c1052e4a9dda21e950b22.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'metrics.flagstat.1d563895d90c1052e4a9dda21e950b22.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools flagstat \
  alignment/LCC9_CTRL_MED1/LCC9_CTRL_MED1.sorted.dup.bam \
  > alignment/LCC9_CTRL_MED1/LCC9_CTRL_MED1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/LCC9_E2_WCE/LCC9_E2_WCE.sorted.dup.bam \
  > alignment/LCC9_E2_WCE/LCC9_E2_WCE.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/LCC9_E2_MED1/LCC9_E2_MED1.sorted.dup.bam \
  > alignment/LCC9_E2_MED1/LCC9_E2_MED1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/LCC9_CTRL_WCE/LCC9_CTRL_WCE.sorted.dup.bam \
  > alignment/LCC9_CTRL_WCE/LCC9_CTRL_WCE.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/LCC9_CTRL_SMC1/LCC9_CTRL_SMC1.sorted.dup.bam \
  > alignment/LCC9_CTRL_SMC1/LCC9_CTRL_SMC1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/LCC9_E2_SMC1/LCC9_E2_SMC1.sorted.dup.bam \
  > alignment/LCC9_E2_SMC1/LCC9_E2_SMC1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/LCC9_E2_NIPBL/LCC9_E2_NIPBL.sorted.dup.bam \
  > alignment/LCC9_E2_NIPBL/LCC9_E2_NIPBL.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/LCC9_CTRL_NIPBL/LCC9_CTRL_NIPBL.sorted.dup.bam \
  > alignment/LCC9_CTRL_NIPBL/LCC9_CTRL_NIPBL.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/LCC9_CTRL_POL2/LCC9_CTRL_POL2.sorted.dup.bam \
  > alignment/LCC9_CTRL_POL2/LCC9_CTRL_POL2.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/LCC9_E2_POL2/LCC9_E2_POL2.sorted.dup.bam \
  > alignment/LCC9_E2_POL2/LCC9_E2_POL2.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/LCC9_CTRL_ERA/LCC9_CTRL_ERA.sorted.dup.bam \
  > alignment/LCC9_CTRL_ERA/LCC9_CTRL_ERA.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/LCC9_E2_ERA/LCC9_E2_ERA.sorted.dup.bam \
  > alignment/LCC9_E2_ERA/LCC9_E2_ERA.sorted.dup.bam.flagstat
metrics.flagstat.1d563895d90c1052e4a9dda21e950b22.mugqic.done
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
JOB_DONE=job_output/metrics/metrics_report.6605e351a453683dfe553bcf66706675.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'metrics_report.6605e351a453683dfe553bcf66706675.mugqic.done'
module load mugqic/pandoc/1.15.2 && \
for sample in LCC9_CTRL_MED1 LCC9_E2_WCE LCC9_E2_MED1 LCC9_CTRL_WCE LCC9_CTRL_SMC1 LCC9_E2_SMC1 LCC9_E2_NIPBL LCC9_CTRL_NIPBL LCC9_CTRL_POL2 LCC9_E2_POL2 LCC9_CTRL_ERA LCC9_E2_ERA
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

metrics_report.6605e351a453683dfe553bcf66706675.mugqic.done
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
# JOB: homer_make_tag_directory_1_JOB_ID: homer_make_tag_directory.LCC9_CTRL_MED1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.LCC9_CTRL_MED1
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.LCC9_CTRL_MED1.f6e8b5eea5a0f1ce9b94eca1a241127a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.LCC9_CTRL_MED1.f6e8b5eea5a0f1ce9b94eca1a241127a.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/LCC9_CTRL_MED1 \
  alignment/LCC9_CTRL_MED1/LCC9_CTRL_MED1.sorted.dup.bam \
  -checkGC -genome GRCh38
homer_make_tag_directory.LCC9_CTRL_MED1.f6e8b5eea5a0f1ce9b94eca1a241127a.mugqic.done
)
homer_make_tag_directory_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_2_JOB_ID: homer_make_tag_directory.LCC9_E2_WCE
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.LCC9_E2_WCE
JOB_DEPENDENCIES=$picard_mark_duplicates_2_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.LCC9_E2_WCE.42152f372a70de4f07d715cfea5b1c7f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.LCC9_E2_WCE.42152f372a70de4f07d715cfea5b1c7f.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/LCC9_E2_WCE \
  alignment/LCC9_E2_WCE/LCC9_E2_WCE.sorted.dup.bam \
  -checkGC -genome GRCh38
homer_make_tag_directory.LCC9_E2_WCE.42152f372a70de4f07d715cfea5b1c7f.mugqic.done
)
homer_make_tag_directory_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_3_JOB_ID: homer_make_tag_directory.LCC9_E2_MED1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.LCC9_E2_MED1
JOB_DEPENDENCIES=$picard_mark_duplicates_3_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.LCC9_E2_MED1.f169916b1d3dce9f70ad676a9c95ed8d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.LCC9_E2_MED1.f169916b1d3dce9f70ad676a9c95ed8d.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/LCC9_E2_MED1 \
  alignment/LCC9_E2_MED1/LCC9_E2_MED1.sorted.dup.bam \
  -checkGC -genome GRCh38
homer_make_tag_directory.LCC9_E2_MED1.f169916b1d3dce9f70ad676a9c95ed8d.mugqic.done
)
homer_make_tag_directory_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_4_JOB_ID: homer_make_tag_directory.LCC9_CTRL_WCE
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.LCC9_CTRL_WCE
JOB_DEPENDENCIES=$picard_mark_duplicates_4_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.LCC9_CTRL_WCE.795615fca2d06f031b68664401f9fc41.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.LCC9_CTRL_WCE.795615fca2d06f031b68664401f9fc41.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/LCC9_CTRL_WCE \
  alignment/LCC9_CTRL_WCE/LCC9_CTRL_WCE.sorted.dup.bam \
  -checkGC -genome GRCh38
homer_make_tag_directory.LCC9_CTRL_WCE.795615fca2d06f031b68664401f9fc41.mugqic.done
)
homer_make_tag_directory_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_5_JOB_ID: homer_make_tag_directory.LCC9_CTRL_SMC1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.LCC9_CTRL_SMC1
JOB_DEPENDENCIES=$picard_mark_duplicates_5_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.LCC9_CTRL_SMC1.7144aa854397d471dca36f6a045cc76b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.LCC9_CTRL_SMC1.7144aa854397d471dca36f6a045cc76b.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/LCC9_CTRL_SMC1 \
  alignment/LCC9_CTRL_SMC1/LCC9_CTRL_SMC1.sorted.dup.bam \
  -checkGC -genome GRCh38
homer_make_tag_directory.LCC9_CTRL_SMC1.7144aa854397d471dca36f6a045cc76b.mugqic.done
)
homer_make_tag_directory_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_6_JOB_ID: homer_make_tag_directory.LCC9_E2_SMC1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.LCC9_E2_SMC1
JOB_DEPENDENCIES=$picard_mark_duplicates_6_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.LCC9_E2_SMC1.13a0607cefaa72ed50d4ac59245f1eb7.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.LCC9_E2_SMC1.13a0607cefaa72ed50d4ac59245f1eb7.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/LCC9_E2_SMC1 \
  alignment/LCC9_E2_SMC1/LCC9_E2_SMC1.sorted.dup.bam \
  -checkGC -genome GRCh38
homer_make_tag_directory.LCC9_E2_SMC1.13a0607cefaa72ed50d4ac59245f1eb7.mugqic.done
)
homer_make_tag_directory_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_7_JOB_ID: homer_make_tag_directory.LCC9_E2_NIPBL
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.LCC9_E2_NIPBL
JOB_DEPENDENCIES=$picard_mark_duplicates_7_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.LCC9_E2_NIPBL.339223377161c0ff6b3d57536c356a9e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.LCC9_E2_NIPBL.339223377161c0ff6b3d57536c356a9e.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/LCC9_E2_NIPBL \
  alignment/LCC9_E2_NIPBL/LCC9_E2_NIPBL.sorted.dup.bam \
  -checkGC -genome GRCh38
homer_make_tag_directory.LCC9_E2_NIPBL.339223377161c0ff6b3d57536c356a9e.mugqic.done
)
homer_make_tag_directory_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_8_JOB_ID: homer_make_tag_directory.LCC9_CTRL_NIPBL
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.LCC9_CTRL_NIPBL
JOB_DEPENDENCIES=$picard_mark_duplicates_8_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.LCC9_CTRL_NIPBL.4d0925ac0af2e2517794ab2b82a78d00.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.LCC9_CTRL_NIPBL.4d0925ac0af2e2517794ab2b82a78d00.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/LCC9_CTRL_NIPBL \
  alignment/LCC9_CTRL_NIPBL/LCC9_CTRL_NIPBL.sorted.dup.bam \
  -checkGC -genome GRCh38
homer_make_tag_directory.LCC9_CTRL_NIPBL.4d0925ac0af2e2517794ab2b82a78d00.mugqic.done
)
homer_make_tag_directory_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_9_JOB_ID: homer_make_tag_directory.LCC9_CTRL_POL2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.LCC9_CTRL_POL2
JOB_DEPENDENCIES=$picard_mark_duplicates_9_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.LCC9_CTRL_POL2.4450c821494cc6461d133495e3ab8e3c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.LCC9_CTRL_POL2.4450c821494cc6461d133495e3ab8e3c.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/LCC9_CTRL_POL2 \
  alignment/LCC9_CTRL_POL2/LCC9_CTRL_POL2.sorted.dup.bam \
  -checkGC -genome GRCh38
homer_make_tag_directory.LCC9_CTRL_POL2.4450c821494cc6461d133495e3ab8e3c.mugqic.done
)
homer_make_tag_directory_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_10_JOB_ID: homer_make_tag_directory.LCC9_E2_POL2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.LCC9_E2_POL2
JOB_DEPENDENCIES=$picard_mark_duplicates_10_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.LCC9_E2_POL2.ae4f7117d61123e4a43ee38f83d96d16.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.LCC9_E2_POL2.ae4f7117d61123e4a43ee38f83d96d16.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/LCC9_E2_POL2 \
  alignment/LCC9_E2_POL2/LCC9_E2_POL2.sorted.dup.bam \
  -checkGC -genome GRCh38
homer_make_tag_directory.LCC9_E2_POL2.ae4f7117d61123e4a43ee38f83d96d16.mugqic.done
)
homer_make_tag_directory_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_11_JOB_ID: homer_make_tag_directory.LCC9_CTRL_ERA
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.LCC9_CTRL_ERA
JOB_DEPENDENCIES=$picard_mark_duplicates_11_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.LCC9_CTRL_ERA.5de251d1fe3f4fc7930309b4d5591ec9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.LCC9_CTRL_ERA.5de251d1fe3f4fc7930309b4d5591ec9.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/LCC9_CTRL_ERA \
  alignment/LCC9_CTRL_ERA/LCC9_CTRL_ERA.sorted.dup.bam \
  -checkGC -genome GRCh38
homer_make_tag_directory.LCC9_CTRL_ERA.5de251d1fe3f4fc7930309b4d5591ec9.mugqic.done
)
homer_make_tag_directory_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_12_JOB_ID: homer_make_tag_directory.LCC9_E2_ERA
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.LCC9_E2_ERA
JOB_DEPENDENCIES=$picard_mark_duplicates_12_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.LCC9_E2_ERA.1f27794689a35304f9202600f030d248.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.LCC9_E2_ERA.1f27794689a35304f9202600f030d248.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/LCC9_E2_ERA \
  alignment/LCC9_E2_ERA/LCC9_E2_ERA.sorted.dup.bam \
  -checkGC -genome GRCh38
homer_make_tag_directory.LCC9_E2_ERA.1f27794689a35304f9202600f030d248.mugqic.done
)
homer_make_tag_directory_12_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: qc_metrics
#-------------------------------------------------------------------------------
STEP=qc_metrics
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: qc_metrics_1_JOB_ID: qc_plots_R
#-------------------------------------------------------------------------------
JOB_NAME=qc_plots_R
JOB_DEPENDENCIES=$homer_make_tag_directory_1_JOB_ID:$homer_make_tag_directory_2_JOB_ID:$homer_make_tag_directory_3_JOB_ID:$homer_make_tag_directory_4_JOB_ID:$homer_make_tag_directory_5_JOB_ID:$homer_make_tag_directory_6_JOB_ID:$homer_make_tag_directory_7_JOB_ID:$homer_make_tag_directory_8_JOB_ID:$homer_make_tag_directory_9_JOB_ID:$homer_make_tag_directory_10_JOB_ID:$homer_make_tag_directory_11_JOB_ID:$homer_make_tag_directory_12_JOB_ID
JOB_DONE=job_output/qc_metrics/qc_plots_R.f9be7274dc2664202a867f6cfe8005bd.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'qc_plots_R.f9be7274dc2664202a867f6cfe8005bd.mugqic.done'
module load mugqic/mugqic_tools/2.1.5 mugqic/R_Bioconductor/3.2.3_3.2 && \
mkdir -p graphs && \
Rscript $R_TOOLS/chipSeqGenerateQCMetrics.R \
  ../../raw/chip-seq/design.txt \
  /gs/scratch/efournier/sb_cofactor_hr/LCC9/output/chip-pipeline-GRCh38 && \
cp /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/ChipSeq.qc_metrics.md report/ChipSeq.qc_metrics.md && \
for sample in LCC9_CTRL_MED1 LCC9_E2_WCE LCC9_E2_MED1 LCC9_CTRL_WCE LCC9_CTRL_SMC1 LCC9_E2_SMC1 LCC9_E2_NIPBL LCC9_CTRL_NIPBL LCC9_CTRL_POL2 LCC9_E2_POL2 LCC9_CTRL_ERA LCC9_E2_ERA
do
  cp --parents graphs/${sample}_QC_Metrics.ps report/
  convert -rotate 90 graphs/${sample}_QC_Metrics.ps report/graphs/${sample}_QC_Metrics.png
  echo -e "----

![QC Metrics for Sample $sample ([download high-res image](graphs/${sample}_QC_Metrics.ps))](graphs/${sample}_QC_Metrics.png)
" \
  >> report/ChipSeq.qc_metrics.md
done
qc_plots_R.f9be7274dc2664202a867f6cfe8005bd.mugqic.done
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
# JOB: homer_make_ucsc_file_1_JOB_ID: homer_make_ucsc_file.LCC9_CTRL_MED1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.LCC9_CTRL_MED1
JOB_DEPENDENCIES=$homer_make_tag_directory_1_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.LCC9_CTRL_MED1.b45fa96c06b27bbd5fa0cc75942adf56.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.LCC9_CTRL_MED1.b45fa96c06b27bbd5fa0cc75942adf56.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/LCC9_CTRL_MED1 && \
makeUCSCfile \
  tags/LCC9_CTRL_MED1 | \
gzip -1 -c > tracks/LCC9_CTRL_MED1/LCC9_CTRL_MED1.ucsc.bedGraph.gz
homer_make_ucsc_file.LCC9_CTRL_MED1.b45fa96c06b27bbd5fa0cc75942adf56.mugqic.done
)
homer_make_ucsc_file_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_2_JOB_ID: homer_make_ucsc_file.LCC9_E2_WCE
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.LCC9_E2_WCE
JOB_DEPENDENCIES=$homer_make_tag_directory_2_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.LCC9_E2_WCE.35851f830bcf45d0da231bb39191c19b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.LCC9_E2_WCE.35851f830bcf45d0da231bb39191c19b.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/LCC9_E2_WCE && \
makeUCSCfile \
  tags/LCC9_E2_WCE | \
gzip -1 -c > tracks/LCC9_E2_WCE/LCC9_E2_WCE.ucsc.bedGraph.gz
homer_make_ucsc_file.LCC9_E2_WCE.35851f830bcf45d0da231bb39191c19b.mugqic.done
)
homer_make_ucsc_file_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_3_JOB_ID: homer_make_ucsc_file.LCC9_E2_MED1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.LCC9_E2_MED1
JOB_DEPENDENCIES=$homer_make_tag_directory_3_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.LCC9_E2_MED1.8925f2488aa591b69593181c4ddabcb4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.LCC9_E2_MED1.8925f2488aa591b69593181c4ddabcb4.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/LCC9_E2_MED1 && \
makeUCSCfile \
  tags/LCC9_E2_MED1 | \
gzip -1 -c > tracks/LCC9_E2_MED1/LCC9_E2_MED1.ucsc.bedGraph.gz
homer_make_ucsc_file.LCC9_E2_MED1.8925f2488aa591b69593181c4ddabcb4.mugqic.done
)
homer_make_ucsc_file_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_4_JOB_ID: homer_make_ucsc_file.LCC9_CTRL_WCE
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.LCC9_CTRL_WCE
JOB_DEPENDENCIES=$homer_make_tag_directory_4_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.LCC9_CTRL_WCE.cc9538c718aedf4e8cc1b5887d542d32.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.LCC9_CTRL_WCE.cc9538c718aedf4e8cc1b5887d542d32.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/LCC9_CTRL_WCE && \
makeUCSCfile \
  tags/LCC9_CTRL_WCE | \
gzip -1 -c > tracks/LCC9_CTRL_WCE/LCC9_CTRL_WCE.ucsc.bedGraph.gz
homer_make_ucsc_file.LCC9_CTRL_WCE.cc9538c718aedf4e8cc1b5887d542d32.mugqic.done
)
homer_make_ucsc_file_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_5_JOB_ID: homer_make_ucsc_file.LCC9_CTRL_SMC1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.LCC9_CTRL_SMC1
JOB_DEPENDENCIES=$homer_make_tag_directory_5_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.LCC9_CTRL_SMC1.d705eb347996bee5f096eaae97690abd.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.LCC9_CTRL_SMC1.d705eb347996bee5f096eaae97690abd.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/LCC9_CTRL_SMC1 && \
makeUCSCfile \
  tags/LCC9_CTRL_SMC1 | \
gzip -1 -c > tracks/LCC9_CTRL_SMC1/LCC9_CTRL_SMC1.ucsc.bedGraph.gz
homer_make_ucsc_file.LCC9_CTRL_SMC1.d705eb347996bee5f096eaae97690abd.mugqic.done
)
homer_make_ucsc_file_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_6_JOB_ID: homer_make_ucsc_file.LCC9_E2_SMC1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.LCC9_E2_SMC1
JOB_DEPENDENCIES=$homer_make_tag_directory_6_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.LCC9_E2_SMC1.1fac41a9ec39b88653f1ece5b72f1e58.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.LCC9_E2_SMC1.1fac41a9ec39b88653f1ece5b72f1e58.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/LCC9_E2_SMC1 && \
makeUCSCfile \
  tags/LCC9_E2_SMC1 | \
gzip -1 -c > tracks/LCC9_E2_SMC1/LCC9_E2_SMC1.ucsc.bedGraph.gz
homer_make_ucsc_file.LCC9_E2_SMC1.1fac41a9ec39b88653f1ece5b72f1e58.mugqic.done
)
homer_make_ucsc_file_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_7_JOB_ID: homer_make_ucsc_file.LCC9_E2_NIPBL
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.LCC9_E2_NIPBL
JOB_DEPENDENCIES=$homer_make_tag_directory_7_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.LCC9_E2_NIPBL.b7576d83c570d794e7a772b4bf2801d0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.LCC9_E2_NIPBL.b7576d83c570d794e7a772b4bf2801d0.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/LCC9_E2_NIPBL && \
makeUCSCfile \
  tags/LCC9_E2_NIPBL | \
gzip -1 -c > tracks/LCC9_E2_NIPBL/LCC9_E2_NIPBL.ucsc.bedGraph.gz
homer_make_ucsc_file.LCC9_E2_NIPBL.b7576d83c570d794e7a772b4bf2801d0.mugqic.done
)
homer_make_ucsc_file_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_8_JOB_ID: homer_make_ucsc_file.LCC9_CTRL_NIPBL
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.LCC9_CTRL_NIPBL
JOB_DEPENDENCIES=$homer_make_tag_directory_8_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.LCC9_CTRL_NIPBL.221b559a960b0712c82c1125e4720865.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.LCC9_CTRL_NIPBL.221b559a960b0712c82c1125e4720865.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/LCC9_CTRL_NIPBL && \
makeUCSCfile \
  tags/LCC9_CTRL_NIPBL | \
gzip -1 -c > tracks/LCC9_CTRL_NIPBL/LCC9_CTRL_NIPBL.ucsc.bedGraph.gz
homer_make_ucsc_file.LCC9_CTRL_NIPBL.221b559a960b0712c82c1125e4720865.mugqic.done
)
homer_make_ucsc_file_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_9_JOB_ID: homer_make_ucsc_file.LCC9_CTRL_POL2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.LCC9_CTRL_POL2
JOB_DEPENDENCIES=$homer_make_tag_directory_9_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.LCC9_CTRL_POL2.edeb13156473397124162682ac03ed9e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.LCC9_CTRL_POL2.edeb13156473397124162682ac03ed9e.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/LCC9_CTRL_POL2 && \
makeUCSCfile \
  tags/LCC9_CTRL_POL2 | \
gzip -1 -c > tracks/LCC9_CTRL_POL2/LCC9_CTRL_POL2.ucsc.bedGraph.gz
homer_make_ucsc_file.LCC9_CTRL_POL2.edeb13156473397124162682ac03ed9e.mugqic.done
)
homer_make_ucsc_file_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_10_JOB_ID: homer_make_ucsc_file.LCC9_E2_POL2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.LCC9_E2_POL2
JOB_DEPENDENCIES=$homer_make_tag_directory_10_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.LCC9_E2_POL2.bbae975e3dfb0ae40deac68146a745f6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.LCC9_E2_POL2.bbae975e3dfb0ae40deac68146a745f6.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/LCC9_E2_POL2 && \
makeUCSCfile \
  tags/LCC9_E2_POL2 | \
gzip -1 -c > tracks/LCC9_E2_POL2/LCC9_E2_POL2.ucsc.bedGraph.gz
homer_make_ucsc_file.LCC9_E2_POL2.bbae975e3dfb0ae40deac68146a745f6.mugqic.done
)
homer_make_ucsc_file_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_11_JOB_ID: homer_make_ucsc_file.LCC9_CTRL_ERA
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.LCC9_CTRL_ERA
JOB_DEPENDENCIES=$homer_make_tag_directory_11_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.LCC9_CTRL_ERA.6e5ea8a2b6e50596ee9eb290ca58bbe0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.LCC9_CTRL_ERA.6e5ea8a2b6e50596ee9eb290ca58bbe0.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/LCC9_CTRL_ERA && \
makeUCSCfile \
  tags/LCC9_CTRL_ERA | \
gzip -1 -c > tracks/LCC9_CTRL_ERA/LCC9_CTRL_ERA.ucsc.bedGraph.gz
homer_make_ucsc_file.LCC9_CTRL_ERA.6e5ea8a2b6e50596ee9eb290ca58bbe0.mugqic.done
)
homer_make_ucsc_file_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_12_JOB_ID: homer_make_ucsc_file.LCC9_E2_ERA
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.LCC9_E2_ERA
JOB_DEPENDENCIES=$homer_make_tag_directory_12_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.LCC9_E2_ERA.2b8caa789c9711c93bcaec1a1db8ba07.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.LCC9_E2_ERA.2b8caa789c9711c93bcaec1a1db8ba07.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/LCC9_E2_ERA && \
makeUCSCfile \
  tags/LCC9_E2_ERA | \
gzip -1 -c > tracks/LCC9_E2_ERA/LCC9_E2_ERA.ucsc.bedGraph.gz
homer_make_ucsc_file.LCC9_E2_ERA.2b8caa789c9711c93bcaec1a1db8ba07.mugqic.done
)
homer_make_ucsc_file_12_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_13_JOB_ID: homer_make_ucsc_file_report
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file_report
JOB_DEPENDENCIES=$homer_make_ucsc_file_12_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file_report.c9cac04b3bc71f34cbcde042cf5b419d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file_report.c9cac04b3bc71f34cbcde042cf5b419d.mugqic.done'
mkdir -p report && \
zip -r report/tracks.zip tracks/*/*.ucsc.bedGraph.gz && \
cp /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/ChipSeq.homer_make_ucsc_file.md report/
homer_make_ucsc_file_report.c9cac04b3bc71f34cbcde042cf5b419d.mugqic.done
)
homer_make_ucsc_file_13_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_13_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: macs2_callpeak
#-------------------------------------------------------------------------------
STEP=macs2_callpeak
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_1_JOB_ID: macs2_callpeak.LCC9_CTRL_MED1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.LCC9_CTRL_MED1
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_4_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.LCC9_CTRL_MED1.88c0ac1b6b3dd43c0f85938fed4a420c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.LCC9_CTRL_MED1.88c0ac1b6b3dd43c0f85938fed4a420c.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/LCC9_CTRL_MED1 && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2479938032.8 \
  --treatment \
  alignment/LCC9_CTRL_WCE/LCC9_CTRL_WCE.sorted.dup.bam \
  --control \
  alignment/LCC9_CTRL_MED1/LCC9_CTRL_MED1.sorted.dup.bam \
  --name peak_call/LCC9_CTRL_MED1/LCC9_CTRL_MED1 \
  >& peak_call/LCC9_CTRL_MED1/LCC9_CTRL_MED1.diag.macs.out
macs2_callpeak.LCC9_CTRL_MED1.88c0ac1b6b3dd43c0f85938fed4a420c.mugqic.done
)
macs2_callpeak_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_2_JOB_ID: macs2_callpeak.LCC9_CTRL_NIPBL
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.LCC9_CTRL_NIPBL
JOB_DEPENDENCIES=$picard_mark_duplicates_4_JOB_ID:$picard_mark_duplicates_8_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.LCC9_CTRL_NIPBL.5422cc78ae5b70cd30b1aa0a77cab91d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.LCC9_CTRL_NIPBL.5422cc78ae5b70cd30b1aa0a77cab91d.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/LCC9_CTRL_NIPBL && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2479938032.8 \
  --treatment \
  alignment/LCC9_CTRL_WCE/LCC9_CTRL_WCE.sorted.dup.bam \
  --control \
  alignment/LCC9_CTRL_NIPBL/LCC9_CTRL_NIPBL.sorted.dup.bam \
  --name peak_call/LCC9_CTRL_NIPBL/LCC9_CTRL_NIPBL \
  >& peak_call/LCC9_CTRL_NIPBL/LCC9_CTRL_NIPBL.diag.macs.out
macs2_callpeak.LCC9_CTRL_NIPBL.5422cc78ae5b70cd30b1aa0a77cab91d.mugqic.done
)
macs2_callpeak_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_3_JOB_ID: macs2_callpeak.LCC9_CTRL_SMC1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.LCC9_CTRL_SMC1
JOB_DEPENDENCIES=$picard_mark_duplicates_4_JOB_ID:$picard_mark_duplicates_5_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.LCC9_CTRL_SMC1.eeaa970118e2b4a5dd73abe8995cd524.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.LCC9_CTRL_SMC1.eeaa970118e2b4a5dd73abe8995cd524.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/LCC9_CTRL_SMC1 && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2479938032.8 \
  --treatment \
  alignment/LCC9_CTRL_WCE/LCC9_CTRL_WCE.sorted.dup.bam \
  --control \
  alignment/LCC9_CTRL_SMC1/LCC9_CTRL_SMC1.sorted.dup.bam \
  --name peak_call/LCC9_CTRL_SMC1/LCC9_CTRL_SMC1 \
  >& peak_call/LCC9_CTRL_SMC1/LCC9_CTRL_SMC1.diag.macs.out
macs2_callpeak.LCC9_CTRL_SMC1.eeaa970118e2b4a5dd73abe8995cd524.mugqic.done
)
macs2_callpeak_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_4_JOB_ID: macs2_callpeak.LCC9_CTRL_POL2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.LCC9_CTRL_POL2
JOB_DEPENDENCIES=$picard_mark_duplicates_4_JOB_ID:$picard_mark_duplicates_9_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.LCC9_CTRL_POL2.4c047a759e9950f211baa9993016600a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.LCC9_CTRL_POL2.4c047a759e9950f211baa9993016600a.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/LCC9_CTRL_POL2 && \
macs2 callpeak --format BAM --broad --nomodel \
  --gsize 2479938032.8 \
  --treatment \
  alignment/LCC9_CTRL_WCE/LCC9_CTRL_WCE.sorted.dup.bam \
  --control \
  alignment/LCC9_CTRL_POL2/LCC9_CTRL_POL2.sorted.dup.bam \
  --name peak_call/LCC9_CTRL_POL2/LCC9_CTRL_POL2 \
  >& peak_call/LCC9_CTRL_POL2/LCC9_CTRL_POL2.diag.macs.out
macs2_callpeak.LCC9_CTRL_POL2.4c047a759e9950f211baa9993016600a.mugqic.done
)
macs2_callpeak_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_5_JOB_ID: macs2_callpeak.LCC9_CTRL_ERA
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.LCC9_CTRL_ERA
JOB_DEPENDENCIES=$picard_mark_duplicates_4_JOB_ID:$picard_mark_duplicates_11_JOB_ID:$picard_mark_duplicates_12_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.LCC9_CTRL_ERA.3d7bab8b0029ffcd2853a3b3ed040390.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.LCC9_CTRL_ERA.3d7bab8b0029ffcd2853a3b3ed040390.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/LCC9_CTRL_ERA && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2479938032.8 \
  --treatment \
  alignment/LCC9_CTRL_WCE/LCC9_CTRL_WCE.sorted.dup.bam \
  --control \
  alignment/LCC9_CTRL_ERA/LCC9_CTRL_ERA.sorted.dup.bam \
  alignment/LCC9_E2_ERA/LCC9_E2_ERA.sorted.dup.bam \
  --name peak_call/LCC9_CTRL_ERA/LCC9_CTRL_ERA \
  >& peak_call/LCC9_CTRL_ERA/LCC9_CTRL_ERA.diag.macs.out
macs2_callpeak.LCC9_CTRL_ERA.3d7bab8b0029ffcd2853a3b3ed040390.mugqic.done
)
macs2_callpeak_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_6_JOB_ID: macs2_callpeak.LCC9_E2_MED1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.LCC9_E2_MED1
JOB_DEPENDENCIES=$picard_mark_duplicates_2_JOB_ID:$picard_mark_duplicates_3_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.LCC9_E2_MED1.e8db3a9eadada0b4de257eef05cf7678.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.LCC9_E2_MED1.e8db3a9eadada0b4de257eef05cf7678.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/LCC9_E2_MED1 && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2479938032.8 \
  --treatment \
  alignment/LCC9_E2_WCE/LCC9_E2_WCE.sorted.dup.bam \
  --control \
  alignment/LCC9_E2_MED1/LCC9_E2_MED1.sorted.dup.bam \
  --name peak_call/LCC9_E2_MED1/LCC9_E2_MED1 \
  >& peak_call/LCC9_E2_MED1/LCC9_E2_MED1.diag.macs.out
macs2_callpeak.LCC9_E2_MED1.e8db3a9eadada0b4de257eef05cf7678.mugqic.done
)
macs2_callpeak_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_7_JOB_ID: macs2_callpeak.LCC9_E2_NIPBL
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.LCC9_E2_NIPBL
JOB_DEPENDENCIES=$picard_mark_duplicates_2_JOB_ID:$picard_mark_duplicates_7_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.LCC9_E2_NIPBL.90d15505d0266c6af4e3b5da0e4e635a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.LCC9_E2_NIPBL.90d15505d0266c6af4e3b5da0e4e635a.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/LCC9_E2_NIPBL && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2479938032.8 \
  --treatment \
  alignment/LCC9_E2_WCE/LCC9_E2_WCE.sorted.dup.bam \
  --control \
  alignment/LCC9_E2_NIPBL/LCC9_E2_NIPBL.sorted.dup.bam \
  --name peak_call/LCC9_E2_NIPBL/LCC9_E2_NIPBL \
  >& peak_call/LCC9_E2_NIPBL/LCC9_E2_NIPBL.diag.macs.out
macs2_callpeak.LCC9_E2_NIPBL.90d15505d0266c6af4e3b5da0e4e635a.mugqic.done
)
macs2_callpeak_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_8_JOB_ID: macs2_callpeak.LCC9_E2_SMC1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.LCC9_E2_SMC1
JOB_DEPENDENCIES=$picard_mark_duplicates_2_JOB_ID:$picard_mark_duplicates_6_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.LCC9_E2_SMC1.043bd803c2d47371b535408795138baf.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.LCC9_E2_SMC1.043bd803c2d47371b535408795138baf.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/LCC9_E2_SMC1 && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2479938032.8 \
  --treatment \
  alignment/LCC9_E2_WCE/LCC9_E2_WCE.sorted.dup.bam \
  --control \
  alignment/LCC9_E2_SMC1/LCC9_E2_SMC1.sorted.dup.bam \
  --name peak_call/LCC9_E2_SMC1/LCC9_E2_SMC1 \
  >& peak_call/LCC9_E2_SMC1/LCC9_E2_SMC1.diag.macs.out
macs2_callpeak.LCC9_E2_SMC1.043bd803c2d47371b535408795138baf.mugqic.done
)
macs2_callpeak_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_9_JOB_ID: macs2_callpeak.LCC9_E2_POL2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.LCC9_E2_POL2
JOB_DEPENDENCIES=$picard_mark_duplicates_2_JOB_ID:$picard_mark_duplicates_10_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.LCC9_E2_POL2.376d4c0b720c5cfe22e5f4e6e4a484d2.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.LCC9_E2_POL2.376d4c0b720c5cfe22e5f4e6e4a484d2.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/LCC9_E2_POL2 && \
macs2 callpeak --format BAM --broad --nomodel \
  --gsize 2479938032.8 \
  --treatment \
  alignment/LCC9_E2_WCE/LCC9_E2_WCE.sorted.dup.bam \
  --control \
  alignment/LCC9_E2_POL2/LCC9_E2_POL2.sorted.dup.bam \
  --name peak_call/LCC9_E2_POL2/LCC9_E2_POL2 \
  >& peak_call/LCC9_E2_POL2/LCC9_E2_POL2.diag.macs.out
macs2_callpeak.LCC9_E2_POL2.376d4c0b720c5cfe22e5f4e6e4a484d2.mugqic.done
)
macs2_callpeak_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_10_JOB_ID: macs2_callpeak.LCC9_E2_ERA
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.LCC9_E2_ERA
JOB_DEPENDENCIES=$picard_mark_duplicates_2_JOB_ID:$picard_mark_duplicates_12_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.LCC9_E2_ERA.f0fd0163f55a3c7325f70c6745d763b1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.LCC9_E2_ERA.f0fd0163f55a3c7325f70c6745d763b1.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/LCC9_E2_ERA && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2479938032.8 \
  --treatment \
  alignment/LCC9_E2_WCE/LCC9_E2_WCE.sorted.dup.bam \
  --control \
  alignment/LCC9_E2_ERA/LCC9_E2_ERA.sorted.dup.bam \
  --name peak_call/LCC9_E2_ERA/LCC9_E2_ERA \
  >& peak_call/LCC9_E2_ERA/LCC9_E2_ERA.diag.macs.out
macs2_callpeak.LCC9_E2_ERA.f0fd0163f55a3c7325f70c6745d763b1.mugqic.done
)
macs2_callpeak_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_11_JOB_ID: macs2_callpeak_report
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_report
JOB_DEPENDENCIES=$macs2_callpeak_1_JOB_ID:$macs2_callpeak_2_JOB_ID:$macs2_callpeak_3_JOB_ID:$macs2_callpeak_4_JOB_ID:$macs2_callpeak_5_JOB_ID:$macs2_callpeak_6_JOB_ID:$macs2_callpeak_7_JOB_ID:$macs2_callpeak_8_JOB_ID:$macs2_callpeak_9_JOB_ID:$macs2_callpeak_10_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_report.51dee2379e1a758f407c4d2e9210ac58.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_report.51dee2379e1a758f407c4d2e9210ac58.mugqic.done'
mkdir -p report && \
cp /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/ChipSeq.macs2_callpeak.md report/ && \
for contrast in LCC9_CTRL_MED1 LCC9_CTRL_NIPBL LCC9_CTRL_SMC1 LCC9_CTRL_POL2 LCC9_CTRL_ERA LCC9_E2_MED1 LCC9_E2_NIPBL LCC9_E2_SMC1 LCC9_E2_POL2 LCC9_E2_ERA
do
  cp -a --parents peak_call/$contrast/ report/ && \
  echo -e "* [Peak Calls File for Design $contrast](peak_call/$contrast/${contrast}_peaks.xls)" \
  >> report/ChipSeq.macs2_callpeak.md
done
macs2_callpeak_report.51dee2379e1a758f407c4d2e9210ac58.mugqic.done
)
macs2_callpeak_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: homer_annotate_peaks
#-------------------------------------------------------------------------------
STEP=homer_annotate_peaks
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_1_JOB_ID: homer_annotate_peaks.LCC9_CTRL_MED1
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.LCC9_CTRL_MED1
JOB_DEPENDENCIES=$macs2_callpeak_1_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.LCC9_CTRL_MED1.730fc8d0af21f7577dccae291d31bd5b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.LCC9_CTRL_MED1.730fc8d0af21f7577dccae291d31bd5b.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/LCC9_CTRL_MED1/LCC9_CTRL_MED1 && \
annotatePeaks.pl \
  peak_call/LCC9_CTRL_MED1/LCC9_CTRL_MED1_peaks.narrowPeak \
  GRCh38 \
  -gsize GRCh38 \
  -cons -CpG \
  -go annotation/LCC9_CTRL_MED1/LCC9_CTRL_MED1 \
  -genomeOntology annotation/LCC9_CTRL_MED1/LCC9_CTRL_MED1 \
  > annotation/LCC9_CTRL_MED1/LCC9_CTRL_MED1.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/LCC9_CTRL_MED1/LCC9_CTRL_MED1.annotated.csv",
  "annotation/LCC9_CTRL_MED1/LCC9_CTRL_MED1",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.LCC9_CTRL_MED1.730fc8d0af21f7577dccae291d31bd5b.mugqic.done
)
homer_annotate_peaks_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_2_JOB_ID: homer_annotate_peaks.LCC9_CTRL_NIPBL
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.LCC9_CTRL_NIPBL
JOB_DEPENDENCIES=$macs2_callpeak_2_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.LCC9_CTRL_NIPBL.dfcb0278af10e732095efaf2bcbcd66f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.LCC9_CTRL_NIPBL.dfcb0278af10e732095efaf2bcbcd66f.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/LCC9_CTRL_NIPBL/LCC9_CTRL_NIPBL && \
annotatePeaks.pl \
  peak_call/LCC9_CTRL_NIPBL/LCC9_CTRL_NIPBL_peaks.narrowPeak \
  GRCh38 \
  -gsize GRCh38 \
  -cons -CpG \
  -go annotation/LCC9_CTRL_NIPBL/LCC9_CTRL_NIPBL \
  -genomeOntology annotation/LCC9_CTRL_NIPBL/LCC9_CTRL_NIPBL \
  > annotation/LCC9_CTRL_NIPBL/LCC9_CTRL_NIPBL.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/LCC9_CTRL_NIPBL/LCC9_CTRL_NIPBL.annotated.csv",
  "annotation/LCC9_CTRL_NIPBL/LCC9_CTRL_NIPBL",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.LCC9_CTRL_NIPBL.dfcb0278af10e732095efaf2bcbcd66f.mugqic.done
)
homer_annotate_peaks_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_3_JOB_ID: homer_annotate_peaks.LCC9_CTRL_SMC1
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.LCC9_CTRL_SMC1
JOB_DEPENDENCIES=$macs2_callpeak_3_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.LCC9_CTRL_SMC1.b697156e4f4ddb066d26da7a233d2e4c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.LCC9_CTRL_SMC1.b697156e4f4ddb066d26da7a233d2e4c.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/LCC9_CTRL_SMC1/LCC9_CTRL_SMC1 && \
annotatePeaks.pl \
  peak_call/LCC9_CTRL_SMC1/LCC9_CTRL_SMC1_peaks.narrowPeak \
  GRCh38 \
  -gsize GRCh38 \
  -cons -CpG \
  -go annotation/LCC9_CTRL_SMC1/LCC9_CTRL_SMC1 \
  -genomeOntology annotation/LCC9_CTRL_SMC1/LCC9_CTRL_SMC1 \
  > annotation/LCC9_CTRL_SMC1/LCC9_CTRL_SMC1.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/LCC9_CTRL_SMC1/LCC9_CTRL_SMC1.annotated.csv",
  "annotation/LCC9_CTRL_SMC1/LCC9_CTRL_SMC1",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.LCC9_CTRL_SMC1.b697156e4f4ddb066d26da7a233d2e4c.mugqic.done
)
homer_annotate_peaks_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_4_JOB_ID: homer_annotate_peaks.LCC9_CTRL_POL2
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.LCC9_CTRL_POL2
JOB_DEPENDENCIES=$macs2_callpeak_4_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.LCC9_CTRL_POL2.b4c57a95bf40f4ef252304d0b77da95b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.LCC9_CTRL_POL2.b4c57a95bf40f4ef252304d0b77da95b.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/LCC9_CTRL_POL2/LCC9_CTRL_POL2 && \
annotatePeaks.pl \
  peak_call/LCC9_CTRL_POL2/LCC9_CTRL_POL2_peaks.broadPeak \
  GRCh38 \
  -gsize GRCh38 \
  -cons -CpG \
  -go annotation/LCC9_CTRL_POL2/LCC9_CTRL_POL2 \
  -genomeOntology annotation/LCC9_CTRL_POL2/LCC9_CTRL_POL2 \
  > annotation/LCC9_CTRL_POL2/LCC9_CTRL_POL2.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/LCC9_CTRL_POL2/LCC9_CTRL_POL2.annotated.csv",
  "annotation/LCC9_CTRL_POL2/LCC9_CTRL_POL2",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.LCC9_CTRL_POL2.b4c57a95bf40f4ef252304d0b77da95b.mugqic.done
)
homer_annotate_peaks_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_5_JOB_ID: homer_annotate_peaks.LCC9_CTRL_ERA
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.LCC9_CTRL_ERA
JOB_DEPENDENCIES=$macs2_callpeak_5_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.LCC9_CTRL_ERA.19b5b5baf97e42b43ba5c8b038029f30.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.LCC9_CTRL_ERA.19b5b5baf97e42b43ba5c8b038029f30.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/LCC9_CTRL_ERA/LCC9_CTRL_ERA && \
annotatePeaks.pl \
  peak_call/LCC9_CTRL_ERA/LCC9_CTRL_ERA_peaks.narrowPeak \
  GRCh38 \
  -gsize GRCh38 \
  -cons -CpG \
  -go annotation/LCC9_CTRL_ERA/LCC9_CTRL_ERA \
  -genomeOntology annotation/LCC9_CTRL_ERA/LCC9_CTRL_ERA \
  > annotation/LCC9_CTRL_ERA/LCC9_CTRL_ERA.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/LCC9_CTRL_ERA/LCC9_CTRL_ERA.annotated.csv",
  "annotation/LCC9_CTRL_ERA/LCC9_CTRL_ERA",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.LCC9_CTRL_ERA.19b5b5baf97e42b43ba5c8b038029f30.mugqic.done
)
homer_annotate_peaks_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_6_JOB_ID: homer_annotate_peaks.LCC9_E2_MED1
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.LCC9_E2_MED1
JOB_DEPENDENCIES=$macs2_callpeak_6_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.LCC9_E2_MED1.6c0d316b6da6ce00ce44d885e5590e94.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.LCC9_E2_MED1.6c0d316b6da6ce00ce44d885e5590e94.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/LCC9_E2_MED1/LCC9_E2_MED1 && \
annotatePeaks.pl \
  peak_call/LCC9_E2_MED1/LCC9_E2_MED1_peaks.narrowPeak \
  GRCh38 \
  -gsize GRCh38 \
  -cons -CpG \
  -go annotation/LCC9_E2_MED1/LCC9_E2_MED1 \
  -genomeOntology annotation/LCC9_E2_MED1/LCC9_E2_MED1 \
  > annotation/LCC9_E2_MED1/LCC9_E2_MED1.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/LCC9_E2_MED1/LCC9_E2_MED1.annotated.csv",
  "annotation/LCC9_E2_MED1/LCC9_E2_MED1",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.LCC9_E2_MED1.6c0d316b6da6ce00ce44d885e5590e94.mugqic.done
)
homer_annotate_peaks_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_7_JOB_ID: homer_annotate_peaks.LCC9_E2_NIPBL
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.LCC9_E2_NIPBL
JOB_DEPENDENCIES=$macs2_callpeak_7_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.LCC9_E2_NIPBL.06b9e8dcfde1a699453760552c87bd5b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.LCC9_E2_NIPBL.06b9e8dcfde1a699453760552c87bd5b.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/LCC9_E2_NIPBL/LCC9_E2_NIPBL && \
annotatePeaks.pl \
  peak_call/LCC9_E2_NIPBL/LCC9_E2_NIPBL_peaks.narrowPeak \
  GRCh38 \
  -gsize GRCh38 \
  -cons -CpG \
  -go annotation/LCC9_E2_NIPBL/LCC9_E2_NIPBL \
  -genomeOntology annotation/LCC9_E2_NIPBL/LCC9_E2_NIPBL \
  > annotation/LCC9_E2_NIPBL/LCC9_E2_NIPBL.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/LCC9_E2_NIPBL/LCC9_E2_NIPBL.annotated.csv",
  "annotation/LCC9_E2_NIPBL/LCC9_E2_NIPBL",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.LCC9_E2_NIPBL.06b9e8dcfde1a699453760552c87bd5b.mugqic.done
)
homer_annotate_peaks_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_8_JOB_ID: homer_annotate_peaks.LCC9_E2_SMC1
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.LCC9_E2_SMC1
JOB_DEPENDENCIES=$macs2_callpeak_8_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.LCC9_E2_SMC1.155bd51378b84e19eceba72c2201fd06.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.LCC9_E2_SMC1.155bd51378b84e19eceba72c2201fd06.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/LCC9_E2_SMC1/LCC9_E2_SMC1 && \
annotatePeaks.pl \
  peak_call/LCC9_E2_SMC1/LCC9_E2_SMC1_peaks.narrowPeak \
  GRCh38 \
  -gsize GRCh38 \
  -cons -CpG \
  -go annotation/LCC9_E2_SMC1/LCC9_E2_SMC1 \
  -genomeOntology annotation/LCC9_E2_SMC1/LCC9_E2_SMC1 \
  > annotation/LCC9_E2_SMC1/LCC9_E2_SMC1.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/LCC9_E2_SMC1/LCC9_E2_SMC1.annotated.csv",
  "annotation/LCC9_E2_SMC1/LCC9_E2_SMC1",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.LCC9_E2_SMC1.155bd51378b84e19eceba72c2201fd06.mugqic.done
)
homer_annotate_peaks_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_9_JOB_ID: homer_annotate_peaks.LCC9_E2_POL2
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.LCC9_E2_POL2
JOB_DEPENDENCIES=$macs2_callpeak_9_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.LCC9_E2_POL2.32af91ba6a20a1621ce2d25299c41d93.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.LCC9_E2_POL2.32af91ba6a20a1621ce2d25299c41d93.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/LCC9_E2_POL2/LCC9_E2_POL2 && \
annotatePeaks.pl \
  peak_call/LCC9_E2_POL2/LCC9_E2_POL2_peaks.broadPeak \
  GRCh38 \
  -gsize GRCh38 \
  -cons -CpG \
  -go annotation/LCC9_E2_POL2/LCC9_E2_POL2 \
  -genomeOntology annotation/LCC9_E2_POL2/LCC9_E2_POL2 \
  > annotation/LCC9_E2_POL2/LCC9_E2_POL2.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/LCC9_E2_POL2/LCC9_E2_POL2.annotated.csv",
  "annotation/LCC9_E2_POL2/LCC9_E2_POL2",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.LCC9_E2_POL2.32af91ba6a20a1621ce2d25299c41d93.mugqic.done
)
homer_annotate_peaks_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_10_JOB_ID: homer_annotate_peaks.LCC9_E2_ERA
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.LCC9_E2_ERA
JOB_DEPENDENCIES=$macs2_callpeak_10_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.LCC9_E2_ERA.4164952b41b52515a73fcdcdfb1857c3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.LCC9_E2_ERA.4164952b41b52515a73fcdcdfb1857c3.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/LCC9_E2_ERA/LCC9_E2_ERA && \
annotatePeaks.pl \
  peak_call/LCC9_E2_ERA/LCC9_E2_ERA_peaks.narrowPeak \
  GRCh38 \
  -gsize GRCh38 \
  -cons -CpG \
  -go annotation/LCC9_E2_ERA/LCC9_E2_ERA \
  -genomeOntology annotation/LCC9_E2_ERA/LCC9_E2_ERA \
  > annotation/LCC9_E2_ERA/LCC9_E2_ERA.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/LCC9_E2_ERA/LCC9_E2_ERA.annotated.csv",
  "annotation/LCC9_E2_ERA/LCC9_E2_ERA",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.LCC9_E2_ERA.4164952b41b52515a73fcdcdfb1857c3.mugqic.done
)
homer_annotate_peaks_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_11_JOB_ID: homer_annotate_peaks_report
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks_report
JOB_DEPENDENCIES=$homer_annotate_peaks_1_JOB_ID:$homer_annotate_peaks_2_JOB_ID:$homer_annotate_peaks_3_JOB_ID:$homer_annotate_peaks_4_JOB_ID:$homer_annotate_peaks_5_JOB_ID:$homer_annotate_peaks_6_JOB_ID:$homer_annotate_peaks_7_JOB_ID:$homer_annotate_peaks_8_JOB_ID:$homer_annotate_peaks_9_JOB_ID:$homer_annotate_peaks_10_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks_report.3f2d8190eb5f8e267649aebf22cd639b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks_report.3f2d8190eb5f8e267649aebf22cd639b.mugqic.done'
mkdir -p report/annotation/ && \
cp /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/ChipSeq.homer_annotate_peaks.md report/ && \
for contrast in LCC9_CTRL_MED1 LCC9_CTRL_NIPBL LCC9_CTRL_SMC1 LCC9_CTRL_POL2 LCC9_CTRL_ERA LCC9_E2_MED1 LCC9_E2_NIPBL LCC9_E2_SMC1 LCC9_E2_POL2 LCC9_E2_ERA
do
  rsync -avP annotation/$contrast report/annotation/ && \
  echo -e "* [Gene Annotations for Design $contrast](annotation/$contrast/${contrast}.annotated.csv)
* [HOMER Gene Ontology Annotations for Design $contrast](annotation/$contrast/$contrast/geneOntology.html)
* [HOMER Genome Ontology Annotations for Design $contrast](annotation/$contrast/$contrast/GenomeOntology.html)" \
  >> report/ChipSeq.homer_annotate_peaks.md
done
homer_annotate_peaks_report.3f2d8190eb5f8e267649aebf22cd639b.mugqic.done
)
homer_annotate_peaks_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: homer_find_motifs_genome
#-------------------------------------------------------------------------------
STEP=homer_find_motifs_genome
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_1_JOB_ID: homer_find_motifs_genome.LCC9_CTRL_MED1
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.LCC9_CTRL_MED1
JOB_DEPENDENCIES=$macs2_callpeak_1_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.LCC9_CTRL_MED1.15b159be99aa35267e01d90362bfd3c0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.LCC9_CTRL_MED1.15b159be99aa35267e01d90362bfd3c0.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/LCC9_CTRL_MED1/LCC9_CTRL_MED1 && \
findMotifsGenome.pl \
  peak_call/LCC9_CTRL_MED1/LCC9_CTRL_MED1_peaks.narrowPeak \
  GRCh38 \
  annotation/LCC9_CTRL_MED1/LCC9_CTRL_MED1 \
  -preparsedDir annotation/LCC9_CTRL_MED1/LCC9_CTRL_MED1/preparsed \
  -p 4
homer_find_motifs_genome.LCC9_CTRL_MED1.15b159be99aa35267e01d90362bfd3c0.mugqic.done
)
homer_find_motifs_genome_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_2_JOB_ID: homer_find_motifs_genome.LCC9_CTRL_NIPBL
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.LCC9_CTRL_NIPBL
JOB_DEPENDENCIES=$macs2_callpeak_2_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.LCC9_CTRL_NIPBL.4322011430b733c8c69d538ac7c61d80.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.LCC9_CTRL_NIPBL.4322011430b733c8c69d538ac7c61d80.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/LCC9_CTRL_NIPBL/LCC9_CTRL_NIPBL && \
findMotifsGenome.pl \
  peak_call/LCC9_CTRL_NIPBL/LCC9_CTRL_NIPBL_peaks.narrowPeak \
  GRCh38 \
  annotation/LCC9_CTRL_NIPBL/LCC9_CTRL_NIPBL \
  -preparsedDir annotation/LCC9_CTRL_NIPBL/LCC9_CTRL_NIPBL/preparsed \
  -p 4
homer_find_motifs_genome.LCC9_CTRL_NIPBL.4322011430b733c8c69d538ac7c61d80.mugqic.done
)
homer_find_motifs_genome_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_3_JOB_ID: homer_find_motifs_genome.LCC9_CTRL_SMC1
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.LCC9_CTRL_SMC1
JOB_DEPENDENCIES=$macs2_callpeak_3_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.LCC9_CTRL_SMC1.0a95e655c5e07c77ae7697d7abef1b1d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.LCC9_CTRL_SMC1.0a95e655c5e07c77ae7697d7abef1b1d.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/LCC9_CTRL_SMC1/LCC9_CTRL_SMC1 && \
findMotifsGenome.pl \
  peak_call/LCC9_CTRL_SMC1/LCC9_CTRL_SMC1_peaks.narrowPeak \
  GRCh38 \
  annotation/LCC9_CTRL_SMC1/LCC9_CTRL_SMC1 \
  -preparsedDir annotation/LCC9_CTRL_SMC1/LCC9_CTRL_SMC1/preparsed \
  -p 4
homer_find_motifs_genome.LCC9_CTRL_SMC1.0a95e655c5e07c77ae7697d7abef1b1d.mugqic.done
)
homer_find_motifs_genome_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_4_JOB_ID: homer_find_motifs_genome.LCC9_CTRL_ERA
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.LCC9_CTRL_ERA
JOB_DEPENDENCIES=$macs2_callpeak_5_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.LCC9_CTRL_ERA.7df91922af80e49e7681a4c88f76982f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.LCC9_CTRL_ERA.7df91922af80e49e7681a4c88f76982f.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/LCC9_CTRL_ERA/LCC9_CTRL_ERA && \
findMotifsGenome.pl \
  peak_call/LCC9_CTRL_ERA/LCC9_CTRL_ERA_peaks.narrowPeak \
  GRCh38 \
  annotation/LCC9_CTRL_ERA/LCC9_CTRL_ERA \
  -preparsedDir annotation/LCC9_CTRL_ERA/LCC9_CTRL_ERA/preparsed \
  -p 4
homer_find_motifs_genome.LCC9_CTRL_ERA.7df91922af80e49e7681a4c88f76982f.mugqic.done
)
homer_find_motifs_genome_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_5_JOB_ID: homer_find_motifs_genome.LCC9_E2_MED1
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.LCC9_E2_MED1
JOB_DEPENDENCIES=$macs2_callpeak_6_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.LCC9_E2_MED1.4fa5df38241a9482e763489b25676988.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.LCC9_E2_MED1.4fa5df38241a9482e763489b25676988.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/LCC9_E2_MED1/LCC9_E2_MED1 && \
findMotifsGenome.pl \
  peak_call/LCC9_E2_MED1/LCC9_E2_MED1_peaks.narrowPeak \
  GRCh38 \
  annotation/LCC9_E2_MED1/LCC9_E2_MED1 \
  -preparsedDir annotation/LCC9_E2_MED1/LCC9_E2_MED1/preparsed \
  -p 4
homer_find_motifs_genome.LCC9_E2_MED1.4fa5df38241a9482e763489b25676988.mugqic.done
)
homer_find_motifs_genome_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_6_JOB_ID: homer_find_motifs_genome.LCC9_E2_NIPBL
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.LCC9_E2_NIPBL
JOB_DEPENDENCIES=$macs2_callpeak_7_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.LCC9_E2_NIPBL.d0d32773ccd58f19bce4d4a21b385c11.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.LCC9_E2_NIPBL.d0d32773ccd58f19bce4d4a21b385c11.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/LCC9_E2_NIPBL/LCC9_E2_NIPBL && \
findMotifsGenome.pl \
  peak_call/LCC9_E2_NIPBL/LCC9_E2_NIPBL_peaks.narrowPeak \
  GRCh38 \
  annotation/LCC9_E2_NIPBL/LCC9_E2_NIPBL \
  -preparsedDir annotation/LCC9_E2_NIPBL/LCC9_E2_NIPBL/preparsed \
  -p 4
homer_find_motifs_genome.LCC9_E2_NIPBL.d0d32773ccd58f19bce4d4a21b385c11.mugqic.done
)
homer_find_motifs_genome_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_7_JOB_ID: homer_find_motifs_genome.LCC9_E2_SMC1
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.LCC9_E2_SMC1
JOB_DEPENDENCIES=$macs2_callpeak_8_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.LCC9_E2_SMC1.bded3fdc2912f48df38ec2b61e23c958.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.LCC9_E2_SMC1.bded3fdc2912f48df38ec2b61e23c958.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/LCC9_E2_SMC1/LCC9_E2_SMC1 && \
findMotifsGenome.pl \
  peak_call/LCC9_E2_SMC1/LCC9_E2_SMC1_peaks.narrowPeak \
  GRCh38 \
  annotation/LCC9_E2_SMC1/LCC9_E2_SMC1 \
  -preparsedDir annotation/LCC9_E2_SMC1/LCC9_E2_SMC1/preparsed \
  -p 4
homer_find_motifs_genome.LCC9_E2_SMC1.bded3fdc2912f48df38ec2b61e23c958.mugqic.done
)
homer_find_motifs_genome_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_8_JOB_ID: homer_find_motifs_genome.LCC9_E2_ERA
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.LCC9_E2_ERA
JOB_DEPENDENCIES=$macs2_callpeak_10_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.LCC9_E2_ERA.2c7f994c61aec2f57332af5f19d2906d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.LCC9_E2_ERA.2c7f994c61aec2f57332af5f19d2906d.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/LCC9_E2_ERA/LCC9_E2_ERA && \
findMotifsGenome.pl \
  peak_call/LCC9_E2_ERA/LCC9_E2_ERA_peaks.narrowPeak \
  GRCh38 \
  annotation/LCC9_E2_ERA/LCC9_E2_ERA \
  -preparsedDir annotation/LCC9_E2_ERA/LCC9_E2_ERA/preparsed \
  -p 4
homer_find_motifs_genome.LCC9_E2_ERA.2c7f994c61aec2f57332af5f19d2906d.mugqic.done
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
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome_report.d3987dd738c7b0424e748d35420eecff.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome_report.d3987dd738c7b0424e748d35420eecff.mugqic.done'
mkdir -p report/annotation/ && \
cp /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/ChipSeq.homer_find_motifs_genome.md report/ && \
for contrast in LCC9_CTRL_MED1 LCC9_CTRL_NIPBL LCC9_CTRL_SMC1 LCC9_CTRL_ERA LCC9_E2_MED1 LCC9_E2_NIPBL LCC9_E2_SMC1 LCC9_E2_ERA
do
  rsync -avP annotation/$contrast report/annotation/ && \
  echo -e "* [HOMER _De Novo_ Motif Results for Design $contrast](annotation/$contrast/$contrast/homerResults.html)
* [HOMER Known Motif Results for Design $contrast](annotation/$contrast/$contrast/knownResults.html)" \
  >> report/ChipSeq.homer_find_motifs_genome.md
done
homer_find_motifs_genome_report.d3987dd738c7b0424e748d35420eecff.mugqic.done
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
JOB_DEPENDENCIES=$homer_annotate_peaks_1_JOB_ID:$homer_annotate_peaks_2_JOB_ID:$homer_annotate_peaks_3_JOB_ID:$homer_annotate_peaks_4_JOB_ID:$homer_annotate_peaks_5_JOB_ID:$homer_annotate_peaks_6_JOB_ID:$homer_annotate_peaks_7_JOB_ID:$homer_annotate_peaks_8_JOB_ID:$homer_annotate_peaks_9_JOB_ID:$homer_annotate_peaks_10_JOB_ID
JOB_DONE=job_output/annotation_graphs/annotation_graphs.e8a6f2ed079ed952a4fcd24acbf89d0a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'annotation_graphs.e8a6f2ed079ed952a4fcd24acbf89d0a.mugqic.done'
module load mugqic/mugqic_tools/2.1.5 mugqic/R_Bioconductor/3.2.3_3.2 mugqic/pandoc/1.15.2 && \
mkdir -p graphs && \
Rscript $R_TOOLS/chipSeqgenerateAnnotationGraphs.R \
  ../../raw/chip-seq/design.txt \
  /gs/scratch/efournier/sb_cofactor_hr/LCC9/output/chip-pipeline-GRCh38 && \
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
for contrast in LCC9_CTRL_MED1 LCC9_CTRL_NIPBL LCC9_CTRL_SMC1 LCC9_CTRL_ERA LCC9_E2_MED1 LCC9_E2_NIPBL LCC9_E2_SMC1 LCC9_E2_ERA
do
  cp --parents graphs/${contrast}_Misc_Graphs.ps report/
  convert -rotate 90 graphs/${contrast}_Misc_Graphs.ps report/graphs/${contrast}_Misc_Graphs.png
  echo -e "----

![Annotation Statistics for Design $contrast ([download high-res image](graphs/${contrast}_Misc_Graphs.ps))](graphs/${contrast}_Misc_Graphs.png)
" \
  >> report/ChipSeq.annotation_graphs.md
done
annotation_graphs.e8a6f2ed079ed952a4fcd24acbf89d0a.mugqic.done
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
wget "http://mugqic.hpc.mcgill.ca/cgi-bin/pipeline.cgi?hostname=lg-1r17-n04&ip=10.241.129.14&pipeline=ChipSeq&steps=picard_sam_to_fastq,trimmomatic,merge_trimmomatic_stats,bwa_mem_picard_sort_sam,samtools_view_filter,picard_merge_sam_files,picard_mark_duplicates,metrics,homer_make_tag_directory,qc_metrics,homer_make_ucsc_file,macs2_callpeak,homer_annotate_peaks,homer_find_motifs_genome,annotation_graphs&samples=12" --quiet --output-document=/dev/null

