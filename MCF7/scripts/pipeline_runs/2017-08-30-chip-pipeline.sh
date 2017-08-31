#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

#-------------------------------------------------------------------------------
# ChipSeq PBSScheduler Job Submission Bash script
# Version: 2.2.0
# Created on: 2017-08-30T20:32:09
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
#   homer_find_motifs_genome: 11 jobs
#   annotation_graphs: 1 job
#   TOTAL: 138 jobs
#-------------------------------------------------------------------------------

OUTPUT_DIR=/gs/scratch/efournier/sb_cofactor_hr/MCF7/output/chip-pipeline-GRCh38
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
# JOB: picard_sam_to_fastq_1_JOB_ID: picard_sam_to_fastq.1
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.1
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.1.b91cc181f5d4145b5cc1dee5a8f7d164.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.1.b91cc181f5d4145b5cc1dee5a8f7d164.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/SamToFastq.jar \
  VALIDATION_STRINGENCY=LENIENT \
  INPUT=/gs/scratch/efournier/sb_cofactor_hr/MCF7/raw/chip-seq/HI.1943.008.Index_2.MCF7_CTRL_WCE_ChIP8_GB_rep1.bam \
  FASTQ=/gs/scratch/efournier/sb_cofactor_hr/MCF7/raw/chip-seq/HI.1943.008.Index_2.MCF7_CTRL_WCE_ChIP8_GB_rep1.single.fastq.gz
picard_sam_to_fastq.1.b91cc181f5d4145b5cc1dee5a8f7d164.mugqic.done
)
picard_sam_to_fastq_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=3 | grep "[0-9]")
echo "$picard_sam_to_fastq_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_2_JOB_ID: picard_sam_to_fastq.2
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.2
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.2.4a2e67ccfd7f3fdb8d61e542a56eb106.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.2.4a2e67ccfd7f3fdb8d61e542a56eb106.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/SamToFastq.jar \
  VALIDATION_STRINGENCY=LENIENT \
  INPUT=/gs/scratch/efournier/sb_cofactor_hr/MCF7/raw/chip-seq/HI.1943.008.Index_7.MCF7_CTRL_MED1_ChIP8_GB_rep1.bam \
  FASTQ=/gs/scratch/efournier/sb_cofactor_hr/MCF7/raw/chip-seq/HI.1943.008.Index_7.MCF7_CTRL_MED1_ChIP8_GB_rep1.single.fastq.gz
picard_sam_to_fastq.2.4a2e67ccfd7f3fdb8d61e542a56eb106.mugqic.done
)
picard_sam_to_fastq_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=3 | grep "[0-9]")
echo "$picard_sam_to_fastq_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_3_JOB_ID: picard_sam_to_fastq.3
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.3
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.3.079eed2420ac63fc85bbe48927b534d6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.3.079eed2420ac63fc85bbe48927b534d6.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/SamToFastq.jar \
  VALIDATION_STRINGENCY=LENIENT \
  INPUT=/gs/scratch/efournier/sb_cofactor_hr/MCF7/raw/chip-seq/HI.1997.004.Index_1.MCF7_CTRL_NIPBL_ChIP8_GB_rep1.bam \
  FASTQ=/gs/scratch/efournier/sb_cofactor_hr/MCF7/raw/chip-seq/HI.1997.004.Index_1.MCF7_CTRL_NIPBL_ChIP8_GB_rep1.single.fastq.gz
picard_sam_to_fastq.3.079eed2420ac63fc85bbe48927b534d6.mugqic.done
)
picard_sam_to_fastq_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=3 | grep "[0-9]")
echo "$picard_sam_to_fastq_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_4_JOB_ID: picard_sam_to_fastq.4
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.4
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.4.e4e27422d9fc90c3e91e8ace29a3ed0e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.4.e4e27422d9fc90c3e91e8ace29a3ed0e.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/SamToFastq.jar \
  VALIDATION_STRINGENCY=LENIENT \
  INPUT=/gs/scratch/efournier/sb_cofactor_hr/MCF7/raw/chip-seq/HI.1997.004.Index_4.MCF7_CTRL_SMC1_ChIP8_GB_rep1.bam \
  FASTQ=/gs/scratch/efournier/sb_cofactor_hr/MCF7/raw/chip-seq/HI.1997.004.Index_4.MCF7_CTRL_SMC1_ChIP8_GB_rep1.single.fastq.gz
picard_sam_to_fastq.4.e4e27422d9fc90c3e91e8ace29a3ed0e.mugqic.done
)
picard_sam_to_fastq_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=3 | grep "[0-9]")
echo "$picard_sam_to_fastq_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_5_JOB_ID: picard_sam_to_fastq.5
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.5
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.5.a9fe582ff0126f9290c70026819cf4d0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.5.a9fe582ff0126f9290c70026819cf4d0.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/SamToFastq.jar \
  VALIDATION_STRINGENCY=LENIENT \
  INPUT=/gs/scratch/efournier/sb_cofactor_hr/MCF7/raw/chip-seq/HI.1997.005.Index_10.MCF7_CTRL_POL2_ChIP8_GB_rep1.bam \
  FASTQ=/gs/scratch/efournier/sb_cofactor_hr/MCF7/raw/chip-seq/HI.1997.005.Index_10.MCF7_CTRL_POL2_ChIP8_GB_rep1.single.fastq.gz
picard_sam_to_fastq.5.a9fe582ff0126f9290c70026819cf4d0.mugqic.done
)
picard_sam_to_fastq_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=3 | grep "[0-9]")
echo "$picard_sam_to_fastq_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_6_JOB_ID: picard_sam_to_fastq.6
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.6
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.6.01b9ec709a8ddaf74bc02a25e043b100.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.6.01b9ec709a8ddaf74bc02a25e043b100.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/SamToFastq.jar \
  VALIDATION_STRINGENCY=LENIENT \
  INPUT=/gs/scratch/efournier/sb_cofactor_hr/MCF7/raw/chip-seq/HI.1997.005.Index_25.MCF7_CTRL_ERA_ChIP8_GB_rep1.bam \
  FASTQ=/gs/scratch/efournier/sb_cofactor_hr/MCF7/raw/chip-seq/HI.1997.005.Index_25.MCF7_CTRL_ERA_ChIP8_GB_rep1.single.fastq.gz
picard_sam_to_fastq.6.01b9ec709a8ddaf74bc02a25e043b100.mugqic.done
)
picard_sam_to_fastq_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=3 | grep "[0-9]")
echo "$picard_sam_to_fastq_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_7_JOB_ID: picard_sam_to_fastq.13
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.13
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.13.5fcfb2955c3b5f9c8f8e61bd223492c8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.13.5fcfb2955c3b5f9c8f8e61bd223492c8.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/SamToFastq.jar \
  VALIDATION_STRINGENCY=LENIENT \
  INPUT=/gs/scratch/efournier/sb_cofactor_hr/MCF7/raw/chip-seq/HI.1943.008.Index_13.MCF7_E2_WCE_ChIP8_GB_rep1.bam \
  FASTQ=/gs/scratch/efournier/sb_cofactor_hr/MCF7/raw/chip-seq/HI.1943.008.Index_13.MCF7_E2_WCE_ChIP8_GB_rep1.single.fastq.gz
picard_sam_to_fastq.13.5fcfb2955c3b5f9c8f8e61bd223492c8.mugqic.done
)
picard_sam_to_fastq_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=3 | grep "[0-9]")
echo "$picard_sam_to_fastq_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_8_JOB_ID: picard_sam_to_fastq.14
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.14
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.14.286f3839ffd614b2102fe54a154302f8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.14.286f3839ffd614b2102fe54a154302f8.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/SamToFastq.jar \
  VALIDATION_STRINGENCY=LENIENT \
  INPUT=/gs/scratch/efournier/sb_cofactor_hr/MCF7/raw/chip-seq/HI.1943.008.Index_18.MCF7_E2_MED1_ChIP8_GB_rep1.bam \
  FASTQ=/gs/scratch/efournier/sb_cofactor_hr/MCF7/raw/chip-seq/HI.1943.008.Index_18.MCF7_E2_MED1_ChIP8_GB_rep1.single.fastq.gz
picard_sam_to_fastq.14.286f3839ffd614b2102fe54a154302f8.mugqic.done
)
picard_sam_to_fastq_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=3 | grep "[0-9]")
echo "$picard_sam_to_fastq_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_9_JOB_ID: picard_sam_to_fastq.15
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.15
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.15.31a0c81820f31def256129f763d865a8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.15.31a0c81820f31def256129f763d865a8.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/SamToFastq.jar \
  VALIDATION_STRINGENCY=LENIENT \
  INPUT=/gs/scratch/efournier/sb_cofactor_hr/MCF7/raw/chip-seq/HI.1997.004.Index_3.MCF7_E2_NIPBL_ChIP8_GB_rep1.bam \
  FASTQ=/gs/scratch/efournier/sb_cofactor_hr/MCF7/raw/chip-seq/HI.1997.004.Index_3.MCF7_E2_NIPBL_ChIP8_GB_rep1.single.fastq.gz
picard_sam_to_fastq.15.31a0c81820f31def256129f763d865a8.mugqic.done
)
picard_sam_to_fastq_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=3 | grep "[0-9]")
echo "$picard_sam_to_fastq_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_10_JOB_ID: picard_sam_to_fastq.16
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.16
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.16.9e516fab64f0a1f3e5206e08cae6c75e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.16.9e516fab64f0a1f3e5206e08cae6c75e.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/SamToFastq.jar \
  VALIDATION_STRINGENCY=LENIENT \
  INPUT=/gs/scratch/efournier/sb_cofactor_hr/MCF7/raw/chip-seq/HI.1997.004.Index_5.MCF7_E2_SMC1_ChIP8_GB_rep1.bam \
  FASTQ=/gs/scratch/efournier/sb_cofactor_hr/MCF7/raw/chip-seq/HI.1997.004.Index_5.MCF7_E2_SMC1_ChIP8_GB_rep1.single.fastq.gz
picard_sam_to_fastq.16.9e516fab64f0a1f3e5206e08cae6c75e.mugqic.done
)
picard_sam_to_fastq_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=3 | grep "[0-9]")
echo "$picard_sam_to_fastq_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_11_JOB_ID: picard_sam_to_fastq.17
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.17
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.17.d43c98646e775e715ed44791c467d98f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.17.d43c98646e775e715ed44791c467d98f.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/SamToFastq.jar \
  VALIDATION_STRINGENCY=LENIENT \
  INPUT=/gs/scratch/efournier/sb_cofactor_hr/MCF7/raw/chip-seq/HI.1997.005.Index_11.MCF7_E2_POL2_ChIP8_GB_rep1.bam \
  FASTQ=/gs/scratch/efournier/sb_cofactor_hr/MCF7/raw/chip-seq/HI.1997.005.Index_11.MCF7_E2_POL2_ChIP8_GB_rep1.single.fastq.gz
picard_sam_to_fastq.17.d43c98646e775e715ed44791c467d98f.mugqic.done
)
picard_sam_to_fastq_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=3 | grep "[0-9]")
echo "$picard_sam_to_fastq_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_12_JOB_ID: picard_sam_to_fastq.18
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.18
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.18.3020f46c6268c1d49d1a983d35499b32.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.18.3020f46c6268c1d49d1a983d35499b32.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/SamToFastq.jar \
  VALIDATION_STRINGENCY=LENIENT \
  INPUT=/gs/scratch/efournier/sb_cofactor_hr/MCF7/raw/chip-seq/HI.1997.005.Index_21.MCF7_E2_ERA_ChIP8_GB_rep1.bam \
  FASTQ=/gs/scratch/efournier/sb_cofactor_hr/MCF7/raw/chip-seq/HI.1997.005.Index_21.MCF7_E2_ERA_ChIP8_GB_rep1.single.fastq.gz
picard_sam_to_fastq.18.3020f46c6268c1d49d1a983d35499b32.mugqic.done
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
# JOB: trimmomatic_1_JOB_ID: trimmomatic.1
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.1
JOB_DEPENDENCIES=$picard_sam_to_fastq_1_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.1.ab608fa3f594b28740cb774174cd5ad9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.1.ab608fa3f594b28740cb774174cd5ad9.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/MCF7_CTRL_WCE && \
`cat > trim/MCF7_CTRL_WCE/1.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/sb_cofactor_hr/MCF7/raw/chip-seq/HI.1943.008.Index_2.MCF7_CTRL_WCE_ChIP8_GB_rep1.single.fastq.gz \
  trim/MCF7_CTRL_WCE/1.trim.single.fastq.gz \
  ILLUMINACLIP:trim/MCF7_CTRL_WCE/1.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/MCF7_CTRL_WCE/1.trim.log
trimmomatic.1.ab608fa3f594b28740cb774174cd5ad9.mugqic.done
)
trimmomatic_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_2_JOB_ID: trimmomatic.2
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.2
JOB_DEPENDENCIES=$picard_sam_to_fastq_2_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.2.3c7d8d5c3a621eb1827ada58cdf23c07.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.2.3c7d8d5c3a621eb1827ada58cdf23c07.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/MCF7_CTRL_MED1 && \
`cat > trim/MCF7_CTRL_MED1/2.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/sb_cofactor_hr/MCF7/raw/chip-seq/HI.1943.008.Index_7.MCF7_CTRL_MED1_ChIP8_GB_rep1.single.fastq.gz \
  trim/MCF7_CTRL_MED1/2.trim.single.fastq.gz \
  ILLUMINACLIP:trim/MCF7_CTRL_MED1/2.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/MCF7_CTRL_MED1/2.trim.log
trimmomatic.2.3c7d8d5c3a621eb1827ada58cdf23c07.mugqic.done
)
trimmomatic_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_3_JOB_ID: trimmomatic.3
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.3
JOB_DEPENDENCIES=$picard_sam_to_fastq_3_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.3.dad7dddd3d1c53571bfe8799a977fa8c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.3.dad7dddd3d1c53571bfe8799a977fa8c.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/MCF7_CTRL_NIPBL && \
`cat > trim/MCF7_CTRL_NIPBL/3.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/sb_cofactor_hr/MCF7/raw/chip-seq/HI.1997.004.Index_1.MCF7_CTRL_NIPBL_ChIP8_GB_rep1.single.fastq.gz \
  trim/MCF7_CTRL_NIPBL/3.trim.single.fastq.gz \
  ILLUMINACLIP:trim/MCF7_CTRL_NIPBL/3.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/MCF7_CTRL_NIPBL/3.trim.log
trimmomatic.3.dad7dddd3d1c53571bfe8799a977fa8c.mugqic.done
)
trimmomatic_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_4_JOB_ID: trimmomatic.4
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.4
JOB_DEPENDENCIES=$picard_sam_to_fastq_4_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.4.62043272f6a436fd33317e5f00528225.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.4.62043272f6a436fd33317e5f00528225.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/MCF7_CTRL_SMC1 && \
`cat > trim/MCF7_CTRL_SMC1/4.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/sb_cofactor_hr/MCF7/raw/chip-seq/HI.1997.004.Index_4.MCF7_CTRL_SMC1_ChIP8_GB_rep1.single.fastq.gz \
  trim/MCF7_CTRL_SMC1/4.trim.single.fastq.gz \
  ILLUMINACLIP:trim/MCF7_CTRL_SMC1/4.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/MCF7_CTRL_SMC1/4.trim.log
trimmomatic.4.62043272f6a436fd33317e5f00528225.mugqic.done
)
trimmomatic_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_5_JOB_ID: trimmomatic.5
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.5
JOB_DEPENDENCIES=$picard_sam_to_fastq_5_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.5.2b9114a021812760139a11a40c04bebe.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.5.2b9114a021812760139a11a40c04bebe.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/MCF7_CTRL_POL2 && \
`cat > trim/MCF7_CTRL_POL2/5.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/sb_cofactor_hr/MCF7/raw/chip-seq/HI.1997.005.Index_10.MCF7_CTRL_POL2_ChIP8_GB_rep1.single.fastq.gz \
  trim/MCF7_CTRL_POL2/5.trim.single.fastq.gz \
  ILLUMINACLIP:trim/MCF7_CTRL_POL2/5.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/MCF7_CTRL_POL2/5.trim.log
trimmomatic.5.2b9114a021812760139a11a40c04bebe.mugqic.done
)
trimmomatic_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_6_JOB_ID: trimmomatic.6
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.6
JOB_DEPENDENCIES=$picard_sam_to_fastq_6_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.6.15af62755f51e078fa006816d63ef15c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.6.15af62755f51e078fa006816d63ef15c.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/MCF7_CTRL_ERA && \
`cat > trim/MCF7_CTRL_ERA/6.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/sb_cofactor_hr/MCF7/raw/chip-seq/HI.1997.005.Index_25.MCF7_CTRL_ERA_ChIP8_GB_rep1.single.fastq.gz \
  trim/MCF7_CTRL_ERA/6.trim.single.fastq.gz \
  ILLUMINACLIP:trim/MCF7_CTRL_ERA/6.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/MCF7_CTRL_ERA/6.trim.log
trimmomatic.6.15af62755f51e078fa006816d63ef15c.mugqic.done
)
trimmomatic_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_7_JOB_ID: trimmomatic.13
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.13
JOB_DEPENDENCIES=$picard_sam_to_fastq_7_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.13.4b2848579a4def43aaf4f69d4547a8af.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.13.4b2848579a4def43aaf4f69d4547a8af.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/MCF7_E2_WCE && \
`cat > trim/MCF7_E2_WCE/13.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/sb_cofactor_hr/MCF7/raw/chip-seq/HI.1943.008.Index_13.MCF7_E2_WCE_ChIP8_GB_rep1.single.fastq.gz \
  trim/MCF7_E2_WCE/13.trim.single.fastq.gz \
  ILLUMINACLIP:trim/MCF7_E2_WCE/13.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/MCF7_E2_WCE/13.trim.log
trimmomatic.13.4b2848579a4def43aaf4f69d4547a8af.mugqic.done
)
trimmomatic_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_8_JOB_ID: trimmomatic.14
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.14
JOB_DEPENDENCIES=$picard_sam_to_fastq_8_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.14.0dbfc95ec3d6135fa92551b4c566b3be.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.14.0dbfc95ec3d6135fa92551b4c566b3be.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/MCF7_E2_MED1 && \
`cat > trim/MCF7_E2_MED1/14.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/sb_cofactor_hr/MCF7/raw/chip-seq/HI.1943.008.Index_18.MCF7_E2_MED1_ChIP8_GB_rep1.single.fastq.gz \
  trim/MCF7_E2_MED1/14.trim.single.fastq.gz \
  ILLUMINACLIP:trim/MCF7_E2_MED1/14.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/MCF7_E2_MED1/14.trim.log
trimmomatic.14.0dbfc95ec3d6135fa92551b4c566b3be.mugqic.done
)
trimmomatic_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_9_JOB_ID: trimmomatic.15
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.15
JOB_DEPENDENCIES=$picard_sam_to_fastq_9_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.15.204146b8b0c1725341663628fbc9fd2c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.15.204146b8b0c1725341663628fbc9fd2c.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/MCF7_E2_NIPBL && \
`cat > trim/MCF7_E2_NIPBL/15.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/sb_cofactor_hr/MCF7/raw/chip-seq/HI.1997.004.Index_3.MCF7_E2_NIPBL_ChIP8_GB_rep1.single.fastq.gz \
  trim/MCF7_E2_NIPBL/15.trim.single.fastq.gz \
  ILLUMINACLIP:trim/MCF7_E2_NIPBL/15.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/MCF7_E2_NIPBL/15.trim.log
trimmomatic.15.204146b8b0c1725341663628fbc9fd2c.mugqic.done
)
trimmomatic_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_10_JOB_ID: trimmomatic.16
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.16
JOB_DEPENDENCIES=$picard_sam_to_fastq_10_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.16.2652e93cc012c28a5202eaea6e514f65.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.16.2652e93cc012c28a5202eaea6e514f65.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/MCF7_E2_SMC1 && \
`cat > trim/MCF7_E2_SMC1/16.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/sb_cofactor_hr/MCF7/raw/chip-seq/HI.1997.004.Index_5.MCF7_E2_SMC1_ChIP8_GB_rep1.single.fastq.gz \
  trim/MCF7_E2_SMC1/16.trim.single.fastq.gz \
  ILLUMINACLIP:trim/MCF7_E2_SMC1/16.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/MCF7_E2_SMC1/16.trim.log
trimmomatic.16.2652e93cc012c28a5202eaea6e514f65.mugqic.done
)
trimmomatic_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_11_JOB_ID: trimmomatic.17
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.17
JOB_DEPENDENCIES=$picard_sam_to_fastq_11_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.17.347f159fa89a9eb6a037452e98cf2441.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.17.347f159fa89a9eb6a037452e98cf2441.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/MCF7_E2_POL2 && \
`cat > trim/MCF7_E2_POL2/17.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/sb_cofactor_hr/MCF7/raw/chip-seq/HI.1997.005.Index_11.MCF7_E2_POL2_ChIP8_GB_rep1.single.fastq.gz \
  trim/MCF7_E2_POL2/17.trim.single.fastq.gz \
  ILLUMINACLIP:trim/MCF7_E2_POL2/17.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/MCF7_E2_POL2/17.trim.log
trimmomatic.17.347f159fa89a9eb6a037452e98cf2441.mugqic.done
)
trimmomatic_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_12_JOB_ID: trimmomatic.18
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.18
JOB_DEPENDENCIES=$picard_sam_to_fastq_12_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.18.5c01567329b9fd0010dc866cfc4b4498.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.18.5c01567329b9fd0010dc866cfc4b4498.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/MCF7_E2_ERA && \
`cat > trim/MCF7_E2_ERA/18.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/sb_cofactor_hr/MCF7/raw/chip-seq/HI.1997.005.Index_21.MCF7_E2_ERA_ChIP8_GB_rep1.single.fastq.gz \
  trim/MCF7_E2_ERA/18.trim.single.fastq.gz \
  ILLUMINACLIP:trim/MCF7_E2_ERA/18.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/MCF7_E2_ERA/18.trim.log
trimmomatic.18.5c01567329b9fd0010dc866cfc4b4498.mugqic.done
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
JOB_DONE=job_output/merge_trimmomatic_stats/merge_trimmomatic_stats.187c7ceabb3c4c3fc58465262338a827.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'merge_trimmomatic_stats.187c7ceabb3c4c3fc58465262338a827.mugqic.done'
module load mugqic/pandoc/1.15.2 && \
mkdir -p metrics && \
echo 'Sample	Readset	Raw Single Reads #	Surviving Single Reads #	Surviving Single Reads %' > metrics/trimReadsetTable.tsv && \
grep ^Input trim/MCF7_CTRL_WCE/1.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/MCF7_CTRL_WCE	1	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/MCF7_CTRL_MED1/2.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/MCF7_CTRL_MED1	2	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/MCF7_CTRL_NIPBL/3.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/MCF7_CTRL_NIPBL	3	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/MCF7_CTRL_SMC1/4.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/MCF7_CTRL_SMC1	4	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/MCF7_CTRL_POL2/5.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/MCF7_CTRL_POL2	5	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/MCF7_CTRL_ERA/6.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/MCF7_CTRL_ERA	6	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/MCF7_E2_WCE/13.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/MCF7_E2_WCE	13	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/MCF7_E2_MED1/14.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/MCF7_E2_MED1	14	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/MCF7_E2_NIPBL/15.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/MCF7_E2_NIPBL	15	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/MCF7_E2_SMC1/16.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/MCF7_E2_SMC1	16	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/MCF7_E2_POL2/17.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/MCF7_E2_POL2	17	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/MCF7_E2_ERA/18.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/MCF7_E2_ERA	18	\1	\2/' | \
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
merge_trimmomatic_stats.187c7ceabb3c4c3fc58465262338a827.mugqic.done
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
# JOB: bwa_mem_picard_sort_sam_1_JOB_ID: bwa_mem_picard_sort_sam.1
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.1
JOB_DEPENDENCIES=$trimmomatic_1_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.1.c7c38f11b1b71db18f231d5566edaf01.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.1.c7c38f11b1b71db18f231d5566edaf01.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/MCF7_CTRL_WCE/1 && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:1	SM:MCF7_CTRL_WCE	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/MCF7_CTRL_WCE/1.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/MCF7_CTRL_WCE/1/1.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.1.c7c38f11b1b71db18f231d5566edaf01.mugqic.done
)
bwa_mem_picard_sort_sam_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_2_JOB_ID: bwa_mem_picard_sort_sam.2
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.2
JOB_DEPENDENCIES=$trimmomatic_2_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.2.f6d6d1c6cc932a91468a3186b40d5a1b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.2.f6d6d1c6cc932a91468a3186b40d5a1b.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/MCF7_CTRL_MED1/2 && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:2	SM:MCF7_CTRL_MED1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/MCF7_CTRL_MED1/2.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/MCF7_CTRL_MED1/2/2.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.2.f6d6d1c6cc932a91468a3186b40d5a1b.mugqic.done
)
bwa_mem_picard_sort_sam_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_3_JOB_ID: bwa_mem_picard_sort_sam.3
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.3
JOB_DEPENDENCIES=$trimmomatic_3_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.3.0e0749767a915469d8b3e7d761d395bd.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.3.0e0749767a915469d8b3e7d761d395bd.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/MCF7_CTRL_NIPBL/3 && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:3	SM:MCF7_CTRL_NIPBL	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/MCF7_CTRL_NIPBL/3.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/MCF7_CTRL_NIPBL/3/3.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.3.0e0749767a915469d8b3e7d761d395bd.mugqic.done
)
bwa_mem_picard_sort_sam_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_4_JOB_ID: bwa_mem_picard_sort_sam.4
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.4
JOB_DEPENDENCIES=$trimmomatic_4_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.4.0f6747e05a2da95bed67daf711670445.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.4.0f6747e05a2da95bed67daf711670445.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/MCF7_CTRL_SMC1/4 && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:4	SM:MCF7_CTRL_SMC1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/MCF7_CTRL_SMC1/4.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/MCF7_CTRL_SMC1/4/4.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.4.0f6747e05a2da95bed67daf711670445.mugqic.done
)
bwa_mem_picard_sort_sam_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_5_JOB_ID: bwa_mem_picard_sort_sam.5
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.5
JOB_DEPENDENCIES=$trimmomatic_5_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.5.c557be9a45d6c8424981d434ceac4fad.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.5.c557be9a45d6c8424981d434ceac4fad.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/MCF7_CTRL_POL2/5 && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:5	SM:MCF7_CTRL_POL2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/MCF7_CTRL_POL2/5.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/MCF7_CTRL_POL2/5/5.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.5.c557be9a45d6c8424981d434ceac4fad.mugqic.done
)
bwa_mem_picard_sort_sam_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_6_JOB_ID: bwa_mem_picard_sort_sam.6
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.6
JOB_DEPENDENCIES=$trimmomatic_6_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.6.69ee914d9ee17ee3c4b4e7c29c9432f7.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.6.69ee914d9ee17ee3c4b4e7c29c9432f7.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/MCF7_CTRL_ERA/6 && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:6	SM:MCF7_CTRL_ERA	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/MCF7_CTRL_ERA/6.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/MCF7_CTRL_ERA/6/6.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.6.69ee914d9ee17ee3c4b4e7c29c9432f7.mugqic.done
)
bwa_mem_picard_sort_sam_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_7_JOB_ID: bwa_mem_picard_sort_sam.13
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.13
JOB_DEPENDENCIES=$trimmomatic_7_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.13.655aef3c756b73c17a52e26fc8f1704f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.13.655aef3c756b73c17a52e26fc8f1704f.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/MCF7_E2_WCE/13 && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:13	SM:MCF7_E2_WCE	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/MCF7_E2_WCE/13.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/MCF7_E2_WCE/13/13.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.13.655aef3c756b73c17a52e26fc8f1704f.mugqic.done
)
bwa_mem_picard_sort_sam_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_8_JOB_ID: bwa_mem_picard_sort_sam.14
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.14
JOB_DEPENDENCIES=$trimmomatic_8_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.14.480fe6948fa594421b4282e977e6251d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.14.480fe6948fa594421b4282e977e6251d.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/MCF7_E2_MED1/14 && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:14	SM:MCF7_E2_MED1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/MCF7_E2_MED1/14.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/MCF7_E2_MED1/14/14.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.14.480fe6948fa594421b4282e977e6251d.mugqic.done
)
bwa_mem_picard_sort_sam_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_9_JOB_ID: bwa_mem_picard_sort_sam.15
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.15
JOB_DEPENDENCIES=$trimmomatic_9_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.15.77b8fbb1b955090f8c854ed8b1d37b63.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.15.77b8fbb1b955090f8c854ed8b1d37b63.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/MCF7_E2_NIPBL/15 && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:15	SM:MCF7_E2_NIPBL	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/MCF7_E2_NIPBL/15.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/MCF7_E2_NIPBL/15/15.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.15.77b8fbb1b955090f8c854ed8b1d37b63.mugqic.done
)
bwa_mem_picard_sort_sam_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_10_JOB_ID: bwa_mem_picard_sort_sam.16
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.16
JOB_DEPENDENCIES=$trimmomatic_10_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.16.73c29113c230a6428910718114d0d6c0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.16.73c29113c230a6428910718114d0d6c0.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/MCF7_E2_SMC1/16 && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:16	SM:MCF7_E2_SMC1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/MCF7_E2_SMC1/16.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/MCF7_E2_SMC1/16/16.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.16.73c29113c230a6428910718114d0d6c0.mugqic.done
)
bwa_mem_picard_sort_sam_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_11_JOB_ID: bwa_mem_picard_sort_sam.17
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.17
JOB_DEPENDENCIES=$trimmomatic_11_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.17.35faa86aaa23c55c719c8a0f1927ce8f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.17.35faa86aaa23c55c719c8a0f1927ce8f.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/MCF7_E2_POL2/17 && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:17	SM:MCF7_E2_POL2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/MCF7_E2_POL2/17.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/MCF7_E2_POL2/17/17.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.17.35faa86aaa23c55c719c8a0f1927ce8f.mugqic.done
)
bwa_mem_picard_sort_sam_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_12_JOB_ID: bwa_mem_picard_sort_sam.18
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.18
JOB_DEPENDENCIES=$trimmomatic_12_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.18.3ca6ee9a93eb3c1089d077fd4e6d9e42.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.18.3ca6ee9a93eb3c1089d077fd4e6d9e42.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/MCF7_E2_ERA/18 && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:18	SM:MCF7_E2_ERA	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/MCF7_E2_ERA/18.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/MCF7_E2_ERA/18/18.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.18.3ca6ee9a93eb3c1089d077fd4e6d9e42.mugqic.done
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
# JOB: samtools_view_filter_1_JOB_ID: samtools_view_filter.1
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.1
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_1_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.1.6f5b948a6c8c7801aec22d3ec81e9ec3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.1.6f5b948a6c8c7801aec22d3ec81e9ec3.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/MCF7_CTRL_WCE/1/1.sorted.bam \
  > alignment/MCF7_CTRL_WCE/1/1.sorted.filtered.bam
samtools_view_filter.1.6f5b948a6c8c7801aec22d3ec81e9ec3.mugqic.done
)
samtools_view_filter_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_2_JOB_ID: samtools_view_filter.2
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.2
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_2_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.2.5ed95d9c9c90ec95f5789cca1cc53934.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.2.5ed95d9c9c90ec95f5789cca1cc53934.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/MCF7_CTRL_MED1/2/2.sorted.bam \
  > alignment/MCF7_CTRL_MED1/2/2.sorted.filtered.bam
samtools_view_filter.2.5ed95d9c9c90ec95f5789cca1cc53934.mugqic.done
)
samtools_view_filter_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_3_JOB_ID: samtools_view_filter.3
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.3
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_3_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.3.c2bd73254a8d14fee3716daa22f70430.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.3.c2bd73254a8d14fee3716daa22f70430.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/MCF7_CTRL_NIPBL/3/3.sorted.bam \
  > alignment/MCF7_CTRL_NIPBL/3/3.sorted.filtered.bam
samtools_view_filter.3.c2bd73254a8d14fee3716daa22f70430.mugqic.done
)
samtools_view_filter_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_4_JOB_ID: samtools_view_filter.4
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.4
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_4_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.4.51b84f1e612b7912d17b93e0a1ccdd39.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.4.51b84f1e612b7912d17b93e0a1ccdd39.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/MCF7_CTRL_SMC1/4/4.sorted.bam \
  > alignment/MCF7_CTRL_SMC1/4/4.sorted.filtered.bam
samtools_view_filter.4.51b84f1e612b7912d17b93e0a1ccdd39.mugqic.done
)
samtools_view_filter_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_5_JOB_ID: samtools_view_filter.5
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.5
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_5_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.5.ba8dd611083fa674a70d17919fdc8c12.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.5.ba8dd611083fa674a70d17919fdc8c12.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/MCF7_CTRL_POL2/5/5.sorted.bam \
  > alignment/MCF7_CTRL_POL2/5/5.sorted.filtered.bam
samtools_view_filter.5.ba8dd611083fa674a70d17919fdc8c12.mugqic.done
)
samtools_view_filter_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_6_JOB_ID: samtools_view_filter.6
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.6
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_6_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.6.9314ffbeaa0fa8fdbdf57453bed0cc2e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.6.9314ffbeaa0fa8fdbdf57453bed0cc2e.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/MCF7_CTRL_ERA/6/6.sorted.bam \
  > alignment/MCF7_CTRL_ERA/6/6.sorted.filtered.bam
samtools_view_filter.6.9314ffbeaa0fa8fdbdf57453bed0cc2e.mugqic.done
)
samtools_view_filter_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_7_JOB_ID: samtools_view_filter.13
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.13
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_7_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.13.de3cfdcb13a082f03675975b485636a3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.13.de3cfdcb13a082f03675975b485636a3.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/MCF7_E2_WCE/13/13.sorted.bam \
  > alignment/MCF7_E2_WCE/13/13.sorted.filtered.bam
samtools_view_filter.13.de3cfdcb13a082f03675975b485636a3.mugqic.done
)
samtools_view_filter_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_8_JOB_ID: samtools_view_filter.14
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.14
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_8_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.14.c78871416019dc979707633cc9394d37.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.14.c78871416019dc979707633cc9394d37.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/MCF7_E2_MED1/14/14.sorted.bam \
  > alignment/MCF7_E2_MED1/14/14.sorted.filtered.bam
samtools_view_filter.14.c78871416019dc979707633cc9394d37.mugqic.done
)
samtools_view_filter_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_9_JOB_ID: samtools_view_filter.15
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.15
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_9_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.15.d213aa02ef4b2a3b281c51f8e365367e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.15.d213aa02ef4b2a3b281c51f8e365367e.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/MCF7_E2_NIPBL/15/15.sorted.bam \
  > alignment/MCF7_E2_NIPBL/15/15.sorted.filtered.bam
samtools_view_filter.15.d213aa02ef4b2a3b281c51f8e365367e.mugqic.done
)
samtools_view_filter_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_10_JOB_ID: samtools_view_filter.16
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.16
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_10_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.16.1f765f03de8cd54362e01128020e9cd2.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.16.1f765f03de8cd54362e01128020e9cd2.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/MCF7_E2_SMC1/16/16.sorted.bam \
  > alignment/MCF7_E2_SMC1/16/16.sorted.filtered.bam
samtools_view_filter.16.1f765f03de8cd54362e01128020e9cd2.mugqic.done
)
samtools_view_filter_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_11_JOB_ID: samtools_view_filter.17
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.17
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_11_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.17.9b104fa3e4e03a5a60b77244e8798ea1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.17.9b104fa3e4e03a5a60b77244e8798ea1.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/MCF7_E2_POL2/17/17.sorted.bam \
  > alignment/MCF7_E2_POL2/17/17.sorted.filtered.bam
samtools_view_filter.17.9b104fa3e4e03a5a60b77244e8798ea1.mugqic.done
)
samtools_view_filter_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_12_JOB_ID: samtools_view_filter.18
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.18
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_12_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.18.7562d092495d8970d8a7b2e647b4c48b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.18.7562d092495d8970d8a7b2e647b4c48b.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/MCF7_E2_ERA/18/18.sorted.bam \
  > alignment/MCF7_E2_ERA/18/18.sorted.filtered.bam
samtools_view_filter.18.7562d092495d8970d8a7b2e647b4c48b.mugqic.done
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
# JOB: picard_merge_sam_files_1_JOB_ID: symlink_readset_sample_bam.MCF7_CTRL_WCE
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.MCF7_CTRL_WCE
JOB_DEPENDENCIES=$samtools_view_filter_1_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.MCF7_CTRL_WCE.6ec034098fee735ab53b05d619f74e12.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.MCF7_CTRL_WCE.6ec034098fee735ab53b05d619f74e12.mugqic.done'
mkdir -p alignment/MCF7_CTRL_WCE && \
ln -s -f 1/1.sorted.filtered.bam alignment/MCF7_CTRL_WCE/MCF7_CTRL_WCE.merged.bam
symlink_readset_sample_bam.MCF7_CTRL_WCE.6ec034098fee735ab53b05d619f74e12.mugqic.done
)
picard_merge_sam_files_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_2_JOB_ID: symlink_readset_sample_bam.MCF7_CTRL_MED1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.MCF7_CTRL_MED1
JOB_DEPENDENCIES=$samtools_view_filter_2_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.MCF7_CTRL_MED1.b0f01ea20b2980c7ef181c7b49fbb720.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.MCF7_CTRL_MED1.b0f01ea20b2980c7ef181c7b49fbb720.mugqic.done'
mkdir -p alignment/MCF7_CTRL_MED1 && \
ln -s -f 2/2.sorted.filtered.bam alignment/MCF7_CTRL_MED1/MCF7_CTRL_MED1.merged.bam
symlink_readset_sample_bam.MCF7_CTRL_MED1.b0f01ea20b2980c7ef181c7b49fbb720.mugqic.done
)
picard_merge_sam_files_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_3_JOB_ID: symlink_readset_sample_bam.MCF7_CTRL_NIPBL
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.MCF7_CTRL_NIPBL
JOB_DEPENDENCIES=$samtools_view_filter_3_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.MCF7_CTRL_NIPBL.ec22e4f0265edf8db0692b92c63eda2e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.MCF7_CTRL_NIPBL.ec22e4f0265edf8db0692b92c63eda2e.mugqic.done'
mkdir -p alignment/MCF7_CTRL_NIPBL && \
ln -s -f 3/3.sorted.filtered.bam alignment/MCF7_CTRL_NIPBL/MCF7_CTRL_NIPBL.merged.bam
symlink_readset_sample_bam.MCF7_CTRL_NIPBL.ec22e4f0265edf8db0692b92c63eda2e.mugqic.done
)
picard_merge_sam_files_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_4_JOB_ID: symlink_readset_sample_bam.MCF7_CTRL_SMC1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.MCF7_CTRL_SMC1
JOB_DEPENDENCIES=$samtools_view_filter_4_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.MCF7_CTRL_SMC1.e21826f86911b8d6741147179409245b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.MCF7_CTRL_SMC1.e21826f86911b8d6741147179409245b.mugqic.done'
mkdir -p alignment/MCF7_CTRL_SMC1 && \
ln -s -f 4/4.sorted.filtered.bam alignment/MCF7_CTRL_SMC1/MCF7_CTRL_SMC1.merged.bam
symlink_readset_sample_bam.MCF7_CTRL_SMC1.e21826f86911b8d6741147179409245b.mugqic.done
)
picard_merge_sam_files_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_5_JOB_ID: symlink_readset_sample_bam.MCF7_CTRL_POL2
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.MCF7_CTRL_POL2
JOB_DEPENDENCIES=$samtools_view_filter_5_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.MCF7_CTRL_POL2.c1ee945125a8a258bc92f8919681779b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.MCF7_CTRL_POL2.c1ee945125a8a258bc92f8919681779b.mugqic.done'
mkdir -p alignment/MCF7_CTRL_POL2 && \
ln -s -f 5/5.sorted.filtered.bam alignment/MCF7_CTRL_POL2/MCF7_CTRL_POL2.merged.bam
symlink_readset_sample_bam.MCF7_CTRL_POL2.c1ee945125a8a258bc92f8919681779b.mugqic.done
)
picard_merge_sam_files_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_6_JOB_ID: symlink_readset_sample_bam.MCF7_CTRL_ERA
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.MCF7_CTRL_ERA
JOB_DEPENDENCIES=$samtools_view_filter_6_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.MCF7_CTRL_ERA.7d720f5a9e20edb32a89f1989bdfd919.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.MCF7_CTRL_ERA.7d720f5a9e20edb32a89f1989bdfd919.mugqic.done'
mkdir -p alignment/MCF7_CTRL_ERA && \
ln -s -f 6/6.sorted.filtered.bam alignment/MCF7_CTRL_ERA/MCF7_CTRL_ERA.merged.bam
symlink_readset_sample_bam.MCF7_CTRL_ERA.7d720f5a9e20edb32a89f1989bdfd919.mugqic.done
)
picard_merge_sam_files_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_7_JOB_ID: symlink_readset_sample_bam.MCF7_E2_WCE
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.MCF7_E2_WCE
JOB_DEPENDENCIES=$samtools_view_filter_7_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.MCF7_E2_WCE.39f87576ec65016687d6613603a99d89.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.MCF7_E2_WCE.39f87576ec65016687d6613603a99d89.mugqic.done'
mkdir -p alignment/MCF7_E2_WCE && \
ln -s -f 13/13.sorted.filtered.bam alignment/MCF7_E2_WCE/MCF7_E2_WCE.merged.bam
symlink_readset_sample_bam.MCF7_E2_WCE.39f87576ec65016687d6613603a99d89.mugqic.done
)
picard_merge_sam_files_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_8_JOB_ID: symlink_readset_sample_bam.MCF7_E2_MED1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.MCF7_E2_MED1
JOB_DEPENDENCIES=$samtools_view_filter_8_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.MCF7_E2_MED1.25e3e7ec0332bf575abf839f25b86570.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.MCF7_E2_MED1.25e3e7ec0332bf575abf839f25b86570.mugqic.done'
mkdir -p alignment/MCF7_E2_MED1 && \
ln -s -f 14/14.sorted.filtered.bam alignment/MCF7_E2_MED1/MCF7_E2_MED1.merged.bam
symlink_readset_sample_bam.MCF7_E2_MED1.25e3e7ec0332bf575abf839f25b86570.mugqic.done
)
picard_merge_sam_files_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_9_JOB_ID: symlink_readset_sample_bam.MCF7_E2_NIPBL
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.MCF7_E2_NIPBL
JOB_DEPENDENCIES=$samtools_view_filter_9_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.MCF7_E2_NIPBL.971024e92a7346e529213939d300d672.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.MCF7_E2_NIPBL.971024e92a7346e529213939d300d672.mugqic.done'
mkdir -p alignment/MCF7_E2_NIPBL && \
ln -s -f 15/15.sorted.filtered.bam alignment/MCF7_E2_NIPBL/MCF7_E2_NIPBL.merged.bam
symlink_readset_sample_bam.MCF7_E2_NIPBL.971024e92a7346e529213939d300d672.mugqic.done
)
picard_merge_sam_files_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_10_JOB_ID: symlink_readset_sample_bam.MCF7_E2_SMC1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.MCF7_E2_SMC1
JOB_DEPENDENCIES=$samtools_view_filter_10_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.MCF7_E2_SMC1.a6ae8a148118a4bbdf2177daf8278415.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.MCF7_E2_SMC1.a6ae8a148118a4bbdf2177daf8278415.mugqic.done'
mkdir -p alignment/MCF7_E2_SMC1 && \
ln -s -f 16/16.sorted.filtered.bam alignment/MCF7_E2_SMC1/MCF7_E2_SMC1.merged.bam
symlink_readset_sample_bam.MCF7_E2_SMC1.a6ae8a148118a4bbdf2177daf8278415.mugqic.done
)
picard_merge_sam_files_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_11_JOB_ID: symlink_readset_sample_bam.MCF7_E2_POL2
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.MCF7_E2_POL2
JOB_DEPENDENCIES=$samtools_view_filter_11_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.MCF7_E2_POL2.d15ce305cdac16ff711f26332982f201.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.MCF7_E2_POL2.d15ce305cdac16ff711f26332982f201.mugqic.done'
mkdir -p alignment/MCF7_E2_POL2 && \
ln -s -f 17/17.sorted.filtered.bam alignment/MCF7_E2_POL2/MCF7_E2_POL2.merged.bam
symlink_readset_sample_bam.MCF7_E2_POL2.d15ce305cdac16ff711f26332982f201.mugqic.done
)
picard_merge_sam_files_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_12_JOB_ID: symlink_readset_sample_bam.MCF7_E2_ERA
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.MCF7_E2_ERA
JOB_DEPENDENCIES=$samtools_view_filter_12_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.MCF7_E2_ERA.55bacac06a89039ca5ffd30beb3be0ab.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.MCF7_E2_ERA.55bacac06a89039ca5ffd30beb3be0ab.mugqic.done'
mkdir -p alignment/MCF7_E2_ERA && \
ln -s -f 18/18.sorted.filtered.bam alignment/MCF7_E2_ERA/MCF7_E2_ERA.merged.bam
symlink_readset_sample_bam.MCF7_E2_ERA.55bacac06a89039ca5ffd30beb3be0ab.mugqic.done
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
# JOB: picard_mark_duplicates_1_JOB_ID: picard_mark_duplicates.MCF7_CTRL_WCE
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.MCF7_CTRL_WCE
JOB_DEPENDENCIES=$picard_merge_sam_files_1_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.MCF7_CTRL_WCE.aab4eacf52ff9d8edc5b454f2a2f6a79.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.MCF7_CTRL_WCE.aab4eacf52ff9d8edc5b454f2a2f6a79.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/MCF7_CTRL_WCE/MCF7_CTRL_WCE.merged.bam \
  OUTPUT=alignment/MCF7_CTRL_WCE/MCF7_CTRL_WCE.sorted.dup.bam \
  METRICS_FILE=alignment/MCF7_CTRL_WCE/MCF7_CTRL_WCE.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.MCF7_CTRL_WCE.aab4eacf52ff9d8edc5b454f2a2f6a79.mugqic.done
)
picard_mark_duplicates_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_2_JOB_ID: picard_mark_duplicates.MCF7_CTRL_MED1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.MCF7_CTRL_MED1
JOB_DEPENDENCIES=$picard_merge_sam_files_2_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.MCF7_CTRL_MED1.75a5127c07957aec6759ca3a5d03a99d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.MCF7_CTRL_MED1.75a5127c07957aec6759ca3a5d03a99d.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/MCF7_CTRL_MED1/MCF7_CTRL_MED1.merged.bam \
  OUTPUT=alignment/MCF7_CTRL_MED1/MCF7_CTRL_MED1.sorted.dup.bam \
  METRICS_FILE=alignment/MCF7_CTRL_MED1/MCF7_CTRL_MED1.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.MCF7_CTRL_MED1.75a5127c07957aec6759ca3a5d03a99d.mugqic.done
)
picard_mark_duplicates_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_3_JOB_ID: picard_mark_duplicates.MCF7_CTRL_NIPBL
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.MCF7_CTRL_NIPBL
JOB_DEPENDENCIES=$picard_merge_sam_files_3_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.MCF7_CTRL_NIPBL.0168f511fd6cfc65359e3d4ae9200bf2.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.MCF7_CTRL_NIPBL.0168f511fd6cfc65359e3d4ae9200bf2.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/MCF7_CTRL_NIPBL/MCF7_CTRL_NIPBL.merged.bam \
  OUTPUT=alignment/MCF7_CTRL_NIPBL/MCF7_CTRL_NIPBL.sorted.dup.bam \
  METRICS_FILE=alignment/MCF7_CTRL_NIPBL/MCF7_CTRL_NIPBL.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.MCF7_CTRL_NIPBL.0168f511fd6cfc65359e3d4ae9200bf2.mugqic.done
)
picard_mark_duplicates_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_4_JOB_ID: picard_mark_duplicates.MCF7_CTRL_SMC1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.MCF7_CTRL_SMC1
JOB_DEPENDENCIES=$picard_merge_sam_files_4_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.MCF7_CTRL_SMC1.9bca22c9ca8895502f7a6c913fa67655.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.MCF7_CTRL_SMC1.9bca22c9ca8895502f7a6c913fa67655.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/MCF7_CTRL_SMC1/MCF7_CTRL_SMC1.merged.bam \
  OUTPUT=alignment/MCF7_CTRL_SMC1/MCF7_CTRL_SMC1.sorted.dup.bam \
  METRICS_FILE=alignment/MCF7_CTRL_SMC1/MCF7_CTRL_SMC1.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.MCF7_CTRL_SMC1.9bca22c9ca8895502f7a6c913fa67655.mugqic.done
)
picard_mark_duplicates_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_5_JOB_ID: picard_mark_duplicates.MCF7_CTRL_POL2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.MCF7_CTRL_POL2
JOB_DEPENDENCIES=$picard_merge_sam_files_5_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.MCF7_CTRL_POL2.fb3b99d51372b4c5dda8a39a7df703b8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.MCF7_CTRL_POL2.fb3b99d51372b4c5dda8a39a7df703b8.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/MCF7_CTRL_POL2/MCF7_CTRL_POL2.merged.bam \
  OUTPUT=alignment/MCF7_CTRL_POL2/MCF7_CTRL_POL2.sorted.dup.bam \
  METRICS_FILE=alignment/MCF7_CTRL_POL2/MCF7_CTRL_POL2.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.MCF7_CTRL_POL2.fb3b99d51372b4c5dda8a39a7df703b8.mugqic.done
)
picard_mark_duplicates_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_6_JOB_ID: picard_mark_duplicates.MCF7_CTRL_ERA
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.MCF7_CTRL_ERA
JOB_DEPENDENCIES=$picard_merge_sam_files_6_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.MCF7_CTRL_ERA.6f0fcfe225482ddd72598e3be5e5e22d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.MCF7_CTRL_ERA.6f0fcfe225482ddd72598e3be5e5e22d.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/MCF7_CTRL_ERA/MCF7_CTRL_ERA.merged.bam \
  OUTPUT=alignment/MCF7_CTRL_ERA/MCF7_CTRL_ERA.sorted.dup.bam \
  METRICS_FILE=alignment/MCF7_CTRL_ERA/MCF7_CTRL_ERA.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.MCF7_CTRL_ERA.6f0fcfe225482ddd72598e3be5e5e22d.mugqic.done
)
picard_mark_duplicates_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_7_JOB_ID: picard_mark_duplicates.MCF7_E2_WCE
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.MCF7_E2_WCE
JOB_DEPENDENCIES=$picard_merge_sam_files_7_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.MCF7_E2_WCE.90413ebf2e2676722f8180f1c50c85bb.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.MCF7_E2_WCE.90413ebf2e2676722f8180f1c50c85bb.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/MCF7_E2_WCE/MCF7_E2_WCE.merged.bam \
  OUTPUT=alignment/MCF7_E2_WCE/MCF7_E2_WCE.sorted.dup.bam \
  METRICS_FILE=alignment/MCF7_E2_WCE/MCF7_E2_WCE.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.MCF7_E2_WCE.90413ebf2e2676722f8180f1c50c85bb.mugqic.done
)
picard_mark_duplicates_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_8_JOB_ID: picard_mark_duplicates.MCF7_E2_MED1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.MCF7_E2_MED1
JOB_DEPENDENCIES=$picard_merge_sam_files_8_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.MCF7_E2_MED1.26a7e225694ef13acd00098728a4af18.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.MCF7_E2_MED1.26a7e225694ef13acd00098728a4af18.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/MCF7_E2_MED1/MCF7_E2_MED1.merged.bam \
  OUTPUT=alignment/MCF7_E2_MED1/MCF7_E2_MED1.sorted.dup.bam \
  METRICS_FILE=alignment/MCF7_E2_MED1/MCF7_E2_MED1.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.MCF7_E2_MED1.26a7e225694ef13acd00098728a4af18.mugqic.done
)
picard_mark_duplicates_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_9_JOB_ID: picard_mark_duplicates.MCF7_E2_NIPBL
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.MCF7_E2_NIPBL
JOB_DEPENDENCIES=$picard_merge_sam_files_9_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.MCF7_E2_NIPBL.0ddb7dc1985a248d1e09d3bf40e20c6a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.MCF7_E2_NIPBL.0ddb7dc1985a248d1e09d3bf40e20c6a.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/MCF7_E2_NIPBL/MCF7_E2_NIPBL.merged.bam \
  OUTPUT=alignment/MCF7_E2_NIPBL/MCF7_E2_NIPBL.sorted.dup.bam \
  METRICS_FILE=alignment/MCF7_E2_NIPBL/MCF7_E2_NIPBL.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.MCF7_E2_NIPBL.0ddb7dc1985a248d1e09d3bf40e20c6a.mugqic.done
)
picard_mark_duplicates_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_10_JOB_ID: picard_mark_duplicates.MCF7_E2_SMC1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.MCF7_E2_SMC1
JOB_DEPENDENCIES=$picard_merge_sam_files_10_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.MCF7_E2_SMC1.2721be97900e799e73834b71618f22c8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.MCF7_E2_SMC1.2721be97900e799e73834b71618f22c8.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/MCF7_E2_SMC1/MCF7_E2_SMC1.merged.bam \
  OUTPUT=alignment/MCF7_E2_SMC1/MCF7_E2_SMC1.sorted.dup.bam \
  METRICS_FILE=alignment/MCF7_E2_SMC1/MCF7_E2_SMC1.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.MCF7_E2_SMC1.2721be97900e799e73834b71618f22c8.mugqic.done
)
picard_mark_duplicates_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_11_JOB_ID: picard_mark_duplicates.MCF7_E2_POL2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.MCF7_E2_POL2
JOB_DEPENDENCIES=$picard_merge_sam_files_11_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.MCF7_E2_POL2.af61d098d996cb0c8706a0023b0e73c5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.MCF7_E2_POL2.af61d098d996cb0c8706a0023b0e73c5.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/MCF7_E2_POL2/MCF7_E2_POL2.merged.bam \
  OUTPUT=alignment/MCF7_E2_POL2/MCF7_E2_POL2.sorted.dup.bam \
  METRICS_FILE=alignment/MCF7_E2_POL2/MCF7_E2_POL2.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.MCF7_E2_POL2.af61d098d996cb0c8706a0023b0e73c5.mugqic.done
)
picard_mark_duplicates_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_12_JOB_ID: picard_mark_duplicates.MCF7_E2_ERA
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.MCF7_E2_ERA
JOB_DEPENDENCIES=$picard_merge_sam_files_12_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.MCF7_E2_ERA.19545eae74811d4f685d8ecaa1dbb556.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.MCF7_E2_ERA.19545eae74811d4f685d8ecaa1dbb556.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/MCF7_E2_ERA/MCF7_E2_ERA.merged.bam \
  OUTPUT=alignment/MCF7_E2_ERA/MCF7_E2_ERA.sorted.dup.bam \
  METRICS_FILE=alignment/MCF7_E2_ERA/MCF7_E2_ERA.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.MCF7_E2_ERA.19545eae74811d4f685d8ecaa1dbb556.mugqic.done
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
JOB_DONE=job_output/metrics/metrics.flagstat.773a2c89a76211336fbfaaa87eef3ea8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'metrics.flagstat.773a2c89a76211336fbfaaa87eef3ea8.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools flagstat \
  alignment/MCF7_CTRL_WCE/MCF7_CTRL_WCE.sorted.dup.bam \
  > alignment/MCF7_CTRL_WCE/MCF7_CTRL_WCE.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/MCF7_CTRL_MED1/MCF7_CTRL_MED1.sorted.dup.bam \
  > alignment/MCF7_CTRL_MED1/MCF7_CTRL_MED1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/MCF7_CTRL_NIPBL/MCF7_CTRL_NIPBL.sorted.dup.bam \
  > alignment/MCF7_CTRL_NIPBL/MCF7_CTRL_NIPBL.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/MCF7_CTRL_SMC1/MCF7_CTRL_SMC1.sorted.dup.bam \
  > alignment/MCF7_CTRL_SMC1/MCF7_CTRL_SMC1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/MCF7_CTRL_POL2/MCF7_CTRL_POL2.sorted.dup.bam \
  > alignment/MCF7_CTRL_POL2/MCF7_CTRL_POL2.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/MCF7_CTRL_ERA/MCF7_CTRL_ERA.sorted.dup.bam \
  > alignment/MCF7_CTRL_ERA/MCF7_CTRL_ERA.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/MCF7_E2_WCE/MCF7_E2_WCE.sorted.dup.bam \
  > alignment/MCF7_E2_WCE/MCF7_E2_WCE.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/MCF7_E2_MED1/MCF7_E2_MED1.sorted.dup.bam \
  > alignment/MCF7_E2_MED1/MCF7_E2_MED1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/MCF7_E2_NIPBL/MCF7_E2_NIPBL.sorted.dup.bam \
  > alignment/MCF7_E2_NIPBL/MCF7_E2_NIPBL.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/MCF7_E2_SMC1/MCF7_E2_SMC1.sorted.dup.bam \
  > alignment/MCF7_E2_SMC1/MCF7_E2_SMC1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/MCF7_E2_POL2/MCF7_E2_POL2.sorted.dup.bam \
  > alignment/MCF7_E2_POL2/MCF7_E2_POL2.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/MCF7_E2_ERA/MCF7_E2_ERA.sorted.dup.bam \
  > alignment/MCF7_E2_ERA/MCF7_E2_ERA.sorted.dup.bam.flagstat
metrics.flagstat.773a2c89a76211336fbfaaa87eef3ea8.mugqic.done
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
JOB_DONE=job_output/metrics/metrics_report.3b2d6071a98f791a3c276aa9b772adc7.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'metrics_report.3b2d6071a98f791a3c276aa9b772adc7.mugqic.done'
module load mugqic/pandoc/1.15.2 && \
for sample in MCF7_CTRL_WCE MCF7_CTRL_MED1 MCF7_CTRL_NIPBL MCF7_CTRL_SMC1 MCF7_CTRL_POL2 MCF7_CTRL_ERA MCF7_E2_WCE MCF7_E2_MED1 MCF7_E2_NIPBL MCF7_E2_SMC1 MCF7_E2_POL2 MCF7_E2_ERA
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

metrics_report.3b2d6071a98f791a3c276aa9b772adc7.mugqic.done
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
# JOB: homer_make_tag_directory_1_JOB_ID: homer_make_tag_directory.MCF7_CTRL_WCE
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.MCF7_CTRL_WCE
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.MCF7_CTRL_WCE.fc3eb953abe9559b525cc6a6fc785c18.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.MCF7_CTRL_WCE.fc3eb953abe9559b525cc6a6fc785c18.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/MCF7_CTRL_WCE \
  alignment/MCF7_CTRL_WCE/MCF7_CTRL_WCE.sorted.dup.bam \
  -checkGC -genome GRCh38
homer_make_tag_directory.MCF7_CTRL_WCE.fc3eb953abe9559b525cc6a6fc785c18.mugqic.done
)
homer_make_tag_directory_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_2_JOB_ID: homer_make_tag_directory.MCF7_CTRL_MED1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.MCF7_CTRL_MED1
JOB_DEPENDENCIES=$picard_mark_duplicates_2_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.MCF7_CTRL_MED1.cc816384daf987e8116c7ecd62ff3ebd.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.MCF7_CTRL_MED1.cc816384daf987e8116c7ecd62ff3ebd.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/MCF7_CTRL_MED1 \
  alignment/MCF7_CTRL_MED1/MCF7_CTRL_MED1.sorted.dup.bam \
  -checkGC -genome GRCh38
homer_make_tag_directory.MCF7_CTRL_MED1.cc816384daf987e8116c7ecd62ff3ebd.mugqic.done
)
homer_make_tag_directory_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_3_JOB_ID: homer_make_tag_directory.MCF7_CTRL_NIPBL
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.MCF7_CTRL_NIPBL
JOB_DEPENDENCIES=$picard_mark_duplicates_3_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.MCF7_CTRL_NIPBL.964d32100210eaf9df2ae53219e1bc3b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.MCF7_CTRL_NIPBL.964d32100210eaf9df2ae53219e1bc3b.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/MCF7_CTRL_NIPBL \
  alignment/MCF7_CTRL_NIPBL/MCF7_CTRL_NIPBL.sorted.dup.bam \
  -checkGC -genome GRCh38
homer_make_tag_directory.MCF7_CTRL_NIPBL.964d32100210eaf9df2ae53219e1bc3b.mugqic.done
)
homer_make_tag_directory_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_4_JOB_ID: homer_make_tag_directory.MCF7_CTRL_SMC1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.MCF7_CTRL_SMC1
JOB_DEPENDENCIES=$picard_mark_duplicates_4_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.MCF7_CTRL_SMC1.81c1e49c1826274db32d19b077cac421.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.MCF7_CTRL_SMC1.81c1e49c1826274db32d19b077cac421.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/MCF7_CTRL_SMC1 \
  alignment/MCF7_CTRL_SMC1/MCF7_CTRL_SMC1.sorted.dup.bam \
  -checkGC -genome GRCh38
homer_make_tag_directory.MCF7_CTRL_SMC1.81c1e49c1826274db32d19b077cac421.mugqic.done
)
homer_make_tag_directory_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_5_JOB_ID: homer_make_tag_directory.MCF7_CTRL_POL2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.MCF7_CTRL_POL2
JOB_DEPENDENCIES=$picard_mark_duplicates_5_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.MCF7_CTRL_POL2.b8a4213d67590b7987c6a6221ff1792d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.MCF7_CTRL_POL2.b8a4213d67590b7987c6a6221ff1792d.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/MCF7_CTRL_POL2 \
  alignment/MCF7_CTRL_POL2/MCF7_CTRL_POL2.sorted.dup.bam \
  -checkGC -genome GRCh38
homer_make_tag_directory.MCF7_CTRL_POL2.b8a4213d67590b7987c6a6221ff1792d.mugqic.done
)
homer_make_tag_directory_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_6_JOB_ID: homer_make_tag_directory.MCF7_CTRL_ERA
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.MCF7_CTRL_ERA
JOB_DEPENDENCIES=$picard_mark_duplicates_6_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.MCF7_CTRL_ERA.aab9478b3617140a8bfe6242f5e22821.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.MCF7_CTRL_ERA.aab9478b3617140a8bfe6242f5e22821.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/MCF7_CTRL_ERA \
  alignment/MCF7_CTRL_ERA/MCF7_CTRL_ERA.sorted.dup.bam \
  -checkGC -genome GRCh38
homer_make_tag_directory.MCF7_CTRL_ERA.aab9478b3617140a8bfe6242f5e22821.mugqic.done
)
homer_make_tag_directory_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_7_JOB_ID: homer_make_tag_directory.MCF7_E2_WCE
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.MCF7_E2_WCE
JOB_DEPENDENCIES=$picard_mark_duplicates_7_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.MCF7_E2_WCE.1b7791b14b156914bccf67473088b6a9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.MCF7_E2_WCE.1b7791b14b156914bccf67473088b6a9.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/MCF7_E2_WCE \
  alignment/MCF7_E2_WCE/MCF7_E2_WCE.sorted.dup.bam \
  -checkGC -genome GRCh38
homer_make_tag_directory.MCF7_E2_WCE.1b7791b14b156914bccf67473088b6a9.mugqic.done
)
homer_make_tag_directory_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_8_JOB_ID: homer_make_tag_directory.MCF7_E2_MED1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.MCF7_E2_MED1
JOB_DEPENDENCIES=$picard_mark_duplicates_8_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.MCF7_E2_MED1.cc0f423bb3d8477df59061100bbb98a9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.MCF7_E2_MED1.cc0f423bb3d8477df59061100bbb98a9.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/MCF7_E2_MED1 \
  alignment/MCF7_E2_MED1/MCF7_E2_MED1.sorted.dup.bam \
  -checkGC -genome GRCh38
homer_make_tag_directory.MCF7_E2_MED1.cc0f423bb3d8477df59061100bbb98a9.mugqic.done
)
homer_make_tag_directory_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_9_JOB_ID: homer_make_tag_directory.MCF7_E2_NIPBL
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.MCF7_E2_NIPBL
JOB_DEPENDENCIES=$picard_mark_duplicates_9_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.MCF7_E2_NIPBL.317c7fc731d900dabf8959c2ac42bfe0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.MCF7_E2_NIPBL.317c7fc731d900dabf8959c2ac42bfe0.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/MCF7_E2_NIPBL \
  alignment/MCF7_E2_NIPBL/MCF7_E2_NIPBL.sorted.dup.bam \
  -checkGC -genome GRCh38
homer_make_tag_directory.MCF7_E2_NIPBL.317c7fc731d900dabf8959c2ac42bfe0.mugqic.done
)
homer_make_tag_directory_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_10_JOB_ID: homer_make_tag_directory.MCF7_E2_SMC1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.MCF7_E2_SMC1
JOB_DEPENDENCIES=$picard_mark_duplicates_10_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.MCF7_E2_SMC1.ff502566f374a43e07da741f50f1c17b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.MCF7_E2_SMC1.ff502566f374a43e07da741f50f1c17b.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/MCF7_E2_SMC1 \
  alignment/MCF7_E2_SMC1/MCF7_E2_SMC1.sorted.dup.bam \
  -checkGC -genome GRCh38
homer_make_tag_directory.MCF7_E2_SMC1.ff502566f374a43e07da741f50f1c17b.mugqic.done
)
homer_make_tag_directory_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_11_JOB_ID: homer_make_tag_directory.MCF7_E2_POL2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.MCF7_E2_POL2
JOB_DEPENDENCIES=$picard_mark_duplicates_11_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.MCF7_E2_POL2.94970ff58feae1fbc3523aa864e58e89.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.MCF7_E2_POL2.94970ff58feae1fbc3523aa864e58e89.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/MCF7_E2_POL2 \
  alignment/MCF7_E2_POL2/MCF7_E2_POL2.sorted.dup.bam \
  -checkGC -genome GRCh38
homer_make_tag_directory.MCF7_E2_POL2.94970ff58feae1fbc3523aa864e58e89.mugqic.done
)
homer_make_tag_directory_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_12_JOB_ID: homer_make_tag_directory.MCF7_E2_ERA
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.MCF7_E2_ERA
JOB_DEPENDENCIES=$picard_mark_duplicates_12_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.MCF7_E2_ERA.61378f4d9a3205155a28871af0d9d98a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.MCF7_E2_ERA.61378f4d9a3205155a28871af0d9d98a.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/MCF7_E2_ERA \
  alignment/MCF7_E2_ERA/MCF7_E2_ERA.sorted.dup.bam \
  -checkGC -genome GRCh38
homer_make_tag_directory.MCF7_E2_ERA.61378f4d9a3205155a28871af0d9d98a.mugqic.done
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
JOB_DONE=job_output/qc_metrics/qc_plots_R.b956215f4d60f01920649ec9601b133f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'qc_plots_R.b956215f4d60f01920649ec9601b133f.mugqic.done'
module load mugqic/mugqic_tools/2.1.5 mugqic/R_Bioconductor/3.2.3_3.2 && \
mkdir -p graphs && \
Rscript $R_TOOLS/chipSeqGenerateQCMetrics.R \
  ../../raw/chip-seq/design.txt \
  /gs/scratch/efournier/sb_cofactor_hr/MCF7/output/chip-pipeline-GRCh38 && \
cp /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/ChipSeq.qc_metrics.md report/ChipSeq.qc_metrics.md && \
for sample in MCF7_CTRL_WCE MCF7_CTRL_MED1 MCF7_CTRL_NIPBL MCF7_CTRL_SMC1 MCF7_CTRL_POL2 MCF7_CTRL_ERA MCF7_E2_WCE MCF7_E2_MED1 MCF7_E2_NIPBL MCF7_E2_SMC1 MCF7_E2_POL2 MCF7_E2_ERA
do
  cp --parents graphs/${sample}_QC_Metrics.ps report/
  convert -rotate 90 graphs/${sample}_QC_Metrics.ps report/graphs/${sample}_QC_Metrics.png
  echo -e "----

![QC Metrics for Sample $sample ([download high-res image](graphs/${sample}_QC_Metrics.ps))](graphs/${sample}_QC_Metrics.png)
" \
  >> report/ChipSeq.qc_metrics.md
done
qc_plots_R.b956215f4d60f01920649ec9601b133f.mugqic.done
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
# JOB: homer_make_ucsc_file_1_JOB_ID: homer_make_ucsc_file.MCF7_CTRL_WCE
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.MCF7_CTRL_WCE
JOB_DEPENDENCIES=$homer_make_tag_directory_1_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.MCF7_CTRL_WCE.a4525fdcf14c8e9f746da173aa285cab.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.MCF7_CTRL_WCE.a4525fdcf14c8e9f746da173aa285cab.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/MCF7_CTRL_WCE && \
makeUCSCfile \
  tags/MCF7_CTRL_WCE | \
gzip -1 -c > tracks/MCF7_CTRL_WCE/MCF7_CTRL_WCE.ucsc.bedGraph.gz
homer_make_ucsc_file.MCF7_CTRL_WCE.a4525fdcf14c8e9f746da173aa285cab.mugqic.done
)
homer_make_ucsc_file_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_2_JOB_ID: homer_make_ucsc_file.MCF7_CTRL_MED1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.MCF7_CTRL_MED1
JOB_DEPENDENCIES=$homer_make_tag_directory_2_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.MCF7_CTRL_MED1.8d01de3fdc44efd4f5ca19a2a2709c2f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.MCF7_CTRL_MED1.8d01de3fdc44efd4f5ca19a2a2709c2f.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/MCF7_CTRL_MED1 && \
makeUCSCfile \
  tags/MCF7_CTRL_MED1 | \
gzip -1 -c > tracks/MCF7_CTRL_MED1/MCF7_CTRL_MED1.ucsc.bedGraph.gz
homer_make_ucsc_file.MCF7_CTRL_MED1.8d01de3fdc44efd4f5ca19a2a2709c2f.mugqic.done
)
homer_make_ucsc_file_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_3_JOB_ID: homer_make_ucsc_file.MCF7_CTRL_NIPBL
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.MCF7_CTRL_NIPBL
JOB_DEPENDENCIES=$homer_make_tag_directory_3_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.MCF7_CTRL_NIPBL.dcfd6270ca176629ecf962981104721b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.MCF7_CTRL_NIPBL.dcfd6270ca176629ecf962981104721b.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/MCF7_CTRL_NIPBL && \
makeUCSCfile \
  tags/MCF7_CTRL_NIPBL | \
gzip -1 -c > tracks/MCF7_CTRL_NIPBL/MCF7_CTRL_NIPBL.ucsc.bedGraph.gz
homer_make_ucsc_file.MCF7_CTRL_NIPBL.dcfd6270ca176629ecf962981104721b.mugqic.done
)
homer_make_ucsc_file_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_4_JOB_ID: homer_make_ucsc_file.MCF7_CTRL_SMC1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.MCF7_CTRL_SMC1
JOB_DEPENDENCIES=$homer_make_tag_directory_4_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.MCF7_CTRL_SMC1.7ab1056ee2093d29810686d4216c86ee.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.MCF7_CTRL_SMC1.7ab1056ee2093d29810686d4216c86ee.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/MCF7_CTRL_SMC1 && \
makeUCSCfile \
  tags/MCF7_CTRL_SMC1 | \
gzip -1 -c > tracks/MCF7_CTRL_SMC1/MCF7_CTRL_SMC1.ucsc.bedGraph.gz
homer_make_ucsc_file.MCF7_CTRL_SMC1.7ab1056ee2093d29810686d4216c86ee.mugqic.done
)
homer_make_ucsc_file_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_5_JOB_ID: homer_make_ucsc_file.MCF7_CTRL_POL2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.MCF7_CTRL_POL2
JOB_DEPENDENCIES=$homer_make_tag_directory_5_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.MCF7_CTRL_POL2.9e87df2a1a83eaf359f0c228c8341c6c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.MCF7_CTRL_POL2.9e87df2a1a83eaf359f0c228c8341c6c.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/MCF7_CTRL_POL2 && \
makeUCSCfile \
  tags/MCF7_CTRL_POL2 | \
gzip -1 -c > tracks/MCF7_CTRL_POL2/MCF7_CTRL_POL2.ucsc.bedGraph.gz
homer_make_ucsc_file.MCF7_CTRL_POL2.9e87df2a1a83eaf359f0c228c8341c6c.mugqic.done
)
homer_make_ucsc_file_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_6_JOB_ID: homer_make_ucsc_file.MCF7_CTRL_ERA
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.MCF7_CTRL_ERA
JOB_DEPENDENCIES=$homer_make_tag_directory_6_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.MCF7_CTRL_ERA.6b3575203e2bb1ad166c72d126a3a2bf.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.MCF7_CTRL_ERA.6b3575203e2bb1ad166c72d126a3a2bf.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/MCF7_CTRL_ERA && \
makeUCSCfile \
  tags/MCF7_CTRL_ERA | \
gzip -1 -c > tracks/MCF7_CTRL_ERA/MCF7_CTRL_ERA.ucsc.bedGraph.gz
homer_make_ucsc_file.MCF7_CTRL_ERA.6b3575203e2bb1ad166c72d126a3a2bf.mugqic.done
)
homer_make_ucsc_file_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_7_JOB_ID: homer_make_ucsc_file.MCF7_E2_WCE
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.MCF7_E2_WCE
JOB_DEPENDENCIES=$homer_make_tag_directory_7_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.MCF7_E2_WCE.d47ce422c7387cc62675c873216b92d5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.MCF7_E2_WCE.d47ce422c7387cc62675c873216b92d5.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/MCF7_E2_WCE && \
makeUCSCfile \
  tags/MCF7_E2_WCE | \
gzip -1 -c > tracks/MCF7_E2_WCE/MCF7_E2_WCE.ucsc.bedGraph.gz
homer_make_ucsc_file.MCF7_E2_WCE.d47ce422c7387cc62675c873216b92d5.mugqic.done
)
homer_make_ucsc_file_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_8_JOB_ID: homer_make_ucsc_file.MCF7_E2_MED1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.MCF7_E2_MED1
JOB_DEPENDENCIES=$homer_make_tag_directory_8_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.MCF7_E2_MED1.99992c5b5266286ffa9cb5c4419dcd87.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.MCF7_E2_MED1.99992c5b5266286ffa9cb5c4419dcd87.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/MCF7_E2_MED1 && \
makeUCSCfile \
  tags/MCF7_E2_MED1 | \
gzip -1 -c > tracks/MCF7_E2_MED1/MCF7_E2_MED1.ucsc.bedGraph.gz
homer_make_ucsc_file.MCF7_E2_MED1.99992c5b5266286ffa9cb5c4419dcd87.mugqic.done
)
homer_make_ucsc_file_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_9_JOB_ID: homer_make_ucsc_file.MCF7_E2_NIPBL
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.MCF7_E2_NIPBL
JOB_DEPENDENCIES=$homer_make_tag_directory_9_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.MCF7_E2_NIPBL.704cac7e471fd7485fbcb79693b16c78.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.MCF7_E2_NIPBL.704cac7e471fd7485fbcb79693b16c78.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/MCF7_E2_NIPBL && \
makeUCSCfile \
  tags/MCF7_E2_NIPBL | \
gzip -1 -c > tracks/MCF7_E2_NIPBL/MCF7_E2_NIPBL.ucsc.bedGraph.gz
homer_make_ucsc_file.MCF7_E2_NIPBL.704cac7e471fd7485fbcb79693b16c78.mugqic.done
)
homer_make_ucsc_file_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_10_JOB_ID: homer_make_ucsc_file.MCF7_E2_SMC1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.MCF7_E2_SMC1
JOB_DEPENDENCIES=$homer_make_tag_directory_10_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.MCF7_E2_SMC1.c306b49b1108e3ec21b3323565ca4778.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.MCF7_E2_SMC1.c306b49b1108e3ec21b3323565ca4778.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/MCF7_E2_SMC1 && \
makeUCSCfile \
  tags/MCF7_E2_SMC1 | \
gzip -1 -c > tracks/MCF7_E2_SMC1/MCF7_E2_SMC1.ucsc.bedGraph.gz
homer_make_ucsc_file.MCF7_E2_SMC1.c306b49b1108e3ec21b3323565ca4778.mugqic.done
)
homer_make_ucsc_file_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_11_JOB_ID: homer_make_ucsc_file.MCF7_E2_POL2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.MCF7_E2_POL2
JOB_DEPENDENCIES=$homer_make_tag_directory_11_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.MCF7_E2_POL2.c1d10d9b8e424db76a4561e5e7105082.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.MCF7_E2_POL2.c1d10d9b8e424db76a4561e5e7105082.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/MCF7_E2_POL2 && \
makeUCSCfile \
  tags/MCF7_E2_POL2 | \
gzip -1 -c > tracks/MCF7_E2_POL2/MCF7_E2_POL2.ucsc.bedGraph.gz
homer_make_ucsc_file.MCF7_E2_POL2.c1d10d9b8e424db76a4561e5e7105082.mugqic.done
)
homer_make_ucsc_file_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_12_JOB_ID: homer_make_ucsc_file.MCF7_E2_ERA
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.MCF7_E2_ERA
JOB_DEPENDENCIES=$homer_make_tag_directory_12_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.MCF7_E2_ERA.d33f65aca95ea9741e5f4278c79199db.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.MCF7_E2_ERA.d33f65aca95ea9741e5f4278c79199db.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/MCF7_E2_ERA && \
makeUCSCfile \
  tags/MCF7_E2_ERA | \
gzip -1 -c > tracks/MCF7_E2_ERA/MCF7_E2_ERA.ucsc.bedGraph.gz
homer_make_ucsc_file.MCF7_E2_ERA.d33f65aca95ea9741e5f4278c79199db.mugqic.done
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
# JOB: macs2_callpeak_1_JOB_ID: macs2_callpeak.MCF7_CTRL_MED1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.MCF7_CTRL_MED1
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_2_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.MCF7_CTRL_MED1.547325fc21fd4c91e992ec02d8d05b5d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.MCF7_CTRL_MED1.547325fc21fd4c91e992ec02d8d05b5d.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/MCF7_CTRL_MED1 && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2479938032.8 \
  --treatment \
  alignment/MCF7_CTRL_MED1/MCF7_CTRL_MED1.sorted.dup.bam \
  --control \
  alignment/MCF7_CTRL_WCE/MCF7_CTRL_WCE.sorted.dup.bam \
  --name peak_call/MCF7_CTRL_MED1/MCF7_CTRL_MED1 \
  >& peak_call/MCF7_CTRL_MED1/MCF7_CTRL_MED1.diag.macs.out
macs2_callpeak.MCF7_CTRL_MED1.547325fc21fd4c91e992ec02d8d05b5d.mugqic.done
)
macs2_callpeak_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_2_JOB_ID: macs2_callpeak.MCF7_CTRL_NIPBL
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.MCF7_CTRL_NIPBL
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_3_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.MCF7_CTRL_NIPBL.258d2981b98d9d1abff814330fd1b321.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.MCF7_CTRL_NIPBL.258d2981b98d9d1abff814330fd1b321.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/MCF7_CTRL_NIPBL && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2479938032.8 \
  --treatment \
  alignment/MCF7_CTRL_NIPBL/MCF7_CTRL_NIPBL.sorted.dup.bam \
  --control \
  alignment/MCF7_CTRL_WCE/MCF7_CTRL_WCE.sorted.dup.bam \
  --name peak_call/MCF7_CTRL_NIPBL/MCF7_CTRL_NIPBL \
  >& peak_call/MCF7_CTRL_NIPBL/MCF7_CTRL_NIPBL.diag.macs.out
macs2_callpeak.MCF7_CTRL_NIPBL.258d2981b98d9d1abff814330fd1b321.mugqic.done
)
macs2_callpeak_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_3_JOB_ID: macs2_callpeak.MCF7_CTRL_SMC1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.MCF7_CTRL_SMC1
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_4_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.MCF7_CTRL_SMC1.9e7b0a565ad81a4d331aeea43017fe35.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.MCF7_CTRL_SMC1.9e7b0a565ad81a4d331aeea43017fe35.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/MCF7_CTRL_SMC1 && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2479938032.8 \
  --treatment \
  alignment/MCF7_CTRL_SMC1/MCF7_CTRL_SMC1.sorted.dup.bam \
  --control \
  alignment/MCF7_CTRL_WCE/MCF7_CTRL_WCE.sorted.dup.bam \
  --name peak_call/MCF7_CTRL_SMC1/MCF7_CTRL_SMC1 \
  >& peak_call/MCF7_CTRL_SMC1/MCF7_CTRL_SMC1.diag.macs.out
macs2_callpeak.MCF7_CTRL_SMC1.9e7b0a565ad81a4d331aeea43017fe35.mugqic.done
)
macs2_callpeak_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_4_JOB_ID: macs2_callpeak.MCF7_CTRL_POL2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.MCF7_CTRL_POL2
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_5_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.MCF7_CTRL_POL2.47aa9469dbb90ee06cf566a7b51abe92.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.MCF7_CTRL_POL2.47aa9469dbb90ee06cf566a7b51abe92.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/MCF7_CTRL_POL2 && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2479938032.8 \
  --treatment \
  alignment/MCF7_CTRL_POL2/MCF7_CTRL_POL2.sorted.dup.bam \
  --control \
  alignment/MCF7_CTRL_WCE/MCF7_CTRL_WCE.sorted.dup.bam \
  --name peak_call/MCF7_CTRL_POL2/MCF7_CTRL_POL2 \
  >& peak_call/MCF7_CTRL_POL2/MCF7_CTRL_POL2.diag.macs.out
macs2_callpeak.MCF7_CTRL_POL2.47aa9469dbb90ee06cf566a7b51abe92.mugqic.done
)
macs2_callpeak_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_5_JOB_ID: macs2_callpeak.MCF7_CTRL_ERA
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.MCF7_CTRL_ERA
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_6_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.MCF7_CTRL_ERA.9d1af4197f2d66815424ab7cc20c94fb.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.MCF7_CTRL_ERA.9d1af4197f2d66815424ab7cc20c94fb.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/MCF7_CTRL_ERA && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2479938032.8 \
  --treatment \
  alignment/MCF7_CTRL_ERA/MCF7_CTRL_ERA.sorted.dup.bam \
  --control \
  alignment/MCF7_CTRL_WCE/MCF7_CTRL_WCE.sorted.dup.bam \
  --name peak_call/MCF7_CTRL_ERA/MCF7_CTRL_ERA \
  >& peak_call/MCF7_CTRL_ERA/MCF7_CTRL_ERA.diag.macs.out
macs2_callpeak.MCF7_CTRL_ERA.9d1af4197f2d66815424ab7cc20c94fb.mugqic.done
)
macs2_callpeak_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_6_JOB_ID: macs2_callpeak.MCF7_E2_MED1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.MCF7_E2_MED1
JOB_DEPENDENCIES=$picard_mark_duplicates_7_JOB_ID:$picard_mark_duplicates_8_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.MCF7_E2_MED1.a19b1689288008382f1734e4879c9c31.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.MCF7_E2_MED1.a19b1689288008382f1734e4879c9c31.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/MCF7_E2_MED1 && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2479938032.8 \
  --treatment \
  alignment/MCF7_E2_MED1/MCF7_E2_MED1.sorted.dup.bam \
  --control \
  alignment/MCF7_E2_WCE/MCF7_E2_WCE.sorted.dup.bam \
  --name peak_call/MCF7_E2_MED1/MCF7_E2_MED1 \
  >& peak_call/MCF7_E2_MED1/MCF7_E2_MED1.diag.macs.out
macs2_callpeak.MCF7_E2_MED1.a19b1689288008382f1734e4879c9c31.mugqic.done
)
macs2_callpeak_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_7_JOB_ID: macs2_callpeak.MCF7_E2_NIPBL
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.MCF7_E2_NIPBL
JOB_DEPENDENCIES=$picard_mark_duplicates_7_JOB_ID:$picard_mark_duplicates_9_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.MCF7_E2_NIPBL.eba101dbc97d5adf3a7b9a3bbf0429d0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.MCF7_E2_NIPBL.eba101dbc97d5adf3a7b9a3bbf0429d0.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/MCF7_E2_NIPBL && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2479938032.8 \
  --treatment \
  alignment/MCF7_E2_NIPBL/MCF7_E2_NIPBL.sorted.dup.bam \
  --control \
  alignment/MCF7_E2_WCE/MCF7_E2_WCE.sorted.dup.bam \
  --name peak_call/MCF7_E2_NIPBL/MCF7_E2_NIPBL \
  >& peak_call/MCF7_E2_NIPBL/MCF7_E2_NIPBL.diag.macs.out
macs2_callpeak.MCF7_E2_NIPBL.eba101dbc97d5adf3a7b9a3bbf0429d0.mugqic.done
)
macs2_callpeak_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_8_JOB_ID: macs2_callpeak.MCF7_E2_SMC1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.MCF7_E2_SMC1
JOB_DEPENDENCIES=$picard_mark_duplicates_7_JOB_ID:$picard_mark_duplicates_10_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.MCF7_E2_SMC1.d4bd99bd37ecf6ced3aa19857bf0596d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.MCF7_E2_SMC1.d4bd99bd37ecf6ced3aa19857bf0596d.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/MCF7_E2_SMC1 && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2479938032.8 \
  --treatment \
  alignment/MCF7_E2_SMC1/MCF7_E2_SMC1.sorted.dup.bam \
  --control \
  alignment/MCF7_E2_WCE/MCF7_E2_WCE.sorted.dup.bam \
  --name peak_call/MCF7_E2_SMC1/MCF7_E2_SMC1 \
  >& peak_call/MCF7_E2_SMC1/MCF7_E2_SMC1.diag.macs.out
macs2_callpeak.MCF7_E2_SMC1.d4bd99bd37ecf6ced3aa19857bf0596d.mugqic.done
)
macs2_callpeak_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_9_JOB_ID: macs2_callpeak.MCF7_E2_POL2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.MCF7_E2_POL2
JOB_DEPENDENCIES=$picard_mark_duplicates_7_JOB_ID:$picard_mark_duplicates_11_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.MCF7_E2_POL2.8a0506b4443fc5683f2933815838cdfb.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.MCF7_E2_POL2.8a0506b4443fc5683f2933815838cdfb.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/MCF7_E2_POL2 && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2479938032.8 \
  --treatment \
  alignment/MCF7_E2_POL2/MCF7_E2_POL2.sorted.dup.bam \
  --control \
  alignment/MCF7_E2_WCE/MCF7_E2_WCE.sorted.dup.bam \
  --name peak_call/MCF7_E2_POL2/MCF7_E2_POL2 \
  >& peak_call/MCF7_E2_POL2/MCF7_E2_POL2.diag.macs.out
macs2_callpeak.MCF7_E2_POL2.8a0506b4443fc5683f2933815838cdfb.mugqic.done
)
macs2_callpeak_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_10_JOB_ID: macs2_callpeak.MCF7_E2_ERA
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.MCF7_E2_ERA
JOB_DEPENDENCIES=$picard_mark_duplicates_7_JOB_ID:$picard_mark_duplicates_12_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.MCF7_E2_ERA.f7af42c2d929d18d50e2c99b0de77572.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.MCF7_E2_ERA.f7af42c2d929d18d50e2c99b0de77572.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/MCF7_E2_ERA && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2479938032.8 \
  --treatment \
  alignment/MCF7_E2_ERA/MCF7_E2_ERA.sorted.dup.bam \
  --control \
  alignment/MCF7_E2_WCE/MCF7_E2_WCE.sorted.dup.bam \
  --name peak_call/MCF7_E2_ERA/MCF7_E2_ERA \
  >& peak_call/MCF7_E2_ERA/MCF7_E2_ERA.diag.macs.out
macs2_callpeak.MCF7_E2_ERA.f7af42c2d929d18d50e2c99b0de77572.mugqic.done
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
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_report.054813886636805c1349b88e9dbd4ff8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_report.054813886636805c1349b88e9dbd4ff8.mugqic.done'
mkdir -p report && \
cp /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/ChipSeq.macs2_callpeak.md report/ && \
for contrast in MCF7_CTRL_MED1 MCF7_CTRL_NIPBL MCF7_CTRL_SMC1 MCF7_CTRL_POL2 MCF7_CTRL_ERA MCF7_E2_MED1 MCF7_E2_NIPBL MCF7_E2_SMC1 MCF7_E2_POL2 MCF7_E2_ERA
do
  cp -a --parents peak_call/$contrast/ report/ && \
  echo -e "* [Peak Calls File for Design $contrast](peak_call/$contrast/${contrast}_peaks.xls)" \
  >> report/ChipSeq.macs2_callpeak.md
done
macs2_callpeak_report.054813886636805c1349b88e9dbd4ff8.mugqic.done
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
# JOB: homer_annotate_peaks_1_JOB_ID: homer_annotate_peaks.MCF7_CTRL_MED1
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.MCF7_CTRL_MED1
JOB_DEPENDENCIES=$macs2_callpeak_1_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.MCF7_CTRL_MED1.c198550d6ef90d06a73ca3a6523b7a13.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.MCF7_CTRL_MED1.c198550d6ef90d06a73ca3a6523b7a13.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/MCF7_CTRL_MED1/MCF7_CTRL_MED1 && \
annotatePeaks.pl \
  peak_call/MCF7_CTRL_MED1/MCF7_CTRL_MED1_peaks.narrowPeak \
  GRCh38 \
  -gsize GRCh38 \
  -cons -CpG \
  -go annotation/MCF7_CTRL_MED1/MCF7_CTRL_MED1 \
  -genomeOntology annotation/MCF7_CTRL_MED1/MCF7_CTRL_MED1 \
  > annotation/MCF7_CTRL_MED1/MCF7_CTRL_MED1.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/MCF7_CTRL_MED1/MCF7_CTRL_MED1.annotated.csv",
  "annotation/MCF7_CTRL_MED1/MCF7_CTRL_MED1",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.MCF7_CTRL_MED1.c198550d6ef90d06a73ca3a6523b7a13.mugqic.done
)
homer_annotate_peaks_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_2_JOB_ID: homer_annotate_peaks.MCF7_CTRL_NIPBL
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.MCF7_CTRL_NIPBL
JOB_DEPENDENCIES=$macs2_callpeak_2_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.MCF7_CTRL_NIPBL.be8d5950d64f0daf6b2832b2c01d8209.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.MCF7_CTRL_NIPBL.be8d5950d64f0daf6b2832b2c01d8209.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/MCF7_CTRL_NIPBL/MCF7_CTRL_NIPBL && \
annotatePeaks.pl \
  peak_call/MCF7_CTRL_NIPBL/MCF7_CTRL_NIPBL_peaks.narrowPeak \
  GRCh38 \
  -gsize GRCh38 \
  -cons -CpG \
  -go annotation/MCF7_CTRL_NIPBL/MCF7_CTRL_NIPBL \
  -genomeOntology annotation/MCF7_CTRL_NIPBL/MCF7_CTRL_NIPBL \
  > annotation/MCF7_CTRL_NIPBL/MCF7_CTRL_NIPBL.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/MCF7_CTRL_NIPBL/MCF7_CTRL_NIPBL.annotated.csv",
  "annotation/MCF7_CTRL_NIPBL/MCF7_CTRL_NIPBL",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.MCF7_CTRL_NIPBL.be8d5950d64f0daf6b2832b2c01d8209.mugqic.done
)
homer_annotate_peaks_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_3_JOB_ID: homer_annotate_peaks.MCF7_CTRL_SMC1
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.MCF7_CTRL_SMC1
JOB_DEPENDENCIES=$macs2_callpeak_3_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.MCF7_CTRL_SMC1.8e2b4569f431072b76ce7de58dd426c7.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.MCF7_CTRL_SMC1.8e2b4569f431072b76ce7de58dd426c7.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/MCF7_CTRL_SMC1/MCF7_CTRL_SMC1 && \
annotatePeaks.pl \
  peak_call/MCF7_CTRL_SMC1/MCF7_CTRL_SMC1_peaks.narrowPeak \
  GRCh38 \
  -gsize GRCh38 \
  -cons -CpG \
  -go annotation/MCF7_CTRL_SMC1/MCF7_CTRL_SMC1 \
  -genomeOntology annotation/MCF7_CTRL_SMC1/MCF7_CTRL_SMC1 \
  > annotation/MCF7_CTRL_SMC1/MCF7_CTRL_SMC1.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/MCF7_CTRL_SMC1/MCF7_CTRL_SMC1.annotated.csv",
  "annotation/MCF7_CTRL_SMC1/MCF7_CTRL_SMC1",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.MCF7_CTRL_SMC1.8e2b4569f431072b76ce7de58dd426c7.mugqic.done
)
homer_annotate_peaks_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_4_JOB_ID: homer_annotate_peaks.MCF7_CTRL_POL2
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.MCF7_CTRL_POL2
JOB_DEPENDENCIES=$macs2_callpeak_4_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.MCF7_CTRL_POL2.d8399fc035c240d655a0bb85e535f331.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.MCF7_CTRL_POL2.d8399fc035c240d655a0bb85e535f331.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/MCF7_CTRL_POL2/MCF7_CTRL_POL2 && \
annotatePeaks.pl \
  peak_call/MCF7_CTRL_POL2/MCF7_CTRL_POL2_peaks.narrowPeak \
  GRCh38 \
  -gsize GRCh38 \
  -cons -CpG \
  -go annotation/MCF7_CTRL_POL2/MCF7_CTRL_POL2 \
  -genomeOntology annotation/MCF7_CTRL_POL2/MCF7_CTRL_POL2 \
  > annotation/MCF7_CTRL_POL2/MCF7_CTRL_POL2.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/MCF7_CTRL_POL2/MCF7_CTRL_POL2.annotated.csv",
  "annotation/MCF7_CTRL_POL2/MCF7_CTRL_POL2",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.MCF7_CTRL_POL2.d8399fc035c240d655a0bb85e535f331.mugqic.done
)
homer_annotate_peaks_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_5_JOB_ID: homer_annotate_peaks.MCF7_CTRL_ERA
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.MCF7_CTRL_ERA
JOB_DEPENDENCIES=$macs2_callpeak_5_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.MCF7_CTRL_ERA.348e503ba6722f3a6614cb288ff4f5fd.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.MCF7_CTRL_ERA.348e503ba6722f3a6614cb288ff4f5fd.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/MCF7_CTRL_ERA/MCF7_CTRL_ERA && \
annotatePeaks.pl \
  peak_call/MCF7_CTRL_ERA/MCF7_CTRL_ERA_peaks.narrowPeak \
  GRCh38 \
  -gsize GRCh38 \
  -cons -CpG \
  -go annotation/MCF7_CTRL_ERA/MCF7_CTRL_ERA \
  -genomeOntology annotation/MCF7_CTRL_ERA/MCF7_CTRL_ERA \
  > annotation/MCF7_CTRL_ERA/MCF7_CTRL_ERA.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/MCF7_CTRL_ERA/MCF7_CTRL_ERA.annotated.csv",
  "annotation/MCF7_CTRL_ERA/MCF7_CTRL_ERA",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.MCF7_CTRL_ERA.348e503ba6722f3a6614cb288ff4f5fd.mugqic.done
)
homer_annotate_peaks_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_6_JOB_ID: homer_annotate_peaks.MCF7_E2_MED1
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.MCF7_E2_MED1
JOB_DEPENDENCIES=$macs2_callpeak_6_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.MCF7_E2_MED1.d3689e0fbf5b40e0ecb183d7647f77a3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.MCF7_E2_MED1.d3689e0fbf5b40e0ecb183d7647f77a3.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/MCF7_E2_MED1/MCF7_E2_MED1 && \
annotatePeaks.pl \
  peak_call/MCF7_E2_MED1/MCF7_E2_MED1_peaks.narrowPeak \
  GRCh38 \
  -gsize GRCh38 \
  -cons -CpG \
  -go annotation/MCF7_E2_MED1/MCF7_E2_MED1 \
  -genomeOntology annotation/MCF7_E2_MED1/MCF7_E2_MED1 \
  > annotation/MCF7_E2_MED1/MCF7_E2_MED1.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/MCF7_E2_MED1/MCF7_E2_MED1.annotated.csv",
  "annotation/MCF7_E2_MED1/MCF7_E2_MED1",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.MCF7_E2_MED1.d3689e0fbf5b40e0ecb183d7647f77a3.mugqic.done
)
homer_annotate_peaks_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_7_JOB_ID: homer_annotate_peaks.MCF7_E2_NIPBL
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.MCF7_E2_NIPBL
JOB_DEPENDENCIES=$macs2_callpeak_7_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.MCF7_E2_NIPBL.ffbe3d13d8c96071d77388900a905751.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.MCF7_E2_NIPBL.ffbe3d13d8c96071d77388900a905751.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/MCF7_E2_NIPBL/MCF7_E2_NIPBL && \
annotatePeaks.pl \
  peak_call/MCF7_E2_NIPBL/MCF7_E2_NIPBL_peaks.narrowPeak \
  GRCh38 \
  -gsize GRCh38 \
  -cons -CpG \
  -go annotation/MCF7_E2_NIPBL/MCF7_E2_NIPBL \
  -genomeOntology annotation/MCF7_E2_NIPBL/MCF7_E2_NIPBL \
  > annotation/MCF7_E2_NIPBL/MCF7_E2_NIPBL.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/MCF7_E2_NIPBL/MCF7_E2_NIPBL.annotated.csv",
  "annotation/MCF7_E2_NIPBL/MCF7_E2_NIPBL",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.MCF7_E2_NIPBL.ffbe3d13d8c96071d77388900a905751.mugqic.done
)
homer_annotate_peaks_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_8_JOB_ID: homer_annotate_peaks.MCF7_E2_SMC1
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.MCF7_E2_SMC1
JOB_DEPENDENCIES=$macs2_callpeak_8_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.MCF7_E2_SMC1.fbd7a45c2edda761981fdf4f687592d7.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.MCF7_E2_SMC1.fbd7a45c2edda761981fdf4f687592d7.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/MCF7_E2_SMC1/MCF7_E2_SMC1 && \
annotatePeaks.pl \
  peak_call/MCF7_E2_SMC1/MCF7_E2_SMC1_peaks.narrowPeak \
  GRCh38 \
  -gsize GRCh38 \
  -cons -CpG \
  -go annotation/MCF7_E2_SMC1/MCF7_E2_SMC1 \
  -genomeOntology annotation/MCF7_E2_SMC1/MCF7_E2_SMC1 \
  > annotation/MCF7_E2_SMC1/MCF7_E2_SMC1.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/MCF7_E2_SMC1/MCF7_E2_SMC1.annotated.csv",
  "annotation/MCF7_E2_SMC1/MCF7_E2_SMC1",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.MCF7_E2_SMC1.fbd7a45c2edda761981fdf4f687592d7.mugqic.done
)
homer_annotate_peaks_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_9_JOB_ID: homer_annotate_peaks.MCF7_E2_POL2
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.MCF7_E2_POL2
JOB_DEPENDENCIES=$macs2_callpeak_9_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.MCF7_E2_POL2.044479935cdb8f6a388988c02a211177.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.MCF7_E2_POL2.044479935cdb8f6a388988c02a211177.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/MCF7_E2_POL2/MCF7_E2_POL2 && \
annotatePeaks.pl \
  peak_call/MCF7_E2_POL2/MCF7_E2_POL2_peaks.narrowPeak \
  GRCh38 \
  -gsize GRCh38 \
  -cons -CpG \
  -go annotation/MCF7_E2_POL2/MCF7_E2_POL2 \
  -genomeOntology annotation/MCF7_E2_POL2/MCF7_E2_POL2 \
  > annotation/MCF7_E2_POL2/MCF7_E2_POL2.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/MCF7_E2_POL2/MCF7_E2_POL2.annotated.csv",
  "annotation/MCF7_E2_POL2/MCF7_E2_POL2",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.MCF7_E2_POL2.044479935cdb8f6a388988c02a211177.mugqic.done
)
homer_annotate_peaks_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_10_JOB_ID: homer_annotate_peaks.MCF7_E2_ERA
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.MCF7_E2_ERA
JOB_DEPENDENCIES=$macs2_callpeak_10_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.MCF7_E2_ERA.66c3cb5d40445de9ddc51cd32ac7ac5e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.MCF7_E2_ERA.66c3cb5d40445de9ddc51cd32ac7ac5e.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/MCF7_E2_ERA/MCF7_E2_ERA && \
annotatePeaks.pl \
  peak_call/MCF7_E2_ERA/MCF7_E2_ERA_peaks.narrowPeak \
  GRCh38 \
  -gsize GRCh38 \
  -cons -CpG \
  -go annotation/MCF7_E2_ERA/MCF7_E2_ERA \
  -genomeOntology annotation/MCF7_E2_ERA/MCF7_E2_ERA \
  > annotation/MCF7_E2_ERA/MCF7_E2_ERA.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/MCF7_E2_ERA/MCF7_E2_ERA.annotated.csv",
  "annotation/MCF7_E2_ERA/MCF7_E2_ERA",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.MCF7_E2_ERA.66c3cb5d40445de9ddc51cd32ac7ac5e.mugqic.done
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
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks_report.5d5d09f3790487c67cd96b43ba46829f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks_report.5d5d09f3790487c67cd96b43ba46829f.mugqic.done'
mkdir -p report/annotation/ && \
cp /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/ChipSeq.homer_annotate_peaks.md report/ && \
for contrast in MCF7_CTRL_MED1 MCF7_CTRL_NIPBL MCF7_CTRL_SMC1 MCF7_CTRL_POL2 MCF7_CTRL_ERA MCF7_E2_MED1 MCF7_E2_NIPBL MCF7_E2_SMC1 MCF7_E2_POL2 MCF7_E2_ERA
do
  rsync -avP annotation/$contrast report/annotation/ && \
  echo -e "* [Gene Annotations for Design $contrast](annotation/$contrast/${contrast}.annotated.csv)
* [HOMER Gene Ontology Annotations for Design $contrast](annotation/$contrast/$contrast/geneOntology.html)
* [HOMER Genome Ontology Annotations for Design $contrast](annotation/$contrast/$contrast/GenomeOntology.html)" \
  >> report/ChipSeq.homer_annotate_peaks.md
done
homer_annotate_peaks_report.5d5d09f3790487c67cd96b43ba46829f.mugqic.done
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
# JOB: homer_find_motifs_genome_1_JOB_ID: homer_find_motifs_genome.MCF7_CTRL_MED1
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.MCF7_CTRL_MED1
JOB_DEPENDENCIES=$macs2_callpeak_1_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.MCF7_CTRL_MED1.79fe453a98c8b63160ee9193350e9a17.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.MCF7_CTRL_MED1.79fe453a98c8b63160ee9193350e9a17.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/MCF7_CTRL_MED1/MCF7_CTRL_MED1 && \
findMotifsGenome.pl \
  peak_call/MCF7_CTRL_MED1/MCF7_CTRL_MED1_peaks.narrowPeak \
  GRCh38 \
  annotation/MCF7_CTRL_MED1/MCF7_CTRL_MED1 \
  -preparsedDir annotation/MCF7_CTRL_MED1/MCF7_CTRL_MED1/preparsed \
  -p 4
homer_find_motifs_genome.MCF7_CTRL_MED1.79fe453a98c8b63160ee9193350e9a17.mugqic.done
)
homer_find_motifs_genome_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_2_JOB_ID: homer_find_motifs_genome.MCF7_CTRL_NIPBL
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.MCF7_CTRL_NIPBL
JOB_DEPENDENCIES=$macs2_callpeak_2_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.MCF7_CTRL_NIPBL.6ff8fc75dfc494f61b8a22c353eb30ed.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.MCF7_CTRL_NIPBL.6ff8fc75dfc494f61b8a22c353eb30ed.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/MCF7_CTRL_NIPBL/MCF7_CTRL_NIPBL && \
findMotifsGenome.pl \
  peak_call/MCF7_CTRL_NIPBL/MCF7_CTRL_NIPBL_peaks.narrowPeak \
  GRCh38 \
  annotation/MCF7_CTRL_NIPBL/MCF7_CTRL_NIPBL \
  -preparsedDir annotation/MCF7_CTRL_NIPBL/MCF7_CTRL_NIPBL/preparsed \
  -p 4
homer_find_motifs_genome.MCF7_CTRL_NIPBL.6ff8fc75dfc494f61b8a22c353eb30ed.mugqic.done
)
homer_find_motifs_genome_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_3_JOB_ID: homer_find_motifs_genome.MCF7_CTRL_SMC1
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.MCF7_CTRL_SMC1
JOB_DEPENDENCIES=$macs2_callpeak_3_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.MCF7_CTRL_SMC1.18544c9172443386fd9f824b75631798.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.MCF7_CTRL_SMC1.18544c9172443386fd9f824b75631798.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/MCF7_CTRL_SMC1/MCF7_CTRL_SMC1 && \
findMotifsGenome.pl \
  peak_call/MCF7_CTRL_SMC1/MCF7_CTRL_SMC1_peaks.narrowPeak \
  GRCh38 \
  annotation/MCF7_CTRL_SMC1/MCF7_CTRL_SMC1 \
  -preparsedDir annotation/MCF7_CTRL_SMC1/MCF7_CTRL_SMC1/preparsed \
  -p 4
homer_find_motifs_genome.MCF7_CTRL_SMC1.18544c9172443386fd9f824b75631798.mugqic.done
)
homer_find_motifs_genome_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_4_JOB_ID: homer_find_motifs_genome.MCF7_CTRL_POL2
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.MCF7_CTRL_POL2
JOB_DEPENDENCIES=$macs2_callpeak_4_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.MCF7_CTRL_POL2.508e49f7fd43498d0a8f2d7e0617d27e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.MCF7_CTRL_POL2.508e49f7fd43498d0a8f2d7e0617d27e.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/MCF7_CTRL_POL2/MCF7_CTRL_POL2 && \
findMotifsGenome.pl \
  peak_call/MCF7_CTRL_POL2/MCF7_CTRL_POL2_peaks.narrowPeak \
  GRCh38 \
  annotation/MCF7_CTRL_POL2/MCF7_CTRL_POL2 \
  -preparsedDir annotation/MCF7_CTRL_POL2/MCF7_CTRL_POL2/preparsed \
  -p 4
homer_find_motifs_genome.MCF7_CTRL_POL2.508e49f7fd43498d0a8f2d7e0617d27e.mugqic.done
)
homer_find_motifs_genome_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_5_JOB_ID: homer_find_motifs_genome.MCF7_CTRL_ERA
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.MCF7_CTRL_ERA
JOB_DEPENDENCIES=$macs2_callpeak_5_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.MCF7_CTRL_ERA.b165007739f19968688f8733311d81db.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.MCF7_CTRL_ERA.b165007739f19968688f8733311d81db.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/MCF7_CTRL_ERA/MCF7_CTRL_ERA && \
findMotifsGenome.pl \
  peak_call/MCF7_CTRL_ERA/MCF7_CTRL_ERA_peaks.narrowPeak \
  GRCh38 \
  annotation/MCF7_CTRL_ERA/MCF7_CTRL_ERA \
  -preparsedDir annotation/MCF7_CTRL_ERA/MCF7_CTRL_ERA/preparsed \
  -p 4
homer_find_motifs_genome.MCF7_CTRL_ERA.b165007739f19968688f8733311d81db.mugqic.done
)
homer_find_motifs_genome_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_6_JOB_ID: homer_find_motifs_genome.MCF7_E2_MED1
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.MCF7_E2_MED1
JOB_DEPENDENCIES=$macs2_callpeak_6_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.MCF7_E2_MED1.fb85fb5912d9cbb63f67e13ce6b04d77.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.MCF7_E2_MED1.fb85fb5912d9cbb63f67e13ce6b04d77.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/MCF7_E2_MED1/MCF7_E2_MED1 && \
findMotifsGenome.pl \
  peak_call/MCF7_E2_MED1/MCF7_E2_MED1_peaks.narrowPeak \
  GRCh38 \
  annotation/MCF7_E2_MED1/MCF7_E2_MED1 \
  -preparsedDir annotation/MCF7_E2_MED1/MCF7_E2_MED1/preparsed \
  -p 4
homer_find_motifs_genome.MCF7_E2_MED1.fb85fb5912d9cbb63f67e13ce6b04d77.mugqic.done
)
homer_find_motifs_genome_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_7_JOB_ID: homer_find_motifs_genome.MCF7_E2_NIPBL
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.MCF7_E2_NIPBL
JOB_DEPENDENCIES=$macs2_callpeak_7_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.MCF7_E2_NIPBL.e74c184b910831b46be03a81c7fd56a7.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.MCF7_E2_NIPBL.e74c184b910831b46be03a81c7fd56a7.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/MCF7_E2_NIPBL/MCF7_E2_NIPBL && \
findMotifsGenome.pl \
  peak_call/MCF7_E2_NIPBL/MCF7_E2_NIPBL_peaks.narrowPeak \
  GRCh38 \
  annotation/MCF7_E2_NIPBL/MCF7_E2_NIPBL \
  -preparsedDir annotation/MCF7_E2_NIPBL/MCF7_E2_NIPBL/preparsed \
  -p 4
homer_find_motifs_genome.MCF7_E2_NIPBL.e74c184b910831b46be03a81c7fd56a7.mugqic.done
)
homer_find_motifs_genome_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_8_JOB_ID: homer_find_motifs_genome.MCF7_E2_SMC1
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.MCF7_E2_SMC1
JOB_DEPENDENCIES=$macs2_callpeak_8_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.MCF7_E2_SMC1.a8404784c840af1d550dc03e9dbe9670.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.MCF7_E2_SMC1.a8404784c840af1d550dc03e9dbe9670.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/MCF7_E2_SMC1/MCF7_E2_SMC1 && \
findMotifsGenome.pl \
  peak_call/MCF7_E2_SMC1/MCF7_E2_SMC1_peaks.narrowPeak \
  GRCh38 \
  annotation/MCF7_E2_SMC1/MCF7_E2_SMC1 \
  -preparsedDir annotation/MCF7_E2_SMC1/MCF7_E2_SMC1/preparsed \
  -p 4
homer_find_motifs_genome.MCF7_E2_SMC1.a8404784c840af1d550dc03e9dbe9670.mugqic.done
)
homer_find_motifs_genome_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_9_JOB_ID: homer_find_motifs_genome.MCF7_E2_POL2
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.MCF7_E2_POL2
JOB_DEPENDENCIES=$macs2_callpeak_9_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.MCF7_E2_POL2.ff97609a1dc934d959936b5860170d02.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.MCF7_E2_POL2.ff97609a1dc934d959936b5860170d02.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/MCF7_E2_POL2/MCF7_E2_POL2 && \
findMotifsGenome.pl \
  peak_call/MCF7_E2_POL2/MCF7_E2_POL2_peaks.narrowPeak \
  GRCh38 \
  annotation/MCF7_E2_POL2/MCF7_E2_POL2 \
  -preparsedDir annotation/MCF7_E2_POL2/MCF7_E2_POL2/preparsed \
  -p 4
homer_find_motifs_genome.MCF7_E2_POL2.ff97609a1dc934d959936b5860170d02.mugqic.done
)
homer_find_motifs_genome_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_10_JOB_ID: homer_find_motifs_genome.MCF7_E2_ERA
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.MCF7_E2_ERA
JOB_DEPENDENCIES=$macs2_callpeak_10_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.MCF7_E2_ERA.ca1c7f307aa706757783e26920a27d5e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.MCF7_E2_ERA.ca1c7f307aa706757783e26920a27d5e.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/MCF7_E2_ERA/MCF7_E2_ERA && \
findMotifsGenome.pl \
  peak_call/MCF7_E2_ERA/MCF7_E2_ERA_peaks.narrowPeak \
  GRCh38 \
  annotation/MCF7_E2_ERA/MCF7_E2_ERA \
  -preparsedDir annotation/MCF7_E2_ERA/MCF7_E2_ERA/preparsed \
  -p 4
homer_find_motifs_genome.MCF7_E2_ERA.ca1c7f307aa706757783e26920a27d5e.mugqic.done
)
homer_find_motifs_genome_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_11_JOB_ID: homer_find_motifs_genome_report
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome_report
JOB_DEPENDENCIES=$homer_find_motifs_genome_1_JOB_ID:$homer_find_motifs_genome_2_JOB_ID:$homer_find_motifs_genome_3_JOB_ID:$homer_find_motifs_genome_4_JOB_ID:$homer_find_motifs_genome_5_JOB_ID:$homer_find_motifs_genome_6_JOB_ID:$homer_find_motifs_genome_7_JOB_ID:$homer_find_motifs_genome_8_JOB_ID:$homer_find_motifs_genome_9_JOB_ID:$homer_find_motifs_genome_10_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome_report.4c95ebb100926c3e7f2a37c329b43e0d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome_report.4c95ebb100926c3e7f2a37c329b43e0d.mugqic.done'
mkdir -p report/annotation/ && \
cp /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/ChipSeq.homer_find_motifs_genome.md report/ && \
for contrast in MCF7_CTRL_MED1 MCF7_CTRL_NIPBL MCF7_CTRL_SMC1 MCF7_CTRL_POL2 MCF7_CTRL_ERA MCF7_E2_MED1 MCF7_E2_NIPBL MCF7_E2_SMC1 MCF7_E2_POL2 MCF7_E2_ERA
do
  rsync -avP annotation/$contrast report/annotation/ && \
  echo -e "* [HOMER _De Novo_ Motif Results for Design $contrast](annotation/$contrast/$contrast/homerResults.html)
* [HOMER Known Motif Results for Design $contrast](annotation/$contrast/$contrast/knownResults.html)" \
  >> report/ChipSeq.homer_find_motifs_genome.md
done
homer_find_motifs_genome_report.4c95ebb100926c3e7f2a37c329b43e0d.mugqic.done
)
homer_find_motifs_genome_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


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
JOB_DONE=job_output/annotation_graphs/annotation_graphs.4ce6516a4417022c103077722c6ddd23.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'annotation_graphs.4ce6516a4417022c103077722c6ddd23.mugqic.done'
module load mugqic/mugqic_tools/2.1.5 mugqic/R_Bioconductor/3.2.3_3.2 mugqic/pandoc/1.15.2 && \
mkdir -p graphs && \
Rscript $R_TOOLS/chipSeqgenerateAnnotationGraphs.R \
  ../../raw/chip-seq/design.txt \
  /gs/scratch/efournier/sb_cofactor_hr/MCF7/output/chip-pipeline-GRCh38 && \
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
for contrast in MCF7_CTRL_MED1 MCF7_CTRL_NIPBL MCF7_CTRL_SMC1 MCF7_CTRL_POL2 MCF7_CTRL_ERA MCF7_E2_MED1 MCF7_E2_NIPBL MCF7_E2_SMC1 MCF7_E2_POL2 MCF7_E2_ERA
do
  cp --parents graphs/${contrast}_Misc_Graphs.ps report/
  convert -rotate 90 graphs/${contrast}_Misc_Graphs.ps report/graphs/${contrast}_Misc_Graphs.png
  echo -e "----

![Annotation Statistics for Design $contrast ([download high-res image](graphs/${contrast}_Misc_Graphs.ps))](graphs/${contrast}_Misc_Graphs.png)
" \
  >> report/ChipSeq.annotation_graphs.md
done
annotation_graphs.4ce6516a4417022c103077722c6ddd23.mugqic.done
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
wget "http://mugqic.hpc.mcgill.ca/cgi-bin/pipeline.cgi?hostname=lg-1r14-n04&ip=10.241.129.4&pipeline=ChipSeq&steps=picard_sam_to_fastq,trimmomatic,merge_trimmomatic_stats,bwa_mem_picard_sort_sam,samtools_view_filter,picard_merge_sam_files,picard_mark_duplicates,metrics,homer_make_tag_directory,qc_metrics,homer_make_ucsc_file,macs2_callpeak,homer_annotate_peaks,homer_find_motifs_genome,annotation_graphs&samples=12" --quiet --output-document=/dev/null

