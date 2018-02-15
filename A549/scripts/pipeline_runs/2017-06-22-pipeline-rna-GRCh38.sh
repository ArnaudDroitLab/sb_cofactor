#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

#-------------------------------------------------------------------------------
# RnaSeq PBSScheduler Job Submission Bash script
# Version: 2.2.1-beta
# Created on: 2017-06-22T19:51:20
# Steps:
#   picard_sam_to_fastq: 20 jobs
#   trimmomatic: 20 jobs
#   merge_trimmomatic_stats: 1 job
#   star: 42 jobs
#   picard_merge_sam_files: 8 jobs
#   picard_sort_sam: 12 jobs
#   picard_mark_duplicates: 12 jobs
#   picard_rna_metrics: 12 jobs
#   estimate_ribosomal_rna: 20 jobs
#   bam_hard_clip: 12 jobs
#   rnaseqc: 2 jobs
#   wiggle: 48 jobs
#   raw_counts: 12 jobs
#   raw_counts_metrics: 4 jobs
#   cufflinks: 12 jobs
#   cuffmerge: 1 job
#   cuffquant: 12 jobs
#   cuffdiff: 13 jobs
#   cuffnorm: 1 job
#   fpkm_correlation_matrix: 2 jobs
#   gq_seq_utils_exploratory_analysis_rnaseq: 3 jobs
#   differential_expression: 1 job
#   TOTAL: 270 jobs
#-------------------------------------------------------------------------------

OUTPUT_DIR=/gs/scratch/efournier/CofactorHR/A549/output/rna-pipeline-GRCh38
JOB_OUTPUT_DIR=$OUTPUT_DIR/job_output
TIMESTAMP=`date +%FT%H.%M.%S`
JOB_LIST=$JOB_OUTPUT_DIR/RnaSeq_job_list_$TIMESTAMP
mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR


#-------------------------------------------------------------------------------
# STEP: picard_sam_to_fastq
#-------------------------------------------------------------------------------
STEP=picard_sam_to_fastq
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_1_JOB_ID: picard_sam_to_fastq.DEX_nonMamm_1_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.DEX_nonMamm_1_1
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.DEX_nonMamm_1_1.7ff67805d9d671234d424e1e466e44e4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.DEX_nonMamm_1_1.7ff67805d9d671234d424e1e466e44e4.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.002.Index_3.SB_A549-DEX-nonMamm.bam \
 FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.002.Index_3.SB_A549-DEX-nonMamm.pair1.fastq.gz \
  SECOND_END_FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.002.Index_3.SB_A549-DEX-nonMamm.pair2.fastq.gz
picard_sam_to_fastq.DEX_nonMamm_1_1.7ff67805d9d671234d424e1e466e44e4.mugqic.done
)
picard_sam_to_fastq_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m | grep "[0-9]")
echo "$picard_sam_to_fastq_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_2_JOB_ID: picard_sam_to_fastq.DEX_nonMamm_1_2
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.DEX_nonMamm_1_2
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.DEX_nonMamm_1_2.16e16e40980adafb622f84758bb31472.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.DEX_nonMamm_1_2.16e16e40980adafb622f84758bb31472.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.003.Index_3.SB_A549-DEX-nonMamm.bam \
 FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.003.Index_3.SB_A549-DEX-nonMamm.pair1.fastq.gz \
  SECOND_END_FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.003.Index_3.SB_A549-DEX-nonMamm.pair2.fastq.gz
picard_sam_to_fastq.DEX_nonMamm_1_2.16e16e40980adafb622f84758bb31472.mugqic.done
)
picard_sam_to_fastq_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m | grep "[0-9]")
echo "$picard_sam_to_fastq_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_3_JOB_ID: picard_sam_to_fastq.DEX_shMED1_1_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.DEX_shMED1_1_1
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.DEX_shMED1_1_1.1ee9357d5706e29c635a1e50326d2f5e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.DEX_shMED1_1_1.1ee9357d5706e29c635a1e50326d2f5e.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.002.Index_22.SB_A549-DEX-shMED1-2.bam \
 FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.002.Index_22.SB_A549-DEX-shMED1-2.pair1.fastq.gz \
  SECOND_END_FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.002.Index_22.SB_A549-DEX-shMED1-2.pair2.fastq.gz
picard_sam_to_fastq.DEX_shMED1_1_1.1ee9357d5706e29c635a1e50326d2f5e.mugqic.done
)
picard_sam_to_fastq_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m | grep "[0-9]")
echo "$picard_sam_to_fastq_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_4_JOB_ID: picard_sam_to_fastq.DEX_shMED1_1_2
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.DEX_shMED1_1_2
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.DEX_shMED1_1_2.703ccd71503023651b7ad0872d1b32b7.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.DEX_shMED1_1_2.703ccd71503023651b7ad0872d1b32b7.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.003.Index_22.SB_A549-DEX-shMED1-2.bam \
 FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.003.Index_22.SB_A549-DEX-shMED1-2.pair1.fastq.gz \
  SECOND_END_FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.003.Index_22.SB_A549-DEX-shMED1-2.pair2.fastq.gz
picard_sam_to_fastq.DEX_shMED1_1_2.703ccd71503023651b7ad0872d1b32b7.mugqic.done
)
picard_sam_to_fastq_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m | grep "[0-9]")
echo "$picard_sam_to_fastq_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_5_JOB_ID: picard_sam_to_fastq.DEX_shNIPBL_1_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.DEX_shNIPBL_1_1
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.DEX_shNIPBL_1_1.279389fb6e667ddd771f1c444a6c6465.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.DEX_shNIPBL_1_1.279389fb6e667ddd771f1c444a6c6465.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.002.Index_9.SB_A549-DEX-shNIPBL-3.bam \
 FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.002.Index_9.SB_A549-DEX-shNIPBL-3.pair1.fastq.gz \
  SECOND_END_FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.002.Index_9.SB_A549-DEX-shNIPBL-3.pair2.fastq.gz
picard_sam_to_fastq.DEX_shNIPBL_1_1.279389fb6e667ddd771f1c444a6c6465.mugqic.done
)
picard_sam_to_fastq_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m | grep "[0-9]")
echo "$picard_sam_to_fastq_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_6_JOB_ID: picard_sam_to_fastq.DEX_shNIPBL_1_2
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.DEX_shNIPBL_1_2
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.DEX_shNIPBL_1_2.2150c73af98bbaf2800e1e11f30e39db.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.DEX_shNIPBL_1_2.2150c73af98bbaf2800e1e11f30e39db.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.003.Index_9.SB_A549-DEX-shNIPBL-3.bam \
 FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.003.Index_9.SB_A549-DEX-shNIPBL-3.pair1.fastq.gz \
  SECOND_END_FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.003.Index_9.SB_A549-DEX-shNIPBL-3.pair2.fastq.gz
picard_sam_to_fastq.DEX_shNIPBL_1_2.2150c73af98bbaf2800e1e11f30e39db.mugqic.done
)
picard_sam_to_fastq_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m | grep "[0-9]")
echo "$picard_sam_to_fastq_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_7_JOB_ID: picard_sam_to_fastq.DEX_shSMC1A_1_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.DEX_shSMC1A_1_1
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.DEX_shSMC1A_1_1.736a9dd3b39659a3403269a5ba4d22a2.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.DEX_shSMC1A_1_1.736a9dd3b39659a3403269a5ba4d22a2.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.002.Index_27.SB_A549-DEX-shSMC1A-3.bam \
 FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.002.Index_27.SB_A549-DEX-shSMC1A-3.pair1.fastq.gz \
  SECOND_END_FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.002.Index_27.SB_A549-DEX-shSMC1A-3.pair2.fastq.gz
picard_sam_to_fastq.DEX_shSMC1A_1_1.736a9dd3b39659a3403269a5ba4d22a2.mugqic.done
)
picard_sam_to_fastq_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m | grep "[0-9]")
echo "$picard_sam_to_fastq_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_8_JOB_ID: picard_sam_to_fastq.DEX_shSMC1A_1_2
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.DEX_shSMC1A_1_2
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.DEX_shSMC1A_1_2.077c2c39fdf888e6cad52a4ea75ee95c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.DEX_shSMC1A_1_2.077c2c39fdf888e6cad52a4ea75ee95c.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.003.Index_27.SB_A549-DEX-shSMC1A-3.bam \
 FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.003.Index_27.SB_A549-DEX-shSMC1A-3.pair1.fastq.gz \
  SECOND_END_FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.003.Index_27.SB_A549-DEX-shSMC1A-3.pair2.fastq.gz
picard_sam_to_fastq.DEX_shSMC1A_1_2.077c2c39fdf888e6cad52a4ea75ee95c.mugqic.done
)
picard_sam_to_fastq_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m | grep "[0-9]")
echo "$picard_sam_to_fastq_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_9_JOB_ID: picard_sam_to_fastq.ETOH_nonMamm_1_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.ETOH_nonMamm_1_1
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.ETOH_nonMamm_1_1.5c7660349eb857c9f3f86b91818d7d60.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.ETOH_nonMamm_1_1.5c7660349eb857c9f3f86b91818d7d60.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.002.Index_1.SB_A549-nonMamm.bam \
 FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.002.Index_1.SB_A549-nonMamm.pair1.fastq.gz \
  SECOND_END_FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.002.Index_1.SB_A549-nonMamm.pair2.fastq.gz
picard_sam_to_fastq.ETOH_nonMamm_1_1.5c7660349eb857c9f3f86b91818d7d60.mugqic.done
)
picard_sam_to_fastq_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m | grep "[0-9]")
echo "$picard_sam_to_fastq_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_10_JOB_ID: picard_sam_to_fastq.ETOH_nonMamm_1_2
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.ETOH_nonMamm_1_2
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.ETOH_nonMamm_1_2.c569bb9aef8411b501b2fab05ce459b8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.ETOH_nonMamm_1_2.c569bb9aef8411b501b2fab05ce459b8.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.003.Index_1.SB_A549-nonMamm.bam \
 FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.003.Index_1.SB_A549-nonMamm.pair1.fastq.gz \
  SECOND_END_FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.003.Index_1.SB_A549-nonMamm.pair2.fastq.gz
picard_sam_to_fastq.ETOH_nonMamm_1_2.c569bb9aef8411b501b2fab05ce459b8.mugqic.done
)
picard_sam_to_fastq_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m | grep "[0-9]")
echo "$picard_sam_to_fastq_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_11_JOB_ID: picard_sam_to_fastq.ETOH_nonMamm_2_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.ETOH_nonMamm_2_1
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.ETOH_nonMamm_2_1.5bfc9e6f196ab268db3fd1cddf25c56e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.ETOH_nonMamm_2_1.5bfc9e6f196ab268db3fd1cddf25c56e.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2430.006.Index_10.A549_ETOH_nonMAMM.bam \
 FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2430.006.Index_10.A549_ETOH_nonMAMM.pair1.fastq.gz \
  SECOND_END_FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2430.006.Index_10.A549_ETOH_nonMAMM.pair2.fastq.gz
picard_sam_to_fastq.ETOH_nonMamm_2_1.5bfc9e6f196ab268db3fd1cddf25c56e.mugqic.done
)
picard_sam_to_fastq_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m | grep "[0-9]")
echo "$picard_sam_to_fastq_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_12_JOB_ID: picard_sam_to_fastq.ETOH_shMED1_1_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.ETOH_shMED1_1_1
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.ETOH_shMED1_1_1.62548c489b5c31483b90bd3a6d9d8c75.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.ETOH_shMED1_1_1.62548c489b5c31483b90bd3a6d9d8c75.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.002.Index_10.SB_A549-shMED1-2.bam \
 FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.002.Index_10.SB_A549-shMED1-2.pair1.fastq.gz \
  SECOND_END_FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.002.Index_10.SB_A549-shMED1-2.pair2.fastq.gz
picard_sam_to_fastq.ETOH_shMED1_1_1.62548c489b5c31483b90bd3a6d9d8c75.mugqic.done
)
picard_sam_to_fastq_12_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m | grep "[0-9]")
echo "$picard_sam_to_fastq_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_13_JOB_ID: picard_sam_to_fastq.ETOH_shMED1_1_2
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.ETOH_shMED1_1_2
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.ETOH_shMED1_1_2.4dc54688cf0f9d868138d6012771913f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.ETOH_shMED1_1_2.4dc54688cf0f9d868138d6012771913f.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.003.Index_10.SB_A549-shMED1-2.bam \
 FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.003.Index_10.SB_A549-shMED1-2.pair1.fastq.gz \
  SECOND_END_FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.003.Index_10.SB_A549-shMED1-2.pair2.fastq.gz
picard_sam_to_fastq.ETOH_shMED1_1_2.4dc54688cf0f9d868138d6012771913f.mugqic.done
)
picard_sam_to_fastq_13_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m | grep "[0-9]")
echo "$picard_sam_to_fastq_13_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_14_JOB_ID: picard_sam_to_fastq.ETOH_shNIPBL_2_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.ETOH_shNIPBL_2_1
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.ETOH_shNIPBL_2_1.62748a81123807288dafb6e52136d32c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.ETOH_shNIPBL_2_1.62748a81123807288dafb6e52136d32c.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2430.006.Index_11.A549_ETOH_shMED1-2.bam \
 FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2430.006.Index_11.A549_ETOH_shMED1-2.pair1.fastq.gz \
  SECOND_END_FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2430.006.Index_11.A549_ETOH_shMED1-2.pair2.fastq.gz
picard_sam_to_fastq.ETOH_shNIPBL_2_1.62748a81123807288dafb6e52136d32c.mugqic.done
)
picard_sam_to_fastq_14_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m | grep "[0-9]")
echo "$picard_sam_to_fastq_14_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_15_JOB_ID: picard_sam_to_fastq.ETOH_shNIPBL_1_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.ETOH_shNIPBL_1_1
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.ETOH_shNIPBL_1_1.be83c123dad247e5050f76b8f231e7a6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.ETOH_shNIPBL_1_1.be83c123dad247e5050f76b8f231e7a6.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.002.Index_8.SB_A549-shNIPBL-3.bam \
 FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.002.Index_8.SB_A549-shNIPBL-3.pair1.fastq.gz \
  SECOND_END_FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.002.Index_8.SB_A549-shNIPBL-3.pair2.fastq.gz
picard_sam_to_fastq.ETOH_shNIPBL_1_1.be83c123dad247e5050f76b8f231e7a6.mugqic.done
)
picard_sam_to_fastq_15_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m | grep "[0-9]")
echo "$picard_sam_to_fastq_15_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_16_JOB_ID: picard_sam_to_fastq.ETOH_shNIPBL_1_2
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.ETOH_shNIPBL_1_2
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.ETOH_shNIPBL_1_2.8cd824686331576d6436ea7fabfb56e4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.ETOH_shNIPBL_1_2.8cd824686331576d6436ea7fabfb56e4.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.003.Index_8.SB_A549-shNIPBL-3.bam \
 FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.003.Index_8.SB_A549-shNIPBL-3.pair1.fastq.gz \
  SECOND_END_FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.003.Index_8.SB_A549-shNIPBL-3.pair2.fastq.gz
picard_sam_to_fastq.ETOH_shNIPBL_1_2.8cd824686331576d6436ea7fabfb56e4.mugqic.done
)
picard_sam_to_fastq_16_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m | grep "[0-9]")
echo "$picard_sam_to_fastq_16_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_17_JOB_ID: picard_sam_to_fastq.ETOH_shMED1_2_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.ETOH_shMED1_2_1
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.ETOH_shMED1_2_1.08ec08b70e2630d28eff5c81eadb93bb.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.ETOH_shMED1_2_1.08ec08b70e2630d28eff5c81eadb93bb.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2430.006.Index_20.A549_ETOH_shNIPBL-3.bam \
 FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2430.006.Index_20.A549_ETOH_shNIPBL-3.pair1.fastq.gz \
  SECOND_END_FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2430.006.Index_20.A549_ETOH_shNIPBL-3.pair2.fastq.gz
picard_sam_to_fastq.ETOH_shMED1_2_1.08ec08b70e2630d28eff5c81eadb93bb.mugqic.done
)
picard_sam_to_fastq_17_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m | grep "[0-9]")
echo "$picard_sam_to_fastq_17_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_18_JOB_ID: picard_sam_to_fastq.ETOH_shSMC1A_1_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.ETOH_shSMC1A_1_1
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.ETOH_shSMC1A_1_1.a01bd281b5f19e2012094eeba8c94195.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.ETOH_shSMC1A_1_1.a01bd281b5f19e2012094eeba8c94195.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.002.Index_11.SB_A549-shSMC1A-3.bam \
 FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.002.Index_11.SB_A549-shSMC1A-3.pair1.fastq.gz \
  SECOND_END_FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.002.Index_11.SB_A549-shSMC1A-3.pair2.fastq.gz
picard_sam_to_fastq.ETOH_shSMC1A_1_1.a01bd281b5f19e2012094eeba8c94195.mugqic.done
)
picard_sam_to_fastq_18_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m | grep "[0-9]")
echo "$picard_sam_to_fastq_18_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_19_JOB_ID: picard_sam_to_fastq.ETOH_shSMC1A_1_2
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.ETOH_shSMC1A_1_2
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.ETOH_shSMC1A_1_2.ab6a63b72c448d5d700936a9b3fd044f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.ETOH_shSMC1A_1_2.ab6a63b72c448d5d700936a9b3fd044f.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.003.Index_11.SB_A549-shSMC1A-3.bam \
 FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.003.Index_11.SB_A549-shSMC1A-3.pair1.fastq.gz \
  SECOND_END_FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.003.Index_11.SB_A549-shSMC1A-3.pair2.fastq.gz
picard_sam_to_fastq.ETOH_shSMC1A_1_2.ab6a63b72c448d5d700936a9b3fd044f.mugqic.done
)
picard_sam_to_fastq_19_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m | grep "[0-9]")
echo "$picard_sam_to_fastq_19_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_20_JOB_ID: picard_sam_to_fastq.ETOH_shSMC1A_2_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.ETOH_shSMC1A_2_1
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.ETOH_shSMC1A_2_1.4d831af28230ae9de0433b8d73e28242.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.ETOH_shSMC1A_2_1.4d831af28230ae9de0433b8d73e28242.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2430.006.Index_22.A549_ETOH_shSMC1A-2.bam \
 FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2430.006.Index_22.A549_ETOH_shSMC1A-2.pair1.fastq.gz \
  SECOND_END_FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2430.006.Index_22.A549_ETOH_shSMC1A-2.pair2.fastq.gz
picard_sam_to_fastq.ETOH_shSMC1A_2_1.4d831af28230ae9de0433b8d73e28242.mugqic.done
)
picard_sam_to_fastq_20_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m | grep "[0-9]")
echo "$picard_sam_to_fastq_20_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: trimmomatic
#-------------------------------------------------------------------------------
STEP=trimmomatic
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: trimmomatic_1_JOB_ID: trimmomatic.DEX_nonMamm_1_1
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.DEX_nonMamm_1_1
JOB_DEPENDENCIES=$picard_sam_to_fastq_1_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.DEX_nonMamm_1_1.e87a13043021ef940795a5723fc6acb4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.DEX_nonMamm_1_1.e87a13043021ef940795a5723fc6acb4.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/DEX_nonMamm_1 && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR PE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.002.Index_3.SB_A549-DEX-nonMamm.pair1.fastq.gz \
  /gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.002.Index_3.SB_A549-DEX-nonMamm.pair2.fastq.gz \
  trim/DEX_nonMamm_1/DEX_nonMamm_1_1.trim.pair1.fastq.gz \
  trim/DEX_nonMamm_1/DEX_nonMamm_1_1.trim.single1.fastq.gz \
  trim/DEX_nonMamm_1/DEX_nonMamm_1_1.trim.pair2.fastq.gz \
  trim/DEX_nonMamm_1/DEX_nonMamm_1_1.trim.single2.fastq.gz \
  ILLUMINACLIP:/cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.0/bfx/adapters-truseq.fa:2:30:15:8:true \
  TRAILING:30 \
  MINLEN:32 \
  2> trim/DEX_nonMamm_1/DEX_nonMamm_1_1.trim.log
trimmomatic.DEX_nonMamm_1_1.e87a13043021ef940795a5723fc6acb4.mugqic.done
)
trimmomatic_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_2_JOB_ID: trimmomatic.DEX_nonMamm_1_2
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.DEX_nonMamm_1_2
JOB_DEPENDENCIES=$picard_sam_to_fastq_2_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.DEX_nonMamm_1_2.6b18475e9cbad6e1317c53c72997686c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.DEX_nonMamm_1_2.6b18475e9cbad6e1317c53c72997686c.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/DEX_nonMamm_1 && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR PE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.003.Index_3.SB_A549-DEX-nonMamm.pair1.fastq.gz \
  /gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.003.Index_3.SB_A549-DEX-nonMamm.pair2.fastq.gz \
  trim/DEX_nonMamm_1/DEX_nonMamm_1_2.trim.pair1.fastq.gz \
  trim/DEX_nonMamm_1/DEX_nonMamm_1_2.trim.single1.fastq.gz \
  trim/DEX_nonMamm_1/DEX_nonMamm_1_2.trim.pair2.fastq.gz \
  trim/DEX_nonMamm_1/DEX_nonMamm_1_2.trim.single2.fastq.gz \
  ILLUMINACLIP:/cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.0/bfx/adapters-truseq.fa:2:30:15:8:true \
  TRAILING:30 \
  MINLEN:32 \
  2> trim/DEX_nonMamm_1/DEX_nonMamm_1_2.trim.log
trimmomatic.DEX_nonMamm_1_2.6b18475e9cbad6e1317c53c72997686c.mugqic.done
)
trimmomatic_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_3_JOB_ID: trimmomatic.DEX_shMED1_1_1
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.DEX_shMED1_1_1
JOB_DEPENDENCIES=$picard_sam_to_fastq_3_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.DEX_shMED1_1_1.c78e10c0ad286ec04904214714f2c9aa.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.DEX_shMED1_1_1.c78e10c0ad286ec04904214714f2c9aa.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/DEX_shMED1_1 && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR PE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.002.Index_22.SB_A549-DEX-shMED1-2.pair1.fastq.gz \
  /gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.002.Index_22.SB_A549-DEX-shMED1-2.pair2.fastq.gz \
  trim/DEX_shMED1_1/DEX_shMED1_1_1.trim.pair1.fastq.gz \
  trim/DEX_shMED1_1/DEX_shMED1_1_1.trim.single1.fastq.gz \
  trim/DEX_shMED1_1/DEX_shMED1_1_1.trim.pair2.fastq.gz \
  trim/DEX_shMED1_1/DEX_shMED1_1_1.trim.single2.fastq.gz \
  ILLUMINACLIP:/cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.0/bfx/adapters-truseq.fa:2:30:15:8:true \
  TRAILING:30 \
  MINLEN:32 \
  2> trim/DEX_shMED1_1/DEX_shMED1_1_1.trim.log
trimmomatic.DEX_shMED1_1_1.c78e10c0ad286ec04904214714f2c9aa.mugqic.done
)
trimmomatic_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_4_JOB_ID: trimmomatic.DEX_shMED1_1_2
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.DEX_shMED1_1_2
JOB_DEPENDENCIES=$picard_sam_to_fastq_4_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.DEX_shMED1_1_2.9866c70efce412febc1b3029c742a205.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.DEX_shMED1_1_2.9866c70efce412febc1b3029c742a205.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/DEX_shMED1_1 && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR PE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.003.Index_22.SB_A549-DEX-shMED1-2.pair1.fastq.gz \
  /gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.003.Index_22.SB_A549-DEX-shMED1-2.pair2.fastq.gz \
  trim/DEX_shMED1_1/DEX_shMED1_1_2.trim.pair1.fastq.gz \
  trim/DEX_shMED1_1/DEX_shMED1_1_2.trim.single1.fastq.gz \
  trim/DEX_shMED1_1/DEX_shMED1_1_2.trim.pair2.fastq.gz \
  trim/DEX_shMED1_1/DEX_shMED1_1_2.trim.single2.fastq.gz \
  ILLUMINACLIP:/cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.0/bfx/adapters-truseq.fa:2:30:15:8:true \
  TRAILING:30 \
  MINLEN:32 \
  2> trim/DEX_shMED1_1/DEX_shMED1_1_2.trim.log
trimmomatic.DEX_shMED1_1_2.9866c70efce412febc1b3029c742a205.mugqic.done
)
trimmomatic_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_5_JOB_ID: trimmomatic.DEX_shNIPBL_1_1
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.DEX_shNIPBL_1_1
JOB_DEPENDENCIES=$picard_sam_to_fastq_5_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.DEX_shNIPBL_1_1.7262dab6fe675e4d9f3e9cc676db98a9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.DEX_shNIPBL_1_1.7262dab6fe675e4d9f3e9cc676db98a9.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/DEX_shNIPBL_1 && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR PE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.002.Index_9.SB_A549-DEX-shNIPBL-3.pair1.fastq.gz \
  /gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.002.Index_9.SB_A549-DEX-shNIPBL-3.pair2.fastq.gz \
  trim/DEX_shNIPBL_1/DEX_shNIPBL_1_1.trim.pair1.fastq.gz \
  trim/DEX_shNIPBL_1/DEX_shNIPBL_1_1.trim.single1.fastq.gz \
  trim/DEX_shNIPBL_1/DEX_shNIPBL_1_1.trim.pair2.fastq.gz \
  trim/DEX_shNIPBL_1/DEX_shNIPBL_1_1.trim.single2.fastq.gz \
  ILLUMINACLIP:/cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.0/bfx/adapters-truseq.fa:2:30:15:8:true \
  TRAILING:30 \
  MINLEN:32 \
  2> trim/DEX_shNIPBL_1/DEX_shNIPBL_1_1.trim.log
trimmomatic.DEX_shNIPBL_1_1.7262dab6fe675e4d9f3e9cc676db98a9.mugqic.done
)
trimmomatic_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_6_JOB_ID: trimmomatic.DEX_shNIPBL_1_2
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.DEX_shNIPBL_1_2
JOB_DEPENDENCIES=$picard_sam_to_fastq_6_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.DEX_shNIPBL_1_2.bcf4bbdb5ff31014c0664d4fc2b5fe5d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.DEX_shNIPBL_1_2.bcf4bbdb5ff31014c0664d4fc2b5fe5d.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/DEX_shNIPBL_1 && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR PE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.003.Index_9.SB_A549-DEX-shNIPBL-3.pair1.fastq.gz \
  /gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.003.Index_9.SB_A549-DEX-shNIPBL-3.pair2.fastq.gz \
  trim/DEX_shNIPBL_1/DEX_shNIPBL_1_2.trim.pair1.fastq.gz \
  trim/DEX_shNIPBL_1/DEX_shNIPBL_1_2.trim.single1.fastq.gz \
  trim/DEX_shNIPBL_1/DEX_shNIPBL_1_2.trim.pair2.fastq.gz \
  trim/DEX_shNIPBL_1/DEX_shNIPBL_1_2.trim.single2.fastq.gz \
  ILLUMINACLIP:/cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.0/bfx/adapters-truseq.fa:2:30:15:8:true \
  TRAILING:30 \
  MINLEN:32 \
  2> trim/DEX_shNIPBL_1/DEX_shNIPBL_1_2.trim.log
trimmomatic.DEX_shNIPBL_1_2.bcf4bbdb5ff31014c0664d4fc2b5fe5d.mugqic.done
)
trimmomatic_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_7_JOB_ID: trimmomatic.DEX_shSMC1A_1_1
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.DEX_shSMC1A_1_1
JOB_DEPENDENCIES=$picard_sam_to_fastq_7_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.DEX_shSMC1A_1_1.b8343ca7644549b88568a6b544b08728.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.DEX_shSMC1A_1_1.b8343ca7644549b88568a6b544b08728.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/DEX_shSMC1A_3 && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR PE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.002.Index_27.SB_A549-DEX-shSMC1A-3.pair1.fastq.gz \
  /gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.002.Index_27.SB_A549-DEX-shSMC1A-3.pair2.fastq.gz \
  trim/DEX_shSMC1A_3/DEX_shSMC1A_1_1.trim.pair1.fastq.gz \
  trim/DEX_shSMC1A_3/DEX_shSMC1A_1_1.trim.single1.fastq.gz \
  trim/DEX_shSMC1A_3/DEX_shSMC1A_1_1.trim.pair2.fastq.gz \
  trim/DEX_shSMC1A_3/DEX_shSMC1A_1_1.trim.single2.fastq.gz \
  ILLUMINACLIP:/cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.0/bfx/adapters-truseq.fa:2:30:15:8:true \
  TRAILING:30 \
  MINLEN:32 \
  2> trim/DEX_shSMC1A_3/DEX_shSMC1A_1_1.trim.log
trimmomatic.DEX_shSMC1A_1_1.b8343ca7644549b88568a6b544b08728.mugqic.done
)
trimmomatic_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_8_JOB_ID: trimmomatic.DEX_shSMC1A_1_2
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.DEX_shSMC1A_1_2
JOB_DEPENDENCIES=$picard_sam_to_fastq_8_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.DEX_shSMC1A_1_2.a6ac01d7acf6874b107b8e10a5f7defb.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.DEX_shSMC1A_1_2.a6ac01d7acf6874b107b8e10a5f7defb.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/DEX_shSMC1A_3 && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR PE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.003.Index_27.SB_A549-DEX-shSMC1A-3.pair1.fastq.gz \
  /gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.003.Index_27.SB_A549-DEX-shSMC1A-3.pair2.fastq.gz \
  trim/DEX_shSMC1A_3/DEX_shSMC1A_1_2.trim.pair1.fastq.gz \
  trim/DEX_shSMC1A_3/DEX_shSMC1A_1_2.trim.single1.fastq.gz \
  trim/DEX_shSMC1A_3/DEX_shSMC1A_1_2.trim.pair2.fastq.gz \
  trim/DEX_shSMC1A_3/DEX_shSMC1A_1_2.trim.single2.fastq.gz \
  ILLUMINACLIP:/cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.0/bfx/adapters-truseq.fa:2:30:15:8:true \
  TRAILING:30 \
  MINLEN:32 \
  2> trim/DEX_shSMC1A_3/DEX_shSMC1A_1_2.trim.log
trimmomatic.DEX_shSMC1A_1_2.a6ac01d7acf6874b107b8e10a5f7defb.mugqic.done
)
trimmomatic_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_9_JOB_ID: trimmomatic.ETOH_nonMamm_1_1
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.ETOH_nonMamm_1_1
JOB_DEPENDENCIES=$picard_sam_to_fastq_9_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.ETOH_nonMamm_1_1.bdff3f8a605ff92b82e2d857036ce79b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.ETOH_nonMamm_1_1.bdff3f8a605ff92b82e2d857036ce79b.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/ETOH_nonMamm_1 && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR PE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.002.Index_1.SB_A549-nonMamm.pair1.fastq.gz \
  /gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.002.Index_1.SB_A549-nonMamm.pair2.fastq.gz \
  trim/ETOH_nonMamm_1/ETOH_nonMamm_1_1.trim.pair1.fastq.gz \
  trim/ETOH_nonMamm_1/ETOH_nonMamm_1_1.trim.single1.fastq.gz \
  trim/ETOH_nonMamm_1/ETOH_nonMamm_1_1.trim.pair2.fastq.gz \
  trim/ETOH_nonMamm_1/ETOH_nonMamm_1_1.trim.single2.fastq.gz \
  ILLUMINACLIP:/cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.0/bfx/adapters-truseq.fa:2:30:15:8:true \
  TRAILING:30 \
  MINLEN:32 \
  2> trim/ETOH_nonMamm_1/ETOH_nonMamm_1_1.trim.log
trimmomatic.ETOH_nonMamm_1_1.bdff3f8a605ff92b82e2d857036ce79b.mugqic.done
)
trimmomatic_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_10_JOB_ID: trimmomatic.ETOH_nonMamm_1_2
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.ETOH_nonMamm_1_2
JOB_DEPENDENCIES=$picard_sam_to_fastq_10_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.ETOH_nonMamm_1_2.0a414645db5134a57d1e589f5191d619.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.ETOH_nonMamm_1_2.0a414645db5134a57d1e589f5191d619.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/ETOH_nonMamm_1 && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR PE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.003.Index_1.SB_A549-nonMamm.pair1.fastq.gz \
  /gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.003.Index_1.SB_A549-nonMamm.pair2.fastq.gz \
  trim/ETOH_nonMamm_1/ETOH_nonMamm_1_2.trim.pair1.fastq.gz \
  trim/ETOH_nonMamm_1/ETOH_nonMamm_1_2.trim.single1.fastq.gz \
  trim/ETOH_nonMamm_1/ETOH_nonMamm_1_2.trim.pair2.fastq.gz \
  trim/ETOH_nonMamm_1/ETOH_nonMamm_1_2.trim.single2.fastq.gz \
  ILLUMINACLIP:/cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.0/bfx/adapters-truseq.fa:2:30:15:8:true \
  TRAILING:30 \
  MINLEN:32 \
  2> trim/ETOH_nonMamm_1/ETOH_nonMamm_1_2.trim.log
trimmomatic.ETOH_nonMamm_1_2.0a414645db5134a57d1e589f5191d619.mugqic.done
)
trimmomatic_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_11_JOB_ID: trimmomatic.ETOH_nonMamm_2_1
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.ETOH_nonMamm_2_1
JOB_DEPENDENCIES=$picard_sam_to_fastq_11_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.ETOH_nonMamm_2_1.f38858bc51c84e22dca3aa00aeb5bc57.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.ETOH_nonMamm_2_1.f38858bc51c84e22dca3aa00aeb5bc57.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/ETOH_nonMamm_2 && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR PE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2430.006.Index_10.A549_ETOH_nonMAMM.pair1.fastq.gz \
  /gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2430.006.Index_10.A549_ETOH_nonMAMM.pair2.fastq.gz \
  trim/ETOH_nonMamm_2/ETOH_nonMamm_2_1.trim.pair1.fastq.gz \
  trim/ETOH_nonMamm_2/ETOH_nonMamm_2_1.trim.single1.fastq.gz \
  trim/ETOH_nonMamm_2/ETOH_nonMamm_2_1.trim.pair2.fastq.gz \
  trim/ETOH_nonMamm_2/ETOH_nonMamm_2_1.trim.single2.fastq.gz \
  ILLUMINACLIP:/cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.0/bfx/adapters-truseq.fa:2:30:15:8:true \
  TRAILING:30 \
  MINLEN:32 \
  2> trim/ETOH_nonMamm_2/ETOH_nonMamm_2_1.trim.log
trimmomatic.ETOH_nonMamm_2_1.f38858bc51c84e22dca3aa00aeb5bc57.mugqic.done
)
trimmomatic_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_12_JOB_ID: trimmomatic.ETOH_shMED1_1_1
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.ETOH_shMED1_1_1
JOB_DEPENDENCIES=$picard_sam_to_fastq_12_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.ETOH_shMED1_1_1.9cbb60ea1bb776f21208d938b3936841.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.ETOH_shMED1_1_1.9cbb60ea1bb776f21208d938b3936841.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/ETOH_shMED1_1 && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR PE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.002.Index_10.SB_A549-shMED1-2.pair1.fastq.gz \
  /gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.002.Index_10.SB_A549-shMED1-2.pair2.fastq.gz \
  trim/ETOH_shMED1_1/ETOH_shMED1_1_1.trim.pair1.fastq.gz \
  trim/ETOH_shMED1_1/ETOH_shMED1_1_1.trim.single1.fastq.gz \
  trim/ETOH_shMED1_1/ETOH_shMED1_1_1.trim.pair2.fastq.gz \
  trim/ETOH_shMED1_1/ETOH_shMED1_1_1.trim.single2.fastq.gz \
  ILLUMINACLIP:/cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.0/bfx/adapters-truseq.fa:2:30:15:8:true \
  TRAILING:30 \
  MINLEN:32 \
  2> trim/ETOH_shMED1_1/ETOH_shMED1_1_1.trim.log
trimmomatic.ETOH_shMED1_1_1.9cbb60ea1bb776f21208d938b3936841.mugqic.done
)
trimmomatic_12_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_13_JOB_ID: trimmomatic.ETOH_shMED1_1_2
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.ETOH_shMED1_1_2
JOB_DEPENDENCIES=$picard_sam_to_fastq_13_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.ETOH_shMED1_1_2.802125a665090108b967e24e3e85ad9e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.ETOH_shMED1_1_2.802125a665090108b967e24e3e85ad9e.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/ETOH_shMED1_1 && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR PE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.003.Index_10.SB_A549-shMED1-2.pair1.fastq.gz \
  /gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.003.Index_10.SB_A549-shMED1-2.pair2.fastq.gz \
  trim/ETOH_shMED1_1/ETOH_shMED1_1_2.trim.pair1.fastq.gz \
  trim/ETOH_shMED1_1/ETOH_shMED1_1_2.trim.single1.fastq.gz \
  trim/ETOH_shMED1_1/ETOH_shMED1_1_2.trim.pair2.fastq.gz \
  trim/ETOH_shMED1_1/ETOH_shMED1_1_2.trim.single2.fastq.gz \
  ILLUMINACLIP:/cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.0/bfx/adapters-truseq.fa:2:30:15:8:true \
  TRAILING:30 \
  MINLEN:32 \
  2> trim/ETOH_shMED1_1/ETOH_shMED1_1_2.trim.log
trimmomatic.ETOH_shMED1_1_2.802125a665090108b967e24e3e85ad9e.mugqic.done
)
trimmomatic_13_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_13_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_14_JOB_ID: trimmomatic.ETOH_shNIPBL_2_1
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.ETOH_shNIPBL_2_1
JOB_DEPENDENCIES=$picard_sam_to_fastq_14_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.ETOH_shNIPBL_2_1.9ff84e22b4e4f09a996d8cac002748c6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.ETOH_shNIPBL_2_1.9ff84e22b4e4f09a996d8cac002748c6.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/ETOH_shNIPBL_2 && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR PE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2430.006.Index_11.A549_ETOH_shMED1-2.pair1.fastq.gz \
  /gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2430.006.Index_11.A549_ETOH_shMED1-2.pair2.fastq.gz \
  trim/ETOH_shNIPBL_2/ETOH_shNIPBL_2_1.trim.pair1.fastq.gz \
  trim/ETOH_shNIPBL_2/ETOH_shNIPBL_2_1.trim.single1.fastq.gz \
  trim/ETOH_shNIPBL_2/ETOH_shNIPBL_2_1.trim.pair2.fastq.gz \
  trim/ETOH_shNIPBL_2/ETOH_shNIPBL_2_1.trim.single2.fastq.gz \
  ILLUMINACLIP:/cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.0/bfx/adapters-truseq.fa:2:30:15:8:true \
  TRAILING:30 \
  MINLEN:32 \
  2> trim/ETOH_shNIPBL_2/ETOH_shNIPBL_2_1.trim.log
trimmomatic.ETOH_shNIPBL_2_1.9ff84e22b4e4f09a996d8cac002748c6.mugqic.done
)
trimmomatic_14_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_14_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_15_JOB_ID: trimmomatic.ETOH_shNIPBL_1_1
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.ETOH_shNIPBL_1_1
JOB_DEPENDENCIES=$picard_sam_to_fastq_15_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.ETOH_shNIPBL_1_1.5e85ce7ab22bfc868d7145f1e280be25.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.ETOH_shNIPBL_1_1.5e85ce7ab22bfc868d7145f1e280be25.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/ETOH_shNIPBL_1 && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR PE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.002.Index_8.SB_A549-shNIPBL-3.pair1.fastq.gz \
  /gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.002.Index_8.SB_A549-shNIPBL-3.pair2.fastq.gz \
  trim/ETOH_shNIPBL_1/ETOH_shNIPBL_1_1.trim.pair1.fastq.gz \
  trim/ETOH_shNIPBL_1/ETOH_shNIPBL_1_1.trim.single1.fastq.gz \
  trim/ETOH_shNIPBL_1/ETOH_shNIPBL_1_1.trim.pair2.fastq.gz \
  trim/ETOH_shNIPBL_1/ETOH_shNIPBL_1_1.trim.single2.fastq.gz \
  ILLUMINACLIP:/cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.0/bfx/adapters-truseq.fa:2:30:15:8:true \
  TRAILING:30 \
  MINLEN:32 \
  2> trim/ETOH_shNIPBL_1/ETOH_shNIPBL_1_1.trim.log
trimmomatic.ETOH_shNIPBL_1_1.5e85ce7ab22bfc868d7145f1e280be25.mugqic.done
)
trimmomatic_15_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_15_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_16_JOB_ID: trimmomatic.ETOH_shNIPBL_1_2
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.ETOH_shNIPBL_1_2
JOB_DEPENDENCIES=$picard_sam_to_fastq_16_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.ETOH_shNIPBL_1_2.cc01175ada6eaafc476c64897e52e217.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.ETOH_shNIPBL_1_2.cc01175ada6eaafc476c64897e52e217.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/ETOH_shNIPBL_1 && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR PE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.003.Index_8.SB_A549-shNIPBL-3.pair1.fastq.gz \
  /gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.003.Index_8.SB_A549-shNIPBL-3.pair2.fastq.gz \
  trim/ETOH_shNIPBL_1/ETOH_shNIPBL_1_2.trim.pair1.fastq.gz \
  trim/ETOH_shNIPBL_1/ETOH_shNIPBL_1_2.trim.single1.fastq.gz \
  trim/ETOH_shNIPBL_1/ETOH_shNIPBL_1_2.trim.pair2.fastq.gz \
  trim/ETOH_shNIPBL_1/ETOH_shNIPBL_1_2.trim.single2.fastq.gz \
  ILLUMINACLIP:/cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.0/bfx/adapters-truseq.fa:2:30:15:8:true \
  TRAILING:30 \
  MINLEN:32 \
  2> trim/ETOH_shNIPBL_1/ETOH_shNIPBL_1_2.trim.log
trimmomatic.ETOH_shNIPBL_1_2.cc01175ada6eaafc476c64897e52e217.mugqic.done
)
trimmomatic_16_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_16_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_17_JOB_ID: trimmomatic.ETOH_shMED1_2_1
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.ETOH_shMED1_2_1
JOB_DEPENDENCIES=$picard_sam_to_fastq_17_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.ETOH_shMED1_2_1.109b08b015756ca2368e264a05d5d0fc.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.ETOH_shMED1_2_1.109b08b015756ca2368e264a05d5d0fc.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/ETOH_shMED1_2 && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR PE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2430.006.Index_20.A549_ETOH_shNIPBL-3.pair1.fastq.gz \
  /gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2430.006.Index_20.A549_ETOH_shNIPBL-3.pair2.fastq.gz \
  trim/ETOH_shMED1_2/ETOH_shMED1_2_1.trim.pair1.fastq.gz \
  trim/ETOH_shMED1_2/ETOH_shMED1_2_1.trim.single1.fastq.gz \
  trim/ETOH_shMED1_2/ETOH_shMED1_2_1.trim.pair2.fastq.gz \
  trim/ETOH_shMED1_2/ETOH_shMED1_2_1.trim.single2.fastq.gz \
  ILLUMINACLIP:/cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.0/bfx/adapters-truseq.fa:2:30:15:8:true \
  TRAILING:30 \
  MINLEN:32 \
  2> trim/ETOH_shMED1_2/ETOH_shMED1_2_1.trim.log
trimmomatic.ETOH_shMED1_2_1.109b08b015756ca2368e264a05d5d0fc.mugqic.done
)
trimmomatic_17_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_17_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_18_JOB_ID: trimmomatic.ETOH_shSMC1A_1_1
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.ETOH_shSMC1A_1_1
JOB_DEPENDENCIES=$picard_sam_to_fastq_18_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.ETOH_shSMC1A_1_1.099da5358b453d64ae85e8662b8ee149.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.ETOH_shSMC1A_1_1.099da5358b453d64ae85e8662b8ee149.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/ETOH_shSMC1A_3 && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR PE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.002.Index_11.SB_A549-shSMC1A-3.pair1.fastq.gz \
  /gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.002.Index_11.SB_A549-shSMC1A-3.pair2.fastq.gz \
  trim/ETOH_shSMC1A_3/ETOH_shSMC1A_1_1.trim.pair1.fastq.gz \
  trim/ETOH_shSMC1A_3/ETOH_shSMC1A_1_1.trim.single1.fastq.gz \
  trim/ETOH_shSMC1A_3/ETOH_shSMC1A_1_1.trim.pair2.fastq.gz \
  trim/ETOH_shSMC1A_3/ETOH_shSMC1A_1_1.trim.single2.fastq.gz \
  ILLUMINACLIP:/cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.0/bfx/adapters-truseq.fa:2:30:15:8:true \
  TRAILING:30 \
  MINLEN:32 \
  2> trim/ETOH_shSMC1A_3/ETOH_shSMC1A_1_1.trim.log
trimmomatic.ETOH_shSMC1A_1_1.099da5358b453d64ae85e8662b8ee149.mugqic.done
)
trimmomatic_18_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_18_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_19_JOB_ID: trimmomatic.ETOH_shSMC1A_1_2
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.ETOH_shSMC1A_1_2
JOB_DEPENDENCIES=$picard_sam_to_fastq_19_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.ETOH_shSMC1A_1_2.31626f369a43eb65578b9c19f07923b0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.ETOH_shSMC1A_1_2.31626f369a43eb65578b9c19f07923b0.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/ETOH_shSMC1A_3 && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR PE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.003.Index_11.SB_A549-shSMC1A-3.pair1.fastq.gz \
  /gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2031.003.Index_11.SB_A549-shSMC1A-3.pair2.fastq.gz \
  trim/ETOH_shSMC1A_3/ETOH_shSMC1A_1_2.trim.pair1.fastq.gz \
  trim/ETOH_shSMC1A_3/ETOH_shSMC1A_1_2.trim.single1.fastq.gz \
  trim/ETOH_shSMC1A_3/ETOH_shSMC1A_1_2.trim.pair2.fastq.gz \
  trim/ETOH_shSMC1A_3/ETOH_shSMC1A_1_2.trim.single2.fastq.gz \
  ILLUMINACLIP:/cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.0/bfx/adapters-truseq.fa:2:30:15:8:true \
  TRAILING:30 \
  MINLEN:32 \
  2> trim/ETOH_shSMC1A_3/ETOH_shSMC1A_1_2.trim.log
trimmomatic.ETOH_shSMC1A_1_2.31626f369a43eb65578b9c19f07923b0.mugqic.done
)
trimmomatic_19_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_19_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_20_JOB_ID: trimmomatic.ETOH_shSMC1A_2_1
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.ETOH_shSMC1A_2_1
JOB_DEPENDENCIES=$picard_sam_to_fastq_20_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.ETOH_shSMC1A_2_1.c41321bbc0db0fad6d83971fb50445af.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.ETOH_shSMC1A_2_1.c41321bbc0db0fad6d83971fb50445af.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/ETOH_shSMC1A_2 && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR PE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2430.006.Index_22.A549_ETOH_shSMC1A-2.pair1.fastq.gz \
  /gs/scratch/efournier/CofactorHR/A549/raw/rna-seq/HI.2430.006.Index_22.A549_ETOH_shSMC1A-2.pair2.fastq.gz \
  trim/ETOH_shSMC1A_2/ETOH_shSMC1A_2_1.trim.pair1.fastq.gz \
  trim/ETOH_shSMC1A_2/ETOH_shSMC1A_2_1.trim.single1.fastq.gz \
  trim/ETOH_shSMC1A_2/ETOH_shSMC1A_2_1.trim.pair2.fastq.gz \
  trim/ETOH_shSMC1A_2/ETOH_shSMC1A_2_1.trim.single2.fastq.gz \
  ILLUMINACLIP:/cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.0/bfx/adapters-truseq.fa:2:30:15:8:true \
  TRAILING:30 \
  MINLEN:32 \
  2> trim/ETOH_shSMC1A_2/ETOH_shSMC1A_2_1.trim.log
trimmomatic.ETOH_shSMC1A_2_1.c41321bbc0db0fad6d83971fb50445af.mugqic.done
)
trimmomatic_20_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_20_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: merge_trimmomatic_stats
#-------------------------------------------------------------------------------
STEP=merge_trimmomatic_stats
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: merge_trimmomatic_stats_1_JOB_ID: merge_trimmomatic_stats
#-------------------------------------------------------------------------------
JOB_NAME=merge_trimmomatic_stats
JOB_DEPENDENCIES=$trimmomatic_1_JOB_ID:$trimmomatic_2_JOB_ID:$trimmomatic_3_JOB_ID:$trimmomatic_4_JOB_ID:$trimmomatic_5_JOB_ID:$trimmomatic_6_JOB_ID:$trimmomatic_7_JOB_ID:$trimmomatic_8_JOB_ID:$trimmomatic_9_JOB_ID:$trimmomatic_10_JOB_ID:$trimmomatic_11_JOB_ID:$trimmomatic_12_JOB_ID:$trimmomatic_13_JOB_ID:$trimmomatic_14_JOB_ID:$trimmomatic_15_JOB_ID:$trimmomatic_16_JOB_ID:$trimmomatic_17_JOB_ID:$trimmomatic_18_JOB_ID:$trimmomatic_19_JOB_ID:$trimmomatic_20_JOB_ID
JOB_DONE=job_output/merge_trimmomatic_stats/merge_trimmomatic_stats.bbdf186a57c819ba050bbc6cce65fc2f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'merge_trimmomatic_stats.bbdf186a57c819ba050bbc6cce65fc2f.mugqic.done'
module load mugqic/pandoc/1.15.2 && \
mkdir -p metrics && \
echo 'Sample	Readset	Raw Paired Reads #	Surviving Paired Reads #	Surviving Paired Reads %' > metrics/trimReadsetTable.tsv && \
grep ^Input trim/DEX_nonMamm_1/DEX_nonMamm_1_1.trim.log | \
perl -pe 's/^Input Read Pairs: (\d+).*Both Surviving: (\d+).*Forward Only Surviving: (\d+).*$/DEX_nonMamm_1	DEX_nonMamm_1_1	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/DEX_nonMamm_1/DEX_nonMamm_1_2.trim.log | \
perl -pe 's/^Input Read Pairs: (\d+).*Both Surviving: (\d+).*Forward Only Surviving: (\d+).*$/DEX_nonMamm_1	DEX_nonMamm_1_2	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/DEX_shMED1_1/DEX_shMED1_1_1.trim.log | \
perl -pe 's/^Input Read Pairs: (\d+).*Both Surviving: (\d+).*Forward Only Surviving: (\d+).*$/DEX_shMED1_1	DEX_shMED1_1_1	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/DEX_shMED1_1/DEX_shMED1_1_2.trim.log | \
perl -pe 's/^Input Read Pairs: (\d+).*Both Surviving: (\d+).*Forward Only Surviving: (\d+).*$/DEX_shMED1_1	DEX_shMED1_1_2	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/DEX_shNIPBL_1/DEX_shNIPBL_1_1.trim.log | \
perl -pe 's/^Input Read Pairs: (\d+).*Both Surviving: (\d+).*Forward Only Surviving: (\d+).*$/DEX_shNIPBL_1	DEX_shNIPBL_1_1	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/DEX_shNIPBL_1/DEX_shNIPBL_1_2.trim.log | \
perl -pe 's/^Input Read Pairs: (\d+).*Both Surviving: (\d+).*Forward Only Surviving: (\d+).*$/DEX_shNIPBL_1	DEX_shNIPBL_1_2	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/DEX_shSMC1A_3/DEX_shSMC1A_1_1.trim.log | \
perl -pe 's/^Input Read Pairs: (\d+).*Both Surviving: (\d+).*Forward Only Surviving: (\d+).*$/DEX_shSMC1A_3	DEX_shSMC1A_1_1	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/DEX_shSMC1A_3/DEX_shSMC1A_1_2.trim.log | \
perl -pe 's/^Input Read Pairs: (\d+).*Both Surviving: (\d+).*Forward Only Surviving: (\d+).*$/DEX_shSMC1A_3	DEX_shSMC1A_1_2	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/ETOH_nonMamm_1/ETOH_nonMamm_1_1.trim.log | \
perl -pe 's/^Input Read Pairs: (\d+).*Both Surviving: (\d+).*Forward Only Surviving: (\d+).*$/ETOH_nonMamm_1	ETOH_nonMamm_1_1	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/ETOH_nonMamm_1/ETOH_nonMamm_1_2.trim.log | \
perl -pe 's/^Input Read Pairs: (\d+).*Both Surviving: (\d+).*Forward Only Surviving: (\d+).*$/ETOH_nonMamm_1	ETOH_nonMamm_1_2	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/ETOH_nonMamm_2/ETOH_nonMamm_2_1.trim.log | \
perl -pe 's/^Input Read Pairs: (\d+).*Both Surviving: (\d+).*Forward Only Surviving: (\d+).*$/ETOH_nonMamm_2	ETOH_nonMamm_2_1	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/ETOH_shMED1_1/ETOH_shMED1_1_1.trim.log | \
perl -pe 's/^Input Read Pairs: (\d+).*Both Surviving: (\d+).*Forward Only Surviving: (\d+).*$/ETOH_shMED1_1	ETOH_shMED1_1_1	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/ETOH_shMED1_1/ETOH_shMED1_1_2.trim.log | \
perl -pe 's/^Input Read Pairs: (\d+).*Both Surviving: (\d+).*Forward Only Surviving: (\d+).*$/ETOH_shMED1_1	ETOH_shMED1_1_2	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/ETOH_shNIPBL_2/ETOH_shNIPBL_2_1.trim.log | \
perl -pe 's/^Input Read Pairs: (\d+).*Both Surviving: (\d+).*Forward Only Surviving: (\d+).*$/ETOH_shNIPBL_2	ETOH_shNIPBL_2_1	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/ETOH_shNIPBL_1/ETOH_shNIPBL_1_1.trim.log | \
perl -pe 's/^Input Read Pairs: (\d+).*Both Surviving: (\d+).*Forward Only Surviving: (\d+).*$/ETOH_shNIPBL_1	ETOH_shNIPBL_1_1	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/ETOH_shNIPBL_1/ETOH_shNIPBL_1_2.trim.log | \
perl -pe 's/^Input Read Pairs: (\d+).*Both Surviving: (\d+).*Forward Only Surviving: (\d+).*$/ETOH_shNIPBL_1	ETOH_shNIPBL_1_2	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/ETOH_shMED1_2/ETOH_shMED1_2_1.trim.log | \
perl -pe 's/^Input Read Pairs: (\d+).*Both Surviving: (\d+).*Forward Only Surviving: (\d+).*$/ETOH_shMED1_2	ETOH_shMED1_2_1	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/ETOH_shSMC1A_3/ETOH_shSMC1A_1_1.trim.log | \
perl -pe 's/^Input Read Pairs: (\d+).*Both Surviving: (\d+).*Forward Only Surviving: (\d+).*$/ETOH_shSMC1A_3	ETOH_shSMC1A_1_1	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/ETOH_shSMC1A_3/ETOH_shSMC1A_1_2.trim.log | \
perl -pe 's/^Input Read Pairs: (\d+).*Both Surviving: (\d+).*Forward Only Surviving: (\d+).*$/ETOH_shSMC1A_3	ETOH_shSMC1A_1_2	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/ETOH_shSMC1A_2/ETOH_shSMC1A_2_1.trim.log | \
perl -pe 's/^Input Read Pairs: (\d+).*Both Surviving: (\d+).*Forward Only Surviving: (\d+).*$/ETOH_shSMC1A_2	ETOH_shSMC1A_2_1	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
cut -f1,3- metrics/trimReadsetTable.tsv | awk -F"	" '{OFS="	"; if (NR==1) {if ($2=="Raw Paired Reads #") {paired=1};print "Sample", "Raw Reads #", "Surviving Reads #", "Surviving %"} else {if (paired) {$2=$2*2; $3=$3*2}; raw[$1]+=$2; surviving[$1]+=$3}}END{for (sample in raw){print sample, raw[sample], surviving[sample], surviving[sample] / raw[sample] * 100}}' \
  > metrics/trimSampleTable.tsv && \
mkdir -p report && \
cp metrics/trimReadsetTable.tsv metrics/trimSampleTable.tsv report/ && \
trim_readset_table_md=`LC_NUMERIC=en_CA awk -F "	" '{OFS="|"; if (NR == 1) {$1 = $1; print $0; print "-----|-----|-----:|-----:|-----:"} else {print $1, $2, sprintf("%\47d", $3), sprintf("%\47d", $4), sprintf("%.1f", $5)}}' metrics/trimReadsetTable.tsv` && \
pandoc \
  /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.3.0/bfx/report/Illumina.merge_trimmomatic_stats.md \
  --template /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.3.0/bfx/report/Illumina.merge_trimmomatic_stats.md \
  --variable trailing_min_quality=30 \
  --variable min_length=32 \
  --variable read_type=Paired \
  --variable trim_readset_table="$trim_readset_table_md" \
  --to markdown \
  > report/Illumina.merge_trimmomatic_stats.md
merge_trimmomatic_stats.bbdf186a57c819ba050bbc6cce65fc2f.mugqic.done
)
merge_trimmomatic_stats_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$merge_trimmomatic_stats_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: star
#-------------------------------------------------------------------------------
STEP=star
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: star_1_JOB_ID: star_align.1.DEX_nonMamm_1_1
#-------------------------------------------------------------------------------
JOB_NAME=star_align.1.DEX_nonMamm_1_1
JOB_DEPENDENCIES=$trimmomatic_1_JOB_ID
JOB_DONE=job_output/star/star_align.1.DEX_nonMamm_1_1.b2d0578f95f5c90cd80c2914c829009b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.1.DEX_nonMamm_1_1.b2d0578f95f5c90cd80c2914c829009b.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment_1stPass/DEX_nonMamm_1/DEX_nonMamm_1_1 && \
STAR --runMode alignReads \
  --genomeDir /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/star_index/Ensembl87.sjdbOverhang99 \
  --readFilesIn \
    trim/DEX_nonMamm_1/DEX_nonMamm_1_1.trim.pair1.fastq.gz \
    trim/DEX_nonMamm_1/DEX_nonMamm_1_1.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix alignment_1stPass/DEX_nonMamm_1/DEX_nonMamm_1_1/ \
  --outSAMattrRGline ID:"DEX_nonMamm_1_1" 	PL:"ILLUMINA" 			SM:"DEX_nonMamm_1" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 70000000000 \
  --limitIObufferSize 1000000000
star_align.1.DEX_nonMamm_1_1.b2d0578f95f5c90cd80c2914c829009b.mugqic.done
)
star_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=5:00:0 -q metaq -l nodes=1:ppn=16 -l pmem=5700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$star_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: star_2_JOB_ID: star_align.1.DEX_nonMamm_1_2
#-------------------------------------------------------------------------------
JOB_NAME=star_align.1.DEX_nonMamm_1_2
JOB_DEPENDENCIES=$trimmomatic_2_JOB_ID
JOB_DONE=job_output/star/star_align.1.DEX_nonMamm_1_2.8a404c9d3b7f13b7701981e348b4d264.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.1.DEX_nonMamm_1_2.8a404c9d3b7f13b7701981e348b4d264.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment_1stPass/DEX_nonMamm_1/DEX_nonMamm_1_2 && \
STAR --runMode alignReads \
  --genomeDir /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/star_index/Ensembl87.sjdbOverhang99 \
  --readFilesIn \
    trim/DEX_nonMamm_1/DEX_nonMamm_1_2.trim.pair1.fastq.gz \
    trim/DEX_nonMamm_1/DEX_nonMamm_1_2.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix alignment_1stPass/DEX_nonMamm_1/DEX_nonMamm_1_2/ \
  --outSAMattrRGline ID:"DEX_nonMamm_1_2" 	PL:"ILLUMINA" 			SM:"DEX_nonMamm_1" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 70000000000 \
  --limitIObufferSize 1000000000
star_align.1.DEX_nonMamm_1_2.8a404c9d3b7f13b7701981e348b4d264.mugqic.done
)
star_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=5:00:0 -q metaq -l nodes=1:ppn=16 -l pmem=5700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$star_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: star_3_JOB_ID: star_align.1.DEX_shMED1_1_1
#-------------------------------------------------------------------------------
JOB_NAME=star_align.1.DEX_shMED1_1_1
JOB_DEPENDENCIES=$trimmomatic_3_JOB_ID
JOB_DONE=job_output/star/star_align.1.DEX_shMED1_1_1.7fb3e8933955a68d16038720c9a5e266.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.1.DEX_shMED1_1_1.7fb3e8933955a68d16038720c9a5e266.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment_1stPass/DEX_shMED1_1/DEX_shMED1_1_1 && \
STAR --runMode alignReads \
  --genomeDir /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/star_index/Ensembl87.sjdbOverhang99 \
  --readFilesIn \
    trim/DEX_shMED1_1/DEX_shMED1_1_1.trim.pair1.fastq.gz \
    trim/DEX_shMED1_1/DEX_shMED1_1_1.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix alignment_1stPass/DEX_shMED1_1/DEX_shMED1_1_1/ \
  --outSAMattrRGline ID:"DEX_shMED1_1_1" 	PL:"ILLUMINA" 			SM:"DEX_shMED1_1" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 70000000000 \
  --limitIObufferSize 1000000000
star_align.1.DEX_shMED1_1_1.7fb3e8933955a68d16038720c9a5e266.mugqic.done
)
star_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=5:00:0 -q metaq -l nodes=1:ppn=16 -l pmem=5700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$star_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: star_4_JOB_ID: star_align.1.DEX_shMED1_1_2
#-------------------------------------------------------------------------------
JOB_NAME=star_align.1.DEX_shMED1_1_2
JOB_DEPENDENCIES=$trimmomatic_4_JOB_ID
JOB_DONE=job_output/star/star_align.1.DEX_shMED1_1_2.24813d3bab959a0546a87211546cf7ce.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.1.DEX_shMED1_1_2.24813d3bab959a0546a87211546cf7ce.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment_1stPass/DEX_shMED1_1/DEX_shMED1_1_2 && \
STAR --runMode alignReads \
  --genomeDir /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/star_index/Ensembl87.sjdbOverhang99 \
  --readFilesIn \
    trim/DEX_shMED1_1/DEX_shMED1_1_2.trim.pair1.fastq.gz \
    trim/DEX_shMED1_1/DEX_shMED1_1_2.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix alignment_1stPass/DEX_shMED1_1/DEX_shMED1_1_2/ \
  --outSAMattrRGline ID:"DEX_shMED1_1_2" 	PL:"ILLUMINA" 			SM:"DEX_shMED1_1" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 70000000000 \
  --limitIObufferSize 1000000000
star_align.1.DEX_shMED1_1_2.24813d3bab959a0546a87211546cf7ce.mugqic.done
)
star_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=5:00:0 -q metaq -l nodes=1:ppn=16 -l pmem=5700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$star_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: star_5_JOB_ID: star_align.1.DEX_shNIPBL_1_1
#-------------------------------------------------------------------------------
JOB_NAME=star_align.1.DEX_shNIPBL_1_1
JOB_DEPENDENCIES=$trimmomatic_5_JOB_ID
JOB_DONE=job_output/star/star_align.1.DEX_shNIPBL_1_1.83fca6fe8730600eb5383971129ce90e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.1.DEX_shNIPBL_1_1.83fca6fe8730600eb5383971129ce90e.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment_1stPass/DEX_shNIPBL_1/DEX_shNIPBL_1_1 && \
STAR --runMode alignReads \
  --genomeDir /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/star_index/Ensembl87.sjdbOverhang99 \
  --readFilesIn \
    trim/DEX_shNIPBL_1/DEX_shNIPBL_1_1.trim.pair1.fastq.gz \
    trim/DEX_shNIPBL_1/DEX_shNIPBL_1_1.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix alignment_1stPass/DEX_shNIPBL_1/DEX_shNIPBL_1_1/ \
  --outSAMattrRGline ID:"DEX_shNIPBL_1_1" 	PL:"ILLUMINA" 			SM:"DEX_shNIPBL_1" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 70000000000 \
  --limitIObufferSize 1000000000
star_align.1.DEX_shNIPBL_1_1.83fca6fe8730600eb5383971129ce90e.mugqic.done
)
star_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=5:00:0 -q metaq -l nodes=1:ppn=16 -l pmem=5700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$star_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: star_6_JOB_ID: star_align.1.DEX_shNIPBL_1_2
#-------------------------------------------------------------------------------
JOB_NAME=star_align.1.DEX_shNIPBL_1_2
JOB_DEPENDENCIES=$trimmomatic_6_JOB_ID
JOB_DONE=job_output/star/star_align.1.DEX_shNIPBL_1_2.28749aa75885db709b8f25fef85c4120.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.1.DEX_shNIPBL_1_2.28749aa75885db709b8f25fef85c4120.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment_1stPass/DEX_shNIPBL_1/DEX_shNIPBL_1_2 && \
STAR --runMode alignReads \
  --genomeDir /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/star_index/Ensembl87.sjdbOverhang99 \
  --readFilesIn \
    trim/DEX_shNIPBL_1/DEX_shNIPBL_1_2.trim.pair1.fastq.gz \
    trim/DEX_shNIPBL_1/DEX_shNIPBL_1_2.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix alignment_1stPass/DEX_shNIPBL_1/DEX_shNIPBL_1_2/ \
  --outSAMattrRGline ID:"DEX_shNIPBL_1_2" 	PL:"ILLUMINA" 			SM:"DEX_shNIPBL_1" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 70000000000 \
  --limitIObufferSize 1000000000
star_align.1.DEX_shNIPBL_1_2.28749aa75885db709b8f25fef85c4120.mugqic.done
)
star_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=5:00:0 -q metaq -l nodes=1:ppn=16 -l pmem=5700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$star_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: star_7_JOB_ID: star_align.1.DEX_shSMC1A_1_1
#-------------------------------------------------------------------------------
JOB_NAME=star_align.1.DEX_shSMC1A_1_1
JOB_DEPENDENCIES=$trimmomatic_7_JOB_ID
JOB_DONE=job_output/star/star_align.1.DEX_shSMC1A_1_1.9fb18858c0f46adcb494073a6bbf3e96.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.1.DEX_shSMC1A_1_1.9fb18858c0f46adcb494073a6bbf3e96.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment_1stPass/DEX_shSMC1A_3/DEX_shSMC1A_1_1 && \
STAR --runMode alignReads \
  --genomeDir /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/star_index/Ensembl87.sjdbOverhang99 \
  --readFilesIn \
    trim/DEX_shSMC1A_3/DEX_shSMC1A_1_1.trim.pair1.fastq.gz \
    trim/DEX_shSMC1A_3/DEX_shSMC1A_1_1.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix alignment_1stPass/DEX_shSMC1A_3/DEX_shSMC1A_1_1/ \
  --outSAMattrRGline ID:"DEX_shSMC1A_1_1" 	PL:"ILLUMINA" 			SM:"DEX_shSMC1A_3" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 70000000000 \
  --limitIObufferSize 1000000000
star_align.1.DEX_shSMC1A_1_1.9fb18858c0f46adcb494073a6bbf3e96.mugqic.done
)
star_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=5:00:0 -q metaq -l nodes=1:ppn=16 -l pmem=5700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$star_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: star_8_JOB_ID: star_align.1.DEX_shSMC1A_1_2
#-------------------------------------------------------------------------------
JOB_NAME=star_align.1.DEX_shSMC1A_1_2
JOB_DEPENDENCIES=$trimmomatic_8_JOB_ID
JOB_DONE=job_output/star/star_align.1.DEX_shSMC1A_1_2.09982c5a24d34a2d77c2ee007b6ae9de.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.1.DEX_shSMC1A_1_2.09982c5a24d34a2d77c2ee007b6ae9de.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment_1stPass/DEX_shSMC1A_3/DEX_shSMC1A_1_2 && \
STAR --runMode alignReads \
  --genomeDir /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/star_index/Ensembl87.sjdbOverhang99 \
  --readFilesIn \
    trim/DEX_shSMC1A_3/DEX_shSMC1A_1_2.trim.pair1.fastq.gz \
    trim/DEX_shSMC1A_3/DEX_shSMC1A_1_2.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix alignment_1stPass/DEX_shSMC1A_3/DEX_shSMC1A_1_2/ \
  --outSAMattrRGline ID:"DEX_shSMC1A_1_2" 	PL:"ILLUMINA" 			SM:"DEX_shSMC1A_3" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 70000000000 \
  --limitIObufferSize 1000000000
star_align.1.DEX_shSMC1A_1_2.09982c5a24d34a2d77c2ee007b6ae9de.mugqic.done
)
star_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=5:00:0 -q metaq -l nodes=1:ppn=16 -l pmem=5700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$star_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: star_9_JOB_ID: star_align.1.ETOH_nonMamm_1_1
#-------------------------------------------------------------------------------
JOB_NAME=star_align.1.ETOH_nonMamm_1_1
JOB_DEPENDENCIES=$trimmomatic_9_JOB_ID
JOB_DONE=job_output/star/star_align.1.ETOH_nonMamm_1_1.dfb9917057b15d33cfa44464308cd588.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.1.ETOH_nonMamm_1_1.dfb9917057b15d33cfa44464308cd588.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment_1stPass/ETOH_nonMamm_1/ETOH_nonMamm_1_1 && \
STAR --runMode alignReads \
  --genomeDir /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/star_index/Ensembl87.sjdbOverhang99 \
  --readFilesIn \
    trim/ETOH_nonMamm_1/ETOH_nonMamm_1_1.trim.pair1.fastq.gz \
    trim/ETOH_nonMamm_1/ETOH_nonMamm_1_1.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix alignment_1stPass/ETOH_nonMamm_1/ETOH_nonMamm_1_1/ \
  --outSAMattrRGline ID:"ETOH_nonMamm_1_1" 	PL:"ILLUMINA" 			SM:"ETOH_nonMamm_1" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 70000000000 \
  --limitIObufferSize 1000000000
star_align.1.ETOH_nonMamm_1_1.dfb9917057b15d33cfa44464308cd588.mugqic.done
)
star_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=5:00:0 -q metaq -l nodes=1:ppn=16 -l pmem=5700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$star_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: star_10_JOB_ID: star_align.1.ETOH_nonMamm_1_2
#-------------------------------------------------------------------------------
JOB_NAME=star_align.1.ETOH_nonMamm_1_2
JOB_DEPENDENCIES=$trimmomatic_10_JOB_ID
JOB_DONE=job_output/star/star_align.1.ETOH_nonMamm_1_2.d0869b1268377483dfd4d096045e834b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.1.ETOH_nonMamm_1_2.d0869b1268377483dfd4d096045e834b.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment_1stPass/ETOH_nonMamm_1/ETOH_nonMamm_1_2 && \
STAR --runMode alignReads \
  --genomeDir /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/star_index/Ensembl87.sjdbOverhang99 \
  --readFilesIn \
    trim/ETOH_nonMamm_1/ETOH_nonMamm_1_2.trim.pair1.fastq.gz \
    trim/ETOH_nonMamm_1/ETOH_nonMamm_1_2.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix alignment_1stPass/ETOH_nonMamm_1/ETOH_nonMamm_1_2/ \
  --outSAMattrRGline ID:"ETOH_nonMamm_1_2" 	PL:"ILLUMINA" 			SM:"ETOH_nonMamm_1" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 70000000000 \
  --limitIObufferSize 1000000000
star_align.1.ETOH_nonMamm_1_2.d0869b1268377483dfd4d096045e834b.mugqic.done
)
star_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=5:00:0 -q metaq -l nodes=1:ppn=16 -l pmem=5700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$star_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: star_11_JOB_ID: star_align.1.ETOH_nonMamm_2_1
#-------------------------------------------------------------------------------
JOB_NAME=star_align.1.ETOH_nonMamm_2_1
JOB_DEPENDENCIES=$trimmomatic_11_JOB_ID
JOB_DONE=job_output/star/star_align.1.ETOH_nonMamm_2_1.988c4693b2dfdbbee7720e77af9fa83d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.1.ETOH_nonMamm_2_1.988c4693b2dfdbbee7720e77af9fa83d.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment_1stPass/ETOH_nonMamm_2/ETOH_nonMamm_2_1 && \
STAR --runMode alignReads \
  --genomeDir /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/star_index/Ensembl87.sjdbOverhang99 \
  --readFilesIn \
    trim/ETOH_nonMamm_2/ETOH_nonMamm_2_1.trim.pair1.fastq.gz \
    trim/ETOH_nonMamm_2/ETOH_nonMamm_2_1.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix alignment_1stPass/ETOH_nonMamm_2/ETOH_nonMamm_2_1/ \
  --outSAMattrRGline ID:"ETOH_nonMamm_2_1" 	PL:"ILLUMINA" 			SM:"ETOH_nonMamm_2" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 70000000000 \
  --limitIObufferSize 1000000000
star_align.1.ETOH_nonMamm_2_1.988c4693b2dfdbbee7720e77af9fa83d.mugqic.done
)
star_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=5:00:0 -q metaq -l nodes=1:ppn=16 -l pmem=5700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$star_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: star_12_JOB_ID: star_align.1.ETOH_shMED1_1_1
#-------------------------------------------------------------------------------
JOB_NAME=star_align.1.ETOH_shMED1_1_1
JOB_DEPENDENCIES=$trimmomatic_12_JOB_ID
JOB_DONE=job_output/star/star_align.1.ETOH_shMED1_1_1.30de45bc7f681cf8c0ea707a4881d244.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.1.ETOH_shMED1_1_1.30de45bc7f681cf8c0ea707a4881d244.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment_1stPass/ETOH_shMED1_1/ETOH_shMED1_1_1 && \
STAR --runMode alignReads \
  --genomeDir /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/star_index/Ensembl87.sjdbOverhang99 \
  --readFilesIn \
    trim/ETOH_shMED1_1/ETOH_shMED1_1_1.trim.pair1.fastq.gz \
    trim/ETOH_shMED1_1/ETOH_shMED1_1_1.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix alignment_1stPass/ETOH_shMED1_1/ETOH_shMED1_1_1/ \
  --outSAMattrRGline ID:"ETOH_shMED1_1_1" 	PL:"ILLUMINA" 			SM:"ETOH_shMED1_1" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 70000000000 \
  --limitIObufferSize 1000000000
star_align.1.ETOH_shMED1_1_1.30de45bc7f681cf8c0ea707a4881d244.mugqic.done
)
star_12_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=5:00:0 -q metaq -l nodes=1:ppn=16 -l pmem=5700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$star_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: star_13_JOB_ID: star_align.1.ETOH_shMED1_1_2
#-------------------------------------------------------------------------------
JOB_NAME=star_align.1.ETOH_shMED1_1_2
JOB_DEPENDENCIES=$trimmomatic_13_JOB_ID
JOB_DONE=job_output/star/star_align.1.ETOH_shMED1_1_2.5dcc3eb5eaf4ede3b3323a3e436b4f0d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.1.ETOH_shMED1_1_2.5dcc3eb5eaf4ede3b3323a3e436b4f0d.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment_1stPass/ETOH_shMED1_1/ETOH_shMED1_1_2 && \
STAR --runMode alignReads \
  --genomeDir /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/star_index/Ensembl87.sjdbOverhang99 \
  --readFilesIn \
    trim/ETOH_shMED1_1/ETOH_shMED1_1_2.trim.pair1.fastq.gz \
    trim/ETOH_shMED1_1/ETOH_shMED1_1_2.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix alignment_1stPass/ETOH_shMED1_1/ETOH_shMED1_1_2/ \
  --outSAMattrRGline ID:"ETOH_shMED1_1_2" 	PL:"ILLUMINA" 			SM:"ETOH_shMED1_1" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 70000000000 \
  --limitIObufferSize 1000000000
star_align.1.ETOH_shMED1_1_2.5dcc3eb5eaf4ede3b3323a3e436b4f0d.mugqic.done
)
star_13_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=5:00:0 -q metaq -l nodes=1:ppn=16 -l pmem=5700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$star_13_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: star_14_JOB_ID: star_align.1.ETOH_shNIPBL_2_1
#-------------------------------------------------------------------------------
JOB_NAME=star_align.1.ETOH_shNIPBL_2_1
JOB_DEPENDENCIES=$trimmomatic_14_JOB_ID
JOB_DONE=job_output/star/star_align.1.ETOH_shNIPBL_2_1.8420a62b40834f5e590a6fdbc6717315.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.1.ETOH_shNIPBL_2_1.8420a62b40834f5e590a6fdbc6717315.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment_1stPass/ETOH_shNIPBL_2/ETOH_shNIPBL_2_1 && \
STAR --runMode alignReads \
  --genomeDir /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/star_index/Ensembl87.sjdbOverhang99 \
  --readFilesIn \
    trim/ETOH_shNIPBL_2/ETOH_shNIPBL_2_1.trim.pair1.fastq.gz \
    trim/ETOH_shNIPBL_2/ETOH_shNIPBL_2_1.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix alignment_1stPass/ETOH_shNIPBL_2/ETOH_shNIPBL_2_1/ \
  --outSAMattrRGline ID:"ETOH_shNIPBL_2_1" 	PL:"ILLUMINA" 			SM:"ETOH_shNIPBL_2" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 70000000000 \
  --limitIObufferSize 1000000000
star_align.1.ETOH_shNIPBL_2_1.8420a62b40834f5e590a6fdbc6717315.mugqic.done
)
star_14_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=5:00:0 -q metaq -l nodes=1:ppn=16 -l pmem=5700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$star_14_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: star_15_JOB_ID: star_align.1.ETOH_shNIPBL_1_1
#-------------------------------------------------------------------------------
JOB_NAME=star_align.1.ETOH_shNIPBL_1_1
JOB_DEPENDENCIES=$trimmomatic_15_JOB_ID
JOB_DONE=job_output/star/star_align.1.ETOH_shNIPBL_1_1.ea1b33dad2ef3713aef07a1db1823d2a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.1.ETOH_shNIPBL_1_1.ea1b33dad2ef3713aef07a1db1823d2a.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment_1stPass/ETOH_shNIPBL_1/ETOH_shNIPBL_1_1 && \
STAR --runMode alignReads \
  --genomeDir /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/star_index/Ensembl87.sjdbOverhang99 \
  --readFilesIn \
    trim/ETOH_shNIPBL_1/ETOH_shNIPBL_1_1.trim.pair1.fastq.gz \
    trim/ETOH_shNIPBL_1/ETOH_shNIPBL_1_1.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix alignment_1stPass/ETOH_shNIPBL_1/ETOH_shNIPBL_1_1/ \
  --outSAMattrRGline ID:"ETOH_shNIPBL_1_1" 	PL:"ILLUMINA" 			SM:"ETOH_shNIPBL_1" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 70000000000 \
  --limitIObufferSize 1000000000
star_align.1.ETOH_shNIPBL_1_1.ea1b33dad2ef3713aef07a1db1823d2a.mugqic.done
)
star_15_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=5:00:0 -q metaq -l nodes=1:ppn=16 -l pmem=5700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$star_15_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: star_16_JOB_ID: star_align.1.ETOH_shNIPBL_1_2
#-------------------------------------------------------------------------------
JOB_NAME=star_align.1.ETOH_shNIPBL_1_2
JOB_DEPENDENCIES=$trimmomatic_16_JOB_ID
JOB_DONE=job_output/star/star_align.1.ETOH_shNIPBL_1_2.1d43e4543580bf1e7dcead38b7410842.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.1.ETOH_shNIPBL_1_2.1d43e4543580bf1e7dcead38b7410842.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment_1stPass/ETOH_shNIPBL_1/ETOH_shNIPBL_1_2 && \
STAR --runMode alignReads \
  --genomeDir /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/star_index/Ensembl87.sjdbOverhang99 \
  --readFilesIn \
    trim/ETOH_shNIPBL_1/ETOH_shNIPBL_1_2.trim.pair1.fastq.gz \
    trim/ETOH_shNIPBL_1/ETOH_shNIPBL_1_2.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix alignment_1stPass/ETOH_shNIPBL_1/ETOH_shNIPBL_1_2/ \
  --outSAMattrRGline ID:"ETOH_shNIPBL_1_2" 	PL:"ILLUMINA" 			SM:"ETOH_shNIPBL_1" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 70000000000 \
  --limitIObufferSize 1000000000
star_align.1.ETOH_shNIPBL_1_2.1d43e4543580bf1e7dcead38b7410842.mugqic.done
)
star_16_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=5:00:0 -q metaq -l nodes=1:ppn=16 -l pmem=5700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$star_16_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: star_17_JOB_ID: star_align.1.ETOH_shMED1_2_1
#-------------------------------------------------------------------------------
JOB_NAME=star_align.1.ETOH_shMED1_2_1
JOB_DEPENDENCIES=$trimmomatic_17_JOB_ID
JOB_DONE=job_output/star/star_align.1.ETOH_shMED1_2_1.aafdc333d270bdaa9993d4ddad739a9a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.1.ETOH_shMED1_2_1.aafdc333d270bdaa9993d4ddad739a9a.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment_1stPass/ETOH_shMED1_2/ETOH_shMED1_2_1 && \
STAR --runMode alignReads \
  --genomeDir /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/star_index/Ensembl87.sjdbOverhang99 \
  --readFilesIn \
    trim/ETOH_shMED1_2/ETOH_shMED1_2_1.trim.pair1.fastq.gz \
    trim/ETOH_shMED1_2/ETOH_shMED1_2_1.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix alignment_1stPass/ETOH_shMED1_2/ETOH_shMED1_2_1/ \
  --outSAMattrRGline ID:"ETOH_shMED1_2_1" 	PL:"ILLUMINA" 			SM:"ETOH_shMED1_2" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 70000000000 \
  --limitIObufferSize 1000000000
star_align.1.ETOH_shMED1_2_1.aafdc333d270bdaa9993d4ddad739a9a.mugqic.done
)
star_17_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=5:00:0 -q metaq -l nodes=1:ppn=16 -l pmem=5700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$star_17_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: star_18_JOB_ID: star_align.1.ETOH_shSMC1A_1_1
#-------------------------------------------------------------------------------
JOB_NAME=star_align.1.ETOH_shSMC1A_1_1
JOB_DEPENDENCIES=$trimmomatic_18_JOB_ID
JOB_DONE=job_output/star/star_align.1.ETOH_shSMC1A_1_1.f509d607efe01d9f811573a030cc1a7e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.1.ETOH_shSMC1A_1_1.f509d607efe01d9f811573a030cc1a7e.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment_1stPass/ETOH_shSMC1A_3/ETOH_shSMC1A_1_1 && \
STAR --runMode alignReads \
  --genomeDir /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/star_index/Ensembl87.sjdbOverhang99 \
  --readFilesIn \
    trim/ETOH_shSMC1A_3/ETOH_shSMC1A_1_1.trim.pair1.fastq.gz \
    trim/ETOH_shSMC1A_3/ETOH_shSMC1A_1_1.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix alignment_1stPass/ETOH_shSMC1A_3/ETOH_shSMC1A_1_1/ \
  --outSAMattrRGline ID:"ETOH_shSMC1A_1_1" 	PL:"ILLUMINA" 			SM:"ETOH_shSMC1A_3" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 70000000000 \
  --limitIObufferSize 1000000000
star_align.1.ETOH_shSMC1A_1_1.f509d607efe01d9f811573a030cc1a7e.mugqic.done
)
star_18_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=5:00:0 -q metaq -l nodes=1:ppn=16 -l pmem=5700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$star_18_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: star_19_JOB_ID: star_align.1.ETOH_shSMC1A_1_2
#-------------------------------------------------------------------------------
JOB_NAME=star_align.1.ETOH_shSMC1A_1_2
JOB_DEPENDENCIES=$trimmomatic_19_JOB_ID
JOB_DONE=job_output/star/star_align.1.ETOH_shSMC1A_1_2.32d0e97557262ef78c1f3ecf216ee899.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.1.ETOH_shSMC1A_1_2.32d0e97557262ef78c1f3ecf216ee899.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment_1stPass/ETOH_shSMC1A_3/ETOH_shSMC1A_1_2 && \
STAR --runMode alignReads \
  --genomeDir /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/star_index/Ensembl87.sjdbOverhang99 \
  --readFilesIn \
    trim/ETOH_shSMC1A_3/ETOH_shSMC1A_1_2.trim.pair1.fastq.gz \
    trim/ETOH_shSMC1A_3/ETOH_shSMC1A_1_2.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix alignment_1stPass/ETOH_shSMC1A_3/ETOH_shSMC1A_1_2/ \
  --outSAMattrRGline ID:"ETOH_shSMC1A_1_2" 	PL:"ILLUMINA" 			SM:"ETOH_shSMC1A_3" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 70000000000 \
  --limitIObufferSize 1000000000
star_align.1.ETOH_shSMC1A_1_2.32d0e97557262ef78c1f3ecf216ee899.mugqic.done
)
star_19_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=5:00:0 -q metaq -l nodes=1:ppn=16 -l pmem=5700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$star_19_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: star_20_JOB_ID: star_align.1.ETOH_shSMC1A_2_1
#-------------------------------------------------------------------------------
JOB_NAME=star_align.1.ETOH_shSMC1A_2_1
JOB_DEPENDENCIES=$trimmomatic_20_JOB_ID
JOB_DONE=job_output/star/star_align.1.ETOH_shSMC1A_2_1.786fb41602229cb940dfea714691dcb9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.1.ETOH_shSMC1A_2_1.786fb41602229cb940dfea714691dcb9.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment_1stPass/ETOH_shSMC1A_2/ETOH_shSMC1A_2_1 && \
STAR --runMode alignReads \
  --genomeDir /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/star_index/Ensembl87.sjdbOverhang99 \
  --readFilesIn \
    trim/ETOH_shSMC1A_2/ETOH_shSMC1A_2_1.trim.pair1.fastq.gz \
    trim/ETOH_shSMC1A_2/ETOH_shSMC1A_2_1.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix alignment_1stPass/ETOH_shSMC1A_2/ETOH_shSMC1A_2_1/ \
  --outSAMattrRGline ID:"ETOH_shSMC1A_2_1" 	PL:"ILLUMINA" 			SM:"ETOH_shSMC1A_2" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 70000000000 \
  --limitIObufferSize 1000000000
star_align.1.ETOH_shSMC1A_2_1.786fb41602229cb940dfea714691dcb9.mugqic.done
)
star_20_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=5:00:0 -q metaq -l nodes=1:ppn=16 -l pmem=5700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$star_20_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: star_21_JOB_ID: star_index.AllSamples
#-------------------------------------------------------------------------------
JOB_NAME=star_index.AllSamples
JOB_DEPENDENCIES=$star_1_JOB_ID:$star_2_JOB_ID:$star_3_JOB_ID:$star_4_JOB_ID:$star_5_JOB_ID:$star_6_JOB_ID:$star_7_JOB_ID:$star_8_JOB_ID:$star_9_JOB_ID:$star_10_JOB_ID:$star_11_JOB_ID:$star_12_JOB_ID:$star_13_JOB_ID:$star_14_JOB_ID:$star_15_JOB_ID:$star_16_JOB_ID:$star_17_JOB_ID:$star_18_JOB_ID:$star_19_JOB_ID:$star_20_JOB_ID
JOB_DONE=job_output/star/star_index.AllSamples.f5d2cb2187c87f55ce5d41b2ae760d12.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_index.AllSamples.f5d2cb2187c87f55ce5d41b2ae760d12.mugqic.done'
module load mugqic/star/2.5.1b && \
cat \
  alignment_1stPass/DEX_nonMamm_1/DEX_nonMamm_1_1/SJ.out.tab \
  alignment_1stPass/DEX_nonMamm_1/DEX_nonMamm_1_2/SJ.out.tab \
  alignment_1stPass/DEX_shMED1_1/DEX_shMED1_1_1/SJ.out.tab \
  alignment_1stPass/DEX_shMED1_1/DEX_shMED1_1_2/SJ.out.tab \
  alignment_1stPass/DEX_shNIPBL_1/DEX_shNIPBL_1_1/SJ.out.tab \
  alignment_1stPass/DEX_shNIPBL_1/DEX_shNIPBL_1_2/SJ.out.tab \
  alignment_1stPass/DEX_shSMC1A_3/DEX_shSMC1A_1_1/SJ.out.tab \
  alignment_1stPass/DEX_shSMC1A_3/DEX_shSMC1A_1_2/SJ.out.tab \
  alignment_1stPass/ETOH_nonMamm_1/ETOH_nonMamm_1_1/SJ.out.tab \
  alignment_1stPass/ETOH_nonMamm_1/ETOH_nonMamm_1_2/SJ.out.tab \
  alignment_1stPass/ETOH_nonMamm_2/ETOH_nonMamm_2_1/SJ.out.tab \
  alignment_1stPass/ETOH_shMED1_1/ETOH_shMED1_1_1/SJ.out.tab \
  alignment_1stPass/ETOH_shMED1_1/ETOH_shMED1_1_2/SJ.out.tab \
  alignment_1stPass/ETOH_shNIPBL_2/ETOH_shNIPBL_2_1/SJ.out.tab \
  alignment_1stPass/ETOH_shNIPBL_1/ETOH_shNIPBL_1_1/SJ.out.tab \
  alignment_1stPass/ETOH_shNIPBL_1/ETOH_shNIPBL_1_2/SJ.out.tab \
  alignment_1stPass/ETOH_shMED1_2/ETOH_shMED1_2_1/SJ.out.tab \
  alignment_1stPass/ETOH_shSMC1A_3/ETOH_shSMC1A_1_1/SJ.out.tab \
  alignment_1stPass/ETOH_shSMC1A_3/ETOH_shSMC1A_1_2/SJ.out.tab \
  alignment_1stPass/ETOH_shSMC1A_2/ETOH_shSMC1A_2_1/SJ.out.tab | \
awk 'BEGIN {OFS="	"; strChar[0]="."; strChar[1]="+"; strChar[2]="-"} {if($5>0){print $1,$2,$3,strChar[$4]}}' | sort -k1,1h -k2,2n > alignment_1stPass/AllSamples.SJ.out.tab && \
mkdir -p reference.Merged && \
STAR --runMode genomeGenerate \
  --genomeDir reference.Merged \
  --genomeFastaFiles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --runThreadN 12 \
  --limitGenomeGenerateRAM 50000000000 \
  --sjdbFileChrStartEnd alignment_1stPass/AllSamples.SJ.out.tab \
  --limitIObufferSize 1000000000 \
  --sjdbOverhang 99
star_index.AllSamples.f5d2cb2187c87f55ce5d41b2ae760d12.mugqic.done
)
star_21_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q lm -l nodes=1:ppn=12 -l pmem=5700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$star_21_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: star_22_JOB_ID: star_align.2.DEX_nonMamm_1_1
#-------------------------------------------------------------------------------
JOB_NAME=star_align.2.DEX_nonMamm_1_1
JOB_DEPENDENCIES=$trimmomatic_1_JOB_ID:$star_21_JOB_ID
JOB_DONE=job_output/star/star_align.2.DEX_nonMamm_1_1.a8428b57c67e96105c734879fbc94e00.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.2.DEX_nonMamm_1_1.a8428b57c67e96105c734879fbc94e00.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment/DEX_nonMamm_1/DEX_nonMamm_1_1 && \
STAR --runMode alignReads \
  --genomeDir reference.Merged \
  --readFilesIn \
    trim/DEX_nonMamm_1/DEX_nonMamm_1_1.trim.pair1.fastq.gz \
    trim/DEX_nonMamm_1/DEX_nonMamm_1_1.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix alignment/DEX_nonMamm_1/DEX_nonMamm_1_1/ \
  --outSAMattrRGline ID:"DEX_nonMamm_1_1" 	PL:"ILLUMINA" 			SM:"DEX_nonMamm_1" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 70000000000 \
  --limitBAMsortRAM 70000000000 \
  --limitIObufferSize 1000000000 \
  --outWigType wiggle read1_5p --outWigStrand Stranded --outWigReferencesPrefix chr \
  --chimSegmentMin 21
star_align.2.DEX_nonMamm_1_1.a8428b57c67e96105c734879fbc94e00.mugqic.done
)
star_22_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=5:00:0 -q metaq -l nodes=1:ppn=16 -l pmem=5700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$star_22_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: star_23_JOB_ID: star_align.2.DEX_nonMamm_1_2
#-------------------------------------------------------------------------------
JOB_NAME=star_align.2.DEX_nonMamm_1_2
JOB_DEPENDENCIES=$trimmomatic_2_JOB_ID:$star_21_JOB_ID
JOB_DONE=job_output/star/star_align.2.DEX_nonMamm_1_2.078dc1ebdb1f53d51bdd5dbfdf2ef69b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.2.DEX_nonMamm_1_2.078dc1ebdb1f53d51bdd5dbfdf2ef69b.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment/DEX_nonMamm_1/DEX_nonMamm_1_2 && \
STAR --runMode alignReads \
  --genomeDir reference.Merged \
  --readFilesIn \
    trim/DEX_nonMamm_1/DEX_nonMamm_1_2.trim.pair1.fastq.gz \
    trim/DEX_nonMamm_1/DEX_nonMamm_1_2.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix alignment/DEX_nonMamm_1/DEX_nonMamm_1_2/ \
  --outSAMattrRGline ID:"DEX_nonMamm_1_2" 	PL:"ILLUMINA" 			SM:"DEX_nonMamm_1" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 70000000000 \
  --limitBAMsortRAM 70000000000 \
  --limitIObufferSize 1000000000 \
  --outWigType wiggle read1_5p --outWigStrand Stranded --outWigReferencesPrefix chr \
  --chimSegmentMin 21
star_align.2.DEX_nonMamm_1_2.078dc1ebdb1f53d51bdd5dbfdf2ef69b.mugqic.done
)
star_23_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=5:00:0 -q metaq -l nodes=1:ppn=16 -l pmem=5700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$star_23_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: star_24_JOB_ID: star_align.2.DEX_shMED1_1_1
#-------------------------------------------------------------------------------
JOB_NAME=star_align.2.DEX_shMED1_1_1
JOB_DEPENDENCIES=$trimmomatic_3_JOB_ID:$star_21_JOB_ID
JOB_DONE=job_output/star/star_align.2.DEX_shMED1_1_1.99dac859adf98ddc9f9bafe036a598f1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.2.DEX_shMED1_1_1.99dac859adf98ddc9f9bafe036a598f1.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment/DEX_shMED1_1/DEX_shMED1_1_1 && \
STAR --runMode alignReads \
  --genomeDir reference.Merged \
  --readFilesIn \
    trim/DEX_shMED1_1/DEX_shMED1_1_1.trim.pair1.fastq.gz \
    trim/DEX_shMED1_1/DEX_shMED1_1_1.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix alignment/DEX_shMED1_1/DEX_shMED1_1_1/ \
  --outSAMattrRGline ID:"DEX_shMED1_1_1" 	PL:"ILLUMINA" 			SM:"DEX_shMED1_1" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 70000000000 \
  --limitBAMsortRAM 70000000000 \
  --limitIObufferSize 1000000000 \
  --outWigType wiggle read1_5p --outWigStrand Stranded --outWigReferencesPrefix chr \
  --chimSegmentMin 21
star_align.2.DEX_shMED1_1_1.99dac859adf98ddc9f9bafe036a598f1.mugqic.done
)
star_24_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=5:00:0 -q metaq -l nodes=1:ppn=16 -l pmem=5700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$star_24_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: star_25_JOB_ID: star_align.2.DEX_shMED1_1_2
#-------------------------------------------------------------------------------
JOB_NAME=star_align.2.DEX_shMED1_1_2
JOB_DEPENDENCIES=$trimmomatic_4_JOB_ID:$star_21_JOB_ID
JOB_DONE=job_output/star/star_align.2.DEX_shMED1_1_2.0d60a57062f39bf64976051234b46e5c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.2.DEX_shMED1_1_2.0d60a57062f39bf64976051234b46e5c.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment/DEX_shMED1_1/DEX_shMED1_1_2 && \
STAR --runMode alignReads \
  --genomeDir reference.Merged \
  --readFilesIn \
    trim/DEX_shMED1_1/DEX_shMED1_1_2.trim.pair1.fastq.gz \
    trim/DEX_shMED1_1/DEX_shMED1_1_2.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix alignment/DEX_shMED1_1/DEX_shMED1_1_2/ \
  --outSAMattrRGline ID:"DEX_shMED1_1_2" 	PL:"ILLUMINA" 			SM:"DEX_shMED1_1" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 70000000000 \
  --limitBAMsortRAM 70000000000 \
  --limitIObufferSize 1000000000 \
  --outWigType wiggle read1_5p --outWigStrand Stranded --outWigReferencesPrefix chr \
  --chimSegmentMin 21
star_align.2.DEX_shMED1_1_2.0d60a57062f39bf64976051234b46e5c.mugqic.done
)
star_25_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=5:00:0 -q metaq -l nodes=1:ppn=16 -l pmem=5700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$star_25_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: star_26_JOB_ID: star_align.2.DEX_shNIPBL_1_1
#-------------------------------------------------------------------------------
JOB_NAME=star_align.2.DEX_shNIPBL_1_1
JOB_DEPENDENCIES=$trimmomatic_5_JOB_ID:$star_21_JOB_ID
JOB_DONE=job_output/star/star_align.2.DEX_shNIPBL_1_1.2c70a5660f8a036cedd5aad391cb03de.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.2.DEX_shNIPBL_1_1.2c70a5660f8a036cedd5aad391cb03de.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment/DEX_shNIPBL_1/DEX_shNIPBL_1_1 && \
STAR --runMode alignReads \
  --genomeDir reference.Merged \
  --readFilesIn \
    trim/DEX_shNIPBL_1/DEX_shNIPBL_1_1.trim.pair1.fastq.gz \
    trim/DEX_shNIPBL_1/DEX_shNIPBL_1_1.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix alignment/DEX_shNIPBL_1/DEX_shNIPBL_1_1/ \
  --outSAMattrRGline ID:"DEX_shNIPBL_1_1" 	PL:"ILLUMINA" 			SM:"DEX_shNIPBL_1" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 70000000000 \
  --limitBAMsortRAM 70000000000 \
  --limitIObufferSize 1000000000 \
  --outWigType wiggle read1_5p --outWigStrand Stranded --outWigReferencesPrefix chr \
  --chimSegmentMin 21
star_align.2.DEX_shNIPBL_1_1.2c70a5660f8a036cedd5aad391cb03de.mugqic.done
)
star_26_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=5:00:0 -q metaq -l nodes=1:ppn=16 -l pmem=5700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$star_26_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: star_27_JOB_ID: star_align.2.DEX_shNIPBL_1_2
#-------------------------------------------------------------------------------
JOB_NAME=star_align.2.DEX_shNIPBL_1_2
JOB_DEPENDENCIES=$trimmomatic_6_JOB_ID:$star_21_JOB_ID
JOB_DONE=job_output/star/star_align.2.DEX_shNIPBL_1_2.b6bc907f8a5cc78a59e39d48cc2a373c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.2.DEX_shNIPBL_1_2.b6bc907f8a5cc78a59e39d48cc2a373c.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment/DEX_shNIPBL_1/DEX_shNIPBL_1_2 && \
STAR --runMode alignReads \
  --genomeDir reference.Merged \
  --readFilesIn \
    trim/DEX_shNIPBL_1/DEX_shNIPBL_1_2.trim.pair1.fastq.gz \
    trim/DEX_shNIPBL_1/DEX_shNIPBL_1_2.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix alignment/DEX_shNIPBL_1/DEX_shNIPBL_1_2/ \
  --outSAMattrRGline ID:"DEX_shNIPBL_1_2" 	PL:"ILLUMINA" 			SM:"DEX_shNIPBL_1" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 70000000000 \
  --limitBAMsortRAM 70000000000 \
  --limitIObufferSize 1000000000 \
  --outWigType wiggle read1_5p --outWigStrand Stranded --outWigReferencesPrefix chr \
  --chimSegmentMin 21
star_align.2.DEX_shNIPBL_1_2.b6bc907f8a5cc78a59e39d48cc2a373c.mugqic.done
)
star_27_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=5:00:0 -q metaq -l nodes=1:ppn=16 -l pmem=5700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$star_27_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: star_28_JOB_ID: star_align.2.DEX_shSMC1A_1_1
#-------------------------------------------------------------------------------
JOB_NAME=star_align.2.DEX_shSMC1A_1_1
JOB_DEPENDENCIES=$trimmomatic_7_JOB_ID:$star_21_JOB_ID
JOB_DONE=job_output/star/star_align.2.DEX_shSMC1A_1_1.cc8923e6a581f03d61f0527d1cf2fb8e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.2.DEX_shSMC1A_1_1.cc8923e6a581f03d61f0527d1cf2fb8e.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment/DEX_shSMC1A_3/DEX_shSMC1A_1_1 && \
STAR --runMode alignReads \
  --genomeDir reference.Merged \
  --readFilesIn \
    trim/DEX_shSMC1A_3/DEX_shSMC1A_1_1.trim.pair1.fastq.gz \
    trim/DEX_shSMC1A_3/DEX_shSMC1A_1_1.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix alignment/DEX_shSMC1A_3/DEX_shSMC1A_1_1/ \
  --outSAMattrRGline ID:"DEX_shSMC1A_1_1" 	PL:"ILLUMINA" 			SM:"DEX_shSMC1A_3" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 70000000000 \
  --limitBAMsortRAM 70000000000 \
  --limitIObufferSize 1000000000 \
  --outWigType wiggle read1_5p --outWigStrand Stranded --outWigReferencesPrefix chr \
  --chimSegmentMin 21
star_align.2.DEX_shSMC1A_1_1.cc8923e6a581f03d61f0527d1cf2fb8e.mugqic.done
)
star_28_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=5:00:0 -q metaq -l nodes=1:ppn=16 -l pmem=5700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$star_28_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: star_29_JOB_ID: star_align.2.DEX_shSMC1A_1_2
#-------------------------------------------------------------------------------
JOB_NAME=star_align.2.DEX_shSMC1A_1_2
JOB_DEPENDENCIES=$trimmomatic_8_JOB_ID:$star_21_JOB_ID
JOB_DONE=job_output/star/star_align.2.DEX_shSMC1A_1_2.ab541441c9bcd707a9677a28c54d5bdc.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.2.DEX_shSMC1A_1_2.ab541441c9bcd707a9677a28c54d5bdc.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment/DEX_shSMC1A_3/DEX_shSMC1A_1_2 && \
STAR --runMode alignReads \
  --genomeDir reference.Merged \
  --readFilesIn \
    trim/DEX_shSMC1A_3/DEX_shSMC1A_1_2.trim.pair1.fastq.gz \
    trim/DEX_shSMC1A_3/DEX_shSMC1A_1_2.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix alignment/DEX_shSMC1A_3/DEX_shSMC1A_1_2/ \
  --outSAMattrRGline ID:"DEX_shSMC1A_1_2" 	PL:"ILLUMINA" 			SM:"DEX_shSMC1A_3" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 70000000000 \
  --limitBAMsortRAM 70000000000 \
  --limitIObufferSize 1000000000 \
  --outWigType wiggle read1_5p --outWigStrand Stranded --outWigReferencesPrefix chr \
  --chimSegmentMin 21
star_align.2.DEX_shSMC1A_1_2.ab541441c9bcd707a9677a28c54d5bdc.mugqic.done
)
star_29_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=5:00:0 -q metaq -l nodes=1:ppn=16 -l pmem=5700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$star_29_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: star_30_JOB_ID: star_align.2.ETOH_nonMamm_1_1
#-------------------------------------------------------------------------------
JOB_NAME=star_align.2.ETOH_nonMamm_1_1
JOB_DEPENDENCIES=$trimmomatic_9_JOB_ID:$star_21_JOB_ID
JOB_DONE=job_output/star/star_align.2.ETOH_nonMamm_1_1.76fd2c0d2aaef0e91e28924854077c43.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.2.ETOH_nonMamm_1_1.76fd2c0d2aaef0e91e28924854077c43.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment/ETOH_nonMamm_1/ETOH_nonMamm_1_1 && \
STAR --runMode alignReads \
  --genomeDir reference.Merged \
  --readFilesIn \
    trim/ETOH_nonMamm_1/ETOH_nonMamm_1_1.trim.pair1.fastq.gz \
    trim/ETOH_nonMamm_1/ETOH_nonMamm_1_1.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix alignment/ETOH_nonMamm_1/ETOH_nonMamm_1_1/ \
  --outSAMattrRGline ID:"ETOH_nonMamm_1_1" 	PL:"ILLUMINA" 			SM:"ETOH_nonMamm_1" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 70000000000 \
  --limitBAMsortRAM 70000000000 \
  --limitIObufferSize 1000000000 \
  --outWigType wiggle read1_5p --outWigStrand Stranded --outWigReferencesPrefix chr \
  --chimSegmentMin 21
star_align.2.ETOH_nonMamm_1_1.76fd2c0d2aaef0e91e28924854077c43.mugqic.done
)
star_30_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=5:00:0 -q metaq -l nodes=1:ppn=16 -l pmem=5700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$star_30_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: star_31_JOB_ID: star_align.2.ETOH_nonMamm_1_2
#-------------------------------------------------------------------------------
JOB_NAME=star_align.2.ETOH_nonMamm_1_2
JOB_DEPENDENCIES=$trimmomatic_10_JOB_ID:$star_21_JOB_ID
JOB_DONE=job_output/star/star_align.2.ETOH_nonMamm_1_2.2b87b5bcc2ec2919f754fa5cfb718c80.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.2.ETOH_nonMamm_1_2.2b87b5bcc2ec2919f754fa5cfb718c80.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment/ETOH_nonMamm_1/ETOH_nonMamm_1_2 && \
STAR --runMode alignReads \
  --genomeDir reference.Merged \
  --readFilesIn \
    trim/ETOH_nonMamm_1/ETOH_nonMamm_1_2.trim.pair1.fastq.gz \
    trim/ETOH_nonMamm_1/ETOH_nonMamm_1_2.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix alignment/ETOH_nonMamm_1/ETOH_nonMamm_1_2/ \
  --outSAMattrRGline ID:"ETOH_nonMamm_1_2" 	PL:"ILLUMINA" 			SM:"ETOH_nonMamm_1" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 70000000000 \
  --limitBAMsortRAM 70000000000 \
  --limitIObufferSize 1000000000 \
  --outWigType wiggle read1_5p --outWigStrand Stranded --outWigReferencesPrefix chr \
  --chimSegmentMin 21
star_align.2.ETOH_nonMamm_1_2.2b87b5bcc2ec2919f754fa5cfb718c80.mugqic.done
)
star_31_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=5:00:0 -q metaq -l nodes=1:ppn=16 -l pmem=5700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$star_31_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: star_32_JOB_ID: star_align.2.ETOH_nonMamm_2_1
#-------------------------------------------------------------------------------
JOB_NAME=star_align.2.ETOH_nonMamm_2_1
JOB_DEPENDENCIES=$trimmomatic_11_JOB_ID:$star_21_JOB_ID
JOB_DONE=job_output/star/star_align.2.ETOH_nonMamm_2_1.45db5a08d8c039e2cd2680d823bd4665.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.2.ETOH_nonMamm_2_1.45db5a08d8c039e2cd2680d823bd4665.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment/ETOH_nonMamm_2/ETOH_nonMamm_2_1 && \
STAR --runMode alignReads \
  --genomeDir reference.Merged \
  --readFilesIn \
    trim/ETOH_nonMamm_2/ETOH_nonMamm_2_1.trim.pair1.fastq.gz \
    trim/ETOH_nonMamm_2/ETOH_nonMamm_2_1.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix alignment/ETOH_nonMamm_2/ETOH_nonMamm_2_1/ \
  --outSAMattrRGline ID:"ETOH_nonMamm_2_1" 	PL:"ILLUMINA" 			SM:"ETOH_nonMamm_2" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 70000000000 \
  --limitBAMsortRAM 70000000000 \
  --limitIObufferSize 1000000000 \
  --outWigType wiggle read1_5p --outWigStrand Stranded --outWigReferencesPrefix chr \
  --chimSegmentMin 21 && \
ln -s -f ETOH_nonMamm_2_1/Aligned.sortedByCoord.out.bam alignment/ETOH_nonMamm_2/ETOH_nonMamm_2.sorted.bam
star_align.2.ETOH_nonMamm_2_1.45db5a08d8c039e2cd2680d823bd4665.mugqic.done
)
star_32_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=5:00:0 -q metaq -l nodes=1:ppn=16 -l pmem=5700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$star_32_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: star_33_JOB_ID: star_align.2.ETOH_shMED1_1_1
#-------------------------------------------------------------------------------
JOB_NAME=star_align.2.ETOH_shMED1_1_1
JOB_DEPENDENCIES=$trimmomatic_12_JOB_ID:$star_21_JOB_ID
JOB_DONE=job_output/star/star_align.2.ETOH_shMED1_1_1.d0fac7df6eff9680d1aaae833ca4ee0d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.2.ETOH_shMED1_1_1.d0fac7df6eff9680d1aaae833ca4ee0d.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment/ETOH_shMED1_1/ETOH_shMED1_1_1 && \
STAR --runMode alignReads \
  --genomeDir reference.Merged \
  --readFilesIn \
    trim/ETOH_shMED1_1/ETOH_shMED1_1_1.trim.pair1.fastq.gz \
    trim/ETOH_shMED1_1/ETOH_shMED1_1_1.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix alignment/ETOH_shMED1_1/ETOH_shMED1_1_1/ \
  --outSAMattrRGline ID:"ETOH_shMED1_1_1" 	PL:"ILLUMINA" 			SM:"ETOH_shMED1_1" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 70000000000 \
  --limitBAMsortRAM 70000000000 \
  --limitIObufferSize 1000000000 \
  --outWigType wiggle read1_5p --outWigStrand Stranded --outWigReferencesPrefix chr \
  --chimSegmentMin 21
star_align.2.ETOH_shMED1_1_1.d0fac7df6eff9680d1aaae833ca4ee0d.mugqic.done
)
star_33_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=5:00:0 -q metaq -l nodes=1:ppn=16 -l pmem=5700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$star_33_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: star_34_JOB_ID: star_align.2.ETOH_shMED1_1_2
#-------------------------------------------------------------------------------
JOB_NAME=star_align.2.ETOH_shMED1_1_2
JOB_DEPENDENCIES=$trimmomatic_13_JOB_ID:$star_21_JOB_ID
JOB_DONE=job_output/star/star_align.2.ETOH_shMED1_1_2.c68660f777194cdcf0ab7e88246d5f13.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.2.ETOH_shMED1_1_2.c68660f777194cdcf0ab7e88246d5f13.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment/ETOH_shMED1_1/ETOH_shMED1_1_2 && \
STAR --runMode alignReads \
  --genomeDir reference.Merged \
  --readFilesIn \
    trim/ETOH_shMED1_1/ETOH_shMED1_1_2.trim.pair1.fastq.gz \
    trim/ETOH_shMED1_1/ETOH_shMED1_1_2.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix alignment/ETOH_shMED1_1/ETOH_shMED1_1_2/ \
  --outSAMattrRGline ID:"ETOH_shMED1_1_2" 	PL:"ILLUMINA" 			SM:"ETOH_shMED1_1" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 70000000000 \
  --limitBAMsortRAM 70000000000 \
  --limitIObufferSize 1000000000 \
  --outWigType wiggle read1_5p --outWigStrand Stranded --outWigReferencesPrefix chr \
  --chimSegmentMin 21
star_align.2.ETOH_shMED1_1_2.c68660f777194cdcf0ab7e88246d5f13.mugqic.done
)
star_34_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=5:00:0 -q metaq -l nodes=1:ppn=16 -l pmem=5700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$star_34_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: star_35_JOB_ID: star_align.2.ETOH_shNIPBL_2_1
#-------------------------------------------------------------------------------
JOB_NAME=star_align.2.ETOH_shNIPBL_2_1
JOB_DEPENDENCIES=$trimmomatic_14_JOB_ID:$star_21_JOB_ID
JOB_DONE=job_output/star/star_align.2.ETOH_shNIPBL_2_1.789840126fea58a8952462edcd554989.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.2.ETOH_shNIPBL_2_1.789840126fea58a8952462edcd554989.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment/ETOH_shNIPBL_2/ETOH_shNIPBL_2_1 && \
STAR --runMode alignReads \
  --genomeDir reference.Merged \
  --readFilesIn \
    trim/ETOH_shNIPBL_2/ETOH_shNIPBL_2_1.trim.pair1.fastq.gz \
    trim/ETOH_shNIPBL_2/ETOH_shNIPBL_2_1.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix alignment/ETOH_shNIPBL_2/ETOH_shNIPBL_2_1/ \
  --outSAMattrRGline ID:"ETOH_shNIPBL_2_1" 	PL:"ILLUMINA" 			SM:"ETOH_shNIPBL_2" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 70000000000 \
  --limitBAMsortRAM 70000000000 \
  --limitIObufferSize 1000000000 \
  --outWigType wiggle read1_5p --outWigStrand Stranded --outWigReferencesPrefix chr \
  --chimSegmentMin 21 && \
ln -s -f ETOH_shNIPBL_2_1/Aligned.sortedByCoord.out.bam alignment/ETOH_shNIPBL_2/ETOH_shNIPBL_2.sorted.bam
star_align.2.ETOH_shNIPBL_2_1.789840126fea58a8952462edcd554989.mugqic.done
)
star_35_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=5:00:0 -q metaq -l nodes=1:ppn=16 -l pmem=5700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$star_35_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: star_36_JOB_ID: star_align.2.ETOH_shNIPBL_1_1
#-------------------------------------------------------------------------------
JOB_NAME=star_align.2.ETOH_shNIPBL_1_1
JOB_DEPENDENCIES=$trimmomatic_15_JOB_ID:$star_21_JOB_ID
JOB_DONE=job_output/star/star_align.2.ETOH_shNIPBL_1_1.78732629bb9506acaabd13bf23330d79.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.2.ETOH_shNIPBL_1_1.78732629bb9506acaabd13bf23330d79.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1_1 && \
STAR --runMode alignReads \
  --genomeDir reference.Merged \
  --readFilesIn \
    trim/ETOH_shNIPBL_1/ETOH_shNIPBL_1_1.trim.pair1.fastq.gz \
    trim/ETOH_shNIPBL_1/ETOH_shNIPBL_1_1.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1_1/ \
  --outSAMattrRGline ID:"ETOH_shNIPBL_1_1" 	PL:"ILLUMINA" 			SM:"ETOH_shNIPBL_1" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 70000000000 \
  --limitBAMsortRAM 70000000000 \
  --limitIObufferSize 1000000000 \
  --outWigType wiggle read1_5p --outWigStrand Stranded --outWigReferencesPrefix chr \
  --chimSegmentMin 21
star_align.2.ETOH_shNIPBL_1_1.78732629bb9506acaabd13bf23330d79.mugqic.done
)
star_36_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=5:00:0 -q metaq -l nodes=1:ppn=16 -l pmem=5700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$star_36_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: star_37_JOB_ID: star_align.2.ETOH_shNIPBL_1_2
#-------------------------------------------------------------------------------
JOB_NAME=star_align.2.ETOH_shNIPBL_1_2
JOB_DEPENDENCIES=$trimmomatic_16_JOB_ID:$star_21_JOB_ID
JOB_DONE=job_output/star/star_align.2.ETOH_shNIPBL_1_2.fada2e707804a51df7ffc2cd0e2ac024.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.2.ETOH_shNIPBL_1_2.fada2e707804a51df7ffc2cd0e2ac024.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1_2 && \
STAR --runMode alignReads \
  --genomeDir reference.Merged \
  --readFilesIn \
    trim/ETOH_shNIPBL_1/ETOH_shNIPBL_1_2.trim.pair1.fastq.gz \
    trim/ETOH_shNIPBL_1/ETOH_shNIPBL_1_2.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1_2/ \
  --outSAMattrRGline ID:"ETOH_shNIPBL_1_2" 	PL:"ILLUMINA" 			SM:"ETOH_shNIPBL_1" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 70000000000 \
  --limitBAMsortRAM 70000000000 \
  --limitIObufferSize 1000000000 \
  --outWigType wiggle read1_5p --outWigStrand Stranded --outWigReferencesPrefix chr \
  --chimSegmentMin 21
star_align.2.ETOH_shNIPBL_1_2.fada2e707804a51df7ffc2cd0e2ac024.mugqic.done
)
star_37_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=5:00:0 -q metaq -l nodes=1:ppn=16 -l pmem=5700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$star_37_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: star_38_JOB_ID: star_align.2.ETOH_shMED1_2_1
#-------------------------------------------------------------------------------
JOB_NAME=star_align.2.ETOH_shMED1_2_1
JOB_DEPENDENCIES=$trimmomatic_17_JOB_ID:$star_21_JOB_ID
JOB_DONE=job_output/star/star_align.2.ETOH_shMED1_2_1.555b3cc2f96f1920e9e01013c5edde9f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.2.ETOH_shMED1_2_1.555b3cc2f96f1920e9e01013c5edde9f.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment/ETOH_shMED1_2/ETOH_shMED1_2_1 && \
STAR --runMode alignReads \
  --genomeDir reference.Merged \
  --readFilesIn \
    trim/ETOH_shMED1_2/ETOH_shMED1_2_1.trim.pair1.fastq.gz \
    trim/ETOH_shMED1_2/ETOH_shMED1_2_1.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix alignment/ETOH_shMED1_2/ETOH_shMED1_2_1/ \
  --outSAMattrRGline ID:"ETOH_shMED1_2_1" 	PL:"ILLUMINA" 			SM:"ETOH_shMED1_2" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 70000000000 \
  --limitBAMsortRAM 70000000000 \
  --limitIObufferSize 1000000000 \
  --outWigType wiggle read1_5p --outWigStrand Stranded --outWigReferencesPrefix chr \
  --chimSegmentMin 21 && \
ln -s -f ETOH_shMED1_2_1/Aligned.sortedByCoord.out.bam alignment/ETOH_shMED1_2/ETOH_shMED1_2.sorted.bam
star_align.2.ETOH_shMED1_2_1.555b3cc2f96f1920e9e01013c5edde9f.mugqic.done
)
star_38_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=5:00:0 -q metaq -l nodes=1:ppn=16 -l pmem=5700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$star_38_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: star_39_JOB_ID: star_align.2.ETOH_shSMC1A_1_1
#-------------------------------------------------------------------------------
JOB_NAME=star_align.2.ETOH_shSMC1A_1_1
JOB_DEPENDENCIES=$trimmomatic_18_JOB_ID:$star_21_JOB_ID
JOB_DONE=job_output/star/star_align.2.ETOH_shSMC1A_1_1.54c28883eb80a6d168a19b164c332220.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.2.ETOH_shSMC1A_1_1.54c28883eb80a6d168a19b164c332220.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment/ETOH_shSMC1A_3/ETOH_shSMC1A_1_1 && \
STAR --runMode alignReads \
  --genomeDir reference.Merged \
  --readFilesIn \
    trim/ETOH_shSMC1A_3/ETOH_shSMC1A_1_1.trim.pair1.fastq.gz \
    trim/ETOH_shSMC1A_3/ETOH_shSMC1A_1_1.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix alignment/ETOH_shSMC1A_3/ETOH_shSMC1A_1_1/ \
  --outSAMattrRGline ID:"ETOH_shSMC1A_1_1" 	PL:"ILLUMINA" 			SM:"ETOH_shSMC1A_3" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 70000000000 \
  --limitBAMsortRAM 70000000000 \
  --limitIObufferSize 1000000000 \
  --outWigType wiggle read1_5p --outWigStrand Stranded --outWigReferencesPrefix chr \
  --chimSegmentMin 21
star_align.2.ETOH_shSMC1A_1_1.54c28883eb80a6d168a19b164c332220.mugqic.done
)
star_39_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=5:00:0 -q metaq -l nodes=1:ppn=16 -l pmem=5700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$star_39_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: star_40_JOB_ID: star_align.2.ETOH_shSMC1A_1_2
#-------------------------------------------------------------------------------
JOB_NAME=star_align.2.ETOH_shSMC1A_1_2
JOB_DEPENDENCIES=$trimmomatic_19_JOB_ID:$star_21_JOB_ID
JOB_DONE=job_output/star/star_align.2.ETOH_shSMC1A_1_2.2bbf81c39e6da8dd0a1404551d804655.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.2.ETOH_shSMC1A_1_2.2bbf81c39e6da8dd0a1404551d804655.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment/ETOH_shSMC1A_3/ETOH_shSMC1A_1_2 && \
STAR --runMode alignReads \
  --genomeDir reference.Merged \
  --readFilesIn \
    trim/ETOH_shSMC1A_3/ETOH_shSMC1A_1_2.trim.pair1.fastq.gz \
    trim/ETOH_shSMC1A_3/ETOH_shSMC1A_1_2.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix alignment/ETOH_shSMC1A_3/ETOH_shSMC1A_1_2/ \
  --outSAMattrRGline ID:"ETOH_shSMC1A_1_2" 	PL:"ILLUMINA" 			SM:"ETOH_shSMC1A_3" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 70000000000 \
  --limitBAMsortRAM 70000000000 \
  --limitIObufferSize 1000000000 \
  --outWigType wiggle read1_5p --outWigStrand Stranded --outWigReferencesPrefix chr \
  --chimSegmentMin 21
star_align.2.ETOH_shSMC1A_1_2.2bbf81c39e6da8dd0a1404551d804655.mugqic.done
)
star_40_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=5:00:0 -q metaq -l nodes=1:ppn=16 -l pmem=5700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$star_40_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: star_41_JOB_ID: star_align.2.ETOH_shSMC1A_2_1
#-------------------------------------------------------------------------------
JOB_NAME=star_align.2.ETOH_shSMC1A_2_1
JOB_DEPENDENCIES=$trimmomatic_20_JOB_ID:$star_21_JOB_ID
JOB_DONE=job_output/star/star_align.2.ETOH_shSMC1A_2_1.d8236853f8a30e60488c02a8a704206f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.2.ETOH_shSMC1A_2_1.d8236853f8a30e60488c02a8a704206f.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment/ETOH_shSMC1A_2/ETOH_shSMC1A_2_1 && \
STAR --runMode alignReads \
  --genomeDir reference.Merged \
  --readFilesIn \
    trim/ETOH_shSMC1A_2/ETOH_shSMC1A_2_1.trim.pair1.fastq.gz \
    trim/ETOH_shSMC1A_2/ETOH_shSMC1A_2_1.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix alignment/ETOH_shSMC1A_2/ETOH_shSMC1A_2_1/ \
  --outSAMattrRGline ID:"ETOH_shSMC1A_2_1" 	PL:"ILLUMINA" 			SM:"ETOH_shSMC1A_2" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 70000000000 \
  --limitBAMsortRAM 70000000000 \
  --limitIObufferSize 1000000000 \
  --outWigType wiggle read1_5p --outWigStrand Stranded --outWigReferencesPrefix chr \
  --chimSegmentMin 21 && \
ln -s -f ETOH_shSMC1A_2_1/Aligned.sortedByCoord.out.bam alignment/ETOH_shSMC1A_2/ETOH_shSMC1A_2.sorted.bam
star_align.2.ETOH_shSMC1A_2_1.d8236853f8a30e60488c02a8a704206f.mugqic.done
)
star_41_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=5:00:0 -q metaq -l nodes=1:ppn=16 -l pmem=5700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$star_41_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: star_42_JOB_ID: star_report
#-------------------------------------------------------------------------------
JOB_NAME=star_report
JOB_DEPENDENCIES=$star_22_JOB_ID:$star_23_JOB_ID:$star_24_JOB_ID:$star_25_JOB_ID:$star_26_JOB_ID:$star_27_JOB_ID:$star_28_JOB_ID:$star_29_JOB_ID:$star_30_JOB_ID:$star_31_JOB_ID:$star_32_JOB_ID:$star_33_JOB_ID:$star_34_JOB_ID:$star_35_JOB_ID:$star_36_JOB_ID:$star_37_JOB_ID:$star_38_JOB_ID:$star_39_JOB_ID:$star_40_JOB_ID:$star_41_JOB_ID
JOB_DONE=job_output/star/star_report.9089faa6a2c55b3e6016d51622a0e925.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_report.9089faa6a2c55b3e6016d51622a0e925.mugqic.done'
module load mugqic/pandoc/1.15.2 && \
mkdir -p report && \
pandoc --to=markdown \
  --template /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.3.0/bfx/report/RnaSeq.star.md \
  --variable scientific_name="Homo_sapiens" \
  --variable assembly="GRCh38" \
  /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.3.0/bfx/report/RnaSeq.star.md \
  > report/RnaSeq.star.md
star_report.9089faa6a2c55b3e6016d51622a0e925.mugqic.done
)
star_42_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$star_42_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: picard_merge_sam_files
#-------------------------------------------------------------------------------
STEP=picard_merge_sam_files
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_1_JOB_ID: picard_merge_sam_files.DEX_nonMamm_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_merge_sam_files.DEX_nonMamm_1
JOB_DEPENDENCIES=$star_22_JOB_ID:$star_23_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/picard_merge_sam_files.DEX_nonMamm_1.110d90cc3e49bae83b80599af692ca9f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_merge_sam_files.DEX_nonMamm_1.110d90cc3e49bae83b80599af692ca9f.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/picard.jar MergeSamFiles \
 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/DEX_nonMamm_1/DEX_nonMamm_1_1/Aligned.sortedByCoord.out.bam \
  INPUT=alignment/DEX_nonMamm_1/DEX_nonMamm_1_2/Aligned.sortedByCoord.out.bam \
 OUTPUT=alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.bam \
 MAX_RECORDS_IN_RAM=5750000
picard_merge_sam_files.DEX_nonMamm_1.110d90cc3e49bae83b80599af692ca9f.mugqic.done
)
picard_merge_sam_files_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_2_JOB_ID: picard_merge_sam_files.DEX_shMED1_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_merge_sam_files.DEX_shMED1_1
JOB_DEPENDENCIES=$star_24_JOB_ID:$star_25_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/picard_merge_sam_files.DEX_shMED1_1.5c9715357d07f9deb857629383177995.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_merge_sam_files.DEX_shMED1_1.5c9715357d07f9deb857629383177995.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/picard.jar MergeSamFiles \
 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/DEX_shMED1_1/DEX_shMED1_1_1/Aligned.sortedByCoord.out.bam \
  INPUT=alignment/DEX_shMED1_1/DEX_shMED1_1_2/Aligned.sortedByCoord.out.bam \
 OUTPUT=alignment/DEX_shMED1_1/DEX_shMED1_1.sorted.bam \
 MAX_RECORDS_IN_RAM=5750000
picard_merge_sam_files.DEX_shMED1_1.5c9715357d07f9deb857629383177995.mugqic.done
)
picard_merge_sam_files_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_3_JOB_ID: picard_merge_sam_files.DEX_shNIPBL_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_merge_sam_files.DEX_shNIPBL_1
JOB_DEPENDENCIES=$star_26_JOB_ID:$star_27_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/picard_merge_sam_files.DEX_shNIPBL_1.3136dfb0df5c40c78daf374e63525e48.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_merge_sam_files.DEX_shNIPBL_1.3136dfb0df5c40c78daf374e63525e48.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/picard.jar MergeSamFiles \
 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/DEX_shNIPBL_1/DEX_shNIPBL_1_1/Aligned.sortedByCoord.out.bam \
  INPUT=alignment/DEX_shNIPBL_1/DEX_shNIPBL_1_2/Aligned.sortedByCoord.out.bam \
 OUTPUT=alignment/DEX_shNIPBL_1/DEX_shNIPBL_1.sorted.bam \
 MAX_RECORDS_IN_RAM=5750000
picard_merge_sam_files.DEX_shNIPBL_1.3136dfb0df5c40c78daf374e63525e48.mugqic.done
)
picard_merge_sam_files_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_4_JOB_ID: picard_merge_sam_files.DEX_shSMC1A_3
#-------------------------------------------------------------------------------
JOB_NAME=picard_merge_sam_files.DEX_shSMC1A_3
JOB_DEPENDENCIES=$star_28_JOB_ID:$star_29_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/picard_merge_sam_files.DEX_shSMC1A_3.4ded1f864fc8d96ed4fa7927a40c2778.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_merge_sam_files.DEX_shSMC1A_3.4ded1f864fc8d96ed4fa7927a40c2778.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/picard.jar MergeSamFiles \
 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/DEX_shSMC1A_3/DEX_shSMC1A_1_1/Aligned.sortedByCoord.out.bam \
  INPUT=alignment/DEX_shSMC1A_3/DEX_shSMC1A_1_2/Aligned.sortedByCoord.out.bam \
 OUTPUT=alignment/DEX_shSMC1A_3/DEX_shSMC1A_3.sorted.bam \
 MAX_RECORDS_IN_RAM=5750000
picard_merge_sam_files.DEX_shSMC1A_3.4ded1f864fc8d96ed4fa7927a40c2778.mugqic.done
)
picard_merge_sam_files_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_5_JOB_ID: picard_merge_sam_files.ETOH_nonMamm_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_merge_sam_files.ETOH_nonMamm_1
JOB_DEPENDENCIES=$star_30_JOB_ID:$star_31_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/picard_merge_sam_files.ETOH_nonMamm_1.4a83e0372e395d637177e9fe8d5a79f6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_merge_sam_files.ETOH_nonMamm_1.4a83e0372e395d637177e9fe8d5a79f6.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/picard.jar MergeSamFiles \
 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/ETOH_nonMamm_1/ETOH_nonMamm_1_1/Aligned.sortedByCoord.out.bam \
  INPUT=alignment/ETOH_nonMamm_1/ETOH_nonMamm_1_2/Aligned.sortedByCoord.out.bam \
 OUTPUT=alignment/ETOH_nonMamm_1/ETOH_nonMamm_1.sorted.bam \
 MAX_RECORDS_IN_RAM=5750000
picard_merge_sam_files.ETOH_nonMamm_1.4a83e0372e395d637177e9fe8d5a79f6.mugqic.done
)
picard_merge_sam_files_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_6_JOB_ID: picard_merge_sam_files.ETOH_shMED1_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_merge_sam_files.ETOH_shMED1_1
JOB_DEPENDENCIES=$star_33_JOB_ID:$star_34_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/picard_merge_sam_files.ETOH_shMED1_1.2e6754c50ea035c0003b7297224f7500.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_merge_sam_files.ETOH_shMED1_1.2e6754c50ea035c0003b7297224f7500.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/picard.jar MergeSamFiles \
 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/ETOH_shMED1_1/ETOH_shMED1_1_1/Aligned.sortedByCoord.out.bam \
  INPUT=alignment/ETOH_shMED1_1/ETOH_shMED1_1_2/Aligned.sortedByCoord.out.bam \
 OUTPUT=alignment/ETOH_shMED1_1/ETOH_shMED1_1.sorted.bam \
 MAX_RECORDS_IN_RAM=5750000
picard_merge_sam_files.ETOH_shMED1_1.2e6754c50ea035c0003b7297224f7500.mugqic.done
)
picard_merge_sam_files_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_7_JOB_ID: picard_merge_sam_files.ETOH_shNIPBL_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_merge_sam_files.ETOH_shNIPBL_1
JOB_DEPENDENCIES=$star_36_JOB_ID:$star_37_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/picard_merge_sam_files.ETOH_shNIPBL_1.5beba5714224b37f5f8a75de5b8edb3c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_merge_sam_files.ETOH_shNIPBL_1.5beba5714224b37f5f8a75de5b8edb3c.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/picard.jar MergeSamFiles \
 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1_1/Aligned.sortedByCoord.out.bam \
  INPUT=alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1_2/Aligned.sortedByCoord.out.bam \
 OUTPUT=alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1.sorted.bam \
 MAX_RECORDS_IN_RAM=5750000
picard_merge_sam_files.ETOH_shNIPBL_1.5beba5714224b37f5f8a75de5b8edb3c.mugqic.done
)
picard_merge_sam_files_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_8_JOB_ID: picard_merge_sam_files.ETOH_shSMC1A_3
#-------------------------------------------------------------------------------
JOB_NAME=picard_merge_sam_files.ETOH_shSMC1A_3
JOB_DEPENDENCIES=$star_39_JOB_ID:$star_40_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/picard_merge_sam_files.ETOH_shSMC1A_3.438e6cb7ecf7940949b453f405880eec.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_merge_sam_files.ETOH_shSMC1A_3.438e6cb7ecf7940949b453f405880eec.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/picard.jar MergeSamFiles \
 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/ETOH_shSMC1A_3/ETOH_shSMC1A_1_1/Aligned.sortedByCoord.out.bam \
  INPUT=alignment/ETOH_shSMC1A_3/ETOH_shSMC1A_1_2/Aligned.sortedByCoord.out.bam \
 OUTPUT=alignment/ETOH_shSMC1A_3/ETOH_shSMC1A_3.sorted.bam \
 MAX_RECORDS_IN_RAM=5750000
picard_merge_sam_files.ETOH_shSMC1A_3.438e6cb7ecf7940949b453f405880eec.mugqic.done
)
picard_merge_sam_files_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: picard_sort_sam
#-------------------------------------------------------------------------------
STEP=picard_sort_sam
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: picard_sort_sam_1_JOB_ID: picard_sort_sam.DEX_nonMamm_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_sort_sam.DEX_nonMamm_1
JOB_DEPENDENCIES=$picard_merge_sam_files_1_JOB_ID
JOB_DONE=job_output/picard_sort_sam/picard_sort_sam.DEX_nonMamm_1.ffc2de65d347f3ba2d5a1b007f284b13.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sort_sam.DEX_nonMamm_1.ffc2de65d347f3ba2d5a1b007f284b13.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
 java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.bam \
 OUTPUT=alignment/DEX_nonMamm_1/DEX_nonMamm_1.QueryNameSorted.bam \
 SORT_ORDER=queryname \
 MAX_RECORDS_IN_RAM=5750000
picard_sort_sam.DEX_nonMamm_1.ffc2de65d347f3ba2d5a1b007f284b13.mugqic.done
)
picard_sort_sam_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_sort_sam_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sort_sam_2_JOB_ID: picard_sort_sam.DEX_shMED1_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_sort_sam.DEX_shMED1_1
JOB_DEPENDENCIES=$picard_merge_sam_files_2_JOB_ID
JOB_DONE=job_output/picard_sort_sam/picard_sort_sam.DEX_shMED1_1.4dbb315428759eebe7b1470c81bfed28.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sort_sam.DEX_shMED1_1.4dbb315428759eebe7b1470c81bfed28.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
 java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/DEX_shMED1_1/DEX_shMED1_1.sorted.bam \
 OUTPUT=alignment/DEX_shMED1_1/DEX_shMED1_1.QueryNameSorted.bam \
 SORT_ORDER=queryname \
 MAX_RECORDS_IN_RAM=5750000
picard_sort_sam.DEX_shMED1_1.4dbb315428759eebe7b1470c81bfed28.mugqic.done
)
picard_sort_sam_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_sort_sam_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sort_sam_3_JOB_ID: picard_sort_sam.DEX_shNIPBL_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_sort_sam.DEX_shNIPBL_1
JOB_DEPENDENCIES=$picard_merge_sam_files_3_JOB_ID
JOB_DONE=job_output/picard_sort_sam/picard_sort_sam.DEX_shNIPBL_1.7816c7eb0d0adc356c0ede473cde9417.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sort_sam.DEX_shNIPBL_1.7816c7eb0d0adc356c0ede473cde9417.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
 java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/DEX_shNIPBL_1/DEX_shNIPBL_1.sorted.bam \
 OUTPUT=alignment/DEX_shNIPBL_1/DEX_shNIPBL_1.QueryNameSorted.bam \
 SORT_ORDER=queryname \
 MAX_RECORDS_IN_RAM=5750000
picard_sort_sam.DEX_shNIPBL_1.7816c7eb0d0adc356c0ede473cde9417.mugqic.done
)
picard_sort_sam_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_sort_sam_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sort_sam_4_JOB_ID: picard_sort_sam.DEX_shSMC1A_3
#-------------------------------------------------------------------------------
JOB_NAME=picard_sort_sam.DEX_shSMC1A_3
JOB_DEPENDENCIES=$picard_merge_sam_files_4_JOB_ID
JOB_DONE=job_output/picard_sort_sam/picard_sort_sam.DEX_shSMC1A_3.6ea3ecaf5cc23b08c69b58d9bb03c63d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sort_sam.DEX_shSMC1A_3.6ea3ecaf5cc23b08c69b58d9bb03c63d.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
 java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/DEX_shSMC1A_3/DEX_shSMC1A_3.sorted.bam \
 OUTPUT=alignment/DEX_shSMC1A_3/DEX_shSMC1A_3.QueryNameSorted.bam \
 SORT_ORDER=queryname \
 MAX_RECORDS_IN_RAM=5750000
picard_sort_sam.DEX_shSMC1A_3.6ea3ecaf5cc23b08c69b58d9bb03c63d.mugqic.done
)
picard_sort_sam_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_sort_sam_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sort_sam_5_JOB_ID: picard_sort_sam.ETOH_nonMamm_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_sort_sam.ETOH_nonMamm_1
JOB_DEPENDENCIES=$picard_merge_sam_files_5_JOB_ID
JOB_DONE=job_output/picard_sort_sam/picard_sort_sam.ETOH_nonMamm_1.b2257f9cca60707c328577a114afa112.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sort_sam.ETOH_nonMamm_1.b2257f9cca60707c328577a114afa112.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
 java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/ETOH_nonMamm_1/ETOH_nonMamm_1.sorted.bam \
 OUTPUT=alignment/ETOH_nonMamm_1/ETOH_nonMamm_1.QueryNameSorted.bam \
 SORT_ORDER=queryname \
 MAX_RECORDS_IN_RAM=5750000
picard_sort_sam.ETOH_nonMamm_1.b2257f9cca60707c328577a114afa112.mugqic.done
)
picard_sort_sam_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_sort_sam_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sort_sam_6_JOB_ID: picard_sort_sam.ETOH_nonMamm_2
#-------------------------------------------------------------------------------
JOB_NAME=picard_sort_sam.ETOH_nonMamm_2
JOB_DEPENDENCIES=$star_32_JOB_ID
JOB_DONE=job_output/picard_sort_sam/picard_sort_sam.ETOH_nonMamm_2.c05c982ccf09f55bb35e289e8ccc145f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sort_sam.ETOH_nonMamm_2.c05c982ccf09f55bb35e289e8ccc145f.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
 java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/ETOH_nonMamm_2/ETOH_nonMamm_2.sorted.bam \
 OUTPUT=alignment/ETOH_nonMamm_2/ETOH_nonMamm_2.QueryNameSorted.bam \
 SORT_ORDER=queryname \
 MAX_RECORDS_IN_RAM=5750000
picard_sort_sam.ETOH_nonMamm_2.c05c982ccf09f55bb35e289e8ccc145f.mugqic.done
)
picard_sort_sam_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_sort_sam_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sort_sam_7_JOB_ID: picard_sort_sam.ETOH_shMED1_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_sort_sam.ETOH_shMED1_1
JOB_DEPENDENCIES=$picard_merge_sam_files_6_JOB_ID
JOB_DONE=job_output/picard_sort_sam/picard_sort_sam.ETOH_shMED1_1.a10f20c9f842c77774459448016d86b3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sort_sam.ETOH_shMED1_1.a10f20c9f842c77774459448016d86b3.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
 java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/ETOH_shMED1_1/ETOH_shMED1_1.sorted.bam \
 OUTPUT=alignment/ETOH_shMED1_1/ETOH_shMED1_1.QueryNameSorted.bam \
 SORT_ORDER=queryname \
 MAX_RECORDS_IN_RAM=5750000
picard_sort_sam.ETOH_shMED1_1.a10f20c9f842c77774459448016d86b3.mugqic.done
)
picard_sort_sam_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_sort_sam_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sort_sam_8_JOB_ID: picard_sort_sam.ETOH_shNIPBL_2
#-------------------------------------------------------------------------------
JOB_NAME=picard_sort_sam.ETOH_shNIPBL_2
JOB_DEPENDENCIES=$star_35_JOB_ID
JOB_DONE=job_output/picard_sort_sam/picard_sort_sam.ETOH_shNIPBL_2.1ac89b195e3191672931dd1a738d901b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sort_sam.ETOH_shNIPBL_2.1ac89b195e3191672931dd1a738d901b.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
 java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/ETOH_shNIPBL_2/ETOH_shNIPBL_2.sorted.bam \
 OUTPUT=alignment/ETOH_shNIPBL_2/ETOH_shNIPBL_2.QueryNameSorted.bam \
 SORT_ORDER=queryname \
 MAX_RECORDS_IN_RAM=5750000
picard_sort_sam.ETOH_shNIPBL_2.1ac89b195e3191672931dd1a738d901b.mugqic.done
)
picard_sort_sam_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_sort_sam_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sort_sam_9_JOB_ID: picard_sort_sam.ETOH_shNIPBL_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_sort_sam.ETOH_shNIPBL_1
JOB_DEPENDENCIES=$picard_merge_sam_files_7_JOB_ID
JOB_DONE=job_output/picard_sort_sam/picard_sort_sam.ETOH_shNIPBL_1.e15b055cfa6d1711335503e6235549ab.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sort_sam.ETOH_shNIPBL_1.e15b055cfa6d1711335503e6235549ab.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
 java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1.sorted.bam \
 OUTPUT=alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1.QueryNameSorted.bam \
 SORT_ORDER=queryname \
 MAX_RECORDS_IN_RAM=5750000
picard_sort_sam.ETOH_shNIPBL_1.e15b055cfa6d1711335503e6235549ab.mugqic.done
)
picard_sort_sam_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_sort_sam_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sort_sam_10_JOB_ID: picard_sort_sam.ETOH_shMED1_2
#-------------------------------------------------------------------------------
JOB_NAME=picard_sort_sam.ETOH_shMED1_2
JOB_DEPENDENCIES=$star_38_JOB_ID
JOB_DONE=job_output/picard_sort_sam/picard_sort_sam.ETOH_shMED1_2.89f64e4307dea61232e9622492b59bca.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sort_sam.ETOH_shMED1_2.89f64e4307dea61232e9622492b59bca.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
 java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/ETOH_shMED1_2/ETOH_shMED1_2.sorted.bam \
 OUTPUT=alignment/ETOH_shMED1_2/ETOH_shMED1_2.QueryNameSorted.bam \
 SORT_ORDER=queryname \
 MAX_RECORDS_IN_RAM=5750000
picard_sort_sam.ETOH_shMED1_2.89f64e4307dea61232e9622492b59bca.mugqic.done
)
picard_sort_sam_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_sort_sam_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sort_sam_11_JOB_ID: picard_sort_sam.ETOH_shSMC1A_3
#-------------------------------------------------------------------------------
JOB_NAME=picard_sort_sam.ETOH_shSMC1A_3
JOB_DEPENDENCIES=$picard_merge_sam_files_8_JOB_ID
JOB_DONE=job_output/picard_sort_sam/picard_sort_sam.ETOH_shSMC1A_3.ac4d076d15fb32862752240206d4f148.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sort_sam.ETOH_shSMC1A_3.ac4d076d15fb32862752240206d4f148.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
 java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/ETOH_shSMC1A_3/ETOH_shSMC1A_3.sorted.bam \
 OUTPUT=alignment/ETOH_shSMC1A_3/ETOH_shSMC1A_3.QueryNameSorted.bam \
 SORT_ORDER=queryname \
 MAX_RECORDS_IN_RAM=5750000
picard_sort_sam.ETOH_shSMC1A_3.ac4d076d15fb32862752240206d4f148.mugqic.done
)
picard_sort_sam_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_sort_sam_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sort_sam_12_JOB_ID: picard_sort_sam.ETOH_shSMC1A_2
#-------------------------------------------------------------------------------
JOB_NAME=picard_sort_sam.ETOH_shSMC1A_2
JOB_DEPENDENCIES=$star_41_JOB_ID
JOB_DONE=job_output/picard_sort_sam/picard_sort_sam.ETOH_shSMC1A_2.5fc1164616c07b95319c7528218d6335.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sort_sam.ETOH_shSMC1A_2.5fc1164616c07b95319c7528218d6335.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
 java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/ETOH_shSMC1A_2/ETOH_shSMC1A_2.sorted.bam \
 OUTPUT=alignment/ETOH_shSMC1A_2/ETOH_shSMC1A_2.QueryNameSorted.bam \
 SORT_ORDER=queryname \
 MAX_RECORDS_IN_RAM=5750000
picard_sort_sam.ETOH_shSMC1A_2.5fc1164616c07b95319c7528218d6335.mugqic.done
)
picard_sort_sam_12_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_sort_sam_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: picard_mark_duplicates
#-------------------------------------------------------------------------------
STEP=picard_mark_duplicates
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_1_JOB_ID: picard_mark_duplicates.DEX_nonMamm_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.DEX_nonMamm_1
JOB_DEPENDENCIES=$picard_merge_sam_files_1_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.DEX_nonMamm_1.2d108e96fdefd1662a9140bce3fd1ddb.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.DEX_nonMamm_1.2d108e96fdefd1662a9140bce3fd1ddb.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx14G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.bam \
 OUTPUT=alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.bam \
 METRICS_FILE=alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.metrics \
 MAX_RECORDS_IN_RAM=3500000
picard_mark_duplicates.DEX_nonMamm_1.2d108e96fdefd1662a9140bce3fd1ddb.mugqic.done
)
picard_mark_duplicates_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_2_JOB_ID: picard_mark_duplicates.DEX_shMED1_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.DEX_shMED1_1
JOB_DEPENDENCIES=$picard_merge_sam_files_2_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.DEX_shMED1_1.01df24bd54641f72d1ac734d963695d0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.DEX_shMED1_1.01df24bd54641f72d1ac734d963695d0.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx14G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/DEX_shMED1_1/DEX_shMED1_1.sorted.bam \
 OUTPUT=alignment/DEX_shMED1_1/DEX_shMED1_1.sorted.mdup.bam \
 METRICS_FILE=alignment/DEX_shMED1_1/DEX_shMED1_1.sorted.mdup.metrics \
 MAX_RECORDS_IN_RAM=3500000
picard_mark_duplicates.DEX_shMED1_1.01df24bd54641f72d1ac734d963695d0.mugqic.done
)
picard_mark_duplicates_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_3_JOB_ID: picard_mark_duplicates.DEX_shNIPBL_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.DEX_shNIPBL_1
JOB_DEPENDENCIES=$picard_merge_sam_files_3_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.DEX_shNIPBL_1.b718e8e46b4dba9de29494cd57857192.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.DEX_shNIPBL_1.b718e8e46b4dba9de29494cd57857192.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx14G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/DEX_shNIPBL_1/DEX_shNIPBL_1.sorted.bam \
 OUTPUT=alignment/DEX_shNIPBL_1/DEX_shNIPBL_1.sorted.mdup.bam \
 METRICS_FILE=alignment/DEX_shNIPBL_1/DEX_shNIPBL_1.sorted.mdup.metrics \
 MAX_RECORDS_IN_RAM=3500000
picard_mark_duplicates.DEX_shNIPBL_1.b718e8e46b4dba9de29494cd57857192.mugqic.done
)
picard_mark_duplicates_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_4_JOB_ID: picard_mark_duplicates.DEX_shSMC1A_3
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.DEX_shSMC1A_3
JOB_DEPENDENCIES=$picard_merge_sam_files_4_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.DEX_shSMC1A_3.9cdb21af7bd448390708e52e7d0c035d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.DEX_shSMC1A_3.9cdb21af7bd448390708e52e7d0c035d.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx14G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/DEX_shSMC1A_3/DEX_shSMC1A_3.sorted.bam \
 OUTPUT=alignment/DEX_shSMC1A_3/DEX_shSMC1A_3.sorted.mdup.bam \
 METRICS_FILE=alignment/DEX_shSMC1A_3/DEX_shSMC1A_3.sorted.mdup.metrics \
 MAX_RECORDS_IN_RAM=3500000
picard_mark_duplicates.DEX_shSMC1A_3.9cdb21af7bd448390708e52e7d0c035d.mugqic.done
)
picard_mark_duplicates_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_5_JOB_ID: picard_mark_duplicates.ETOH_nonMamm_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.ETOH_nonMamm_1
JOB_DEPENDENCIES=$picard_merge_sam_files_5_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.ETOH_nonMamm_1.10f0c0128fe4c427596bc9aa04d80a9c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.ETOH_nonMamm_1.10f0c0128fe4c427596bc9aa04d80a9c.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx14G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/ETOH_nonMamm_1/ETOH_nonMamm_1.sorted.bam \
 OUTPUT=alignment/ETOH_nonMamm_1/ETOH_nonMamm_1.sorted.mdup.bam \
 METRICS_FILE=alignment/ETOH_nonMamm_1/ETOH_nonMamm_1.sorted.mdup.metrics \
 MAX_RECORDS_IN_RAM=3500000
picard_mark_duplicates.ETOH_nonMamm_1.10f0c0128fe4c427596bc9aa04d80a9c.mugqic.done
)
picard_mark_duplicates_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_6_JOB_ID: picard_mark_duplicates.ETOH_nonMamm_2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.ETOH_nonMamm_2
JOB_DEPENDENCIES=$star_32_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.ETOH_nonMamm_2.5670ec632283c361b7171a4e80d205c5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.ETOH_nonMamm_2.5670ec632283c361b7171a4e80d205c5.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx14G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/ETOH_nonMamm_2/ETOH_nonMamm_2.sorted.bam \
 OUTPUT=alignment/ETOH_nonMamm_2/ETOH_nonMamm_2.sorted.mdup.bam \
 METRICS_FILE=alignment/ETOH_nonMamm_2/ETOH_nonMamm_2.sorted.mdup.metrics \
 MAX_RECORDS_IN_RAM=3500000
picard_mark_duplicates.ETOH_nonMamm_2.5670ec632283c361b7171a4e80d205c5.mugqic.done
)
picard_mark_duplicates_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_7_JOB_ID: picard_mark_duplicates.ETOH_shMED1_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.ETOH_shMED1_1
JOB_DEPENDENCIES=$picard_merge_sam_files_6_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.ETOH_shMED1_1.0b5bc0b37da1bee558f730e4aa01d60c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.ETOH_shMED1_1.0b5bc0b37da1bee558f730e4aa01d60c.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx14G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/ETOH_shMED1_1/ETOH_shMED1_1.sorted.bam \
 OUTPUT=alignment/ETOH_shMED1_1/ETOH_shMED1_1.sorted.mdup.bam \
 METRICS_FILE=alignment/ETOH_shMED1_1/ETOH_shMED1_1.sorted.mdup.metrics \
 MAX_RECORDS_IN_RAM=3500000
picard_mark_duplicates.ETOH_shMED1_1.0b5bc0b37da1bee558f730e4aa01d60c.mugqic.done
)
picard_mark_duplicates_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_8_JOB_ID: picard_mark_duplicates.ETOH_shNIPBL_2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.ETOH_shNIPBL_2
JOB_DEPENDENCIES=$star_35_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.ETOH_shNIPBL_2.b3d3e413a215ac0412da8464d4a32136.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.ETOH_shNIPBL_2.b3d3e413a215ac0412da8464d4a32136.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx14G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/ETOH_shNIPBL_2/ETOH_shNIPBL_2.sorted.bam \
 OUTPUT=alignment/ETOH_shNIPBL_2/ETOH_shNIPBL_2.sorted.mdup.bam \
 METRICS_FILE=alignment/ETOH_shNIPBL_2/ETOH_shNIPBL_2.sorted.mdup.metrics \
 MAX_RECORDS_IN_RAM=3500000
picard_mark_duplicates.ETOH_shNIPBL_2.b3d3e413a215ac0412da8464d4a32136.mugqic.done
)
picard_mark_duplicates_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_9_JOB_ID: picard_mark_duplicates.ETOH_shNIPBL_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.ETOH_shNIPBL_1
JOB_DEPENDENCIES=$picard_merge_sam_files_7_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.ETOH_shNIPBL_1.e2933699572961d661c2290acddd5c06.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.ETOH_shNIPBL_1.e2933699572961d661c2290acddd5c06.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx14G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1.sorted.bam \
 OUTPUT=alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1.sorted.mdup.bam \
 METRICS_FILE=alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1.sorted.mdup.metrics \
 MAX_RECORDS_IN_RAM=3500000
picard_mark_duplicates.ETOH_shNIPBL_1.e2933699572961d661c2290acddd5c06.mugqic.done
)
picard_mark_duplicates_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_10_JOB_ID: picard_mark_duplicates.ETOH_shMED1_2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.ETOH_shMED1_2
JOB_DEPENDENCIES=$star_38_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.ETOH_shMED1_2.50533eb7dcf0e7719a5b8aa42544d451.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.ETOH_shMED1_2.50533eb7dcf0e7719a5b8aa42544d451.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx14G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/ETOH_shMED1_2/ETOH_shMED1_2.sorted.bam \
 OUTPUT=alignment/ETOH_shMED1_2/ETOH_shMED1_2.sorted.mdup.bam \
 METRICS_FILE=alignment/ETOH_shMED1_2/ETOH_shMED1_2.sorted.mdup.metrics \
 MAX_RECORDS_IN_RAM=3500000
picard_mark_duplicates.ETOH_shMED1_2.50533eb7dcf0e7719a5b8aa42544d451.mugqic.done
)
picard_mark_duplicates_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_11_JOB_ID: picard_mark_duplicates.ETOH_shSMC1A_3
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.ETOH_shSMC1A_3
JOB_DEPENDENCIES=$picard_merge_sam_files_8_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.ETOH_shSMC1A_3.8470c9a65936a38db2ded1f78086dabc.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.ETOH_shSMC1A_3.8470c9a65936a38db2ded1f78086dabc.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx14G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/ETOH_shSMC1A_3/ETOH_shSMC1A_3.sorted.bam \
 OUTPUT=alignment/ETOH_shSMC1A_3/ETOH_shSMC1A_3.sorted.mdup.bam \
 METRICS_FILE=alignment/ETOH_shSMC1A_3/ETOH_shSMC1A_3.sorted.mdup.metrics \
 MAX_RECORDS_IN_RAM=3500000
picard_mark_duplicates.ETOH_shSMC1A_3.8470c9a65936a38db2ded1f78086dabc.mugqic.done
)
picard_mark_duplicates_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_12_JOB_ID: picard_mark_duplicates.ETOH_shSMC1A_2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.ETOH_shSMC1A_2
JOB_DEPENDENCIES=$star_41_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.ETOH_shSMC1A_2.311141c563300ce27d4bd99256436506.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.ETOH_shSMC1A_2.311141c563300ce27d4bd99256436506.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx14G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/ETOH_shSMC1A_2/ETOH_shSMC1A_2.sorted.bam \
 OUTPUT=alignment/ETOH_shSMC1A_2/ETOH_shSMC1A_2.sorted.mdup.bam \
 METRICS_FILE=alignment/ETOH_shSMC1A_2/ETOH_shSMC1A_2.sorted.mdup.metrics \
 MAX_RECORDS_IN_RAM=3500000
picard_mark_duplicates.ETOH_shSMC1A_2.311141c563300ce27d4bd99256436506.mugqic.done
)
picard_mark_duplicates_12_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: picard_rna_metrics
#-------------------------------------------------------------------------------
STEP=picard_rna_metrics
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: picard_rna_metrics_1_JOB_ID: picard_rna_metrics.DEX_nonMamm_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_rna_metrics.DEX_nonMamm_1
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID
JOB_DONE=job_output/picard_rna_metrics/picard_rna_metrics.DEX_nonMamm_1.70cd462844f2f6762bd8d78d02a766d9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_rna_metrics.DEX_nonMamm_1.70cd462844f2f6762bd8d78d02a766d9.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 mugqic/R_Bioconductor/3.2.3_3.2 && \
mkdir -p metrics/DEX_nonMamm_1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=4 -Xmx27G -jar $PICARD_HOME/picard.jar CollectMultipleMetrics \
 PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT \
 TMP_DIR=/gs/scratch/$USER \
 REFERENCE_SEQUENCE=/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
 INPUT=alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.bam \
 OUTPUT=metrics/DEX_nonMamm_1/DEX_nonMamm_1 \
 MAX_RECORDS_IN_RAM=5750000 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=4 -Xmx27G -jar $PICARD_HOME/picard.jar CollectRnaSeqMetrics \
 VALIDATION_STRINGENCY=SILENT  \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.bam \
 OUTPUT=metrics/DEX_nonMamm_1/DEX_nonMamm_1.picard_rna_metrics \
 REF_FLAT=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.ref_flat.tsv \
 STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
 MINIMUM_LENGTH=200 \
 REFERENCE_SEQUENCE=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
 MAX_RECORDS_IN_RAM=5750000
picard_rna_metrics.DEX_nonMamm_1.70cd462844f2f6762bd8d78d02a766d9.mugqic.done
)
picard_rna_metrics_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_rna_metrics_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_rna_metrics_2_JOB_ID: picard_rna_metrics.DEX_shMED1_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_rna_metrics.DEX_shMED1_1
JOB_DEPENDENCIES=$picard_mark_duplicates_2_JOB_ID
JOB_DONE=job_output/picard_rna_metrics/picard_rna_metrics.DEX_shMED1_1.201ce1ac6abd995746639b5e81e830dd.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_rna_metrics.DEX_shMED1_1.201ce1ac6abd995746639b5e81e830dd.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 mugqic/R_Bioconductor/3.2.3_3.2 && \
mkdir -p metrics/DEX_shMED1_1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=4 -Xmx27G -jar $PICARD_HOME/picard.jar CollectMultipleMetrics \
 PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT \
 TMP_DIR=/gs/scratch/$USER \
 REFERENCE_SEQUENCE=/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
 INPUT=alignment/DEX_shMED1_1/DEX_shMED1_1.sorted.mdup.bam \
 OUTPUT=metrics/DEX_shMED1_1/DEX_shMED1_1 \
 MAX_RECORDS_IN_RAM=5750000 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=4 -Xmx27G -jar $PICARD_HOME/picard.jar CollectRnaSeqMetrics \
 VALIDATION_STRINGENCY=SILENT  \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/DEX_shMED1_1/DEX_shMED1_1.sorted.mdup.bam \
 OUTPUT=metrics/DEX_shMED1_1/DEX_shMED1_1.picard_rna_metrics \
 REF_FLAT=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.ref_flat.tsv \
 STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
 MINIMUM_LENGTH=200 \
 REFERENCE_SEQUENCE=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
 MAX_RECORDS_IN_RAM=5750000
picard_rna_metrics.DEX_shMED1_1.201ce1ac6abd995746639b5e81e830dd.mugqic.done
)
picard_rna_metrics_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_rna_metrics_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_rna_metrics_3_JOB_ID: picard_rna_metrics.DEX_shNIPBL_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_rna_metrics.DEX_shNIPBL_1
JOB_DEPENDENCIES=$picard_mark_duplicates_3_JOB_ID
JOB_DONE=job_output/picard_rna_metrics/picard_rna_metrics.DEX_shNIPBL_1.091c2c821518e612fac38f0682e8f552.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_rna_metrics.DEX_shNIPBL_1.091c2c821518e612fac38f0682e8f552.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 mugqic/R_Bioconductor/3.2.3_3.2 && \
mkdir -p metrics/DEX_shNIPBL_1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=4 -Xmx27G -jar $PICARD_HOME/picard.jar CollectMultipleMetrics \
 PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT \
 TMP_DIR=/gs/scratch/$USER \
 REFERENCE_SEQUENCE=/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
 INPUT=alignment/DEX_shNIPBL_1/DEX_shNIPBL_1.sorted.mdup.bam \
 OUTPUT=metrics/DEX_shNIPBL_1/DEX_shNIPBL_1 \
 MAX_RECORDS_IN_RAM=5750000 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=4 -Xmx27G -jar $PICARD_HOME/picard.jar CollectRnaSeqMetrics \
 VALIDATION_STRINGENCY=SILENT  \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/DEX_shNIPBL_1/DEX_shNIPBL_1.sorted.mdup.bam \
 OUTPUT=metrics/DEX_shNIPBL_1/DEX_shNIPBL_1.picard_rna_metrics \
 REF_FLAT=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.ref_flat.tsv \
 STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
 MINIMUM_LENGTH=200 \
 REFERENCE_SEQUENCE=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
 MAX_RECORDS_IN_RAM=5750000
picard_rna_metrics.DEX_shNIPBL_1.091c2c821518e612fac38f0682e8f552.mugqic.done
)
picard_rna_metrics_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_rna_metrics_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_rna_metrics_4_JOB_ID: picard_rna_metrics.DEX_shSMC1A_3
#-------------------------------------------------------------------------------
JOB_NAME=picard_rna_metrics.DEX_shSMC1A_3
JOB_DEPENDENCIES=$picard_mark_duplicates_4_JOB_ID
JOB_DONE=job_output/picard_rna_metrics/picard_rna_metrics.DEX_shSMC1A_3.bca84c36e57e1187f8511392d56836ab.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_rna_metrics.DEX_shSMC1A_3.bca84c36e57e1187f8511392d56836ab.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 mugqic/R_Bioconductor/3.2.3_3.2 && \
mkdir -p metrics/DEX_shSMC1A_3 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=4 -Xmx27G -jar $PICARD_HOME/picard.jar CollectMultipleMetrics \
 PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT \
 TMP_DIR=/gs/scratch/$USER \
 REFERENCE_SEQUENCE=/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
 INPUT=alignment/DEX_shSMC1A_3/DEX_shSMC1A_3.sorted.mdup.bam \
 OUTPUT=metrics/DEX_shSMC1A_3/DEX_shSMC1A_3 \
 MAX_RECORDS_IN_RAM=5750000 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=4 -Xmx27G -jar $PICARD_HOME/picard.jar CollectRnaSeqMetrics \
 VALIDATION_STRINGENCY=SILENT  \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/DEX_shSMC1A_3/DEX_shSMC1A_3.sorted.mdup.bam \
 OUTPUT=metrics/DEX_shSMC1A_3/DEX_shSMC1A_3.picard_rna_metrics \
 REF_FLAT=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.ref_flat.tsv \
 STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
 MINIMUM_LENGTH=200 \
 REFERENCE_SEQUENCE=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
 MAX_RECORDS_IN_RAM=5750000
picard_rna_metrics.DEX_shSMC1A_3.bca84c36e57e1187f8511392d56836ab.mugqic.done
)
picard_rna_metrics_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_rna_metrics_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_rna_metrics_5_JOB_ID: picard_rna_metrics.ETOH_nonMamm_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_rna_metrics.ETOH_nonMamm_1
JOB_DEPENDENCIES=$picard_mark_duplicates_5_JOB_ID
JOB_DONE=job_output/picard_rna_metrics/picard_rna_metrics.ETOH_nonMamm_1.297988a11e8de4db57fd827af95ff8dd.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_rna_metrics.ETOH_nonMamm_1.297988a11e8de4db57fd827af95ff8dd.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 mugqic/R_Bioconductor/3.2.3_3.2 && \
mkdir -p metrics/ETOH_nonMamm_1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=4 -Xmx27G -jar $PICARD_HOME/picard.jar CollectMultipleMetrics \
 PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT \
 TMP_DIR=/gs/scratch/$USER \
 REFERENCE_SEQUENCE=/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
 INPUT=alignment/ETOH_nonMamm_1/ETOH_nonMamm_1.sorted.mdup.bam \
 OUTPUT=metrics/ETOH_nonMamm_1/ETOH_nonMamm_1 \
 MAX_RECORDS_IN_RAM=5750000 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=4 -Xmx27G -jar $PICARD_HOME/picard.jar CollectRnaSeqMetrics \
 VALIDATION_STRINGENCY=SILENT  \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/ETOH_nonMamm_1/ETOH_nonMamm_1.sorted.mdup.bam \
 OUTPUT=metrics/ETOH_nonMamm_1/ETOH_nonMamm_1.picard_rna_metrics \
 REF_FLAT=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.ref_flat.tsv \
 STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
 MINIMUM_LENGTH=200 \
 REFERENCE_SEQUENCE=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
 MAX_RECORDS_IN_RAM=5750000
picard_rna_metrics.ETOH_nonMamm_1.297988a11e8de4db57fd827af95ff8dd.mugqic.done
)
picard_rna_metrics_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_rna_metrics_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_rna_metrics_6_JOB_ID: picard_rna_metrics.ETOH_nonMamm_2
#-------------------------------------------------------------------------------
JOB_NAME=picard_rna_metrics.ETOH_nonMamm_2
JOB_DEPENDENCIES=$picard_mark_duplicates_6_JOB_ID
JOB_DONE=job_output/picard_rna_metrics/picard_rna_metrics.ETOH_nonMamm_2.25cf60f6b2c7f30ae258b69a76e3b18e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_rna_metrics.ETOH_nonMamm_2.25cf60f6b2c7f30ae258b69a76e3b18e.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 mugqic/R_Bioconductor/3.2.3_3.2 && \
mkdir -p metrics/ETOH_nonMamm_2 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=4 -Xmx27G -jar $PICARD_HOME/picard.jar CollectMultipleMetrics \
 PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT \
 TMP_DIR=/gs/scratch/$USER \
 REFERENCE_SEQUENCE=/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
 INPUT=alignment/ETOH_nonMamm_2/ETOH_nonMamm_2.sorted.mdup.bam \
 OUTPUT=metrics/ETOH_nonMamm_2/ETOH_nonMamm_2 \
 MAX_RECORDS_IN_RAM=5750000 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=4 -Xmx27G -jar $PICARD_HOME/picard.jar CollectRnaSeqMetrics \
 VALIDATION_STRINGENCY=SILENT  \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/ETOH_nonMamm_2/ETOH_nonMamm_2.sorted.mdup.bam \
 OUTPUT=metrics/ETOH_nonMamm_2/ETOH_nonMamm_2.picard_rna_metrics \
 REF_FLAT=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.ref_flat.tsv \
 STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
 MINIMUM_LENGTH=200 \
 REFERENCE_SEQUENCE=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
 MAX_RECORDS_IN_RAM=5750000
picard_rna_metrics.ETOH_nonMamm_2.25cf60f6b2c7f30ae258b69a76e3b18e.mugqic.done
)
picard_rna_metrics_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_rna_metrics_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_rna_metrics_7_JOB_ID: picard_rna_metrics.ETOH_shMED1_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_rna_metrics.ETOH_shMED1_1
JOB_DEPENDENCIES=$picard_mark_duplicates_7_JOB_ID
JOB_DONE=job_output/picard_rna_metrics/picard_rna_metrics.ETOH_shMED1_1.1a0d5720584690a2a4b490f7fbb93e5b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_rna_metrics.ETOH_shMED1_1.1a0d5720584690a2a4b490f7fbb93e5b.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 mugqic/R_Bioconductor/3.2.3_3.2 && \
mkdir -p metrics/ETOH_shMED1_1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=4 -Xmx27G -jar $PICARD_HOME/picard.jar CollectMultipleMetrics \
 PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT \
 TMP_DIR=/gs/scratch/$USER \
 REFERENCE_SEQUENCE=/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
 INPUT=alignment/ETOH_shMED1_1/ETOH_shMED1_1.sorted.mdup.bam \
 OUTPUT=metrics/ETOH_shMED1_1/ETOH_shMED1_1 \
 MAX_RECORDS_IN_RAM=5750000 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=4 -Xmx27G -jar $PICARD_HOME/picard.jar CollectRnaSeqMetrics \
 VALIDATION_STRINGENCY=SILENT  \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/ETOH_shMED1_1/ETOH_shMED1_1.sorted.mdup.bam \
 OUTPUT=metrics/ETOH_shMED1_1/ETOH_shMED1_1.picard_rna_metrics \
 REF_FLAT=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.ref_flat.tsv \
 STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
 MINIMUM_LENGTH=200 \
 REFERENCE_SEQUENCE=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
 MAX_RECORDS_IN_RAM=5750000
picard_rna_metrics.ETOH_shMED1_1.1a0d5720584690a2a4b490f7fbb93e5b.mugqic.done
)
picard_rna_metrics_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_rna_metrics_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_rna_metrics_8_JOB_ID: picard_rna_metrics.ETOH_shNIPBL_2
#-------------------------------------------------------------------------------
JOB_NAME=picard_rna_metrics.ETOH_shNIPBL_2
JOB_DEPENDENCIES=$picard_mark_duplicates_8_JOB_ID
JOB_DONE=job_output/picard_rna_metrics/picard_rna_metrics.ETOH_shNIPBL_2.908ed5634f17ccd598ac4869ab272e51.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_rna_metrics.ETOH_shNIPBL_2.908ed5634f17ccd598ac4869ab272e51.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 mugqic/R_Bioconductor/3.2.3_3.2 && \
mkdir -p metrics/ETOH_shNIPBL_2 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=4 -Xmx27G -jar $PICARD_HOME/picard.jar CollectMultipleMetrics \
 PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT \
 TMP_DIR=/gs/scratch/$USER \
 REFERENCE_SEQUENCE=/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
 INPUT=alignment/ETOH_shNIPBL_2/ETOH_shNIPBL_2.sorted.mdup.bam \
 OUTPUT=metrics/ETOH_shNIPBL_2/ETOH_shNIPBL_2 \
 MAX_RECORDS_IN_RAM=5750000 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=4 -Xmx27G -jar $PICARD_HOME/picard.jar CollectRnaSeqMetrics \
 VALIDATION_STRINGENCY=SILENT  \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/ETOH_shNIPBL_2/ETOH_shNIPBL_2.sorted.mdup.bam \
 OUTPUT=metrics/ETOH_shNIPBL_2/ETOH_shNIPBL_2.picard_rna_metrics \
 REF_FLAT=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.ref_flat.tsv \
 STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
 MINIMUM_LENGTH=200 \
 REFERENCE_SEQUENCE=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
 MAX_RECORDS_IN_RAM=5750000
picard_rna_metrics.ETOH_shNIPBL_2.908ed5634f17ccd598ac4869ab272e51.mugqic.done
)
picard_rna_metrics_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_rna_metrics_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_rna_metrics_9_JOB_ID: picard_rna_metrics.ETOH_shNIPBL_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_rna_metrics.ETOH_shNIPBL_1
JOB_DEPENDENCIES=$picard_mark_duplicates_9_JOB_ID
JOB_DONE=job_output/picard_rna_metrics/picard_rna_metrics.ETOH_shNIPBL_1.b1ce01cab910490fecba0c7e17d914e8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_rna_metrics.ETOH_shNIPBL_1.b1ce01cab910490fecba0c7e17d914e8.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 mugqic/R_Bioconductor/3.2.3_3.2 && \
mkdir -p metrics/ETOH_shNIPBL_1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=4 -Xmx27G -jar $PICARD_HOME/picard.jar CollectMultipleMetrics \
 PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT \
 TMP_DIR=/gs/scratch/$USER \
 REFERENCE_SEQUENCE=/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
 INPUT=alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1.sorted.mdup.bam \
 OUTPUT=metrics/ETOH_shNIPBL_1/ETOH_shNIPBL_1 \
 MAX_RECORDS_IN_RAM=5750000 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=4 -Xmx27G -jar $PICARD_HOME/picard.jar CollectRnaSeqMetrics \
 VALIDATION_STRINGENCY=SILENT  \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1.sorted.mdup.bam \
 OUTPUT=metrics/ETOH_shNIPBL_1/ETOH_shNIPBL_1.picard_rna_metrics \
 REF_FLAT=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.ref_flat.tsv \
 STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
 MINIMUM_LENGTH=200 \
 REFERENCE_SEQUENCE=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
 MAX_RECORDS_IN_RAM=5750000
picard_rna_metrics.ETOH_shNIPBL_1.b1ce01cab910490fecba0c7e17d914e8.mugqic.done
)
picard_rna_metrics_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_rna_metrics_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_rna_metrics_10_JOB_ID: picard_rna_metrics.ETOH_shMED1_2
#-------------------------------------------------------------------------------
JOB_NAME=picard_rna_metrics.ETOH_shMED1_2
JOB_DEPENDENCIES=$picard_mark_duplicates_10_JOB_ID
JOB_DONE=job_output/picard_rna_metrics/picard_rna_metrics.ETOH_shMED1_2.045c60ede89b4b3a7346f3e9f4564160.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_rna_metrics.ETOH_shMED1_2.045c60ede89b4b3a7346f3e9f4564160.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 mugqic/R_Bioconductor/3.2.3_3.2 && \
mkdir -p metrics/ETOH_shMED1_2 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=4 -Xmx27G -jar $PICARD_HOME/picard.jar CollectMultipleMetrics \
 PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT \
 TMP_DIR=/gs/scratch/$USER \
 REFERENCE_SEQUENCE=/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
 INPUT=alignment/ETOH_shMED1_2/ETOH_shMED1_2.sorted.mdup.bam \
 OUTPUT=metrics/ETOH_shMED1_2/ETOH_shMED1_2 \
 MAX_RECORDS_IN_RAM=5750000 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=4 -Xmx27G -jar $PICARD_HOME/picard.jar CollectRnaSeqMetrics \
 VALIDATION_STRINGENCY=SILENT  \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/ETOH_shMED1_2/ETOH_shMED1_2.sorted.mdup.bam \
 OUTPUT=metrics/ETOH_shMED1_2/ETOH_shMED1_2.picard_rna_metrics \
 REF_FLAT=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.ref_flat.tsv \
 STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
 MINIMUM_LENGTH=200 \
 REFERENCE_SEQUENCE=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
 MAX_RECORDS_IN_RAM=5750000
picard_rna_metrics.ETOH_shMED1_2.045c60ede89b4b3a7346f3e9f4564160.mugqic.done
)
picard_rna_metrics_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_rna_metrics_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_rna_metrics_11_JOB_ID: picard_rna_metrics.ETOH_shSMC1A_3
#-------------------------------------------------------------------------------
JOB_NAME=picard_rna_metrics.ETOH_shSMC1A_3
JOB_DEPENDENCIES=$picard_mark_duplicates_11_JOB_ID
JOB_DONE=job_output/picard_rna_metrics/picard_rna_metrics.ETOH_shSMC1A_3.b15241ccec600777d9266932aafb94a4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_rna_metrics.ETOH_shSMC1A_3.b15241ccec600777d9266932aafb94a4.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 mugqic/R_Bioconductor/3.2.3_3.2 && \
mkdir -p metrics/ETOH_shSMC1A_3 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=4 -Xmx27G -jar $PICARD_HOME/picard.jar CollectMultipleMetrics \
 PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT \
 TMP_DIR=/gs/scratch/$USER \
 REFERENCE_SEQUENCE=/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
 INPUT=alignment/ETOH_shSMC1A_3/ETOH_shSMC1A_3.sorted.mdup.bam \
 OUTPUT=metrics/ETOH_shSMC1A_3/ETOH_shSMC1A_3 \
 MAX_RECORDS_IN_RAM=5750000 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=4 -Xmx27G -jar $PICARD_HOME/picard.jar CollectRnaSeqMetrics \
 VALIDATION_STRINGENCY=SILENT  \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/ETOH_shSMC1A_3/ETOH_shSMC1A_3.sorted.mdup.bam \
 OUTPUT=metrics/ETOH_shSMC1A_3/ETOH_shSMC1A_3.picard_rna_metrics \
 REF_FLAT=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.ref_flat.tsv \
 STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
 MINIMUM_LENGTH=200 \
 REFERENCE_SEQUENCE=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
 MAX_RECORDS_IN_RAM=5750000
picard_rna_metrics.ETOH_shSMC1A_3.b15241ccec600777d9266932aafb94a4.mugqic.done
)
picard_rna_metrics_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_rna_metrics_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_rna_metrics_12_JOB_ID: picard_rna_metrics.ETOH_shSMC1A_2
#-------------------------------------------------------------------------------
JOB_NAME=picard_rna_metrics.ETOH_shSMC1A_2
JOB_DEPENDENCIES=$picard_mark_duplicates_12_JOB_ID
JOB_DONE=job_output/picard_rna_metrics/picard_rna_metrics.ETOH_shSMC1A_2.2fd61173ffd8231090221f5c3f4760ef.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_rna_metrics.ETOH_shSMC1A_2.2fd61173ffd8231090221f5c3f4760ef.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 mugqic/R_Bioconductor/3.2.3_3.2 && \
mkdir -p metrics/ETOH_shSMC1A_2 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=4 -Xmx27G -jar $PICARD_HOME/picard.jar CollectMultipleMetrics \
 PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT \
 TMP_DIR=/gs/scratch/$USER \
 REFERENCE_SEQUENCE=/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
 INPUT=alignment/ETOH_shSMC1A_2/ETOH_shSMC1A_2.sorted.mdup.bam \
 OUTPUT=metrics/ETOH_shSMC1A_2/ETOH_shSMC1A_2 \
 MAX_RECORDS_IN_RAM=5750000 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=4 -Xmx27G -jar $PICARD_HOME/picard.jar CollectRnaSeqMetrics \
 VALIDATION_STRINGENCY=SILENT  \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/ETOH_shSMC1A_2/ETOH_shSMC1A_2.sorted.mdup.bam \
 OUTPUT=metrics/ETOH_shSMC1A_2/ETOH_shSMC1A_2.picard_rna_metrics \
 REF_FLAT=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.ref_flat.tsv \
 STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
 MINIMUM_LENGTH=200 \
 REFERENCE_SEQUENCE=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
 MAX_RECORDS_IN_RAM=5750000
picard_rna_metrics.ETOH_shSMC1A_2.2fd61173ffd8231090221f5c3f4760ef.mugqic.done
)
picard_rna_metrics_12_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_rna_metrics_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: estimate_ribosomal_rna
#-------------------------------------------------------------------------------
STEP=estimate_ribosomal_rna
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: estimate_ribosomal_rna_1_JOB_ID: bwa_mem_rRNA.DEX_nonMamm_1_1
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_rRNA.DEX_nonMamm_1_1
JOB_DEPENDENCIES=$star_22_JOB_ID
JOB_DONE=job_output/estimate_ribosomal_rna/bwa_mem_rRNA.DEX_nonMamm_1_1.2ac311dc791ce64a9f80306d262820f1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_rRNA.DEX_nonMamm_1_1.2ac311dc791ce64a9f80306d262820f1.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/bvatools/1.6 mugqic/bwa/0.7.12 mugqic/picard/2.0.1 mugqic/mugqic_tools/2.1.7 mugqic/python/2.7.13 && \
mkdir -p alignment/DEX_nonMamm_1/DEX_nonMamm_1_1 metrics/DEX_nonMamm_1/DEX_nonMamm_1_1 && \
java -XX:ParallelGCThreads=1 -Xmx10G -jar $BVATOOLS_JAR \
  bam2fq --mapped ONLY \
  --bam alignment/DEX_nonMamm_1/DEX_nonMamm_1_1/Aligned.sortedByCoord.out.bam    | \
bwa mem  \
  -M -t 10 \
  -R '@RG	ID:DEX_nonMamm_1_1	SM:DEX_nonMamm_1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/rrna_bwa_index/Homo_sapiens.GRCh38.Ensembl87.rrna.fa \
  /dev/stdin | \
 java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx7G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=/dev/stdin \
 OUTPUT=metrics/DEX_nonMamm_1/DEX_nonMamm_1_1/DEX_nonMamm_1_1rRNA.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=1750000 && \
python $PYTHON_TOOLS/rrnaBAMcounter.py \
  -i metrics/DEX_nonMamm_1/DEX_nonMamm_1_1/DEX_nonMamm_1_1rRNA.bam \
  -g $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  -o metrics/DEX_nonMamm_1/DEX_nonMamm_1_1/DEX_nonMamm_1_1rRNA.stats.tsv \
  -t transcript
bwa_mem_rRNA.DEX_nonMamm_1_1.2ac311dc791ce64a9f80306d262820f1.mugqic.done
)
estimate_ribosomal_rna_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$estimate_ribosomal_rna_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: estimate_ribosomal_rna_2_JOB_ID: bwa_mem_rRNA.DEX_nonMamm_1_2
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_rRNA.DEX_nonMamm_1_2
JOB_DEPENDENCIES=$star_23_JOB_ID
JOB_DONE=job_output/estimate_ribosomal_rna/bwa_mem_rRNA.DEX_nonMamm_1_2.a98990d058a052bdff3eb682908116cf.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_rRNA.DEX_nonMamm_1_2.a98990d058a052bdff3eb682908116cf.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/bvatools/1.6 mugqic/bwa/0.7.12 mugqic/picard/2.0.1 mugqic/mugqic_tools/2.1.7 mugqic/python/2.7.13 && \
mkdir -p alignment/DEX_nonMamm_1/DEX_nonMamm_1_2 metrics/DEX_nonMamm_1/DEX_nonMamm_1_2 && \
java -XX:ParallelGCThreads=1 -Xmx10G -jar $BVATOOLS_JAR \
  bam2fq --mapped ONLY \
  --bam alignment/DEX_nonMamm_1/DEX_nonMamm_1_2/Aligned.sortedByCoord.out.bam    | \
bwa mem  \
  -M -t 10 \
  -R '@RG	ID:DEX_nonMamm_1_2	SM:DEX_nonMamm_1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/rrna_bwa_index/Homo_sapiens.GRCh38.Ensembl87.rrna.fa \
  /dev/stdin | \
 java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx7G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=/dev/stdin \
 OUTPUT=metrics/DEX_nonMamm_1/DEX_nonMamm_1_2/DEX_nonMamm_1_2rRNA.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=1750000 && \
python $PYTHON_TOOLS/rrnaBAMcounter.py \
  -i metrics/DEX_nonMamm_1/DEX_nonMamm_1_2/DEX_nonMamm_1_2rRNA.bam \
  -g $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  -o metrics/DEX_nonMamm_1/DEX_nonMamm_1_2/DEX_nonMamm_1_2rRNA.stats.tsv \
  -t transcript
bwa_mem_rRNA.DEX_nonMamm_1_2.a98990d058a052bdff3eb682908116cf.mugqic.done
)
estimate_ribosomal_rna_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$estimate_ribosomal_rna_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: estimate_ribosomal_rna_3_JOB_ID: bwa_mem_rRNA.DEX_shMED1_1_1
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_rRNA.DEX_shMED1_1_1
JOB_DEPENDENCIES=$star_24_JOB_ID
JOB_DONE=job_output/estimate_ribosomal_rna/bwa_mem_rRNA.DEX_shMED1_1_1.d8a3f3bc689268875e60e16c96fbf21f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_rRNA.DEX_shMED1_1_1.d8a3f3bc689268875e60e16c96fbf21f.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/bvatools/1.6 mugqic/bwa/0.7.12 mugqic/picard/2.0.1 mugqic/mugqic_tools/2.1.7 mugqic/python/2.7.13 && \
mkdir -p alignment/DEX_shMED1_1/DEX_shMED1_1_1 metrics/DEX_shMED1_1/DEX_shMED1_1_1 && \
java -XX:ParallelGCThreads=1 -Xmx10G -jar $BVATOOLS_JAR \
  bam2fq --mapped ONLY \
  --bam alignment/DEX_shMED1_1/DEX_shMED1_1_1/Aligned.sortedByCoord.out.bam    | \
bwa mem  \
  -M -t 10 \
  -R '@RG	ID:DEX_shMED1_1_1	SM:DEX_shMED1_1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/rrna_bwa_index/Homo_sapiens.GRCh38.Ensembl87.rrna.fa \
  /dev/stdin | \
 java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx7G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=/dev/stdin \
 OUTPUT=metrics/DEX_shMED1_1/DEX_shMED1_1_1/DEX_shMED1_1_1rRNA.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=1750000 && \
python $PYTHON_TOOLS/rrnaBAMcounter.py \
  -i metrics/DEX_shMED1_1/DEX_shMED1_1_1/DEX_shMED1_1_1rRNA.bam \
  -g $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  -o metrics/DEX_shMED1_1/DEX_shMED1_1_1/DEX_shMED1_1_1rRNA.stats.tsv \
  -t transcript
bwa_mem_rRNA.DEX_shMED1_1_1.d8a3f3bc689268875e60e16c96fbf21f.mugqic.done
)
estimate_ribosomal_rna_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$estimate_ribosomal_rna_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: estimate_ribosomal_rna_4_JOB_ID: bwa_mem_rRNA.DEX_shMED1_1_2
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_rRNA.DEX_shMED1_1_2
JOB_DEPENDENCIES=$star_25_JOB_ID
JOB_DONE=job_output/estimate_ribosomal_rna/bwa_mem_rRNA.DEX_shMED1_1_2.e780fc13b6a00384a56adc70c4bb91f4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_rRNA.DEX_shMED1_1_2.e780fc13b6a00384a56adc70c4bb91f4.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/bvatools/1.6 mugqic/bwa/0.7.12 mugqic/picard/2.0.1 mugqic/mugqic_tools/2.1.7 mugqic/python/2.7.13 && \
mkdir -p alignment/DEX_shMED1_1/DEX_shMED1_1_2 metrics/DEX_shMED1_1/DEX_shMED1_1_2 && \
java -XX:ParallelGCThreads=1 -Xmx10G -jar $BVATOOLS_JAR \
  bam2fq --mapped ONLY \
  --bam alignment/DEX_shMED1_1/DEX_shMED1_1_2/Aligned.sortedByCoord.out.bam    | \
bwa mem  \
  -M -t 10 \
  -R '@RG	ID:DEX_shMED1_1_2	SM:DEX_shMED1_1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/rrna_bwa_index/Homo_sapiens.GRCh38.Ensembl87.rrna.fa \
  /dev/stdin | \
 java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx7G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=/dev/stdin \
 OUTPUT=metrics/DEX_shMED1_1/DEX_shMED1_1_2/DEX_shMED1_1_2rRNA.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=1750000 && \
python $PYTHON_TOOLS/rrnaBAMcounter.py \
  -i metrics/DEX_shMED1_1/DEX_shMED1_1_2/DEX_shMED1_1_2rRNA.bam \
  -g $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  -o metrics/DEX_shMED1_1/DEX_shMED1_1_2/DEX_shMED1_1_2rRNA.stats.tsv \
  -t transcript
bwa_mem_rRNA.DEX_shMED1_1_2.e780fc13b6a00384a56adc70c4bb91f4.mugqic.done
)
estimate_ribosomal_rna_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$estimate_ribosomal_rna_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: estimate_ribosomal_rna_5_JOB_ID: bwa_mem_rRNA.DEX_shNIPBL_1_1
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_rRNA.DEX_shNIPBL_1_1
JOB_DEPENDENCIES=$star_26_JOB_ID
JOB_DONE=job_output/estimate_ribosomal_rna/bwa_mem_rRNA.DEX_shNIPBL_1_1.3c6fcea94292f66b1ad053a2a0d5ee26.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_rRNA.DEX_shNIPBL_1_1.3c6fcea94292f66b1ad053a2a0d5ee26.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/bvatools/1.6 mugqic/bwa/0.7.12 mugqic/picard/2.0.1 mugqic/mugqic_tools/2.1.7 mugqic/python/2.7.13 && \
mkdir -p alignment/DEX_shNIPBL_1/DEX_shNIPBL_1_1 metrics/DEX_shNIPBL_1/DEX_shNIPBL_1_1 && \
java -XX:ParallelGCThreads=1 -Xmx10G -jar $BVATOOLS_JAR \
  bam2fq --mapped ONLY \
  --bam alignment/DEX_shNIPBL_1/DEX_shNIPBL_1_1/Aligned.sortedByCoord.out.bam    | \
bwa mem  \
  -M -t 10 \
  -R '@RG	ID:DEX_shNIPBL_1_1	SM:DEX_shNIPBL_1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/rrna_bwa_index/Homo_sapiens.GRCh38.Ensembl87.rrna.fa \
  /dev/stdin | \
 java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx7G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=/dev/stdin \
 OUTPUT=metrics/DEX_shNIPBL_1/DEX_shNIPBL_1_1/DEX_shNIPBL_1_1rRNA.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=1750000 && \
python $PYTHON_TOOLS/rrnaBAMcounter.py \
  -i metrics/DEX_shNIPBL_1/DEX_shNIPBL_1_1/DEX_shNIPBL_1_1rRNA.bam \
  -g $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  -o metrics/DEX_shNIPBL_1/DEX_shNIPBL_1_1/DEX_shNIPBL_1_1rRNA.stats.tsv \
  -t transcript
bwa_mem_rRNA.DEX_shNIPBL_1_1.3c6fcea94292f66b1ad053a2a0d5ee26.mugqic.done
)
estimate_ribosomal_rna_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$estimate_ribosomal_rna_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: estimate_ribosomal_rna_6_JOB_ID: bwa_mem_rRNA.DEX_shNIPBL_1_2
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_rRNA.DEX_shNIPBL_1_2
JOB_DEPENDENCIES=$star_27_JOB_ID
JOB_DONE=job_output/estimate_ribosomal_rna/bwa_mem_rRNA.DEX_shNIPBL_1_2.aa7b0f7e778e73e7a846be0a60acbf5e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_rRNA.DEX_shNIPBL_1_2.aa7b0f7e778e73e7a846be0a60acbf5e.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/bvatools/1.6 mugqic/bwa/0.7.12 mugqic/picard/2.0.1 mugqic/mugqic_tools/2.1.7 mugqic/python/2.7.13 && \
mkdir -p alignment/DEX_shNIPBL_1/DEX_shNIPBL_1_2 metrics/DEX_shNIPBL_1/DEX_shNIPBL_1_2 && \
java -XX:ParallelGCThreads=1 -Xmx10G -jar $BVATOOLS_JAR \
  bam2fq --mapped ONLY \
  --bam alignment/DEX_shNIPBL_1/DEX_shNIPBL_1_2/Aligned.sortedByCoord.out.bam    | \
bwa mem  \
  -M -t 10 \
  -R '@RG	ID:DEX_shNIPBL_1_2	SM:DEX_shNIPBL_1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/rrna_bwa_index/Homo_sapiens.GRCh38.Ensembl87.rrna.fa \
  /dev/stdin | \
 java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx7G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=/dev/stdin \
 OUTPUT=metrics/DEX_shNIPBL_1/DEX_shNIPBL_1_2/DEX_shNIPBL_1_2rRNA.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=1750000 && \
python $PYTHON_TOOLS/rrnaBAMcounter.py \
  -i metrics/DEX_shNIPBL_1/DEX_shNIPBL_1_2/DEX_shNIPBL_1_2rRNA.bam \
  -g $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  -o metrics/DEX_shNIPBL_1/DEX_shNIPBL_1_2/DEX_shNIPBL_1_2rRNA.stats.tsv \
  -t transcript
bwa_mem_rRNA.DEX_shNIPBL_1_2.aa7b0f7e778e73e7a846be0a60acbf5e.mugqic.done
)
estimate_ribosomal_rna_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$estimate_ribosomal_rna_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: estimate_ribosomal_rna_7_JOB_ID: bwa_mem_rRNA.DEX_shSMC1A_1_1
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_rRNA.DEX_shSMC1A_1_1
JOB_DEPENDENCIES=$star_28_JOB_ID
JOB_DONE=job_output/estimate_ribosomal_rna/bwa_mem_rRNA.DEX_shSMC1A_1_1.b3b66d00f2280c00c4975670336d4b78.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_rRNA.DEX_shSMC1A_1_1.b3b66d00f2280c00c4975670336d4b78.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/bvatools/1.6 mugqic/bwa/0.7.12 mugqic/picard/2.0.1 mugqic/mugqic_tools/2.1.7 mugqic/python/2.7.13 && \
mkdir -p alignment/DEX_shSMC1A_3/DEX_shSMC1A_1_1 metrics/DEX_shSMC1A_3/DEX_shSMC1A_1_1 && \
java -XX:ParallelGCThreads=1 -Xmx10G -jar $BVATOOLS_JAR \
  bam2fq --mapped ONLY \
  --bam alignment/DEX_shSMC1A_3/DEX_shSMC1A_1_1/Aligned.sortedByCoord.out.bam    | \
bwa mem  \
  -M -t 10 \
  -R '@RG	ID:DEX_shSMC1A_1_1	SM:DEX_shSMC1A_3	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/rrna_bwa_index/Homo_sapiens.GRCh38.Ensembl87.rrna.fa \
  /dev/stdin | \
 java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx7G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=/dev/stdin \
 OUTPUT=metrics/DEX_shSMC1A_3/DEX_shSMC1A_1_1/DEX_shSMC1A_1_1rRNA.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=1750000 && \
python $PYTHON_TOOLS/rrnaBAMcounter.py \
  -i metrics/DEX_shSMC1A_3/DEX_shSMC1A_1_1/DEX_shSMC1A_1_1rRNA.bam \
  -g $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  -o metrics/DEX_shSMC1A_3/DEX_shSMC1A_1_1/DEX_shSMC1A_1_1rRNA.stats.tsv \
  -t transcript
bwa_mem_rRNA.DEX_shSMC1A_1_1.b3b66d00f2280c00c4975670336d4b78.mugqic.done
)
estimate_ribosomal_rna_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$estimate_ribosomal_rna_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: estimate_ribosomal_rna_8_JOB_ID: bwa_mem_rRNA.DEX_shSMC1A_1_2
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_rRNA.DEX_shSMC1A_1_2
JOB_DEPENDENCIES=$star_29_JOB_ID
JOB_DONE=job_output/estimate_ribosomal_rna/bwa_mem_rRNA.DEX_shSMC1A_1_2.a3bc4e82477d6b335c132d99ad9e6644.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_rRNA.DEX_shSMC1A_1_2.a3bc4e82477d6b335c132d99ad9e6644.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/bvatools/1.6 mugqic/bwa/0.7.12 mugqic/picard/2.0.1 mugqic/mugqic_tools/2.1.7 mugqic/python/2.7.13 && \
mkdir -p alignment/DEX_shSMC1A_3/DEX_shSMC1A_1_2 metrics/DEX_shSMC1A_3/DEX_shSMC1A_1_2 && \
java -XX:ParallelGCThreads=1 -Xmx10G -jar $BVATOOLS_JAR \
  bam2fq --mapped ONLY \
  --bam alignment/DEX_shSMC1A_3/DEX_shSMC1A_1_2/Aligned.sortedByCoord.out.bam    | \
bwa mem  \
  -M -t 10 \
  -R '@RG	ID:DEX_shSMC1A_1_2	SM:DEX_shSMC1A_3	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/rrna_bwa_index/Homo_sapiens.GRCh38.Ensembl87.rrna.fa \
  /dev/stdin | \
 java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx7G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=/dev/stdin \
 OUTPUT=metrics/DEX_shSMC1A_3/DEX_shSMC1A_1_2/DEX_shSMC1A_1_2rRNA.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=1750000 && \
python $PYTHON_TOOLS/rrnaBAMcounter.py \
  -i metrics/DEX_shSMC1A_3/DEX_shSMC1A_1_2/DEX_shSMC1A_1_2rRNA.bam \
  -g $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  -o metrics/DEX_shSMC1A_3/DEX_shSMC1A_1_2/DEX_shSMC1A_1_2rRNA.stats.tsv \
  -t transcript
bwa_mem_rRNA.DEX_shSMC1A_1_2.a3bc4e82477d6b335c132d99ad9e6644.mugqic.done
)
estimate_ribosomal_rna_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$estimate_ribosomal_rna_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: estimate_ribosomal_rna_9_JOB_ID: bwa_mem_rRNA.ETOH_nonMamm_1_1
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_rRNA.ETOH_nonMamm_1_1
JOB_DEPENDENCIES=$star_30_JOB_ID
JOB_DONE=job_output/estimate_ribosomal_rna/bwa_mem_rRNA.ETOH_nonMamm_1_1.742624ad6cbfbd7cfcb99121303a3bec.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_rRNA.ETOH_nonMamm_1_1.742624ad6cbfbd7cfcb99121303a3bec.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/bvatools/1.6 mugqic/bwa/0.7.12 mugqic/picard/2.0.1 mugqic/mugqic_tools/2.1.7 mugqic/python/2.7.13 && \
mkdir -p alignment/ETOH_nonMamm_1/ETOH_nonMamm_1_1 metrics/ETOH_nonMamm_1/ETOH_nonMamm_1_1 && \
java -XX:ParallelGCThreads=1 -Xmx10G -jar $BVATOOLS_JAR \
  bam2fq --mapped ONLY \
  --bam alignment/ETOH_nonMamm_1/ETOH_nonMamm_1_1/Aligned.sortedByCoord.out.bam    | \
bwa mem  \
  -M -t 10 \
  -R '@RG	ID:ETOH_nonMamm_1_1	SM:ETOH_nonMamm_1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/rrna_bwa_index/Homo_sapiens.GRCh38.Ensembl87.rrna.fa \
  /dev/stdin | \
 java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx7G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=/dev/stdin \
 OUTPUT=metrics/ETOH_nonMamm_1/ETOH_nonMamm_1_1/ETOH_nonMamm_1_1rRNA.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=1750000 && \
python $PYTHON_TOOLS/rrnaBAMcounter.py \
  -i metrics/ETOH_nonMamm_1/ETOH_nonMamm_1_1/ETOH_nonMamm_1_1rRNA.bam \
  -g $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  -o metrics/ETOH_nonMamm_1/ETOH_nonMamm_1_1/ETOH_nonMamm_1_1rRNA.stats.tsv \
  -t transcript
bwa_mem_rRNA.ETOH_nonMamm_1_1.742624ad6cbfbd7cfcb99121303a3bec.mugqic.done
)
estimate_ribosomal_rna_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$estimate_ribosomal_rna_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: estimate_ribosomal_rna_10_JOB_ID: bwa_mem_rRNA.ETOH_nonMamm_1_2
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_rRNA.ETOH_nonMamm_1_2
JOB_DEPENDENCIES=$star_31_JOB_ID
JOB_DONE=job_output/estimate_ribosomal_rna/bwa_mem_rRNA.ETOH_nonMamm_1_2.723e6d271bc6be4bf4073feb12a9488a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_rRNA.ETOH_nonMamm_1_2.723e6d271bc6be4bf4073feb12a9488a.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/bvatools/1.6 mugqic/bwa/0.7.12 mugqic/picard/2.0.1 mugqic/mugqic_tools/2.1.7 mugqic/python/2.7.13 && \
mkdir -p alignment/ETOH_nonMamm_1/ETOH_nonMamm_1_2 metrics/ETOH_nonMamm_1/ETOH_nonMamm_1_2 && \
java -XX:ParallelGCThreads=1 -Xmx10G -jar $BVATOOLS_JAR \
  bam2fq --mapped ONLY \
  --bam alignment/ETOH_nonMamm_1/ETOH_nonMamm_1_2/Aligned.sortedByCoord.out.bam    | \
bwa mem  \
  -M -t 10 \
  -R '@RG	ID:ETOH_nonMamm_1_2	SM:ETOH_nonMamm_1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/rrna_bwa_index/Homo_sapiens.GRCh38.Ensembl87.rrna.fa \
  /dev/stdin | \
 java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx7G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=/dev/stdin \
 OUTPUT=metrics/ETOH_nonMamm_1/ETOH_nonMamm_1_2/ETOH_nonMamm_1_2rRNA.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=1750000 && \
python $PYTHON_TOOLS/rrnaBAMcounter.py \
  -i metrics/ETOH_nonMamm_1/ETOH_nonMamm_1_2/ETOH_nonMamm_1_2rRNA.bam \
  -g $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  -o metrics/ETOH_nonMamm_1/ETOH_nonMamm_1_2/ETOH_nonMamm_1_2rRNA.stats.tsv \
  -t transcript
bwa_mem_rRNA.ETOH_nonMamm_1_2.723e6d271bc6be4bf4073feb12a9488a.mugqic.done
)
estimate_ribosomal_rna_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$estimate_ribosomal_rna_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: estimate_ribosomal_rna_11_JOB_ID: bwa_mem_rRNA.ETOH_nonMamm_2_1
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_rRNA.ETOH_nonMamm_2_1
JOB_DEPENDENCIES=$star_32_JOB_ID
JOB_DONE=job_output/estimate_ribosomal_rna/bwa_mem_rRNA.ETOH_nonMamm_2_1.b853e457aab20e5309642a321c1b500a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_rRNA.ETOH_nonMamm_2_1.b853e457aab20e5309642a321c1b500a.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/bvatools/1.6 mugqic/bwa/0.7.12 mugqic/picard/2.0.1 mugqic/mugqic_tools/2.1.7 mugqic/python/2.7.13 && \
mkdir -p alignment/ETOH_nonMamm_2/ETOH_nonMamm_2_1 metrics/ETOH_nonMamm_2/ETOH_nonMamm_2_1 && \
java -XX:ParallelGCThreads=1 -Xmx10G -jar $BVATOOLS_JAR \
  bam2fq --mapped ONLY \
  --bam alignment/ETOH_nonMamm_2/ETOH_nonMamm_2_1/Aligned.sortedByCoord.out.bam    | \
bwa mem  \
  -M -t 10 \
  -R '@RG	ID:ETOH_nonMamm_2_1	SM:ETOH_nonMamm_2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/rrna_bwa_index/Homo_sapiens.GRCh38.Ensembl87.rrna.fa \
  /dev/stdin | \
 java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx7G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=/dev/stdin \
 OUTPUT=metrics/ETOH_nonMamm_2/ETOH_nonMamm_2_1/ETOH_nonMamm_2_1rRNA.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=1750000 && \
python $PYTHON_TOOLS/rrnaBAMcounter.py \
  -i metrics/ETOH_nonMamm_2/ETOH_nonMamm_2_1/ETOH_nonMamm_2_1rRNA.bam \
  -g $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  -o metrics/ETOH_nonMamm_2/ETOH_nonMamm_2_1/ETOH_nonMamm_2_1rRNA.stats.tsv \
  -t transcript
bwa_mem_rRNA.ETOH_nonMamm_2_1.b853e457aab20e5309642a321c1b500a.mugqic.done
)
estimate_ribosomal_rna_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$estimate_ribosomal_rna_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: estimate_ribosomal_rna_12_JOB_ID: bwa_mem_rRNA.ETOH_shMED1_1_1
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_rRNA.ETOH_shMED1_1_1
JOB_DEPENDENCIES=$star_33_JOB_ID
JOB_DONE=job_output/estimate_ribosomal_rna/bwa_mem_rRNA.ETOH_shMED1_1_1.f66db57c5c608290acc5df35c5cb462c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_rRNA.ETOH_shMED1_1_1.f66db57c5c608290acc5df35c5cb462c.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/bvatools/1.6 mugqic/bwa/0.7.12 mugqic/picard/2.0.1 mugqic/mugqic_tools/2.1.7 mugqic/python/2.7.13 && \
mkdir -p alignment/ETOH_shMED1_1/ETOH_shMED1_1_1 metrics/ETOH_shMED1_1/ETOH_shMED1_1_1 && \
java -XX:ParallelGCThreads=1 -Xmx10G -jar $BVATOOLS_JAR \
  bam2fq --mapped ONLY \
  --bam alignment/ETOH_shMED1_1/ETOH_shMED1_1_1/Aligned.sortedByCoord.out.bam    | \
bwa mem  \
  -M -t 10 \
  -R '@RG	ID:ETOH_shMED1_1_1	SM:ETOH_shMED1_1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/rrna_bwa_index/Homo_sapiens.GRCh38.Ensembl87.rrna.fa \
  /dev/stdin | \
 java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx7G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=/dev/stdin \
 OUTPUT=metrics/ETOH_shMED1_1/ETOH_shMED1_1_1/ETOH_shMED1_1_1rRNA.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=1750000 && \
python $PYTHON_TOOLS/rrnaBAMcounter.py \
  -i metrics/ETOH_shMED1_1/ETOH_shMED1_1_1/ETOH_shMED1_1_1rRNA.bam \
  -g $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  -o metrics/ETOH_shMED1_1/ETOH_shMED1_1_1/ETOH_shMED1_1_1rRNA.stats.tsv \
  -t transcript
bwa_mem_rRNA.ETOH_shMED1_1_1.f66db57c5c608290acc5df35c5cb462c.mugqic.done
)
estimate_ribosomal_rna_12_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$estimate_ribosomal_rna_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: estimate_ribosomal_rna_13_JOB_ID: bwa_mem_rRNA.ETOH_shMED1_1_2
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_rRNA.ETOH_shMED1_1_2
JOB_DEPENDENCIES=$star_34_JOB_ID
JOB_DONE=job_output/estimate_ribosomal_rna/bwa_mem_rRNA.ETOH_shMED1_1_2.97e82006d0194b8f9ec2a053f74ec8e4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_rRNA.ETOH_shMED1_1_2.97e82006d0194b8f9ec2a053f74ec8e4.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/bvatools/1.6 mugqic/bwa/0.7.12 mugqic/picard/2.0.1 mugqic/mugqic_tools/2.1.7 mugqic/python/2.7.13 && \
mkdir -p alignment/ETOH_shMED1_1/ETOH_shMED1_1_2 metrics/ETOH_shMED1_1/ETOH_shMED1_1_2 && \
java -XX:ParallelGCThreads=1 -Xmx10G -jar $BVATOOLS_JAR \
  bam2fq --mapped ONLY \
  --bam alignment/ETOH_shMED1_1/ETOH_shMED1_1_2/Aligned.sortedByCoord.out.bam    | \
bwa mem  \
  -M -t 10 \
  -R '@RG	ID:ETOH_shMED1_1_2	SM:ETOH_shMED1_1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/rrna_bwa_index/Homo_sapiens.GRCh38.Ensembl87.rrna.fa \
  /dev/stdin | \
 java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx7G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=/dev/stdin \
 OUTPUT=metrics/ETOH_shMED1_1/ETOH_shMED1_1_2/ETOH_shMED1_1_2rRNA.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=1750000 && \
python $PYTHON_TOOLS/rrnaBAMcounter.py \
  -i metrics/ETOH_shMED1_1/ETOH_shMED1_1_2/ETOH_shMED1_1_2rRNA.bam \
  -g $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  -o metrics/ETOH_shMED1_1/ETOH_shMED1_1_2/ETOH_shMED1_1_2rRNA.stats.tsv \
  -t transcript
bwa_mem_rRNA.ETOH_shMED1_1_2.97e82006d0194b8f9ec2a053f74ec8e4.mugqic.done
)
estimate_ribosomal_rna_13_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$estimate_ribosomal_rna_13_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: estimate_ribosomal_rna_14_JOB_ID: bwa_mem_rRNA.ETOH_shNIPBL_2_1
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_rRNA.ETOH_shNIPBL_2_1
JOB_DEPENDENCIES=$star_35_JOB_ID
JOB_DONE=job_output/estimate_ribosomal_rna/bwa_mem_rRNA.ETOH_shNIPBL_2_1.6254e238c641643ab8dced232fdd9c06.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_rRNA.ETOH_shNIPBL_2_1.6254e238c641643ab8dced232fdd9c06.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/bvatools/1.6 mugqic/bwa/0.7.12 mugqic/picard/2.0.1 mugqic/mugqic_tools/2.1.7 mugqic/python/2.7.13 && \
mkdir -p alignment/ETOH_shNIPBL_2/ETOH_shNIPBL_2_1 metrics/ETOH_shNIPBL_2/ETOH_shNIPBL_2_1 && \
java -XX:ParallelGCThreads=1 -Xmx10G -jar $BVATOOLS_JAR \
  bam2fq --mapped ONLY \
  --bam alignment/ETOH_shNIPBL_2/ETOH_shNIPBL_2_1/Aligned.sortedByCoord.out.bam    | \
bwa mem  \
  -M -t 10 \
  -R '@RG	ID:ETOH_shNIPBL_2_1	SM:ETOH_shNIPBL_2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/rrna_bwa_index/Homo_sapiens.GRCh38.Ensembl87.rrna.fa \
  /dev/stdin | \
 java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx7G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=/dev/stdin \
 OUTPUT=metrics/ETOH_shNIPBL_2/ETOH_shNIPBL_2_1/ETOH_shNIPBL_2_1rRNA.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=1750000 && \
python $PYTHON_TOOLS/rrnaBAMcounter.py \
  -i metrics/ETOH_shNIPBL_2/ETOH_shNIPBL_2_1/ETOH_shNIPBL_2_1rRNA.bam \
  -g $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  -o metrics/ETOH_shNIPBL_2/ETOH_shNIPBL_2_1/ETOH_shNIPBL_2_1rRNA.stats.tsv \
  -t transcript
bwa_mem_rRNA.ETOH_shNIPBL_2_1.6254e238c641643ab8dced232fdd9c06.mugqic.done
)
estimate_ribosomal_rna_14_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$estimate_ribosomal_rna_14_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: estimate_ribosomal_rna_15_JOB_ID: bwa_mem_rRNA.ETOH_shNIPBL_1_1
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_rRNA.ETOH_shNIPBL_1_1
JOB_DEPENDENCIES=$star_36_JOB_ID
JOB_DONE=job_output/estimate_ribosomal_rna/bwa_mem_rRNA.ETOH_shNIPBL_1_1.8cf1d0f75d5a531391cad9ee3782f2c8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_rRNA.ETOH_shNIPBL_1_1.8cf1d0f75d5a531391cad9ee3782f2c8.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/bvatools/1.6 mugqic/bwa/0.7.12 mugqic/picard/2.0.1 mugqic/mugqic_tools/2.1.7 mugqic/python/2.7.13 && \
mkdir -p alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1_1 metrics/ETOH_shNIPBL_1/ETOH_shNIPBL_1_1 && \
java -XX:ParallelGCThreads=1 -Xmx10G -jar $BVATOOLS_JAR \
  bam2fq --mapped ONLY \
  --bam alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1_1/Aligned.sortedByCoord.out.bam    | \
bwa mem  \
  -M -t 10 \
  -R '@RG	ID:ETOH_shNIPBL_1_1	SM:ETOH_shNIPBL_1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/rrna_bwa_index/Homo_sapiens.GRCh38.Ensembl87.rrna.fa \
  /dev/stdin | \
 java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx7G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=/dev/stdin \
 OUTPUT=metrics/ETOH_shNIPBL_1/ETOH_shNIPBL_1_1/ETOH_shNIPBL_1_1rRNA.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=1750000 && \
python $PYTHON_TOOLS/rrnaBAMcounter.py \
  -i metrics/ETOH_shNIPBL_1/ETOH_shNIPBL_1_1/ETOH_shNIPBL_1_1rRNA.bam \
  -g $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  -o metrics/ETOH_shNIPBL_1/ETOH_shNIPBL_1_1/ETOH_shNIPBL_1_1rRNA.stats.tsv \
  -t transcript
bwa_mem_rRNA.ETOH_shNIPBL_1_1.8cf1d0f75d5a531391cad9ee3782f2c8.mugqic.done
)
estimate_ribosomal_rna_15_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$estimate_ribosomal_rna_15_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: estimate_ribosomal_rna_16_JOB_ID: bwa_mem_rRNA.ETOH_shNIPBL_1_2
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_rRNA.ETOH_shNIPBL_1_2
JOB_DEPENDENCIES=$star_37_JOB_ID
JOB_DONE=job_output/estimate_ribosomal_rna/bwa_mem_rRNA.ETOH_shNIPBL_1_2.44e5b80c4c40fb1cee899221602e72a6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_rRNA.ETOH_shNIPBL_1_2.44e5b80c4c40fb1cee899221602e72a6.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/bvatools/1.6 mugqic/bwa/0.7.12 mugqic/picard/2.0.1 mugqic/mugqic_tools/2.1.7 mugqic/python/2.7.13 && \
mkdir -p alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1_2 metrics/ETOH_shNIPBL_1/ETOH_shNIPBL_1_2 && \
java -XX:ParallelGCThreads=1 -Xmx10G -jar $BVATOOLS_JAR \
  bam2fq --mapped ONLY \
  --bam alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1_2/Aligned.sortedByCoord.out.bam    | \
bwa mem  \
  -M -t 10 \
  -R '@RG	ID:ETOH_shNIPBL_1_2	SM:ETOH_shNIPBL_1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/rrna_bwa_index/Homo_sapiens.GRCh38.Ensembl87.rrna.fa \
  /dev/stdin | \
 java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx7G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=/dev/stdin \
 OUTPUT=metrics/ETOH_shNIPBL_1/ETOH_shNIPBL_1_2/ETOH_shNIPBL_1_2rRNA.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=1750000 && \
python $PYTHON_TOOLS/rrnaBAMcounter.py \
  -i metrics/ETOH_shNIPBL_1/ETOH_shNIPBL_1_2/ETOH_shNIPBL_1_2rRNA.bam \
  -g $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  -o metrics/ETOH_shNIPBL_1/ETOH_shNIPBL_1_2/ETOH_shNIPBL_1_2rRNA.stats.tsv \
  -t transcript
bwa_mem_rRNA.ETOH_shNIPBL_1_2.44e5b80c4c40fb1cee899221602e72a6.mugqic.done
)
estimate_ribosomal_rna_16_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$estimate_ribosomal_rna_16_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: estimate_ribosomal_rna_17_JOB_ID: bwa_mem_rRNA.ETOH_shMED1_2_1
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_rRNA.ETOH_shMED1_2_1
JOB_DEPENDENCIES=$star_38_JOB_ID
JOB_DONE=job_output/estimate_ribosomal_rna/bwa_mem_rRNA.ETOH_shMED1_2_1.1829b5aeafe7fde6db8842b803354b72.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_rRNA.ETOH_shMED1_2_1.1829b5aeafe7fde6db8842b803354b72.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/bvatools/1.6 mugqic/bwa/0.7.12 mugqic/picard/2.0.1 mugqic/mugqic_tools/2.1.7 mugqic/python/2.7.13 && \
mkdir -p alignment/ETOH_shMED1_2/ETOH_shMED1_2_1 metrics/ETOH_shMED1_2/ETOH_shMED1_2_1 && \
java -XX:ParallelGCThreads=1 -Xmx10G -jar $BVATOOLS_JAR \
  bam2fq --mapped ONLY \
  --bam alignment/ETOH_shMED1_2/ETOH_shMED1_2_1/Aligned.sortedByCoord.out.bam    | \
bwa mem  \
  -M -t 10 \
  -R '@RG	ID:ETOH_shMED1_2_1	SM:ETOH_shMED1_2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/rrna_bwa_index/Homo_sapiens.GRCh38.Ensembl87.rrna.fa \
  /dev/stdin | \
 java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx7G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=/dev/stdin \
 OUTPUT=metrics/ETOH_shMED1_2/ETOH_shMED1_2_1/ETOH_shMED1_2_1rRNA.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=1750000 && \
python $PYTHON_TOOLS/rrnaBAMcounter.py \
  -i metrics/ETOH_shMED1_2/ETOH_shMED1_2_1/ETOH_shMED1_2_1rRNA.bam \
  -g $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  -o metrics/ETOH_shMED1_2/ETOH_shMED1_2_1/ETOH_shMED1_2_1rRNA.stats.tsv \
  -t transcript
bwa_mem_rRNA.ETOH_shMED1_2_1.1829b5aeafe7fde6db8842b803354b72.mugqic.done
)
estimate_ribosomal_rna_17_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$estimate_ribosomal_rna_17_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: estimate_ribosomal_rna_18_JOB_ID: bwa_mem_rRNA.ETOH_shSMC1A_1_1
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_rRNA.ETOH_shSMC1A_1_1
JOB_DEPENDENCIES=$star_39_JOB_ID
JOB_DONE=job_output/estimate_ribosomal_rna/bwa_mem_rRNA.ETOH_shSMC1A_1_1.c1f5a2661cfe215e6b33b491bb6c7cad.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_rRNA.ETOH_shSMC1A_1_1.c1f5a2661cfe215e6b33b491bb6c7cad.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/bvatools/1.6 mugqic/bwa/0.7.12 mugqic/picard/2.0.1 mugqic/mugqic_tools/2.1.7 mugqic/python/2.7.13 && \
mkdir -p alignment/ETOH_shSMC1A_3/ETOH_shSMC1A_1_1 metrics/ETOH_shSMC1A_3/ETOH_shSMC1A_1_1 && \
java -XX:ParallelGCThreads=1 -Xmx10G -jar $BVATOOLS_JAR \
  bam2fq --mapped ONLY \
  --bam alignment/ETOH_shSMC1A_3/ETOH_shSMC1A_1_1/Aligned.sortedByCoord.out.bam    | \
bwa mem  \
  -M -t 10 \
  -R '@RG	ID:ETOH_shSMC1A_1_1	SM:ETOH_shSMC1A_3	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/rrna_bwa_index/Homo_sapiens.GRCh38.Ensembl87.rrna.fa \
  /dev/stdin | \
 java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx7G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=/dev/stdin \
 OUTPUT=metrics/ETOH_shSMC1A_3/ETOH_shSMC1A_1_1/ETOH_shSMC1A_1_1rRNA.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=1750000 && \
python $PYTHON_TOOLS/rrnaBAMcounter.py \
  -i metrics/ETOH_shSMC1A_3/ETOH_shSMC1A_1_1/ETOH_shSMC1A_1_1rRNA.bam \
  -g $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  -o metrics/ETOH_shSMC1A_3/ETOH_shSMC1A_1_1/ETOH_shSMC1A_1_1rRNA.stats.tsv \
  -t transcript
bwa_mem_rRNA.ETOH_shSMC1A_1_1.c1f5a2661cfe215e6b33b491bb6c7cad.mugqic.done
)
estimate_ribosomal_rna_18_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$estimate_ribosomal_rna_18_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: estimate_ribosomal_rna_19_JOB_ID: bwa_mem_rRNA.ETOH_shSMC1A_1_2
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_rRNA.ETOH_shSMC1A_1_2
JOB_DEPENDENCIES=$star_40_JOB_ID
JOB_DONE=job_output/estimate_ribosomal_rna/bwa_mem_rRNA.ETOH_shSMC1A_1_2.ca6318b38a76493d0019920ba27bdb91.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_rRNA.ETOH_shSMC1A_1_2.ca6318b38a76493d0019920ba27bdb91.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/bvatools/1.6 mugqic/bwa/0.7.12 mugqic/picard/2.0.1 mugqic/mugqic_tools/2.1.7 mugqic/python/2.7.13 && \
mkdir -p alignment/ETOH_shSMC1A_3/ETOH_shSMC1A_1_2 metrics/ETOH_shSMC1A_3/ETOH_shSMC1A_1_2 && \
java -XX:ParallelGCThreads=1 -Xmx10G -jar $BVATOOLS_JAR \
  bam2fq --mapped ONLY \
  --bam alignment/ETOH_shSMC1A_3/ETOH_shSMC1A_1_2/Aligned.sortedByCoord.out.bam    | \
bwa mem  \
  -M -t 10 \
  -R '@RG	ID:ETOH_shSMC1A_1_2	SM:ETOH_shSMC1A_3	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/rrna_bwa_index/Homo_sapiens.GRCh38.Ensembl87.rrna.fa \
  /dev/stdin | \
 java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx7G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=/dev/stdin \
 OUTPUT=metrics/ETOH_shSMC1A_3/ETOH_shSMC1A_1_2/ETOH_shSMC1A_1_2rRNA.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=1750000 && \
python $PYTHON_TOOLS/rrnaBAMcounter.py \
  -i metrics/ETOH_shSMC1A_3/ETOH_shSMC1A_1_2/ETOH_shSMC1A_1_2rRNA.bam \
  -g $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  -o metrics/ETOH_shSMC1A_3/ETOH_shSMC1A_1_2/ETOH_shSMC1A_1_2rRNA.stats.tsv \
  -t transcript
bwa_mem_rRNA.ETOH_shSMC1A_1_2.ca6318b38a76493d0019920ba27bdb91.mugqic.done
)
estimate_ribosomal_rna_19_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$estimate_ribosomal_rna_19_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: estimate_ribosomal_rna_20_JOB_ID: bwa_mem_rRNA.ETOH_shSMC1A_2_1
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_rRNA.ETOH_shSMC1A_2_1
JOB_DEPENDENCIES=$star_41_JOB_ID
JOB_DONE=job_output/estimate_ribosomal_rna/bwa_mem_rRNA.ETOH_shSMC1A_2_1.54cbff1a30cd9462602a36fc929ab07f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_rRNA.ETOH_shSMC1A_2_1.54cbff1a30cd9462602a36fc929ab07f.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/bvatools/1.6 mugqic/bwa/0.7.12 mugqic/picard/2.0.1 mugqic/mugqic_tools/2.1.7 mugqic/python/2.7.13 && \
mkdir -p alignment/ETOH_shSMC1A_2/ETOH_shSMC1A_2_1 metrics/ETOH_shSMC1A_2/ETOH_shSMC1A_2_1 && \
java -XX:ParallelGCThreads=1 -Xmx10G -jar $BVATOOLS_JAR \
  bam2fq --mapped ONLY \
  --bam alignment/ETOH_shSMC1A_2/ETOH_shSMC1A_2_1/Aligned.sortedByCoord.out.bam    | \
bwa mem  \
  -M -t 10 \
  -R '@RG	ID:ETOH_shSMC1A_2_1	SM:ETOH_shSMC1A_2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/rrna_bwa_index/Homo_sapiens.GRCh38.Ensembl87.rrna.fa \
  /dev/stdin | \
 java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx7G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=/dev/stdin \
 OUTPUT=metrics/ETOH_shSMC1A_2/ETOH_shSMC1A_2_1/ETOH_shSMC1A_2_1rRNA.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=1750000 && \
python $PYTHON_TOOLS/rrnaBAMcounter.py \
  -i metrics/ETOH_shSMC1A_2/ETOH_shSMC1A_2_1/ETOH_shSMC1A_2_1rRNA.bam \
  -g $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  -o metrics/ETOH_shSMC1A_2/ETOH_shSMC1A_2_1/ETOH_shSMC1A_2_1rRNA.stats.tsv \
  -t transcript
bwa_mem_rRNA.ETOH_shSMC1A_2_1.54cbff1a30cd9462602a36fc929ab07f.mugqic.done
)
estimate_ribosomal_rna_20_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$estimate_ribosomal_rna_20_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: bam_hard_clip
#-------------------------------------------------------------------------------
STEP=bam_hard_clip
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: bam_hard_clip_1_JOB_ID: tuxedo_hard_clip.DEX_nonMamm_1
#-------------------------------------------------------------------------------
JOB_NAME=tuxedo_hard_clip.DEX_nonMamm_1
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID
JOB_DONE=job_output/bam_hard_clip/tuxedo_hard_clip.DEX_nonMamm_1.a4a742da1196170007d2c1912ccf229e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'tuxedo_hard_clip.DEX_nonMamm_1.a4a742da1196170007d2c1912ccf229e.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -h \
  alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.bam | \
awk 'BEGIN {OFS="\t"} {if (substr($1,1,1)=="@") {print;next}; split($6,C,/[0-9]*/); split($6,L,/[SMDIN]/); if (C[2]=="S") {$10=substr($10,L[1]+1); $11=substr($11,L[1]+1)}; if (C[length(C)]=="S") {L1=length($10)-L[length(L)-1]; $10=substr($10,1,L1); $11=substr($11,1,L1); }; gsub(/[0-9]*S/,"",$6); print}'  | \
samtools view -hbS \
  - \
  > alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.hardClip.bam
tuxedo_hard_clip.DEX_nonMamm_1.a4a742da1196170007d2c1912ccf229e.mugqic.done
)
bam_hard_clip_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bam_hard_clip_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bam_hard_clip_2_JOB_ID: tuxedo_hard_clip.DEX_shMED1_1
#-------------------------------------------------------------------------------
JOB_NAME=tuxedo_hard_clip.DEX_shMED1_1
JOB_DEPENDENCIES=$picard_mark_duplicates_2_JOB_ID
JOB_DONE=job_output/bam_hard_clip/tuxedo_hard_clip.DEX_shMED1_1.3932f8cabb2245e8d875d3ad4490cc35.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'tuxedo_hard_clip.DEX_shMED1_1.3932f8cabb2245e8d875d3ad4490cc35.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -h \
  alignment/DEX_shMED1_1/DEX_shMED1_1.sorted.mdup.bam | \
awk 'BEGIN {OFS="\t"} {if (substr($1,1,1)=="@") {print;next}; split($6,C,/[0-9]*/); split($6,L,/[SMDIN]/); if (C[2]=="S") {$10=substr($10,L[1]+1); $11=substr($11,L[1]+1)}; if (C[length(C)]=="S") {L1=length($10)-L[length(L)-1]; $10=substr($10,1,L1); $11=substr($11,1,L1); }; gsub(/[0-9]*S/,"",$6); print}'  | \
samtools view -hbS \
  - \
  > alignment/DEX_shMED1_1/DEX_shMED1_1.sorted.mdup.hardClip.bam
tuxedo_hard_clip.DEX_shMED1_1.3932f8cabb2245e8d875d3ad4490cc35.mugqic.done
)
bam_hard_clip_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bam_hard_clip_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bam_hard_clip_3_JOB_ID: tuxedo_hard_clip.DEX_shNIPBL_1
#-------------------------------------------------------------------------------
JOB_NAME=tuxedo_hard_clip.DEX_shNIPBL_1
JOB_DEPENDENCIES=$picard_mark_duplicates_3_JOB_ID
JOB_DONE=job_output/bam_hard_clip/tuxedo_hard_clip.DEX_shNIPBL_1.12f152a379c2418291ab64861f8fe5be.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'tuxedo_hard_clip.DEX_shNIPBL_1.12f152a379c2418291ab64861f8fe5be.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -h \
  alignment/DEX_shNIPBL_1/DEX_shNIPBL_1.sorted.mdup.bam | \
awk 'BEGIN {OFS="\t"} {if (substr($1,1,1)=="@") {print;next}; split($6,C,/[0-9]*/); split($6,L,/[SMDIN]/); if (C[2]=="S") {$10=substr($10,L[1]+1); $11=substr($11,L[1]+1)}; if (C[length(C)]=="S") {L1=length($10)-L[length(L)-1]; $10=substr($10,1,L1); $11=substr($11,1,L1); }; gsub(/[0-9]*S/,"",$6); print}'  | \
samtools view -hbS \
  - \
  > alignment/DEX_shNIPBL_1/DEX_shNIPBL_1.sorted.mdup.hardClip.bam
tuxedo_hard_clip.DEX_shNIPBL_1.12f152a379c2418291ab64861f8fe5be.mugqic.done
)
bam_hard_clip_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bam_hard_clip_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bam_hard_clip_4_JOB_ID: tuxedo_hard_clip.DEX_shSMC1A_3
#-------------------------------------------------------------------------------
JOB_NAME=tuxedo_hard_clip.DEX_shSMC1A_3
JOB_DEPENDENCIES=$picard_mark_duplicates_4_JOB_ID
JOB_DONE=job_output/bam_hard_clip/tuxedo_hard_clip.DEX_shSMC1A_3.7c941d5d3b726a681f3fd77c05be9d93.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'tuxedo_hard_clip.DEX_shSMC1A_3.7c941d5d3b726a681f3fd77c05be9d93.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -h \
  alignment/DEX_shSMC1A_3/DEX_shSMC1A_3.sorted.mdup.bam | \
awk 'BEGIN {OFS="\t"} {if (substr($1,1,1)=="@") {print;next}; split($6,C,/[0-9]*/); split($6,L,/[SMDIN]/); if (C[2]=="S") {$10=substr($10,L[1]+1); $11=substr($11,L[1]+1)}; if (C[length(C)]=="S") {L1=length($10)-L[length(L)-1]; $10=substr($10,1,L1); $11=substr($11,1,L1); }; gsub(/[0-9]*S/,"",$6); print}'  | \
samtools view -hbS \
  - \
  > alignment/DEX_shSMC1A_3/DEX_shSMC1A_3.sorted.mdup.hardClip.bam
tuxedo_hard_clip.DEX_shSMC1A_3.7c941d5d3b726a681f3fd77c05be9d93.mugqic.done
)
bam_hard_clip_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bam_hard_clip_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bam_hard_clip_5_JOB_ID: tuxedo_hard_clip.ETOH_nonMamm_1
#-------------------------------------------------------------------------------
JOB_NAME=tuxedo_hard_clip.ETOH_nonMamm_1
JOB_DEPENDENCIES=$picard_mark_duplicates_5_JOB_ID
JOB_DONE=job_output/bam_hard_clip/tuxedo_hard_clip.ETOH_nonMamm_1.9652261dd7328b5e0fd9faf0a777a5d0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'tuxedo_hard_clip.ETOH_nonMamm_1.9652261dd7328b5e0fd9faf0a777a5d0.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -h \
  alignment/ETOH_nonMamm_1/ETOH_nonMamm_1.sorted.mdup.bam | \
awk 'BEGIN {OFS="\t"} {if (substr($1,1,1)=="@") {print;next}; split($6,C,/[0-9]*/); split($6,L,/[SMDIN]/); if (C[2]=="S") {$10=substr($10,L[1]+1); $11=substr($11,L[1]+1)}; if (C[length(C)]=="S") {L1=length($10)-L[length(L)-1]; $10=substr($10,1,L1); $11=substr($11,1,L1); }; gsub(/[0-9]*S/,"",$6); print}'  | \
samtools view -hbS \
  - \
  > alignment/ETOH_nonMamm_1/ETOH_nonMamm_1.sorted.mdup.hardClip.bam
tuxedo_hard_clip.ETOH_nonMamm_1.9652261dd7328b5e0fd9faf0a777a5d0.mugqic.done
)
bam_hard_clip_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bam_hard_clip_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bam_hard_clip_6_JOB_ID: tuxedo_hard_clip.ETOH_nonMamm_2
#-------------------------------------------------------------------------------
JOB_NAME=tuxedo_hard_clip.ETOH_nonMamm_2
JOB_DEPENDENCIES=$picard_mark_duplicates_6_JOB_ID
JOB_DONE=job_output/bam_hard_clip/tuxedo_hard_clip.ETOH_nonMamm_2.6327dedc2fe339cc44f545c449ebde2f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'tuxedo_hard_clip.ETOH_nonMamm_2.6327dedc2fe339cc44f545c449ebde2f.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -h \
  alignment/ETOH_nonMamm_2/ETOH_nonMamm_2.sorted.mdup.bam | \
awk 'BEGIN {OFS="\t"} {if (substr($1,1,1)=="@") {print;next}; split($6,C,/[0-9]*/); split($6,L,/[SMDIN]/); if (C[2]=="S") {$10=substr($10,L[1]+1); $11=substr($11,L[1]+1)}; if (C[length(C)]=="S") {L1=length($10)-L[length(L)-1]; $10=substr($10,1,L1); $11=substr($11,1,L1); }; gsub(/[0-9]*S/,"",$6); print}'  | \
samtools view -hbS \
  - \
  > alignment/ETOH_nonMamm_2/ETOH_nonMamm_2.sorted.mdup.hardClip.bam
tuxedo_hard_clip.ETOH_nonMamm_2.6327dedc2fe339cc44f545c449ebde2f.mugqic.done
)
bam_hard_clip_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bam_hard_clip_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bam_hard_clip_7_JOB_ID: tuxedo_hard_clip.ETOH_shMED1_1
#-------------------------------------------------------------------------------
JOB_NAME=tuxedo_hard_clip.ETOH_shMED1_1
JOB_DEPENDENCIES=$picard_mark_duplicates_7_JOB_ID
JOB_DONE=job_output/bam_hard_clip/tuxedo_hard_clip.ETOH_shMED1_1.67840f50fcd002400d1307172012d86f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'tuxedo_hard_clip.ETOH_shMED1_1.67840f50fcd002400d1307172012d86f.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -h \
  alignment/ETOH_shMED1_1/ETOH_shMED1_1.sorted.mdup.bam | \
awk 'BEGIN {OFS="\t"} {if (substr($1,1,1)=="@") {print;next}; split($6,C,/[0-9]*/); split($6,L,/[SMDIN]/); if (C[2]=="S") {$10=substr($10,L[1]+1); $11=substr($11,L[1]+1)}; if (C[length(C)]=="S") {L1=length($10)-L[length(L)-1]; $10=substr($10,1,L1); $11=substr($11,1,L1); }; gsub(/[0-9]*S/,"",$6); print}'  | \
samtools view -hbS \
  - \
  > alignment/ETOH_shMED1_1/ETOH_shMED1_1.sorted.mdup.hardClip.bam
tuxedo_hard_clip.ETOH_shMED1_1.67840f50fcd002400d1307172012d86f.mugqic.done
)
bam_hard_clip_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bam_hard_clip_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bam_hard_clip_8_JOB_ID: tuxedo_hard_clip.ETOH_shNIPBL_2
#-------------------------------------------------------------------------------
JOB_NAME=tuxedo_hard_clip.ETOH_shNIPBL_2
JOB_DEPENDENCIES=$picard_mark_duplicates_8_JOB_ID
JOB_DONE=job_output/bam_hard_clip/tuxedo_hard_clip.ETOH_shNIPBL_2.2c64e5cb11b9bf47baa46a8fbb9216fb.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'tuxedo_hard_clip.ETOH_shNIPBL_2.2c64e5cb11b9bf47baa46a8fbb9216fb.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -h \
  alignment/ETOH_shNIPBL_2/ETOH_shNIPBL_2.sorted.mdup.bam | \
awk 'BEGIN {OFS="\t"} {if (substr($1,1,1)=="@") {print;next}; split($6,C,/[0-9]*/); split($6,L,/[SMDIN]/); if (C[2]=="S") {$10=substr($10,L[1]+1); $11=substr($11,L[1]+1)}; if (C[length(C)]=="S") {L1=length($10)-L[length(L)-1]; $10=substr($10,1,L1); $11=substr($11,1,L1); }; gsub(/[0-9]*S/,"",$6); print}'  | \
samtools view -hbS \
  - \
  > alignment/ETOH_shNIPBL_2/ETOH_shNIPBL_2.sorted.mdup.hardClip.bam
tuxedo_hard_clip.ETOH_shNIPBL_2.2c64e5cb11b9bf47baa46a8fbb9216fb.mugqic.done
)
bam_hard_clip_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bam_hard_clip_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bam_hard_clip_9_JOB_ID: tuxedo_hard_clip.ETOH_shNIPBL_1
#-------------------------------------------------------------------------------
JOB_NAME=tuxedo_hard_clip.ETOH_shNIPBL_1
JOB_DEPENDENCIES=$picard_mark_duplicates_9_JOB_ID
JOB_DONE=job_output/bam_hard_clip/tuxedo_hard_clip.ETOH_shNIPBL_1.da7d2e013568a776ff5bfd711b6691db.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'tuxedo_hard_clip.ETOH_shNIPBL_1.da7d2e013568a776ff5bfd711b6691db.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -h \
  alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1.sorted.mdup.bam | \
awk 'BEGIN {OFS="\t"} {if (substr($1,1,1)=="@") {print;next}; split($6,C,/[0-9]*/); split($6,L,/[SMDIN]/); if (C[2]=="S") {$10=substr($10,L[1]+1); $11=substr($11,L[1]+1)}; if (C[length(C)]=="S") {L1=length($10)-L[length(L)-1]; $10=substr($10,1,L1); $11=substr($11,1,L1); }; gsub(/[0-9]*S/,"",$6); print}'  | \
samtools view -hbS \
  - \
  > alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1.sorted.mdup.hardClip.bam
tuxedo_hard_clip.ETOH_shNIPBL_1.da7d2e013568a776ff5bfd711b6691db.mugqic.done
)
bam_hard_clip_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bam_hard_clip_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bam_hard_clip_10_JOB_ID: tuxedo_hard_clip.ETOH_shMED1_2
#-------------------------------------------------------------------------------
JOB_NAME=tuxedo_hard_clip.ETOH_shMED1_2
JOB_DEPENDENCIES=$picard_mark_duplicates_10_JOB_ID
JOB_DONE=job_output/bam_hard_clip/tuxedo_hard_clip.ETOH_shMED1_2.463d08357a7ce47d87c23483906f8898.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'tuxedo_hard_clip.ETOH_shMED1_2.463d08357a7ce47d87c23483906f8898.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -h \
  alignment/ETOH_shMED1_2/ETOH_shMED1_2.sorted.mdup.bam | \
awk 'BEGIN {OFS="\t"} {if (substr($1,1,1)=="@") {print;next}; split($6,C,/[0-9]*/); split($6,L,/[SMDIN]/); if (C[2]=="S") {$10=substr($10,L[1]+1); $11=substr($11,L[1]+1)}; if (C[length(C)]=="S") {L1=length($10)-L[length(L)-1]; $10=substr($10,1,L1); $11=substr($11,1,L1); }; gsub(/[0-9]*S/,"",$6); print}'  | \
samtools view -hbS \
  - \
  > alignment/ETOH_shMED1_2/ETOH_shMED1_2.sorted.mdup.hardClip.bam
tuxedo_hard_clip.ETOH_shMED1_2.463d08357a7ce47d87c23483906f8898.mugqic.done
)
bam_hard_clip_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bam_hard_clip_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bam_hard_clip_11_JOB_ID: tuxedo_hard_clip.ETOH_shSMC1A_3
#-------------------------------------------------------------------------------
JOB_NAME=tuxedo_hard_clip.ETOH_shSMC1A_3
JOB_DEPENDENCIES=$picard_mark_duplicates_11_JOB_ID
JOB_DONE=job_output/bam_hard_clip/tuxedo_hard_clip.ETOH_shSMC1A_3.22dac7bf0933fe6e1b36cbcc19810983.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'tuxedo_hard_clip.ETOH_shSMC1A_3.22dac7bf0933fe6e1b36cbcc19810983.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -h \
  alignment/ETOH_shSMC1A_3/ETOH_shSMC1A_3.sorted.mdup.bam | \
awk 'BEGIN {OFS="\t"} {if (substr($1,1,1)=="@") {print;next}; split($6,C,/[0-9]*/); split($6,L,/[SMDIN]/); if (C[2]=="S") {$10=substr($10,L[1]+1); $11=substr($11,L[1]+1)}; if (C[length(C)]=="S") {L1=length($10)-L[length(L)-1]; $10=substr($10,1,L1); $11=substr($11,1,L1); }; gsub(/[0-9]*S/,"",$6); print}'  | \
samtools view -hbS \
  - \
  > alignment/ETOH_shSMC1A_3/ETOH_shSMC1A_3.sorted.mdup.hardClip.bam
tuxedo_hard_clip.ETOH_shSMC1A_3.22dac7bf0933fe6e1b36cbcc19810983.mugqic.done
)
bam_hard_clip_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bam_hard_clip_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bam_hard_clip_12_JOB_ID: tuxedo_hard_clip.ETOH_shSMC1A_2
#-------------------------------------------------------------------------------
JOB_NAME=tuxedo_hard_clip.ETOH_shSMC1A_2
JOB_DEPENDENCIES=$picard_mark_duplicates_12_JOB_ID
JOB_DONE=job_output/bam_hard_clip/tuxedo_hard_clip.ETOH_shSMC1A_2.f6b3fd9015dd19836b071f83ee47208c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'tuxedo_hard_clip.ETOH_shSMC1A_2.f6b3fd9015dd19836b071f83ee47208c.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -h \
  alignment/ETOH_shSMC1A_2/ETOH_shSMC1A_2.sorted.mdup.bam | \
awk 'BEGIN {OFS="\t"} {if (substr($1,1,1)=="@") {print;next}; split($6,C,/[0-9]*/); split($6,L,/[SMDIN]/); if (C[2]=="S") {$10=substr($10,L[1]+1); $11=substr($11,L[1]+1)}; if (C[length(C)]=="S") {L1=length($10)-L[length(L)-1]; $10=substr($10,1,L1); $11=substr($11,1,L1); }; gsub(/[0-9]*S/,"",$6); print}'  | \
samtools view -hbS \
  - \
  > alignment/ETOH_shSMC1A_2/ETOH_shSMC1A_2.sorted.mdup.hardClip.bam
tuxedo_hard_clip.ETOH_shSMC1A_2.f6b3fd9015dd19836b071f83ee47208c.mugqic.done
)
bam_hard_clip_12_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bam_hard_clip_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: rnaseqc
#-------------------------------------------------------------------------------
STEP=rnaseqc
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: rnaseqc_1_JOB_ID: rnaseqc
#-------------------------------------------------------------------------------
JOB_NAME=rnaseqc
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_2_JOB_ID:$picard_mark_duplicates_3_JOB_ID:$picard_mark_duplicates_4_JOB_ID:$picard_mark_duplicates_5_JOB_ID:$picard_mark_duplicates_6_JOB_ID:$picard_mark_duplicates_7_JOB_ID:$picard_mark_duplicates_8_JOB_ID:$picard_mark_duplicates_9_JOB_ID:$picard_mark_duplicates_10_JOB_ID:$picard_mark_duplicates_11_JOB_ID:$picard_mark_duplicates_12_JOB_ID
JOB_DONE=job_output/rnaseqc/rnaseqc.a24cbfb49d415a98809304bc6ce2c13d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'rnaseqc.a24cbfb49d415a98809304bc6ce2c13d.mugqic.done'
module load mugqic/java/openjdk-jdk1.7.0_60 mugqic/bwa/0.7.12 mugqic/rnaseqc/1.1.8 && \
mkdir -p metrics/rnaseqRep && \
echo "Sample	BamFile	Note
DEX_nonMamm_1	alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.bam	RNAseq
DEX_shMED1_1	alignment/DEX_shMED1_1/DEX_shMED1_1.sorted.mdup.bam	RNAseq
DEX_shNIPBL_1	alignment/DEX_shNIPBL_1/DEX_shNIPBL_1.sorted.mdup.bam	RNAseq
DEX_shSMC1A_3	alignment/DEX_shSMC1A_3/DEX_shSMC1A_3.sorted.mdup.bam	RNAseq
ETOH_nonMamm_1	alignment/ETOH_nonMamm_1/ETOH_nonMamm_1.sorted.mdup.bam	RNAseq
ETOH_nonMamm_2	alignment/ETOH_nonMamm_2/ETOH_nonMamm_2.sorted.mdup.bam	RNAseq
ETOH_shMED1_1	alignment/ETOH_shMED1_1/ETOH_shMED1_1.sorted.mdup.bam	RNAseq
ETOH_shNIPBL_2	alignment/ETOH_shNIPBL_2/ETOH_shNIPBL_2.sorted.mdup.bam	RNAseq
ETOH_shNIPBL_1	alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1.sorted.mdup.bam	RNAseq
ETOH_shMED1_2	alignment/ETOH_shMED1_2/ETOH_shMED1_2.sorted.mdup.bam	RNAseq
ETOH_shSMC1A_3	alignment/ETOH_shSMC1A_3/ETOH_shSMC1A_3.sorted.mdup.bam	RNAseq
ETOH_shSMC1A_2	alignment/ETOH_shSMC1A_2/ETOH_shSMC1A_2.sorted.mdup.bam	RNAseq" \
  > alignment/rnaseqc.samples.txt && \
touch dummy_rRNA.fa && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=4 -Xmx27G -jar $RNASEQC_JAR \
  -n 1000 \
  -o metrics/rnaseqRep \
  -r /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  -s alignment/rnaseqc.samples.txt \
  -t /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.transcript_id.gtf \
  -ttype 2\
  -BWArRNA dummy_rRNA.fa && \
zip -r metrics/rnaseqRep.zip metrics/rnaseqRep
rnaseqc.a24cbfb49d415a98809304bc6ce2c13d.mugqic.done
)
rnaseqc_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:00 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$rnaseqc_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: rnaseqc_2_JOB_ID: rnaseqc_report
#-------------------------------------------------------------------------------
JOB_NAME=rnaseqc_report
JOB_DEPENDENCIES=$rnaseqc_1_JOB_ID
JOB_DONE=job_output/rnaseqc/rnaseqc_report.b8996c235432adcc10558e9a1e05cafc.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'rnaseqc_report.b8996c235432adcc10558e9a1e05cafc.mugqic.done'
module load mugqic/python/2.7.13 mugqic/pandoc/1.15.2 && \
mkdir -p report && \
cp metrics/rnaseqRep.zip report/reportRNAseqQC.zip && \
python -c 'import csv; csv_in = csv.DictReader(open("metrics/rnaseqRep/metrics.tsv"), delimiter="	")
print "	".join(["Sample", "Aligned Reads", "Alternative Alignments", "%", "rRNA Reads", "Coverage", "Exonic Rate", "Genes"])
print "\n".join(["	".join([
    line["Sample"],
    line["Mapped"],
    line["Alternative Aligments"],
    str(float(line["Alternative Aligments"]) / float(line["Mapped"]) * 100),
    line["rRNA"],
    line["Mean Per Base Cov."],
    line["Exonic Rate"],
    line["Genes Detected"]
]) for line in csv_in])' \
  > report/trimAlignmentTable.tsv.tmp && \
if [[ -f metrics/trimSampleTable.tsv ]]
then
  awk -F"	" 'FNR==NR{raw_reads[$1]=$2; surviving_reads[$1]=$3; surviving_pct[$1]=$4; next}{OFS="	"; if ($2=="Aligned Reads"){surviving_pct[$1]="%"; aligned_pct="%"; rrna_pct="%"} else {aligned_pct=($2 / surviving_reads[$1] * 100); rrna_pct=($5 / surviving_reads[$1] * 100)}; printf $1"	"raw_reads[$1]"	"surviving_reads[$1]"	"surviving_pct[$1]"	"$2"	"aligned_pct"	"$3"	"$4"	"$5"	"rrna_pct; for (i = 6; i<= NF; i++) {printf "	"$i}; print ""}' \
  metrics/trimSampleTable.tsv \
  report/trimAlignmentTable.tsv.tmp \
  > report/trimAlignmentTable.tsv
else
  cp report/trimAlignmentTable.tsv.tmp report/trimAlignmentTable.tsv
fi && \
rm report/trimAlignmentTable.tsv.tmp && \
trim_alignment_table_md=`if [[ -f metrics/trimSampleTable.tsv ]] ; then cut -f1-13 report/trimAlignmentTable.tsv | LC_NUMERIC=en_CA awk -F "	" '{OFS="|"; if (NR == 1) {$1 = $1; print $0; print "-----|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:"} else {print $1, sprintf("%\47d", $2), sprintf("%\47d", $3), sprintf("%.1f", $4), sprintf("%\47d", $5), sprintf("%.1f", $6), sprintf("%\47d", $7), sprintf("%.1f", $8), sprintf("%\47d", $9), sprintf("%.1f", $10), sprintf("%.2f", $11), sprintf("%.2f", $12), sprintf("%\47d", $13)}}' ; else cat report/trimAlignmentTable.tsv | LC_NUMERIC=en_CA awk -F "	" '{OFS="|"; if (NR == 1) {$1 = $1; print $0; print "-----|-----:|-----:|-----:|-----:|-----:|-----:|-----:"} else {print $1, sprintf("%\47d", $2), sprintf("%\47d", $3), sprintf("%.1f", $4), sprintf("%\47d", $5), sprintf("%.2f", $6), sprintf("%.2f", $7), $8}}' ; fi`
pandoc \
  /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.3.0/bfx/report/RnaSeq.rnaseqc.md \
  --template /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.3.0/bfx/report/RnaSeq.rnaseqc.md \
  --variable trim_alignment_table="$trim_alignment_table_md" \
  --to markdown \
  > report/RnaSeq.rnaseqc.md
rnaseqc_report.b8996c235432adcc10558e9a1e05cafc.mugqic.done
)
rnaseqc_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$rnaseqc_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: wiggle
#-------------------------------------------------------------------------------
STEP=wiggle
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: wiggle_1_JOB_ID: wiggle.DEX_nonMamm_1.forward_strandspec
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.DEX_nonMamm_1.forward_strandspec
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.DEX_nonMamm_1.forward_strandspec.6f29228edfd56eef45024b68fc34dc1e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.DEX_nonMamm_1.forward_strandspec.6f29228edfd56eef45024b68fc34dc1e.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
samtools view -bh -F 256 -f 81 \
  alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.bam \
  > alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.tmp1.forward.bam && \
samtools view -bh -F 256 -f 161 \
  alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.bam \
  > alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.tmp2.forward.bam && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/picard.jar MergeSamFiles \
 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.tmp1.forward.bam \
  INPUT=alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.tmp2.forward.bam \
 OUTPUT=alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.forward.bam \
 MAX_RECORDS_IN_RAM=5750000 && \
rm alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.tmp1.forward.bam alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.tmp2.forward.bam
wiggle.DEX_nonMamm_1.forward_strandspec.6f29228edfd56eef45024b68fc34dc1e.mugqic.done
)
wiggle_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_2_JOB_ID: wiggle.DEX_nonMamm_1.reverse_strandspec
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.DEX_nonMamm_1.reverse_strandspec
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.DEX_nonMamm_1.reverse_strandspec.9914a7e9fc06aaa45404b6d54cf3e3d2.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.DEX_nonMamm_1.reverse_strandspec.9914a7e9fc06aaa45404b6d54cf3e3d2.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p tracks/DEX_nonMamm_1 tracks/bigWig && \
samtools view -bh -F 256 -f 97 \
  alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.bam \
  > alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.tmp1.reverse.bam && \
samtools view -bh -F 256 -f 145 \
  alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.bam \
  > alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.tmp2.reverse.bam && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/picard.jar MergeSamFiles \
 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.tmp1.reverse.bam \
  INPUT=alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.tmp2.reverse.bam \
 OUTPUT=alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.reverse.bam \
 MAX_RECORDS_IN_RAM=5750000 && \
rm alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.tmp1.reverse.bam alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.tmp2.reverse.bam
wiggle.DEX_nonMamm_1.reverse_strandspec.9914a7e9fc06aaa45404b6d54cf3e3d2.mugqic.done
)
wiggle_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_3_JOB_ID: wiggle.DEX_nonMamm_1.forward
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.DEX_nonMamm_1.forward
JOB_DEPENDENCIES=$wiggle_1_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.DEX_nonMamm_1.forward.1386bffdb19d2363a053af6c41921b1a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.DEX_nonMamm_1.forward.1386bffdb19d2363a053af6c41921b1a.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/bedtools/2.25.0 mugqic/ucsc/v326 && \
mkdir -p tracks/DEX_nonMamm_1 tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 81  alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.forward.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.forward.bam \
  -g /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  > tracks/DEX_nonMamm_1/DEX_nonMamm_1.forward.bedGraph && \
sort -k1,1 -k2,2n tracks/DEX_nonMamm_1/DEX_nonMamm_1.forward.bedGraph > tracks/DEX_nonMamm_1/DEX_nonMamm_1.forward.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/DEX_nonMamm_1/DEX_nonMamm_1.forward.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  tracks/bigWig/DEX_nonMamm_1.forward.bw
wiggle.DEX_nonMamm_1.forward.1386bffdb19d2363a053af6c41921b1a.mugqic.done
)
wiggle_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_4_JOB_ID: wiggle.DEX_nonMamm_1.reverse
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.DEX_nonMamm_1.reverse
JOB_DEPENDENCIES=$wiggle_2_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.DEX_nonMamm_1.reverse.184c4a6c18a7c315b7e9951ce2625b62.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.DEX_nonMamm_1.reverse.184c4a6c18a7c315b7e9951ce2625b62.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/bedtools/2.25.0 mugqic/ucsc/v326 && \
mkdir -p tracks/DEX_nonMamm_1 tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 97  alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.reverse.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.reverse.bam \
  -g /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  > tracks/DEX_nonMamm_1/DEX_nonMamm_1.reverse.bedGraph && \
sort -k1,1 -k2,2n tracks/DEX_nonMamm_1/DEX_nonMamm_1.reverse.bedGraph > tracks/DEX_nonMamm_1/DEX_nonMamm_1.reverse.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/DEX_nonMamm_1/DEX_nonMamm_1.reverse.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  tracks/bigWig/DEX_nonMamm_1.reverse.bw
wiggle.DEX_nonMamm_1.reverse.184c4a6c18a7c315b7e9951ce2625b62.mugqic.done
)
wiggle_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_5_JOB_ID: wiggle.DEX_shMED1_1.forward_strandspec
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.DEX_shMED1_1.forward_strandspec
JOB_DEPENDENCIES=$picard_mark_duplicates_2_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.DEX_shMED1_1.forward_strandspec.63f93290254ed3b7f2acd4f894e81867.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.DEX_shMED1_1.forward_strandspec.63f93290254ed3b7f2acd4f894e81867.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
samtools view -bh -F 256 -f 81 \
  alignment/DEX_shMED1_1/DEX_shMED1_1.sorted.mdup.bam \
  > alignment/DEX_shMED1_1/DEX_shMED1_1.sorted.mdup.tmp1.forward.bam && \
samtools view -bh -F 256 -f 161 \
  alignment/DEX_shMED1_1/DEX_shMED1_1.sorted.mdup.bam \
  > alignment/DEX_shMED1_1/DEX_shMED1_1.sorted.mdup.tmp2.forward.bam && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/picard.jar MergeSamFiles \
 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/DEX_shMED1_1/DEX_shMED1_1.sorted.mdup.tmp1.forward.bam \
  INPUT=alignment/DEX_shMED1_1/DEX_shMED1_1.sorted.mdup.tmp2.forward.bam \
 OUTPUT=alignment/DEX_shMED1_1/DEX_shMED1_1.sorted.mdup.forward.bam \
 MAX_RECORDS_IN_RAM=5750000 && \
rm alignment/DEX_shMED1_1/DEX_shMED1_1.sorted.mdup.tmp1.forward.bam alignment/DEX_shMED1_1/DEX_shMED1_1.sorted.mdup.tmp2.forward.bam
wiggle.DEX_shMED1_1.forward_strandspec.63f93290254ed3b7f2acd4f894e81867.mugqic.done
)
wiggle_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_6_JOB_ID: wiggle.DEX_shMED1_1.reverse_strandspec
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.DEX_shMED1_1.reverse_strandspec
JOB_DEPENDENCIES=$picard_mark_duplicates_2_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.DEX_shMED1_1.reverse_strandspec.0e590e70e920965b2c078ae9666a9778.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.DEX_shMED1_1.reverse_strandspec.0e590e70e920965b2c078ae9666a9778.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p tracks/DEX_shMED1_1 tracks/bigWig && \
samtools view -bh -F 256 -f 97 \
  alignment/DEX_shMED1_1/DEX_shMED1_1.sorted.mdup.bam \
  > alignment/DEX_shMED1_1/DEX_shMED1_1.sorted.mdup.tmp1.reverse.bam && \
samtools view -bh -F 256 -f 145 \
  alignment/DEX_shMED1_1/DEX_shMED1_1.sorted.mdup.bam \
  > alignment/DEX_shMED1_1/DEX_shMED1_1.sorted.mdup.tmp2.reverse.bam && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/picard.jar MergeSamFiles \
 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/DEX_shMED1_1/DEX_shMED1_1.sorted.mdup.tmp1.reverse.bam \
  INPUT=alignment/DEX_shMED1_1/DEX_shMED1_1.sorted.mdup.tmp2.reverse.bam \
 OUTPUT=alignment/DEX_shMED1_1/DEX_shMED1_1.sorted.mdup.reverse.bam \
 MAX_RECORDS_IN_RAM=5750000 && \
rm alignment/DEX_shMED1_1/DEX_shMED1_1.sorted.mdup.tmp1.reverse.bam alignment/DEX_shMED1_1/DEX_shMED1_1.sorted.mdup.tmp2.reverse.bam
wiggle.DEX_shMED1_1.reverse_strandspec.0e590e70e920965b2c078ae9666a9778.mugqic.done
)
wiggle_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_7_JOB_ID: wiggle.DEX_shMED1_1.forward
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.DEX_shMED1_1.forward
JOB_DEPENDENCIES=$wiggle_5_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.DEX_shMED1_1.forward.8d163db4ec5e73f76d075e0ce8857989.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.DEX_shMED1_1.forward.8d163db4ec5e73f76d075e0ce8857989.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/bedtools/2.25.0 mugqic/ucsc/v326 && \
mkdir -p tracks/DEX_shMED1_1 tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 81  alignment/DEX_shMED1_1/DEX_shMED1_1.sorted.mdup.forward.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/DEX_shMED1_1/DEX_shMED1_1.sorted.mdup.forward.bam \
  -g /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  > tracks/DEX_shMED1_1/DEX_shMED1_1.forward.bedGraph && \
sort -k1,1 -k2,2n tracks/DEX_shMED1_1/DEX_shMED1_1.forward.bedGraph > tracks/DEX_shMED1_1/DEX_shMED1_1.forward.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/DEX_shMED1_1/DEX_shMED1_1.forward.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  tracks/bigWig/DEX_shMED1_1.forward.bw
wiggle.DEX_shMED1_1.forward.8d163db4ec5e73f76d075e0ce8857989.mugqic.done
)
wiggle_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_8_JOB_ID: wiggle.DEX_shMED1_1.reverse
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.DEX_shMED1_1.reverse
JOB_DEPENDENCIES=$wiggle_6_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.DEX_shMED1_1.reverse.038efe4af703552bb2e15d356e144716.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.DEX_shMED1_1.reverse.038efe4af703552bb2e15d356e144716.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/bedtools/2.25.0 mugqic/ucsc/v326 && \
mkdir -p tracks/DEX_shMED1_1 tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 97  alignment/DEX_shMED1_1/DEX_shMED1_1.sorted.mdup.reverse.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/DEX_shMED1_1/DEX_shMED1_1.sorted.mdup.reverse.bam \
  -g /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  > tracks/DEX_shMED1_1/DEX_shMED1_1.reverse.bedGraph && \
sort -k1,1 -k2,2n tracks/DEX_shMED1_1/DEX_shMED1_1.reverse.bedGraph > tracks/DEX_shMED1_1/DEX_shMED1_1.reverse.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/DEX_shMED1_1/DEX_shMED1_1.reverse.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  tracks/bigWig/DEX_shMED1_1.reverse.bw
wiggle.DEX_shMED1_1.reverse.038efe4af703552bb2e15d356e144716.mugqic.done
)
wiggle_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_9_JOB_ID: wiggle.DEX_shNIPBL_1.forward_strandspec
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.DEX_shNIPBL_1.forward_strandspec
JOB_DEPENDENCIES=$picard_mark_duplicates_3_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.DEX_shNIPBL_1.forward_strandspec.a3f5ddba53344029c55b6dbe95ee86d8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.DEX_shNIPBL_1.forward_strandspec.a3f5ddba53344029c55b6dbe95ee86d8.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
samtools view -bh -F 256 -f 81 \
  alignment/DEX_shNIPBL_1/DEX_shNIPBL_1.sorted.mdup.bam \
  > alignment/DEX_shNIPBL_1/DEX_shNIPBL_1.sorted.mdup.tmp1.forward.bam && \
samtools view -bh -F 256 -f 161 \
  alignment/DEX_shNIPBL_1/DEX_shNIPBL_1.sorted.mdup.bam \
  > alignment/DEX_shNIPBL_1/DEX_shNIPBL_1.sorted.mdup.tmp2.forward.bam && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/picard.jar MergeSamFiles \
 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/DEX_shNIPBL_1/DEX_shNIPBL_1.sorted.mdup.tmp1.forward.bam \
  INPUT=alignment/DEX_shNIPBL_1/DEX_shNIPBL_1.sorted.mdup.tmp2.forward.bam \
 OUTPUT=alignment/DEX_shNIPBL_1/DEX_shNIPBL_1.sorted.mdup.forward.bam \
 MAX_RECORDS_IN_RAM=5750000 && \
rm alignment/DEX_shNIPBL_1/DEX_shNIPBL_1.sorted.mdup.tmp1.forward.bam alignment/DEX_shNIPBL_1/DEX_shNIPBL_1.sorted.mdup.tmp2.forward.bam
wiggle.DEX_shNIPBL_1.forward_strandspec.a3f5ddba53344029c55b6dbe95ee86d8.mugqic.done
)
wiggle_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_10_JOB_ID: wiggle.DEX_shNIPBL_1.reverse_strandspec
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.DEX_shNIPBL_1.reverse_strandspec
JOB_DEPENDENCIES=$picard_mark_duplicates_3_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.DEX_shNIPBL_1.reverse_strandspec.a4008734beff3a15a510550f39b1f9cb.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.DEX_shNIPBL_1.reverse_strandspec.a4008734beff3a15a510550f39b1f9cb.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p tracks/DEX_shNIPBL_1 tracks/bigWig && \
samtools view -bh -F 256 -f 97 \
  alignment/DEX_shNIPBL_1/DEX_shNIPBL_1.sorted.mdup.bam \
  > alignment/DEX_shNIPBL_1/DEX_shNIPBL_1.sorted.mdup.tmp1.reverse.bam && \
samtools view -bh -F 256 -f 145 \
  alignment/DEX_shNIPBL_1/DEX_shNIPBL_1.sorted.mdup.bam \
  > alignment/DEX_shNIPBL_1/DEX_shNIPBL_1.sorted.mdup.tmp2.reverse.bam && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/picard.jar MergeSamFiles \
 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/DEX_shNIPBL_1/DEX_shNIPBL_1.sorted.mdup.tmp1.reverse.bam \
  INPUT=alignment/DEX_shNIPBL_1/DEX_shNIPBL_1.sorted.mdup.tmp2.reverse.bam \
 OUTPUT=alignment/DEX_shNIPBL_1/DEX_shNIPBL_1.sorted.mdup.reverse.bam \
 MAX_RECORDS_IN_RAM=5750000 && \
rm alignment/DEX_shNIPBL_1/DEX_shNIPBL_1.sorted.mdup.tmp1.reverse.bam alignment/DEX_shNIPBL_1/DEX_shNIPBL_1.sorted.mdup.tmp2.reverse.bam
wiggle.DEX_shNIPBL_1.reverse_strandspec.a4008734beff3a15a510550f39b1f9cb.mugqic.done
)
wiggle_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_11_JOB_ID: wiggle.DEX_shNIPBL_1.forward
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.DEX_shNIPBL_1.forward
JOB_DEPENDENCIES=$wiggle_9_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.DEX_shNIPBL_1.forward.f1d24fce0af2acce39dd6be2cb2c5aa4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.DEX_shNIPBL_1.forward.f1d24fce0af2acce39dd6be2cb2c5aa4.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/bedtools/2.25.0 mugqic/ucsc/v326 && \
mkdir -p tracks/DEX_shNIPBL_1 tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 81  alignment/DEX_shNIPBL_1/DEX_shNIPBL_1.sorted.mdup.forward.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/DEX_shNIPBL_1/DEX_shNIPBL_1.sorted.mdup.forward.bam \
  -g /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  > tracks/DEX_shNIPBL_1/DEX_shNIPBL_1.forward.bedGraph && \
sort -k1,1 -k2,2n tracks/DEX_shNIPBL_1/DEX_shNIPBL_1.forward.bedGraph > tracks/DEX_shNIPBL_1/DEX_shNIPBL_1.forward.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/DEX_shNIPBL_1/DEX_shNIPBL_1.forward.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  tracks/bigWig/DEX_shNIPBL_1.forward.bw
wiggle.DEX_shNIPBL_1.forward.f1d24fce0af2acce39dd6be2cb2c5aa4.mugqic.done
)
wiggle_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_12_JOB_ID: wiggle.DEX_shNIPBL_1.reverse
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.DEX_shNIPBL_1.reverse
JOB_DEPENDENCIES=$wiggle_10_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.DEX_shNIPBL_1.reverse.f9207521e1f7f9db50d1f4c842dbe0ae.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.DEX_shNIPBL_1.reverse.f9207521e1f7f9db50d1f4c842dbe0ae.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/bedtools/2.25.0 mugqic/ucsc/v326 && \
mkdir -p tracks/DEX_shNIPBL_1 tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 97  alignment/DEX_shNIPBL_1/DEX_shNIPBL_1.sorted.mdup.reverse.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/DEX_shNIPBL_1/DEX_shNIPBL_1.sorted.mdup.reverse.bam \
  -g /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  > tracks/DEX_shNIPBL_1/DEX_shNIPBL_1.reverse.bedGraph && \
sort -k1,1 -k2,2n tracks/DEX_shNIPBL_1/DEX_shNIPBL_1.reverse.bedGraph > tracks/DEX_shNIPBL_1/DEX_shNIPBL_1.reverse.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/DEX_shNIPBL_1/DEX_shNIPBL_1.reverse.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  tracks/bigWig/DEX_shNIPBL_1.reverse.bw
wiggle.DEX_shNIPBL_1.reverse.f9207521e1f7f9db50d1f4c842dbe0ae.mugqic.done
)
wiggle_12_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_13_JOB_ID: wiggle.DEX_shSMC1A_3.forward_strandspec
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.DEX_shSMC1A_3.forward_strandspec
JOB_DEPENDENCIES=$picard_mark_duplicates_4_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.DEX_shSMC1A_3.forward_strandspec.3cb8472038afa769e7a68890fc8f5a2b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.DEX_shSMC1A_3.forward_strandspec.3cb8472038afa769e7a68890fc8f5a2b.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
samtools view -bh -F 256 -f 81 \
  alignment/DEX_shSMC1A_3/DEX_shSMC1A_3.sorted.mdup.bam \
  > alignment/DEX_shSMC1A_3/DEX_shSMC1A_3.sorted.mdup.tmp1.forward.bam && \
samtools view -bh -F 256 -f 161 \
  alignment/DEX_shSMC1A_3/DEX_shSMC1A_3.sorted.mdup.bam \
  > alignment/DEX_shSMC1A_3/DEX_shSMC1A_3.sorted.mdup.tmp2.forward.bam && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/picard.jar MergeSamFiles \
 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/DEX_shSMC1A_3/DEX_shSMC1A_3.sorted.mdup.tmp1.forward.bam \
  INPUT=alignment/DEX_shSMC1A_3/DEX_shSMC1A_3.sorted.mdup.tmp2.forward.bam \
 OUTPUT=alignment/DEX_shSMC1A_3/DEX_shSMC1A_3.sorted.mdup.forward.bam \
 MAX_RECORDS_IN_RAM=5750000 && \
rm alignment/DEX_shSMC1A_3/DEX_shSMC1A_3.sorted.mdup.tmp1.forward.bam alignment/DEX_shSMC1A_3/DEX_shSMC1A_3.sorted.mdup.tmp2.forward.bam
wiggle.DEX_shSMC1A_3.forward_strandspec.3cb8472038afa769e7a68890fc8f5a2b.mugqic.done
)
wiggle_13_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_13_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_14_JOB_ID: wiggle.DEX_shSMC1A_3.reverse_strandspec
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.DEX_shSMC1A_3.reverse_strandspec
JOB_DEPENDENCIES=$picard_mark_duplicates_4_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.DEX_shSMC1A_3.reverse_strandspec.985f8ac0fb00dbd412920cbbfe16181d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.DEX_shSMC1A_3.reverse_strandspec.985f8ac0fb00dbd412920cbbfe16181d.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p tracks/DEX_shSMC1A_3 tracks/bigWig && \
samtools view -bh -F 256 -f 97 \
  alignment/DEX_shSMC1A_3/DEX_shSMC1A_3.sorted.mdup.bam \
  > alignment/DEX_shSMC1A_3/DEX_shSMC1A_3.sorted.mdup.tmp1.reverse.bam && \
samtools view -bh -F 256 -f 145 \
  alignment/DEX_shSMC1A_3/DEX_shSMC1A_3.sorted.mdup.bam \
  > alignment/DEX_shSMC1A_3/DEX_shSMC1A_3.sorted.mdup.tmp2.reverse.bam && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/picard.jar MergeSamFiles \
 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/DEX_shSMC1A_3/DEX_shSMC1A_3.sorted.mdup.tmp1.reverse.bam \
  INPUT=alignment/DEX_shSMC1A_3/DEX_shSMC1A_3.sorted.mdup.tmp2.reverse.bam \
 OUTPUT=alignment/DEX_shSMC1A_3/DEX_shSMC1A_3.sorted.mdup.reverse.bam \
 MAX_RECORDS_IN_RAM=5750000 && \
rm alignment/DEX_shSMC1A_3/DEX_shSMC1A_3.sorted.mdup.tmp1.reverse.bam alignment/DEX_shSMC1A_3/DEX_shSMC1A_3.sorted.mdup.tmp2.reverse.bam
wiggle.DEX_shSMC1A_3.reverse_strandspec.985f8ac0fb00dbd412920cbbfe16181d.mugqic.done
)
wiggle_14_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_14_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_15_JOB_ID: wiggle.DEX_shSMC1A_3.forward
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.DEX_shSMC1A_3.forward
JOB_DEPENDENCIES=$wiggle_13_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.DEX_shSMC1A_3.forward.c8a5254924956b8eef43055c6c7ba905.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.DEX_shSMC1A_3.forward.c8a5254924956b8eef43055c6c7ba905.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/bedtools/2.25.0 mugqic/ucsc/v326 && \
mkdir -p tracks/DEX_shSMC1A_3 tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 81  alignment/DEX_shSMC1A_3/DEX_shSMC1A_3.sorted.mdup.forward.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/DEX_shSMC1A_3/DEX_shSMC1A_3.sorted.mdup.forward.bam \
  -g /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  > tracks/DEX_shSMC1A_3/DEX_shSMC1A_3.forward.bedGraph && \
sort -k1,1 -k2,2n tracks/DEX_shSMC1A_3/DEX_shSMC1A_3.forward.bedGraph > tracks/DEX_shSMC1A_3/DEX_shSMC1A_3.forward.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/DEX_shSMC1A_3/DEX_shSMC1A_3.forward.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  tracks/bigWig/DEX_shSMC1A_3.forward.bw
wiggle.DEX_shSMC1A_3.forward.c8a5254924956b8eef43055c6c7ba905.mugqic.done
)
wiggle_15_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_15_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_16_JOB_ID: wiggle.DEX_shSMC1A_3.reverse
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.DEX_shSMC1A_3.reverse
JOB_DEPENDENCIES=$wiggle_14_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.DEX_shSMC1A_3.reverse.1613c70d295fdbc0110e79aded9aa6c2.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.DEX_shSMC1A_3.reverse.1613c70d295fdbc0110e79aded9aa6c2.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/bedtools/2.25.0 mugqic/ucsc/v326 && \
mkdir -p tracks/DEX_shSMC1A_3 tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 97  alignment/DEX_shSMC1A_3/DEX_shSMC1A_3.sorted.mdup.reverse.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/DEX_shSMC1A_3/DEX_shSMC1A_3.sorted.mdup.reverse.bam \
  -g /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  > tracks/DEX_shSMC1A_3/DEX_shSMC1A_3.reverse.bedGraph && \
sort -k1,1 -k2,2n tracks/DEX_shSMC1A_3/DEX_shSMC1A_3.reverse.bedGraph > tracks/DEX_shSMC1A_3/DEX_shSMC1A_3.reverse.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/DEX_shSMC1A_3/DEX_shSMC1A_3.reverse.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  tracks/bigWig/DEX_shSMC1A_3.reverse.bw
wiggle.DEX_shSMC1A_3.reverse.1613c70d295fdbc0110e79aded9aa6c2.mugqic.done
)
wiggle_16_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_16_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_17_JOB_ID: wiggle.ETOH_nonMamm_1.forward_strandspec
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.ETOH_nonMamm_1.forward_strandspec
JOB_DEPENDENCIES=$picard_mark_duplicates_5_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.ETOH_nonMamm_1.forward_strandspec.d10ec11faa70f0e740d30e14b6fd051a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.ETOH_nonMamm_1.forward_strandspec.d10ec11faa70f0e740d30e14b6fd051a.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
samtools view -bh -F 256 -f 81 \
  alignment/ETOH_nonMamm_1/ETOH_nonMamm_1.sorted.mdup.bam \
  > alignment/ETOH_nonMamm_1/ETOH_nonMamm_1.sorted.mdup.tmp1.forward.bam && \
samtools view -bh -F 256 -f 161 \
  alignment/ETOH_nonMamm_1/ETOH_nonMamm_1.sorted.mdup.bam \
  > alignment/ETOH_nonMamm_1/ETOH_nonMamm_1.sorted.mdup.tmp2.forward.bam && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/picard.jar MergeSamFiles \
 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/ETOH_nonMamm_1/ETOH_nonMamm_1.sorted.mdup.tmp1.forward.bam \
  INPUT=alignment/ETOH_nonMamm_1/ETOH_nonMamm_1.sorted.mdup.tmp2.forward.bam \
 OUTPUT=alignment/ETOH_nonMamm_1/ETOH_nonMamm_1.sorted.mdup.forward.bam \
 MAX_RECORDS_IN_RAM=5750000 && \
rm alignment/ETOH_nonMamm_1/ETOH_nonMamm_1.sorted.mdup.tmp1.forward.bam alignment/ETOH_nonMamm_1/ETOH_nonMamm_1.sorted.mdup.tmp2.forward.bam
wiggle.ETOH_nonMamm_1.forward_strandspec.d10ec11faa70f0e740d30e14b6fd051a.mugqic.done
)
wiggle_17_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_17_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_18_JOB_ID: wiggle.ETOH_nonMamm_1.reverse_strandspec
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.ETOH_nonMamm_1.reverse_strandspec
JOB_DEPENDENCIES=$picard_mark_duplicates_5_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.ETOH_nonMamm_1.reverse_strandspec.40f6cbb7974f8b80a2db4bfaa2914958.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.ETOH_nonMamm_1.reverse_strandspec.40f6cbb7974f8b80a2db4bfaa2914958.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p tracks/ETOH_nonMamm_1 tracks/bigWig && \
samtools view -bh -F 256 -f 97 \
  alignment/ETOH_nonMamm_1/ETOH_nonMamm_1.sorted.mdup.bam \
  > alignment/ETOH_nonMamm_1/ETOH_nonMamm_1.sorted.mdup.tmp1.reverse.bam && \
samtools view -bh -F 256 -f 145 \
  alignment/ETOH_nonMamm_1/ETOH_nonMamm_1.sorted.mdup.bam \
  > alignment/ETOH_nonMamm_1/ETOH_nonMamm_1.sorted.mdup.tmp2.reverse.bam && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/picard.jar MergeSamFiles \
 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/ETOH_nonMamm_1/ETOH_nonMamm_1.sorted.mdup.tmp1.reverse.bam \
  INPUT=alignment/ETOH_nonMamm_1/ETOH_nonMamm_1.sorted.mdup.tmp2.reverse.bam \
 OUTPUT=alignment/ETOH_nonMamm_1/ETOH_nonMamm_1.sorted.mdup.reverse.bam \
 MAX_RECORDS_IN_RAM=5750000 && \
rm alignment/ETOH_nonMamm_1/ETOH_nonMamm_1.sorted.mdup.tmp1.reverse.bam alignment/ETOH_nonMamm_1/ETOH_nonMamm_1.sorted.mdup.tmp2.reverse.bam
wiggle.ETOH_nonMamm_1.reverse_strandspec.40f6cbb7974f8b80a2db4bfaa2914958.mugqic.done
)
wiggle_18_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_18_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_19_JOB_ID: wiggle.ETOH_nonMamm_1.forward
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.ETOH_nonMamm_1.forward
JOB_DEPENDENCIES=$wiggle_17_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.ETOH_nonMamm_1.forward.7eca418dd992609f7ee6518294989225.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.ETOH_nonMamm_1.forward.7eca418dd992609f7ee6518294989225.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/bedtools/2.25.0 mugqic/ucsc/v326 && \
mkdir -p tracks/ETOH_nonMamm_1 tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 81  alignment/ETOH_nonMamm_1/ETOH_nonMamm_1.sorted.mdup.forward.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/ETOH_nonMamm_1/ETOH_nonMamm_1.sorted.mdup.forward.bam \
  -g /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  > tracks/ETOH_nonMamm_1/ETOH_nonMamm_1.forward.bedGraph && \
sort -k1,1 -k2,2n tracks/ETOH_nonMamm_1/ETOH_nonMamm_1.forward.bedGraph > tracks/ETOH_nonMamm_1/ETOH_nonMamm_1.forward.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/ETOH_nonMamm_1/ETOH_nonMamm_1.forward.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  tracks/bigWig/ETOH_nonMamm_1.forward.bw
wiggle.ETOH_nonMamm_1.forward.7eca418dd992609f7ee6518294989225.mugqic.done
)
wiggle_19_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_19_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_20_JOB_ID: wiggle.ETOH_nonMamm_1.reverse
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.ETOH_nonMamm_1.reverse
JOB_DEPENDENCIES=$wiggle_18_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.ETOH_nonMamm_1.reverse.acc3fceeaed9285d3f2a789a05658060.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.ETOH_nonMamm_1.reverse.acc3fceeaed9285d3f2a789a05658060.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/bedtools/2.25.0 mugqic/ucsc/v326 && \
mkdir -p tracks/ETOH_nonMamm_1 tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 97  alignment/ETOH_nonMamm_1/ETOH_nonMamm_1.sorted.mdup.reverse.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/ETOH_nonMamm_1/ETOH_nonMamm_1.sorted.mdup.reverse.bam \
  -g /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  > tracks/ETOH_nonMamm_1/ETOH_nonMamm_1.reverse.bedGraph && \
sort -k1,1 -k2,2n tracks/ETOH_nonMamm_1/ETOH_nonMamm_1.reverse.bedGraph > tracks/ETOH_nonMamm_1/ETOH_nonMamm_1.reverse.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/ETOH_nonMamm_1/ETOH_nonMamm_1.reverse.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  tracks/bigWig/ETOH_nonMamm_1.reverse.bw
wiggle.ETOH_nonMamm_1.reverse.acc3fceeaed9285d3f2a789a05658060.mugqic.done
)
wiggle_20_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_20_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_21_JOB_ID: wiggle.ETOH_nonMamm_2.forward_strandspec
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.ETOH_nonMamm_2.forward_strandspec
JOB_DEPENDENCIES=$picard_mark_duplicates_6_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.ETOH_nonMamm_2.forward_strandspec.1cdc628ca483c1172bc6ebd655ef3855.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.ETOH_nonMamm_2.forward_strandspec.1cdc628ca483c1172bc6ebd655ef3855.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
samtools view -bh -F 256 -f 81 \
  alignment/ETOH_nonMamm_2/ETOH_nonMamm_2.sorted.mdup.bam \
  > alignment/ETOH_nonMamm_2/ETOH_nonMamm_2.sorted.mdup.tmp1.forward.bam && \
samtools view -bh -F 256 -f 161 \
  alignment/ETOH_nonMamm_2/ETOH_nonMamm_2.sorted.mdup.bam \
  > alignment/ETOH_nonMamm_2/ETOH_nonMamm_2.sorted.mdup.tmp2.forward.bam && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/picard.jar MergeSamFiles \
 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/ETOH_nonMamm_2/ETOH_nonMamm_2.sorted.mdup.tmp1.forward.bam \
  INPUT=alignment/ETOH_nonMamm_2/ETOH_nonMamm_2.sorted.mdup.tmp2.forward.bam \
 OUTPUT=alignment/ETOH_nonMamm_2/ETOH_nonMamm_2.sorted.mdup.forward.bam \
 MAX_RECORDS_IN_RAM=5750000 && \
rm alignment/ETOH_nonMamm_2/ETOH_nonMamm_2.sorted.mdup.tmp1.forward.bam alignment/ETOH_nonMamm_2/ETOH_nonMamm_2.sorted.mdup.tmp2.forward.bam
wiggle.ETOH_nonMamm_2.forward_strandspec.1cdc628ca483c1172bc6ebd655ef3855.mugqic.done
)
wiggle_21_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_21_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_22_JOB_ID: wiggle.ETOH_nonMamm_2.reverse_strandspec
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.ETOH_nonMamm_2.reverse_strandspec
JOB_DEPENDENCIES=$picard_mark_duplicates_6_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.ETOH_nonMamm_2.reverse_strandspec.8dfe20db2df287c8478cd4762a40fc85.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.ETOH_nonMamm_2.reverse_strandspec.8dfe20db2df287c8478cd4762a40fc85.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p tracks/ETOH_nonMamm_2 tracks/bigWig && \
samtools view -bh -F 256 -f 97 \
  alignment/ETOH_nonMamm_2/ETOH_nonMamm_2.sorted.mdup.bam \
  > alignment/ETOH_nonMamm_2/ETOH_nonMamm_2.sorted.mdup.tmp1.reverse.bam && \
samtools view -bh -F 256 -f 145 \
  alignment/ETOH_nonMamm_2/ETOH_nonMamm_2.sorted.mdup.bam \
  > alignment/ETOH_nonMamm_2/ETOH_nonMamm_2.sorted.mdup.tmp2.reverse.bam && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/picard.jar MergeSamFiles \
 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/ETOH_nonMamm_2/ETOH_nonMamm_2.sorted.mdup.tmp1.reverse.bam \
  INPUT=alignment/ETOH_nonMamm_2/ETOH_nonMamm_2.sorted.mdup.tmp2.reverse.bam \
 OUTPUT=alignment/ETOH_nonMamm_2/ETOH_nonMamm_2.sorted.mdup.reverse.bam \
 MAX_RECORDS_IN_RAM=5750000 && \
rm alignment/ETOH_nonMamm_2/ETOH_nonMamm_2.sorted.mdup.tmp1.reverse.bam alignment/ETOH_nonMamm_2/ETOH_nonMamm_2.sorted.mdup.tmp2.reverse.bam
wiggle.ETOH_nonMamm_2.reverse_strandspec.8dfe20db2df287c8478cd4762a40fc85.mugqic.done
)
wiggle_22_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_22_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_23_JOB_ID: wiggle.ETOH_nonMamm_2.forward
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.ETOH_nonMamm_2.forward
JOB_DEPENDENCIES=$wiggle_21_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.ETOH_nonMamm_2.forward.da5488a63cb49a5495bf61d7ceef2f75.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.ETOH_nonMamm_2.forward.da5488a63cb49a5495bf61d7ceef2f75.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/bedtools/2.25.0 mugqic/ucsc/v326 && \
mkdir -p tracks/ETOH_nonMamm_2 tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 81  alignment/ETOH_nonMamm_2/ETOH_nonMamm_2.sorted.mdup.forward.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/ETOH_nonMamm_2/ETOH_nonMamm_2.sorted.mdup.forward.bam \
  -g /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  > tracks/ETOH_nonMamm_2/ETOH_nonMamm_2.forward.bedGraph && \
sort -k1,1 -k2,2n tracks/ETOH_nonMamm_2/ETOH_nonMamm_2.forward.bedGraph > tracks/ETOH_nonMamm_2/ETOH_nonMamm_2.forward.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/ETOH_nonMamm_2/ETOH_nonMamm_2.forward.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  tracks/bigWig/ETOH_nonMamm_2.forward.bw
wiggle.ETOH_nonMamm_2.forward.da5488a63cb49a5495bf61d7ceef2f75.mugqic.done
)
wiggle_23_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_23_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_24_JOB_ID: wiggle.ETOH_nonMamm_2.reverse
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.ETOH_nonMamm_2.reverse
JOB_DEPENDENCIES=$wiggle_22_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.ETOH_nonMamm_2.reverse.a8af2dedc708d40aa65689da9282307e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.ETOH_nonMamm_2.reverse.a8af2dedc708d40aa65689da9282307e.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/bedtools/2.25.0 mugqic/ucsc/v326 && \
mkdir -p tracks/ETOH_nonMamm_2 tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 97  alignment/ETOH_nonMamm_2/ETOH_nonMamm_2.sorted.mdup.reverse.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/ETOH_nonMamm_2/ETOH_nonMamm_2.sorted.mdup.reverse.bam \
  -g /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  > tracks/ETOH_nonMamm_2/ETOH_nonMamm_2.reverse.bedGraph && \
sort -k1,1 -k2,2n tracks/ETOH_nonMamm_2/ETOH_nonMamm_2.reverse.bedGraph > tracks/ETOH_nonMamm_2/ETOH_nonMamm_2.reverse.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/ETOH_nonMamm_2/ETOH_nonMamm_2.reverse.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  tracks/bigWig/ETOH_nonMamm_2.reverse.bw
wiggle.ETOH_nonMamm_2.reverse.a8af2dedc708d40aa65689da9282307e.mugqic.done
)
wiggle_24_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_24_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_25_JOB_ID: wiggle.ETOH_shMED1_1.forward_strandspec
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.ETOH_shMED1_1.forward_strandspec
JOB_DEPENDENCIES=$picard_mark_duplicates_7_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.ETOH_shMED1_1.forward_strandspec.9124ed26457f74dbac45bce4c8deede4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.ETOH_shMED1_1.forward_strandspec.9124ed26457f74dbac45bce4c8deede4.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
samtools view -bh -F 256 -f 81 \
  alignment/ETOH_shMED1_1/ETOH_shMED1_1.sorted.mdup.bam \
  > alignment/ETOH_shMED1_1/ETOH_shMED1_1.sorted.mdup.tmp1.forward.bam && \
samtools view -bh -F 256 -f 161 \
  alignment/ETOH_shMED1_1/ETOH_shMED1_1.sorted.mdup.bam \
  > alignment/ETOH_shMED1_1/ETOH_shMED1_1.sorted.mdup.tmp2.forward.bam && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/picard.jar MergeSamFiles \
 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/ETOH_shMED1_1/ETOH_shMED1_1.sorted.mdup.tmp1.forward.bam \
  INPUT=alignment/ETOH_shMED1_1/ETOH_shMED1_1.sorted.mdup.tmp2.forward.bam \
 OUTPUT=alignment/ETOH_shMED1_1/ETOH_shMED1_1.sorted.mdup.forward.bam \
 MAX_RECORDS_IN_RAM=5750000 && \
rm alignment/ETOH_shMED1_1/ETOH_shMED1_1.sorted.mdup.tmp1.forward.bam alignment/ETOH_shMED1_1/ETOH_shMED1_1.sorted.mdup.tmp2.forward.bam
wiggle.ETOH_shMED1_1.forward_strandspec.9124ed26457f74dbac45bce4c8deede4.mugqic.done
)
wiggle_25_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_25_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_26_JOB_ID: wiggle.ETOH_shMED1_1.reverse_strandspec
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.ETOH_shMED1_1.reverse_strandspec
JOB_DEPENDENCIES=$picard_mark_duplicates_7_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.ETOH_shMED1_1.reverse_strandspec.9871396af5cb4e29877851d0df076eb1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.ETOH_shMED1_1.reverse_strandspec.9871396af5cb4e29877851d0df076eb1.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p tracks/ETOH_shMED1_1 tracks/bigWig && \
samtools view -bh -F 256 -f 97 \
  alignment/ETOH_shMED1_1/ETOH_shMED1_1.sorted.mdup.bam \
  > alignment/ETOH_shMED1_1/ETOH_shMED1_1.sorted.mdup.tmp1.reverse.bam && \
samtools view -bh -F 256 -f 145 \
  alignment/ETOH_shMED1_1/ETOH_shMED1_1.sorted.mdup.bam \
  > alignment/ETOH_shMED1_1/ETOH_shMED1_1.sorted.mdup.tmp2.reverse.bam && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/picard.jar MergeSamFiles \
 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/ETOH_shMED1_1/ETOH_shMED1_1.sorted.mdup.tmp1.reverse.bam \
  INPUT=alignment/ETOH_shMED1_1/ETOH_shMED1_1.sorted.mdup.tmp2.reverse.bam \
 OUTPUT=alignment/ETOH_shMED1_1/ETOH_shMED1_1.sorted.mdup.reverse.bam \
 MAX_RECORDS_IN_RAM=5750000 && \
rm alignment/ETOH_shMED1_1/ETOH_shMED1_1.sorted.mdup.tmp1.reverse.bam alignment/ETOH_shMED1_1/ETOH_shMED1_1.sorted.mdup.tmp2.reverse.bam
wiggle.ETOH_shMED1_1.reverse_strandspec.9871396af5cb4e29877851d0df076eb1.mugqic.done
)
wiggle_26_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_26_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_27_JOB_ID: wiggle.ETOH_shMED1_1.forward
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.ETOH_shMED1_1.forward
JOB_DEPENDENCIES=$wiggle_25_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.ETOH_shMED1_1.forward.acb9f3623484a7099cfac9ef93384e63.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.ETOH_shMED1_1.forward.acb9f3623484a7099cfac9ef93384e63.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/bedtools/2.25.0 mugqic/ucsc/v326 && \
mkdir -p tracks/ETOH_shMED1_1 tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 81  alignment/ETOH_shMED1_1/ETOH_shMED1_1.sorted.mdup.forward.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/ETOH_shMED1_1/ETOH_shMED1_1.sorted.mdup.forward.bam \
  -g /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  > tracks/ETOH_shMED1_1/ETOH_shMED1_1.forward.bedGraph && \
sort -k1,1 -k2,2n tracks/ETOH_shMED1_1/ETOH_shMED1_1.forward.bedGraph > tracks/ETOH_shMED1_1/ETOH_shMED1_1.forward.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/ETOH_shMED1_1/ETOH_shMED1_1.forward.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  tracks/bigWig/ETOH_shMED1_1.forward.bw
wiggle.ETOH_shMED1_1.forward.acb9f3623484a7099cfac9ef93384e63.mugqic.done
)
wiggle_27_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_27_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_28_JOB_ID: wiggle.ETOH_shMED1_1.reverse
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.ETOH_shMED1_1.reverse
JOB_DEPENDENCIES=$wiggle_26_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.ETOH_shMED1_1.reverse.a63258d9c84a165392584991e2febab2.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.ETOH_shMED1_1.reverse.a63258d9c84a165392584991e2febab2.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/bedtools/2.25.0 mugqic/ucsc/v326 && \
mkdir -p tracks/ETOH_shMED1_1 tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 97  alignment/ETOH_shMED1_1/ETOH_shMED1_1.sorted.mdup.reverse.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/ETOH_shMED1_1/ETOH_shMED1_1.sorted.mdup.reverse.bam \
  -g /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  > tracks/ETOH_shMED1_1/ETOH_shMED1_1.reverse.bedGraph && \
sort -k1,1 -k2,2n tracks/ETOH_shMED1_1/ETOH_shMED1_1.reverse.bedGraph > tracks/ETOH_shMED1_1/ETOH_shMED1_1.reverse.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/ETOH_shMED1_1/ETOH_shMED1_1.reverse.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  tracks/bigWig/ETOH_shMED1_1.reverse.bw
wiggle.ETOH_shMED1_1.reverse.a63258d9c84a165392584991e2febab2.mugqic.done
)
wiggle_28_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_28_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_29_JOB_ID: wiggle.ETOH_shNIPBL_2.forward_strandspec
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.ETOH_shNIPBL_2.forward_strandspec
JOB_DEPENDENCIES=$picard_mark_duplicates_8_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.ETOH_shNIPBL_2.forward_strandspec.79771964b5e12feae1b94b0245615780.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.ETOH_shNIPBL_2.forward_strandspec.79771964b5e12feae1b94b0245615780.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
samtools view -bh -F 256 -f 81 \
  alignment/ETOH_shNIPBL_2/ETOH_shNIPBL_2.sorted.mdup.bam \
  > alignment/ETOH_shNIPBL_2/ETOH_shNIPBL_2.sorted.mdup.tmp1.forward.bam && \
samtools view -bh -F 256 -f 161 \
  alignment/ETOH_shNIPBL_2/ETOH_shNIPBL_2.sorted.mdup.bam \
  > alignment/ETOH_shNIPBL_2/ETOH_shNIPBL_2.sorted.mdup.tmp2.forward.bam && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/picard.jar MergeSamFiles \
 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/ETOH_shNIPBL_2/ETOH_shNIPBL_2.sorted.mdup.tmp1.forward.bam \
  INPUT=alignment/ETOH_shNIPBL_2/ETOH_shNIPBL_2.sorted.mdup.tmp2.forward.bam \
 OUTPUT=alignment/ETOH_shNIPBL_2/ETOH_shNIPBL_2.sorted.mdup.forward.bam \
 MAX_RECORDS_IN_RAM=5750000 && \
rm alignment/ETOH_shNIPBL_2/ETOH_shNIPBL_2.sorted.mdup.tmp1.forward.bam alignment/ETOH_shNIPBL_2/ETOH_shNIPBL_2.sorted.mdup.tmp2.forward.bam
wiggle.ETOH_shNIPBL_2.forward_strandspec.79771964b5e12feae1b94b0245615780.mugqic.done
)
wiggle_29_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_29_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_30_JOB_ID: wiggle.ETOH_shNIPBL_2.reverse_strandspec
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.ETOH_shNIPBL_2.reverse_strandspec
JOB_DEPENDENCIES=$picard_mark_duplicates_8_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.ETOH_shNIPBL_2.reverse_strandspec.6cf4bbc6f696bd00de74656006142cb6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.ETOH_shNIPBL_2.reverse_strandspec.6cf4bbc6f696bd00de74656006142cb6.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p tracks/ETOH_shNIPBL_2 tracks/bigWig && \
samtools view -bh -F 256 -f 97 \
  alignment/ETOH_shNIPBL_2/ETOH_shNIPBL_2.sorted.mdup.bam \
  > alignment/ETOH_shNIPBL_2/ETOH_shNIPBL_2.sorted.mdup.tmp1.reverse.bam && \
samtools view -bh -F 256 -f 145 \
  alignment/ETOH_shNIPBL_2/ETOH_shNIPBL_2.sorted.mdup.bam \
  > alignment/ETOH_shNIPBL_2/ETOH_shNIPBL_2.sorted.mdup.tmp2.reverse.bam && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/picard.jar MergeSamFiles \
 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/ETOH_shNIPBL_2/ETOH_shNIPBL_2.sorted.mdup.tmp1.reverse.bam \
  INPUT=alignment/ETOH_shNIPBL_2/ETOH_shNIPBL_2.sorted.mdup.tmp2.reverse.bam \
 OUTPUT=alignment/ETOH_shNIPBL_2/ETOH_shNIPBL_2.sorted.mdup.reverse.bam \
 MAX_RECORDS_IN_RAM=5750000 && \
rm alignment/ETOH_shNIPBL_2/ETOH_shNIPBL_2.sorted.mdup.tmp1.reverse.bam alignment/ETOH_shNIPBL_2/ETOH_shNIPBL_2.sorted.mdup.tmp2.reverse.bam
wiggle.ETOH_shNIPBL_2.reverse_strandspec.6cf4bbc6f696bd00de74656006142cb6.mugqic.done
)
wiggle_30_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_30_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_31_JOB_ID: wiggle.ETOH_shNIPBL_2.forward
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.ETOH_shNIPBL_2.forward
JOB_DEPENDENCIES=$wiggle_29_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.ETOH_shNIPBL_2.forward.18113da464a3d4afde761a4a61c61a70.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.ETOH_shNIPBL_2.forward.18113da464a3d4afde761a4a61c61a70.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/bedtools/2.25.0 mugqic/ucsc/v326 && \
mkdir -p tracks/ETOH_shNIPBL_2 tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 81  alignment/ETOH_shNIPBL_2/ETOH_shNIPBL_2.sorted.mdup.forward.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/ETOH_shNIPBL_2/ETOH_shNIPBL_2.sorted.mdup.forward.bam \
  -g /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  > tracks/ETOH_shNIPBL_2/ETOH_shNIPBL_2.forward.bedGraph && \
sort -k1,1 -k2,2n tracks/ETOH_shNIPBL_2/ETOH_shNIPBL_2.forward.bedGraph > tracks/ETOH_shNIPBL_2/ETOH_shNIPBL_2.forward.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/ETOH_shNIPBL_2/ETOH_shNIPBL_2.forward.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  tracks/bigWig/ETOH_shNIPBL_2.forward.bw
wiggle.ETOH_shNIPBL_2.forward.18113da464a3d4afde761a4a61c61a70.mugqic.done
)
wiggle_31_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_31_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_32_JOB_ID: wiggle.ETOH_shNIPBL_2.reverse
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.ETOH_shNIPBL_2.reverse
JOB_DEPENDENCIES=$wiggle_30_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.ETOH_shNIPBL_2.reverse.8204caba8c2cec51821db1cc2e4085c0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.ETOH_shNIPBL_2.reverse.8204caba8c2cec51821db1cc2e4085c0.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/bedtools/2.25.0 mugqic/ucsc/v326 && \
mkdir -p tracks/ETOH_shNIPBL_2 tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 97  alignment/ETOH_shNIPBL_2/ETOH_shNIPBL_2.sorted.mdup.reverse.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/ETOH_shNIPBL_2/ETOH_shNIPBL_2.sorted.mdup.reverse.bam \
  -g /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  > tracks/ETOH_shNIPBL_2/ETOH_shNIPBL_2.reverse.bedGraph && \
sort -k1,1 -k2,2n tracks/ETOH_shNIPBL_2/ETOH_shNIPBL_2.reverse.bedGraph > tracks/ETOH_shNIPBL_2/ETOH_shNIPBL_2.reverse.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/ETOH_shNIPBL_2/ETOH_shNIPBL_2.reverse.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  tracks/bigWig/ETOH_shNIPBL_2.reverse.bw
wiggle.ETOH_shNIPBL_2.reverse.8204caba8c2cec51821db1cc2e4085c0.mugqic.done
)
wiggle_32_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_32_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_33_JOB_ID: wiggle.ETOH_shNIPBL_1.forward_strandspec
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.ETOH_shNIPBL_1.forward_strandspec
JOB_DEPENDENCIES=$picard_mark_duplicates_9_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.ETOH_shNIPBL_1.forward_strandspec.2cb17251e0122be1963b847842e4e499.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.ETOH_shNIPBL_1.forward_strandspec.2cb17251e0122be1963b847842e4e499.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
samtools view -bh -F 256 -f 81 \
  alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1.sorted.mdup.bam \
  > alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1.sorted.mdup.tmp1.forward.bam && \
samtools view -bh -F 256 -f 161 \
  alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1.sorted.mdup.bam \
  > alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1.sorted.mdup.tmp2.forward.bam && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/picard.jar MergeSamFiles \
 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1.sorted.mdup.tmp1.forward.bam \
  INPUT=alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1.sorted.mdup.tmp2.forward.bam \
 OUTPUT=alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1.sorted.mdup.forward.bam \
 MAX_RECORDS_IN_RAM=5750000 && \
rm alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1.sorted.mdup.tmp1.forward.bam alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1.sorted.mdup.tmp2.forward.bam
wiggle.ETOH_shNIPBL_1.forward_strandspec.2cb17251e0122be1963b847842e4e499.mugqic.done
)
wiggle_33_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_33_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_34_JOB_ID: wiggle.ETOH_shNIPBL_1.reverse_strandspec
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.ETOH_shNIPBL_1.reverse_strandspec
JOB_DEPENDENCIES=$picard_mark_duplicates_9_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.ETOH_shNIPBL_1.reverse_strandspec.aa64f17e13b705e43b1ecd55e43a3998.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.ETOH_shNIPBL_1.reverse_strandspec.aa64f17e13b705e43b1ecd55e43a3998.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p tracks/ETOH_shNIPBL_1 tracks/bigWig && \
samtools view -bh -F 256 -f 97 \
  alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1.sorted.mdup.bam \
  > alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1.sorted.mdup.tmp1.reverse.bam && \
samtools view -bh -F 256 -f 145 \
  alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1.sorted.mdup.bam \
  > alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1.sorted.mdup.tmp2.reverse.bam && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/picard.jar MergeSamFiles \
 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1.sorted.mdup.tmp1.reverse.bam \
  INPUT=alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1.sorted.mdup.tmp2.reverse.bam \
 OUTPUT=alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1.sorted.mdup.reverse.bam \
 MAX_RECORDS_IN_RAM=5750000 && \
rm alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1.sorted.mdup.tmp1.reverse.bam alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1.sorted.mdup.tmp2.reverse.bam
wiggle.ETOH_shNIPBL_1.reverse_strandspec.aa64f17e13b705e43b1ecd55e43a3998.mugqic.done
)
wiggle_34_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_34_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_35_JOB_ID: wiggle.ETOH_shNIPBL_1.forward
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.ETOH_shNIPBL_1.forward
JOB_DEPENDENCIES=$wiggle_33_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.ETOH_shNIPBL_1.forward.fd8033acfd78a4742a258c8a7accbb12.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.ETOH_shNIPBL_1.forward.fd8033acfd78a4742a258c8a7accbb12.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/bedtools/2.25.0 mugqic/ucsc/v326 && \
mkdir -p tracks/ETOH_shNIPBL_1 tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 81  alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1.sorted.mdup.forward.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1.sorted.mdup.forward.bam \
  -g /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  > tracks/ETOH_shNIPBL_1/ETOH_shNIPBL_1.forward.bedGraph && \
sort -k1,1 -k2,2n tracks/ETOH_shNIPBL_1/ETOH_shNIPBL_1.forward.bedGraph > tracks/ETOH_shNIPBL_1/ETOH_shNIPBL_1.forward.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/ETOH_shNIPBL_1/ETOH_shNIPBL_1.forward.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  tracks/bigWig/ETOH_shNIPBL_1.forward.bw
wiggle.ETOH_shNIPBL_1.forward.fd8033acfd78a4742a258c8a7accbb12.mugqic.done
)
wiggle_35_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_35_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_36_JOB_ID: wiggle.ETOH_shNIPBL_1.reverse
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.ETOH_shNIPBL_1.reverse
JOB_DEPENDENCIES=$wiggle_34_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.ETOH_shNIPBL_1.reverse.bc683a6de82497ed93d598622616f51f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.ETOH_shNIPBL_1.reverse.bc683a6de82497ed93d598622616f51f.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/bedtools/2.25.0 mugqic/ucsc/v326 && \
mkdir -p tracks/ETOH_shNIPBL_1 tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 97  alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1.sorted.mdup.reverse.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1.sorted.mdup.reverse.bam \
  -g /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  > tracks/ETOH_shNIPBL_1/ETOH_shNIPBL_1.reverse.bedGraph && \
sort -k1,1 -k2,2n tracks/ETOH_shNIPBL_1/ETOH_shNIPBL_1.reverse.bedGraph > tracks/ETOH_shNIPBL_1/ETOH_shNIPBL_1.reverse.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/ETOH_shNIPBL_1/ETOH_shNIPBL_1.reverse.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  tracks/bigWig/ETOH_shNIPBL_1.reverse.bw
wiggle.ETOH_shNIPBL_1.reverse.bc683a6de82497ed93d598622616f51f.mugqic.done
)
wiggle_36_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_36_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_37_JOB_ID: wiggle.ETOH_shMED1_2.forward_strandspec
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.ETOH_shMED1_2.forward_strandspec
JOB_DEPENDENCIES=$picard_mark_duplicates_10_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.ETOH_shMED1_2.forward_strandspec.ca26139895c1fa14011becfa2ccecb2d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.ETOH_shMED1_2.forward_strandspec.ca26139895c1fa14011becfa2ccecb2d.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
samtools view -bh -F 256 -f 81 \
  alignment/ETOH_shMED1_2/ETOH_shMED1_2.sorted.mdup.bam \
  > alignment/ETOH_shMED1_2/ETOH_shMED1_2.sorted.mdup.tmp1.forward.bam && \
samtools view -bh -F 256 -f 161 \
  alignment/ETOH_shMED1_2/ETOH_shMED1_2.sorted.mdup.bam \
  > alignment/ETOH_shMED1_2/ETOH_shMED1_2.sorted.mdup.tmp2.forward.bam && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/picard.jar MergeSamFiles \
 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/ETOH_shMED1_2/ETOH_shMED1_2.sorted.mdup.tmp1.forward.bam \
  INPUT=alignment/ETOH_shMED1_2/ETOH_shMED1_2.sorted.mdup.tmp2.forward.bam \
 OUTPUT=alignment/ETOH_shMED1_2/ETOH_shMED1_2.sorted.mdup.forward.bam \
 MAX_RECORDS_IN_RAM=5750000 && \
rm alignment/ETOH_shMED1_2/ETOH_shMED1_2.sorted.mdup.tmp1.forward.bam alignment/ETOH_shMED1_2/ETOH_shMED1_2.sorted.mdup.tmp2.forward.bam
wiggle.ETOH_shMED1_2.forward_strandspec.ca26139895c1fa14011becfa2ccecb2d.mugqic.done
)
wiggle_37_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_37_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_38_JOB_ID: wiggle.ETOH_shMED1_2.reverse_strandspec
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.ETOH_shMED1_2.reverse_strandspec
JOB_DEPENDENCIES=$picard_mark_duplicates_10_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.ETOH_shMED1_2.reverse_strandspec.f203ee45f2eca684fc312abe3f95ed02.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.ETOH_shMED1_2.reverse_strandspec.f203ee45f2eca684fc312abe3f95ed02.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p tracks/ETOH_shMED1_2 tracks/bigWig && \
samtools view -bh -F 256 -f 97 \
  alignment/ETOH_shMED1_2/ETOH_shMED1_2.sorted.mdup.bam \
  > alignment/ETOH_shMED1_2/ETOH_shMED1_2.sorted.mdup.tmp1.reverse.bam && \
samtools view -bh -F 256 -f 145 \
  alignment/ETOH_shMED1_2/ETOH_shMED1_2.sorted.mdup.bam \
  > alignment/ETOH_shMED1_2/ETOH_shMED1_2.sorted.mdup.tmp2.reverse.bam && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/picard.jar MergeSamFiles \
 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/ETOH_shMED1_2/ETOH_shMED1_2.sorted.mdup.tmp1.reverse.bam \
  INPUT=alignment/ETOH_shMED1_2/ETOH_shMED1_2.sorted.mdup.tmp2.reverse.bam \
 OUTPUT=alignment/ETOH_shMED1_2/ETOH_shMED1_2.sorted.mdup.reverse.bam \
 MAX_RECORDS_IN_RAM=5750000 && \
rm alignment/ETOH_shMED1_2/ETOH_shMED1_2.sorted.mdup.tmp1.reverse.bam alignment/ETOH_shMED1_2/ETOH_shMED1_2.sorted.mdup.tmp2.reverse.bam
wiggle.ETOH_shMED1_2.reverse_strandspec.f203ee45f2eca684fc312abe3f95ed02.mugqic.done
)
wiggle_38_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_38_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_39_JOB_ID: wiggle.ETOH_shMED1_2.forward
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.ETOH_shMED1_2.forward
JOB_DEPENDENCIES=$wiggle_37_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.ETOH_shMED1_2.forward.7202ba429a9572196211dbca34384010.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.ETOH_shMED1_2.forward.7202ba429a9572196211dbca34384010.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/bedtools/2.25.0 mugqic/ucsc/v326 && \
mkdir -p tracks/ETOH_shMED1_2 tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 81  alignment/ETOH_shMED1_2/ETOH_shMED1_2.sorted.mdup.forward.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/ETOH_shMED1_2/ETOH_shMED1_2.sorted.mdup.forward.bam \
  -g /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  > tracks/ETOH_shMED1_2/ETOH_shMED1_2.forward.bedGraph && \
sort -k1,1 -k2,2n tracks/ETOH_shMED1_2/ETOH_shMED1_2.forward.bedGraph > tracks/ETOH_shMED1_2/ETOH_shMED1_2.forward.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/ETOH_shMED1_2/ETOH_shMED1_2.forward.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  tracks/bigWig/ETOH_shMED1_2.forward.bw
wiggle.ETOH_shMED1_2.forward.7202ba429a9572196211dbca34384010.mugqic.done
)
wiggle_39_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_39_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_40_JOB_ID: wiggle.ETOH_shMED1_2.reverse
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.ETOH_shMED1_2.reverse
JOB_DEPENDENCIES=$wiggle_38_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.ETOH_shMED1_2.reverse.56ccade4cb49ff4855a6424c12102c24.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.ETOH_shMED1_2.reverse.56ccade4cb49ff4855a6424c12102c24.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/bedtools/2.25.0 mugqic/ucsc/v326 && \
mkdir -p tracks/ETOH_shMED1_2 tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 97  alignment/ETOH_shMED1_2/ETOH_shMED1_2.sorted.mdup.reverse.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/ETOH_shMED1_2/ETOH_shMED1_2.sorted.mdup.reverse.bam \
  -g /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  > tracks/ETOH_shMED1_2/ETOH_shMED1_2.reverse.bedGraph && \
sort -k1,1 -k2,2n tracks/ETOH_shMED1_2/ETOH_shMED1_2.reverse.bedGraph > tracks/ETOH_shMED1_2/ETOH_shMED1_2.reverse.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/ETOH_shMED1_2/ETOH_shMED1_2.reverse.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  tracks/bigWig/ETOH_shMED1_2.reverse.bw
wiggle.ETOH_shMED1_2.reverse.56ccade4cb49ff4855a6424c12102c24.mugqic.done
)
wiggle_40_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_40_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_41_JOB_ID: wiggle.ETOH_shSMC1A_3.forward_strandspec
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.ETOH_shSMC1A_3.forward_strandspec
JOB_DEPENDENCIES=$picard_mark_duplicates_11_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.ETOH_shSMC1A_3.forward_strandspec.c571bed62c5390a7ab36ccc7684ac07d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.ETOH_shSMC1A_3.forward_strandspec.c571bed62c5390a7ab36ccc7684ac07d.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
samtools view -bh -F 256 -f 81 \
  alignment/ETOH_shSMC1A_3/ETOH_shSMC1A_3.sorted.mdup.bam \
  > alignment/ETOH_shSMC1A_3/ETOH_shSMC1A_3.sorted.mdup.tmp1.forward.bam && \
samtools view -bh -F 256 -f 161 \
  alignment/ETOH_shSMC1A_3/ETOH_shSMC1A_3.sorted.mdup.bam \
  > alignment/ETOH_shSMC1A_3/ETOH_shSMC1A_3.sorted.mdup.tmp2.forward.bam && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/picard.jar MergeSamFiles \
 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/ETOH_shSMC1A_3/ETOH_shSMC1A_3.sorted.mdup.tmp1.forward.bam \
  INPUT=alignment/ETOH_shSMC1A_3/ETOH_shSMC1A_3.sorted.mdup.tmp2.forward.bam \
 OUTPUT=alignment/ETOH_shSMC1A_3/ETOH_shSMC1A_3.sorted.mdup.forward.bam \
 MAX_RECORDS_IN_RAM=5750000 && \
rm alignment/ETOH_shSMC1A_3/ETOH_shSMC1A_3.sorted.mdup.tmp1.forward.bam alignment/ETOH_shSMC1A_3/ETOH_shSMC1A_3.sorted.mdup.tmp2.forward.bam
wiggle.ETOH_shSMC1A_3.forward_strandspec.c571bed62c5390a7ab36ccc7684ac07d.mugqic.done
)
wiggle_41_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_41_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_42_JOB_ID: wiggle.ETOH_shSMC1A_3.reverse_strandspec
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.ETOH_shSMC1A_3.reverse_strandspec
JOB_DEPENDENCIES=$picard_mark_duplicates_11_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.ETOH_shSMC1A_3.reverse_strandspec.9dabbc0ac2e27f9a164e09760f316330.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.ETOH_shSMC1A_3.reverse_strandspec.9dabbc0ac2e27f9a164e09760f316330.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p tracks/ETOH_shSMC1A_3 tracks/bigWig && \
samtools view -bh -F 256 -f 97 \
  alignment/ETOH_shSMC1A_3/ETOH_shSMC1A_3.sorted.mdup.bam \
  > alignment/ETOH_shSMC1A_3/ETOH_shSMC1A_3.sorted.mdup.tmp1.reverse.bam && \
samtools view -bh -F 256 -f 145 \
  alignment/ETOH_shSMC1A_3/ETOH_shSMC1A_3.sorted.mdup.bam \
  > alignment/ETOH_shSMC1A_3/ETOH_shSMC1A_3.sorted.mdup.tmp2.reverse.bam && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/picard.jar MergeSamFiles \
 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/ETOH_shSMC1A_3/ETOH_shSMC1A_3.sorted.mdup.tmp1.reverse.bam \
  INPUT=alignment/ETOH_shSMC1A_3/ETOH_shSMC1A_3.sorted.mdup.tmp2.reverse.bam \
 OUTPUT=alignment/ETOH_shSMC1A_3/ETOH_shSMC1A_3.sorted.mdup.reverse.bam \
 MAX_RECORDS_IN_RAM=5750000 && \
rm alignment/ETOH_shSMC1A_3/ETOH_shSMC1A_3.sorted.mdup.tmp1.reverse.bam alignment/ETOH_shSMC1A_3/ETOH_shSMC1A_3.sorted.mdup.tmp2.reverse.bam
wiggle.ETOH_shSMC1A_3.reverse_strandspec.9dabbc0ac2e27f9a164e09760f316330.mugqic.done
)
wiggle_42_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_42_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_43_JOB_ID: wiggle.ETOH_shSMC1A_3.forward
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.ETOH_shSMC1A_3.forward
JOB_DEPENDENCIES=$wiggle_41_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.ETOH_shSMC1A_3.forward.be7ef23716e9b33e8b74c7ec38f0a96d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.ETOH_shSMC1A_3.forward.be7ef23716e9b33e8b74c7ec38f0a96d.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/bedtools/2.25.0 mugqic/ucsc/v326 && \
mkdir -p tracks/ETOH_shSMC1A_3 tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 81  alignment/ETOH_shSMC1A_3/ETOH_shSMC1A_3.sorted.mdup.forward.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/ETOH_shSMC1A_3/ETOH_shSMC1A_3.sorted.mdup.forward.bam \
  -g /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  > tracks/ETOH_shSMC1A_3/ETOH_shSMC1A_3.forward.bedGraph && \
sort -k1,1 -k2,2n tracks/ETOH_shSMC1A_3/ETOH_shSMC1A_3.forward.bedGraph > tracks/ETOH_shSMC1A_3/ETOH_shSMC1A_3.forward.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/ETOH_shSMC1A_3/ETOH_shSMC1A_3.forward.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  tracks/bigWig/ETOH_shSMC1A_3.forward.bw
wiggle.ETOH_shSMC1A_3.forward.be7ef23716e9b33e8b74c7ec38f0a96d.mugqic.done
)
wiggle_43_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_43_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_44_JOB_ID: wiggle.ETOH_shSMC1A_3.reverse
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.ETOH_shSMC1A_3.reverse
JOB_DEPENDENCIES=$wiggle_42_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.ETOH_shSMC1A_3.reverse.d92bd5db471dd0c799bb041834089add.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.ETOH_shSMC1A_3.reverse.d92bd5db471dd0c799bb041834089add.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/bedtools/2.25.0 mugqic/ucsc/v326 && \
mkdir -p tracks/ETOH_shSMC1A_3 tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 97  alignment/ETOH_shSMC1A_3/ETOH_shSMC1A_3.sorted.mdup.reverse.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/ETOH_shSMC1A_3/ETOH_shSMC1A_3.sorted.mdup.reverse.bam \
  -g /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  > tracks/ETOH_shSMC1A_3/ETOH_shSMC1A_3.reverse.bedGraph && \
sort -k1,1 -k2,2n tracks/ETOH_shSMC1A_3/ETOH_shSMC1A_3.reverse.bedGraph > tracks/ETOH_shSMC1A_3/ETOH_shSMC1A_3.reverse.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/ETOH_shSMC1A_3/ETOH_shSMC1A_3.reverse.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  tracks/bigWig/ETOH_shSMC1A_3.reverse.bw
wiggle.ETOH_shSMC1A_3.reverse.d92bd5db471dd0c799bb041834089add.mugqic.done
)
wiggle_44_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_44_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_45_JOB_ID: wiggle.ETOH_shSMC1A_2.forward_strandspec
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.ETOH_shSMC1A_2.forward_strandspec
JOB_DEPENDENCIES=$picard_mark_duplicates_12_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.ETOH_shSMC1A_2.forward_strandspec.ffaffe0bb0a549d033b5559bc5861232.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.ETOH_shSMC1A_2.forward_strandspec.ffaffe0bb0a549d033b5559bc5861232.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
samtools view -bh -F 256 -f 81 \
  alignment/ETOH_shSMC1A_2/ETOH_shSMC1A_2.sorted.mdup.bam \
  > alignment/ETOH_shSMC1A_2/ETOH_shSMC1A_2.sorted.mdup.tmp1.forward.bam && \
samtools view -bh -F 256 -f 161 \
  alignment/ETOH_shSMC1A_2/ETOH_shSMC1A_2.sorted.mdup.bam \
  > alignment/ETOH_shSMC1A_2/ETOH_shSMC1A_2.sorted.mdup.tmp2.forward.bam && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/picard.jar MergeSamFiles \
 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/ETOH_shSMC1A_2/ETOH_shSMC1A_2.sorted.mdup.tmp1.forward.bam \
  INPUT=alignment/ETOH_shSMC1A_2/ETOH_shSMC1A_2.sorted.mdup.tmp2.forward.bam \
 OUTPUT=alignment/ETOH_shSMC1A_2/ETOH_shSMC1A_2.sorted.mdup.forward.bam \
 MAX_RECORDS_IN_RAM=5750000 && \
rm alignment/ETOH_shSMC1A_2/ETOH_shSMC1A_2.sorted.mdup.tmp1.forward.bam alignment/ETOH_shSMC1A_2/ETOH_shSMC1A_2.sorted.mdup.tmp2.forward.bam
wiggle.ETOH_shSMC1A_2.forward_strandspec.ffaffe0bb0a549d033b5559bc5861232.mugqic.done
)
wiggle_45_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_45_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_46_JOB_ID: wiggle.ETOH_shSMC1A_2.reverse_strandspec
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.ETOH_shSMC1A_2.reverse_strandspec
JOB_DEPENDENCIES=$picard_mark_duplicates_12_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.ETOH_shSMC1A_2.reverse_strandspec.c2712773d44eae76fa881bbadce25704.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.ETOH_shSMC1A_2.reverse_strandspec.c2712773d44eae76fa881bbadce25704.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p tracks/ETOH_shSMC1A_2 tracks/bigWig && \
samtools view -bh -F 256 -f 97 \
  alignment/ETOH_shSMC1A_2/ETOH_shSMC1A_2.sorted.mdup.bam \
  > alignment/ETOH_shSMC1A_2/ETOH_shSMC1A_2.sorted.mdup.tmp1.reverse.bam && \
samtools view -bh -F 256 -f 145 \
  alignment/ETOH_shSMC1A_2/ETOH_shSMC1A_2.sorted.mdup.bam \
  > alignment/ETOH_shSMC1A_2/ETOH_shSMC1A_2.sorted.mdup.tmp2.reverse.bam && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/picard.jar MergeSamFiles \
 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/ETOH_shSMC1A_2/ETOH_shSMC1A_2.sorted.mdup.tmp1.reverse.bam \
  INPUT=alignment/ETOH_shSMC1A_2/ETOH_shSMC1A_2.sorted.mdup.tmp2.reverse.bam \
 OUTPUT=alignment/ETOH_shSMC1A_2/ETOH_shSMC1A_2.sorted.mdup.reverse.bam \
 MAX_RECORDS_IN_RAM=5750000 && \
rm alignment/ETOH_shSMC1A_2/ETOH_shSMC1A_2.sorted.mdup.tmp1.reverse.bam alignment/ETOH_shSMC1A_2/ETOH_shSMC1A_2.sorted.mdup.tmp2.reverse.bam
wiggle.ETOH_shSMC1A_2.reverse_strandspec.c2712773d44eae76fa881bbadce25704.mugqic.done
)
wiggle_46_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_46_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_47_JOB_ID: wiggle.ETOH_shSMC1A_2.forward
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.ETOH_shSMC1A_2.forward
JOB_DEPENDENCIES=$wiggle_45_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.ETOH_shSMC1A_2.forward.3b73ad12341f4cf3007b585096d60fe4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.ETOH_shSMC1A_2.forward.3b73ad12341f4cf3007b585096d60fe4.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/bedtools/2.25.0 mugqic/ucsc/v326 && \
mkdir -p tracks/ETOH_shSMC1A_2 tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 81  alignment/ETOH_shSMC1A_2/ETOH_shSMC1A_2.sorted.mdup.forward.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/ETOH_shSMC1A_2/ETOH_shSMC1A_2.sorted.mdup.forward.bam \
  -g /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  > tracks/ETOH_shSMC1A_2/ETOH_shSMC1A_2.forward.bedGraph && \
sort -k1,1 -k2,2n tracks/ETOH_shSMC1A_2/ETOH_shSMC1A_2.forward.bedGraph > tracks/ETOH_shSMC1A_2/ETOH_shSMC1A_2.forward.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/ETOH_shSMC1A_2/ETOH_shSMC1A_2.forward.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  tracks/bigWig/ETOH_shSMC1A_2.forward.bw
wiggle.ETOH_shSMC1A_2.forward.3b73ad12341f4cf3007b585096d60fe4.mugqic.done
)
wiggle_47_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_47_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_48_JOB_ID: wiggle.ETOH_shSMC1A_2.reverse
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.ETOH_shSMC1A_2.reverse
JOB_DEPENDENCIES=$wiggle_46_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.ETOH_shSMC1A_2.reverse.872e5491eddde51f7601f61d12e9a3e6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.ETOH_shSMC1A_2.reverse.872e5491eddde51f7601f61d12e9a3e6.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/bedtools/2.25.0 mugqic/ucsc/v326 && \
mkdir -p tracks/ETOH_shSMC1A_2 tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 97  alignment/ETOH_shSMC1A_2/ETOH_shSMC1A_2.sorted.mdup.reverse.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/ETOH_shSMC1A_2/ETOH_shSMC1A_2.sorted.mdup.reverse.bam \
  -g /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  > tracks/ETOH_shSMC1A_2/ETOH_shSMC1A_2.reverse.bedGraph && \
sort -k1,1 -k2,2n tracks/ETOH_shSMC1A_2/ETOH_shSMC1A_2.reverse.bedGraph > tracks/ETOH_shSMC1A_2/ETOH_shSMC1A_2.reverse.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/ETOH_shSMC1A_2/ETOH_shSMC1A_2.reverse.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  tracks/bigWig/ETOH_shSMC1A_2.reverse.bw
wiggle.ETOH_shSMC1A_2.reverse.872e5491eddde51f7601f61d12e9a3e6.mugqic.done
)
wiggle_48_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_48_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: raw_counts
#-------------------------------------------------------------------------------
STEP=raw_counts
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: raw_counts_1_JOB_ID: htseq_count.DEX_nonMamm_1
#-------------------------------------------------------------------------------
JOB_NAME=htseq_count.DEX_nonMamm_1
JOB_DEPENDENCIES=$picard_sort_sam_1_JOB_ID
JOB_DONE=job_output/raw_counts/htseq_count.DEX_nonMamm_1.ac90ae914dc076f1901eb2dc3b8a5826.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'htseq_count.DEX_nonMamm_1.ac90ae914dc076f1901eb2dc3b8a5826.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/python/2.7.13 && \
mkdir -p raw_counts && \
samtools view -F 4 \
  alignment/DEX_nonMamm_1/DEX_nonMamm_1.QueryNameSorted.bam | \
htseq-count - \
  -m intersection-nonempty \
  --stranded=reverse \
  --format=sam \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  > raw_counts/DEX_nonMamm_1.readcounts.csv
htseq_count.DEX_nonMamm_1.ac90ae914dc076f1901eb2dc3b8a5826.mugqic.done
)
raw_counts_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$raw_counts_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: raw_counts_2_JOB_ID: htseq_count.DEX_shMED1_1
#-------------------------------------------------------------------------------
JOB_NAME=htseq_count.DEX_shMED1_1
JOB_DEPENDENCIES=$picard_sort_sam_2_JOB_ID
JOB_DONE=job_output/raw_counts/htseq_count.DEX_shMED1_1.1fc9349fcb691fb7fe95dfae80c440cb.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'htseq_count.DEX_shMED1_1.1fc9349fcb691fb7fe95dfae80c440cb.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/python/2.7.13 && \
mkdir -p raw_counts && \
samtools view -F 4 \
  alignment/DEX_shMED1_1/DEX_shMED1_1.QueryNameSorted.bam | \
htseq-count - \
  -m intersection-nonempty \
  --stranded=reverse \
  --format=sam \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  > raw_counts/DEX_shMED1_1.readcounts.csv
htseq_count.DEX_shMED1_1.1fc9349fcb691fb7fe95dfae80c440cb.mugqic.done
)
raw_counts_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$raw_counts_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: raw_counts_3_JOB_ID: htseq_count.DEX_shNIPBL_1
#-------------------------------------------------------------------------------
JOB_NAME=htseq_count.DEX_shNIPBL_1
JOB_DEPENDENCIES=$picard_sort_sam_3_JOB_ID
JOB_DONE=job_output/raw_counts/htseq_count.DEX_shNIPBL_1.e2a573b5338903d695dcfbd704be057c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'htseq_count.DEX_shNIPBL_1.e2a573b5338903d695dcfbd704be057c.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/python/2.7.13 && \
mkdir -p raw_counts && \
samtools view -F 4 \
  alignment/DEX_shNIPBL_1/DEX_shNIPBL_1.QueryNameSorted.bam | \
htseq-count - \
  -m intersection-nonempty \
  --stranded=reverse \
  --format=sam \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  > raw_counts/DEX_shNIPBL_1.readcounts.csv
htseq_count.DEX_shNIPBL_1.e2a573b5338903d695dcfbd704be057c.mugqic.done
)
raw_counts_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$raw_counts_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: raw_counts_4_JOB_ID: htseq_count.DEX_shSMC1A_3
#-------------------------------------------------------------------------------
JOB_NAME=htseq_count.DEX_shSMC1A_3
JOB_DEPENDENCIES=$picard_sort_sam_4_JOB_ID
JOB_DONE=job_output/raw_counts/htseq_count.DEX_shSMC1A_3.ddd82f2e3a62e737f2b61fce52ba8daf.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'htseq_count.DEX_shSMC1A_3.ddd82f2e3a62e737f2b61fce52ba8daf.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/python/2.7.13 && \
mkdir -p raw_counts && \
samtools view -F 4 \
  alignment/DEX_shSMC1A_3/DEX_shSMC1A_3.QueryNameSorted.bam | \
htseq-count - \
  -m intersection-nonempty \
  --stranded=reverse \
  --format=sam \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  > raw_counts/DEX_shSMC1A_3.readcounts.csv
htseq_count.DEX_shSMC1A_3.ddd82f2e3a62e737f2b61fce52ba8daf.mugqic.done
)
raw_counts_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$raw_counts_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: raw_counts_5_JOB_ID: htseq_count.ETOH_nonMamm_1
#-------------------------------------------------------------------------------
JOB_NAME=htseq_count.ETOH_nonMamm_1
JOB_DEPENDENCIES=$picard_sort_sam_5_JOB_ID
JOB_DONE=job_output/raw_counts/htseq_count.ETOH_nonMamm_1.91e4ec2d6216095106e2c9dd9075ef4b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'htseq_count.ETOH_nonMamm_1.91e4ec2d6216095106e2c9dd9075ef4b.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/python/2.7.13 && \
mkdir -p raw_counts && \
samtools view -F 4 \
  alignment/ETOH_nonMamm_1/ETOH_nonMamm_1.QueryNameSorted.bam | \
htseq-count - \
  -m intersection-nonempty \
  --stranded=reverse \
  --format=sam \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  > raw_counts/ETOH_nonMamm_1.readcounts.csv
htseq_count.ETOH_nonMamm_1.91e4ec2d6216095106e2c9dd9075ef4b.mugqic.done
)
raw_counts_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$raw_counts_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: raw_counts_6_JOB_ID: htseq_count.ETOH_nonMamm_2
#-------------------------------------------------------------------------------
JOB_NAME=htseq_count.ETOH_nonMamm_2
JOB_DEPENDENCIES=$picard_sort_sam_6_JOB_ID
JOB_DONE=job_output/raw_counts/htseq_count.ETOH_nonMamm_2.86661f838a28211cd19ef6647176e7af.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'htseq_count.ETOH_nonMamm_2.86661f838a28211cd19ef6647176e7af.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/python/2.7.13 && \
mkdir -p raw_counts && \
samtools view -F 4 \
  alignment/ETOH_nonMamm_2/ETOH_nonMamm_2.QueryNameSorted.bam | \
htseq-count - \
  -m intersection-nonempty \
  --stranded=reverse \
  --format=sam \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  > raw_counts/ETOH_nonMamm_2.readcounts.csv
htseq_count.ETOH_nonMamm_2.86661f838a28211cd19ef6647176e7af.mugqic.done
)
raw_counts_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$raw_counts_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: raw_counts_7_JOB_ID: htseq_count.ETOH_shMED1_1
#-------------------------------------------------------------------------------
JOB_NAME=htseq_count.ETOH_shMED1_1
JOB_DEPENDENCIES=$picard_sort_sam_7_JOB_ID
JOB_DONE=job_output/raw_counts/htseq_count.ETOH_shMED1_1.527736c4bea4038f169f297488da16ec.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'htseq_count.ETOH_shMED1_1.527736c4bea4038f169f297488da16ec.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/python/2.7.13 && \
mkdir -p raw_counts && \
samtools view -F 4 \
  alignment/ETOH_shMED1_1/ETOH_shMED1_1.QueryNameSorted.bam | \
htseq-count - \
  -m intersection-nonempty \
  --stranded=reverse \
  --format=sam \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  > raw_counts/ETOH_shMED1_1.readcounts.csv
htseq_count.ETOH_shMED1_1.527736c4bea4038f169f297488da16ec.mugqic.done
)
raw_counts_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$raw_counts_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: raw_counts_8_JOB_ID: htseq_count.ETOH_shNIPBL_2
#-------------------------------------------------------------------------------
JOB_NAME=htseq_count.ETOH_shNIPBL_2
JOB_DEPENDENCIES=$picard_sort_sam_8_JOB_ID
JOB_DONE=job_output/raw_counts/htseq_count.ETOH_shNIPBL_2.ead30f3e293512fc3cbcd02c5db89f7e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'htseq_count.ETOH_shNIPBL_2.ead30f3e293512fc3cbcd02c5db89f7e.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/python/2.7.13 && \
mkdir -p raw_counts && \
samtools view -F 4 \
  alignment/ETOH_shNIPBL_2/ETOH_shNIPBL_2.QueryNameSorted.bam | \
htseq-count - \
  -m intersection-nonempty \
  --stranded=reverse \
  --format=sam \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  > raw_counts/ETOH_shNIPBL_2.readcounts.csv
htseq_count.ETOH_shNIPBL_2.ead30f3e293512fc3cbcd02c5db89f7e.mugqic.done
)
raw_counts_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$raw_counts_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: raw_counts_9_JOB_ID: htseq_count.ETOH_shNIPBL_1
#-------------------------------------------------------------------------------
JOB_NAME=htseq_count.ETOH_shNIPBL_1
JOB_DEPENDENCIES=$picard_sort_sam_9_JOB_ID
JOB_DONE=job_output/raw_counts/htseq_count.ETOH_shNIPBL_1.06c65df61e0e38abcc01dc89f5908f58.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'htseq_count.ETOH_shNIPBL_1.06c65df61e0e38abcc01dc89f5908f58.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/python/2.7.13 && \
mkdir -p raw_counts && \
samtools view -F 4 \
  alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1.QueryNameSorted.bam | \
htseq-count - \
  -m intersection-nonempty \
  --stranded=reverse \
  --format=sam \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  > raw_counts/ETOH_shNIPBL_1.readcounts.csv
htseq_count.ETOH_shNIPBL_1.06c65df61e0e38abcc01dc89f5908f58.mugqic.done
)
raw_counts_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$raw_counts_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: raw_counts_10_JOB_ID: htseq_count.ETOH_shMED1_2
#-------------------------------------------------------------------------------
JOB_NAME=htseq_count.ETOH_shMED1_2
JOB_DEPENDENCIES=$picard_sort_sam_10_JOB_ID
JOB_DONE=job_output/raw_counts/htseq_count.ETOH_shMED1_2.2d9aea802edc099bb2e11817e4114a28.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'htseq_count.ETOH_shMED1_2.2d9aea802edc099bb2e11817e4114a28.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/python/2.7.13 && \
mkdir -p raw_counts && \
samtools view -F 4 \
  alignment/ETOH_shMED1_2/ETOH_shMED1_2.QueryNameSorted.bam | \
htseq-count - \
  -m intersection-nonempty \
  --stranded=reverse \
  --format=sam \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  > raw_counts/ETOH_shMED1_2.readcounts.csv
htseq_count.ETOH_shMED1_2.2d9aea802edc099bb2e11817e4114a28.mugqic.done
)
raw_counts_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$raw_counts_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: raw_counts_11_JOB_ID: htseq_count.ETOH_shSMC1A_3
#-------------------------------------------------------------------------------
JOB_NAME=htseq_count.ETOH_shSMC1A_3
JOB_DEPENDENCIES=$picard_sort_sam_11_JOB_ID
JOB_DONE=job_output/raw_counts/htseq_count.ETOH_shSMC1A_3.2eef0b3d3487853aef8c533ca70da912.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'htseq_count.ETOH_shSMC1A_3.2eef0b3d3487853aef8c533ca70da912.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/python/2.7.13 && \
mkdir -p raw_counts && \
samtools view -F 4 \
  alignment/ETOH_shSMC1A_3/ETOH_shSMC1A_3.QueryNameSorted.bam | \
htseq-count - \
  -m intersection-nonempty \
  --stranded=reverse \
  --format=sam \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  > raw_counts/ETOH_shSMC1A_3.readcounts.csv
htseq_count.ETOH_shSMC1A_3.2eef0b3d3487853aef8c533ca70da912.mugqic.done
)
raw_counts_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$raw_counts_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: raw_counts_12_JOB_ID: htseq_count.ETOH_shSMC1A_2
#-------------------------------------------------------------------------------
JOB_NAME=htseq_count.ETOH_shSMC1A_2
JOB_DEPENDENCIES=$picard_sort_sam_12_JOB_ID
JOB_DONE=job_output/raw_counts/htseq_count.ETOH_shSMC1A_2.cadae720da8534f0cf17a74f1a0e7ad8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'htseq_count.ETOH_shSMC1A_2.cadae720da8534f0cf17a74f1a0e7ad8.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/python/2.7.13 && \
mkdir -p raw_counts && \
samtools view -F 4 \
  alignment/ETOH_shSMC1A_2/ETOH_shSMC1A_2.QueryNameSorted.bam | \
htseq-count - \
  -m intersection-nonempty \
  --stranded=reverse \
  --format=sam \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  > raw_counts/ETOH_shSMC1A_2.readcounts.csv
htseq_count.ETOH_shSMC1A_2.cadae720da8534f0cf17a74f1a0e7ad8.mugqic.done
)
raw_counts_12_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$raw_counts_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: raw_counts_metrics
#-------------------------------------------------------------------------------
STEP=raw_counts_metrics
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: raw_counts_metrics_1_JOB_ID: metrics.matrix
#-------------------------------------------------------------------------------
JOB_NAME=metrics.matrix
JOB_DEPENDENCIES=$raw_counts_1_JOB_ID:$raw_counts_2_JOB_ID:$raw_counts_3_JOB_ID:$raw_counts_4_JOB_ID:$raw_counts_5_JOB_ID:$raw_counts_6_JOB_ID:$raw_counts_7_JOB_ID:$raw_counts_8_JOB_ID:$raw_counts_9_JOB_ID:$raw_counts_10_JOB_ID:$raw_counts_11_JOB_ID:$raw_counts_12_JOB_ID
JOB_DONE=job_output/raw_counts_metrics/metrics.matrix.719e4b825d3a4bd02939f34b5219d5e9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'metrics.matrix.719e4b825d3a4bd02939f34b5219d5e9.mugqic.done'
module load mugqic/mugqic_tools/2.1.7 && \
mkdir -p DGE && \
gtf2tmpMatrix.awk \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  DGE/tmpMatrix.txt && \
HEAD='Gene\tSymbol' && \
for read_count_file in \
  raw_counts/DEX_nonMamm_1.readcounts.csv \
  raw_counts/DEX_shMED1_1.readcounts.csv \
  raw_counts/DEX_shNIPBL_1.readcounts.csv \
  raw_counts/DEX_shSMC1A_3.readcounts.csv \
  raw_counts/ETOH_nonMamm_1.readcounts.csv \
  raw_counts/ETOH_nonMamm_2.readcounts.csv \
  raw_counts/ETOH_shMED1_1.readcounts.csv \
  raw_counts/ETOH_shNIPBL_2.readcounts.csv \
  raw_counts/ETOH_shNIPBL_1.readcounts.csv \
  raw_counts/ETOH_shMED1_2.readcounts.csv \
  raw_counts/ETOH_shSMC1A_3.readcounts.csv \
  raw_counts/ETOH_shSMC1A_2.readcounts.csv
do
  sort -k1,1 $read_count_file > DGE/tmpSort.txt && \
  join -1 1 -2 1 <(sort -k1,1 DGE/tmpMatrix.txt) DGE/tmpSort.txt > DGE/tmpMatrix.2.txt && \
  mv DGE/tmpMatrix.2.txt DGE/tmpMatrix.txt && \
  na=$(basename $read_count_file | rev | cut -d. -f3- | rev) && \
  HEAD="$HEAD\t$na"
done && \
echo -e $HEAD | cat - DGE/tmpMatrix.txt | tr ' ' '\t' > DGE/rawCountMatrix.csv && \
rm DGE/tmpSort.txt DGE/tmpMatrix.txt
metrics.matrix.719e4b825d3a4bd02939f34b5219d5e9.mugqic.done
)
raw_counts_metrics_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=5:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$raw_counts_metrics_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: raw_counts_metrics_2_JOB_ID: metrics.wigzip
#-------------------------------------------------------------------------------
JOB_NAME=metrics.wigzip
JOB_DEPENDENCIES=$wiggle_3_JOB_ID:$wiggle_4_JOB_ID:$wiggle_7_JOB_ID:$wiggle_8_JOB_ID:$wiggle_11_JOB_ID:$wiggle_12_JOB_ID:$wiggle_15_JOB_ID:$wiggle_16_JOB_ID:$wiggle_19_JOB_ID:$wiggle_20_JOB_ID:$wiggle_23_JOB_ID:$wiggle_24_JOB_ID:$wiggle_27_JOB_ID:$wiggle_28_JOB_ID:$wiggle_31_JOB_ID:$wiggle_32_JOB_ID:$wiggle_35_JOB_ID:$wiggle_36_JOB_ID:$wiggle_39_JOB_ID:$wiggle_40_JOB_ID:$wiggle_43_JOB_ID:$wiggle_44_JOB_ID:$wiggle_47_JOB_ID:$wiggle_48_JOB_ID
JOB_DONE=job_output/raw_counts_metrics/metrics.wigzip.edc4e268c60b94d072db60163e113f9c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'metrics.wigzip.edc4e268c60b94d072db60163e113f9c.mugqic.done'
zip -r tracks.zip tracks/bigWig
metrics.wigzip.edc4e268c60b94d072db60163e113f9c.mugqic.done
)
raw_counts_metrics_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=5:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$raw_counts_metrics_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: raw_counts_metrics_3_JOB_ID: rpkm_saturation
#-------------------------------------------------------------------------------
JOB_NAME=rpkm_saturation
JOB_DEPENDENCIES=$raw_counts_metrics_1_JOB_ID
JOB_DONE=job_output/raw_counts_metrics/rpkm_saturation.af537e86a0c3c293b0b93c99233044a7.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'rpkm_saturation.af537e86a0c3c293b0b93c99233044a7.mugqic.done'
module load mugqic/R_Bioconductor/3.2.3_3.2 mugqic/mugqic_tools/2.1.7 && \
mkdir -p metrics/saturation && \
Rscript $R_TOOLS/rpkmSaturation.R \
  DGE/rawCountMatrix.csv \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.genes.length.tsv \
  raw_counts \
  metrics/saturation \
  11 \
  1 && \
zip -r metrics/saturation.zip metrics/saturation
rpkm_saturation.af537e86a0c3c293b0b93c99233044a7.mugqic.done
)
raw_counts_metrics_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q lm -l nodes=1:ppn=12 -l pmem=5700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$raw_counts_metrics_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: raw_counts_metrics_4_JOB_ID: raw_count_metrics_report
#-------------------------------------------------------------------------------
JOB_NAME=raw_count_metrics_report
JOB_DEPENDENCIES=$rnaseqc_1_JOB_ID:$raw_counts_metrics_2_JOB_ID:$raw_counts_metrics_3_JOB_ID
JOB_DONE=job_output/raw_counts_metrics/raw_count_metrics_report.9c41ac325d0250b7c03db57f66938c5b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'raw_count_metrics_report.9c41ac325d0250b7c03db57f66938c5b.mugqic.done'
module load mugqic/pandoc/1.15.2 && \
mkdir -p report && \
cp metrics/rnaseqRep/corrMatrixSpearman.txt report/corrMatrixSpearman.tsv && \
cp tracks.zip report/ && \
cp metrics/saturation.zip report/ && \
pandoc --to=markdown \
  --template /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.3.0/bfx/report/RnaSeq.raw_counts_metrics.md \
  --variable corr_matrix_spearman_table="`head -16 report/corrMatrixSpearman.tsv | cut -f-16| awk -F"	" '{OFS="	"; if (NR==1) {$0="Vs"$0; print; gsub(/[^	]/, "-"); print} else {printf $1; for (i=2; i<=NF; i++) {printf "	"sprintf("%.2f", $i)}; print ""}}' | sed 's/	/|/g'`" \
  /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.3.0/bfx/report/RnaSeq.raw_counts_metrics.md \
  > report/RnaSeq.raw_counts_metrics.md
raw_count_metrics_report.9c41ac325d0250b7c03db57f66938c5b.mugqic.done
)
raw_counts_metrics_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$raw_counts_metrics_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: cufflinks
#-------------------------------------------------------------------------------
STEP=cufflinks
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: cufflinks_1_JOB_ID: cufflinks.DEX_nonMamm_1
#-------------------------------------------------------------------------------
JOB_NAME=cufflinks.DEX_nonMamm_1
JOB_DEPENDENCIES=$bam_hard_clip_1_JOB_ID
JOB_DONE=job_output/cufflinks/cufflinks.DEX_nonMamm_1.35e48a755d9e0f008a1c7e01145745ed.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cufflinks.DEX_nonMamm_1.35e48a755d9e0f008a1c7e01145745ed.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/DEX_nonMamm_1 && \
cufflinks -q  \
  --GTF-guide /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/DEX_nonMamm_1 \
  --num-threads 12 \
  alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.hardClip.bam
cufflinks.DEX_nonMamm_1.35e48a755d9e0f008a1c7e01145745ed.mugqic.done
)
cufflinks_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cufflinks_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cufflinks_2_JOB_ID: cufflinks.DEX_shMED1_1
#-------------------------------------------------------------------------------
JOB_NAME=cufflinks.DEX_shMED1_1
JOB_DEPENDENCIES=$bam_hard_clip_2_JOB_ID
JOB_DONE=job_output/cufflinks/cufflinks.DEX_shMED1_1.051b3e058056d2df2c4feee3f7dcc303.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cufflinks.DEX_shMED1_1.051b3e058056d2df2c4feee3f7dcc303.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/DEX_shMED1_1 && \
cufflinks -q  \
  --GTF-guide /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/DEX_shMED1_1 \
  --num-threads 12 \
  alignment/DEX_shMED1_1/DEX_shMED1_1.sorted.mdup.hardClip.bam
cufflinks.DEX_shMED1_1.051b3e058056d2df2c4feee3f7dcc303.mugqic.done
)
cufflinks_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cufflinks_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cufflinks_3_JOB_ID: cufflinks.DEX_shNIPBL_1
#-------------------------------------------------------------------------------
JOB_NAME=cufflinks.DEX_shNIPBL_1
JOB_DEPENDENCIES=$bam_hard_clip_3_JOB_ID
JOB_DONE=job_output/cufflinks/cufflinks.DEX_shNIPBL_1.860951d557f37ad6bae7422940eefda0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cufflinks.DEX_shNIPBL_1.860951d557f37ad6bae7422940eefda0.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/DEX_shNIPBL_1 && \
cufflinks -q  \
  --GTF-guide /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/DEX_shNIPBL_1 \
  --num-threads 12 \
  alignment/DEX_shNIPBL_1/DEX_shNIPBL_1.sorted.mdup.hardClip.bam
cufflinks.DEX_shNIPBL_1.860951d557f37ad6bae7422940eefda0.mugqic.done
)
cufflinks_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cufflinks_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cufflinks_4_JOB_ID: cufflinks.DEX_shSMC1A_3
#-------------------------------------------------------------------------------
JOB_NAME=cufflinks.DEX_shSMC1A_3
JOB_DEPENDENCIES=$bam_hard_clip_4_JOB_ID
JOB_DONE=job_output/cufflinks/cufflinks.DEX_shSMC1A_3.1988b74015190e20e799368d55d9043a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cufflinks.DEX_shSMC1A_3.1988b74015190e20e799368d55d9043a.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/DEX_shSMC1A_3 && \
cufflinks -q  \
  --GTF-guide /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/DEX_shSMC1A_3 \
  --num-threads 12 \
  alignment/DEX_shSMC1A_3/DEX_shSMC1A_3.sorted.mdup.hardClip.bam
cufflinks.DEX_shSMC1A_3.1988b74015190e20e799368d55d9043a.mugqic.done
)
cufflinks_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cufflinks_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cufflinks_5_JOB_ID: cufflinks.ETOH_nonMamm_1
#-------------------------------------------------------------------------------
JOB_NAME=cufflinks.ETOH_nonMamm_1
JOB_DEPENDENCIES=$bam_hard_clip_5_JOB_ID
JOB_DONE=job_output/cufflinks/cufflinks.ETOH_nonMamm_1.7602cd82db01c4b297befa737de596fb.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cufflinks.ETOH_nonMamm_1.7602cd82db01c4b297befa737de596fb.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/ETOH_nonMamm_1 && \
cufflinks -q  \
  --GTF-guide /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/ETOH_nonMamm_1 \
  --num-threads 12 \
  alignment/ETOH_nonMamm_1/ETOH_nonMamm_1.sorted.mdup.hardClip.bam
cufflinks.ETOH_nonMamm_1.7602cd82db01c4b297befa737de596fb.mugqic.done
)
cufflinks_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cufflinks_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cufflinks_6_JOB_ID: cufflinks.ETOH_nonMamm_2
#-------------------------------------------------------------------------------
JOB_NAME=cufflinks.ETOH_nonMamm_2
JOB_DEPENDENCIES=$bam_hard_clip_6_JOB_ID
JOB_DONE=job_output/cufflinks/cufflinks.ETOH_nonMamm_2.1f67ec4062da0d1363f7fcd8f4b5c8d0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cufflinks.ETOH_nonMamm_2.1f67ec4062da0d1363f7fcd8f4b5c8d0.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/ETOH_nonMamm_2 && \
cufflinks -q  \
  --GTF-guide /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/ETOH_nonMamm_2 \
  --num-threads 12 \
  alignment/ETOH_nonMamm_2/ETOH_nonMamm_2.sorted.mdup.hardClip.bam
cufflinks.ETOH_nonMamm_2.1f67ec4062da0d1363f7fcd8f4b5c8d0.mugqic.done
)
cufflinks_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cufflinks_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cufflinks_7_JOB_ID: cufflinks.ETOH_shMED1_1
#-------------------------------------------------------------------------------
JOB_NAME=cufflinks.ETOH_shMED1_1
JOB_DEPENDENCIES=$bam_hard_clip_7_JOB_ID
JOB_DONE=job_output/cufflinks/cufflinks.ETOH_shMED1_1.45193084a12aa2be3e52feb23f287eb4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cufflinks.ETOH_shMED1_1.45193084a12aa2be3e52feb23f287eb4.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/ETOH_shMED1_1 && \
cufflinks -q  \
  --GTF-guide /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/ETOH_shMED1_1 \
  --num-threads 12 \
  alignment/ETOH_shMED1_1/ETOH_shMED1_1.sorted.mdup.hardClip.bam
cufflinks.ETOH_shMED1_1.45193084a12aa2be3e52feb23f287eb4.mugqic.done
)
cufflinks_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cufflinks_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cufflinks_8_JOB_ID: cufflinks.ETOH_shNIPBL_2
#-------------------------------------------------------------------------------
JOB_NAME=cufflinks.ETOH_shNIPBL_2
JOB_DEPENDENCIES=$bam_hard_clip_8_JOB_ID
JOB_DONE=job_output/cufflinks/cufflinks.ETOH_shNIPBL_2.b1f4b10129ef64c13d3570c38906e763.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cufflinks.ETOH_shNIPBL_2.b1f4b10129ef64c13d3570c38906e763.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/ETOH_shNIPBL_2 && \
cufflinks -q  \
  --GTF-guide /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/ETOH_shNIPBL_2 \
  --num-threads 12 \
  alignment/ETOH_shNIPBL_2/ETOH_shNIPBL_2.sorted.mdup.hardClip.bam
cufflinks.ETOH_shNIPBL_2.b1f4b10129ef64c13d3570c38906e763.mugqic.done
)
cufflinks_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cufflinks_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cufflinks_9_JOB_ID: cufflinks.ETOH_shNIPBL_1
#-------------------------------------------------------------------------------
JOB_NAME=cufflinks.ETOH_shNIPBL_1
JOB_DEPENDENCIES=$bam_hard_clip_9_JOB_ID
JOB_DONE=job_output/cufflinks/cufflinks.ETOH_shNIPBL_1.8b8d312645ddb9b13d736c862bfad84c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cufflinks.ETOH_shNIPBL_1.8b8d312645ddb9b13d736c862bfad84c.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/ETOH_shNIPBL_1 && \
cufflinks -q  \
  --GTF-guide /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/ETOH_shNIPBL_1 \
  --num-threads 12 \
  alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1.sorted.mdup.hardClip.bam
cufflinks.ETOH_shNIPBL_1.8b8d312645ddb9b13d736c862bfad84c.mugqic.done
)
cufflinks_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cufflinks_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cufflinks_10_JOB_ID: cufflinks.ETOH_shMED1_2
#-------------------------------------------------------------------------------
JOB_NAME=cufflinks.ETOH_shMED1_2
JOB_DEPENDENCIES=$bam_hard_clip_10_JOB_ID
JOB_DONE=job_output/cufflinks/cufflinks.ETOH_shMED1_2.8e3c75636314f36b1364b35df6c13997.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cufflinks.ETOH_shMED1_2.8e3c75636314f36b1364b35df6c13997.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/ETOH_shMED1_2 && \
cufflinks -q  \
  --GTF-guide /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/ETOH_shMED1_2 \
  --num-threads 12 \
  alignment/ETOH_shMED1_2/ETOH_shMED1_2.sorted.mdup.hardClip.bam
cufflinks.ETOH_shMED1_2.8e3c75636314f36b1364b35df6c13997.mugqic.done
)
cufflinks_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cufflinks_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cufflinks_11_JOB_ID: cufflinks.ETOH_shSMC1A_3
#-------------------------------------------------------------------------------
JOB_NAME=cufflinks.ETOH_shSMC1A_3
JOB_DEPENDENCIES=$bam_hard_clip_11_JOB_ID
JOB_DONE=job_output/cufflinks/cufflinks.ETOH_shSMC1A_3.caa5ff40aeeec2b92495b9cd2610b220.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cufflinks.ETOH_shSMC1A_3.caa5ff40aeeec2b92495b9cd2610b220.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/ETOH_shSMC1A_3 && \
cufflinks -q  \
  --GTF-guide /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/ETOH_shSMC1A_3 \
  --num-threads 12 \
  alignment/ETOH_shSMC1A_3/ETOH_shSMC1A_3.sorted.mdup.hardClip.bam
cufflinks.ETOH_shSMC1A_3.caa5ff40aeeec2b92495b9cd2610b220.mugqic.done
)
cufflinks_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cufflinks_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cufflinks_12_JOB_ID: cufflinks.ETOH_shSMC1A_2
#-------------------------------------------------------------------------------
JOB_NAME=cufflinks.ETOH_shSMC1A_2
JOB_DEPENDENCIES=$bam_hard_clip_12_JOB_ID
JOB_DONE=job_output/cufflinks/cufflinks.ETOH_shSMC1A_2.7900a44c0c6a43e28cb1cc6d26caeebd.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cufflinks.ETOH_shSMC1A_2.7900a44c0c6a43e28cb1cc6d26caeebd.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/ETOH_shSMC1A_2 && \
cufflinks -q  \
  --GTF-guide /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/ETOH_shSMC1A_2 \
  --num-threads 12 \
  alignment/ETOH_shSMC1A_2/ETOH_shSMC1A_2.sorted.mdup.hardClip.bam
cufflinks.ETOH_shSMC1A_2.7900a44c0c6a43e28cb1cc6d26caeebd.mugqic.done
)
cufflinks_12_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cufflinks_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: cuffmerge
#-------------------------------------------------------------------------------
STEP=cuffmerge
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: cuffmerge_1_JOB_ID: cuffmerge
#-------------------------------------------------------------------------------
JOB_NAME=cuffmerge
JOB_DEPENDENCIES=$cufflinks_1_JOB_ID:$cufflinks_2_JOB_ID:$cufflinks_3_JOB_ID:$cufflinks_4_JOB_ID:$cufflinks_5_JOB_ID:$cufflinks_6_JOB_ID:$cufflinks_7_JOB_ID:$cufflinks_8_JOB_ID:$cufflinks_9_JOB_ID:$cufflinks_10_JOB_ID:$cufflinks_11_JOB_ID:$cufflinks_12_JOB_ID
JOB_DONE=job_output/cuffmerge/cuffmerge.d4f0bf35965424caf7e718ea833e1961.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffmerge.d4f0bf35965424caf7e718ea833e1961.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/AllSamples && \
`cat > cufflinks/cuffmerge.samples.txt << END
cufflinks/DEX_nonMamm_1/transcripts.gtf
cufflinks/DEX_shMED1_1/transcripts.gtf
cufflinks/DEX_shNIPBL_1/transcripts.gtf
cufflinks/DEX_shSMC1A_3/transcripts.gtf
cufflinks/ETOH_nonMamm_1/transcripts.gtf
cufflinks/ETOH_nonMamm_2/transcripts.gtf
cufflinks/ETOH_shMED1_1/transcripts.gtf
cufflinks/ETOH_shNIPBL_2/transcripts.gtf
cufflinks/ETOH_shNIPBL_1/transcripts.gtf
cufflinks/ETOH_shMED1_2/transcripts.gtf
cufflinks/ETOH_shSMC1A_3/transcripts.gtf
cufflinks/ETOH_shSMC1A_2/transcripts.gtf
END

` && \
cuffmerge  \
  --ref-gtf /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  --ref-sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  -o cufflinks/AllSamples \
  --num-threads 12 \
  cufflinks/cuffmerge.samples.txt
cuffmerge.d4f0bf35965424caf7e718ea833e1961.mugqic.done
)
cuffmerge_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cuffmerge_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: cuffquant
#-------------------------------------------------------------------------------
STEP=cuffquant
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: cuffquant_1_JOB_ID: cuffquant.DEX_nonMamm_1
#-------------------------------------------------------------------------------
JOB_NAME=cuffquant.DEX_nonMamm_1
JOB_DEPENDENCIES=$bam_hard_clip_1_JOB_ID:$cuffmerge_1_JOB_ID
JOB_DONE=job_output/cuffquant/cuffquant.DEX_nonMamm_1.dc6af3be68f0ef73335e978d79288a55.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffquant.DEX_nonMamm_1.dc6af3be68f0ef73335e978d79288a55.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/DEX_nonMamm_1 && \
cuffquant -q  \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/DEX_nonMamm_1 \
  --num-threads 12 \
  cufflinks/AllSamples/merged.gtf \
  alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.hardClip.bam
cuffquant.DEX_nonMamm_1.dc6af3be68f0ef73335e978d79288a55.mugqic.done
)
cuffquant_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cuffquant_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cuffquant_2_JOB_ID: cuffquant.DEX_shMED1_1
#-------------------------------------------------------------------------------
JOB_NAME=cuffquant.DEX_shMED1_1
JOB_DEPENDENCIES=$bam_hard_clip_2_JOB_ID:$cuffmerge_1_JOB_ID
JOB_DONE=job_output/cuffquant/cuffquant.DEX_shMED1_1.14082c65efb9e33dcebebbb735372436.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffquant.DEX_shMED1_1.14082c65efb9e33dcebebbb735372436.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/DEX_shMED1_1 && \
cuffquant -q  \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/DEX_shMED1_1 \
  --num-threads 12 \
  cufflinks/AllSamples/merged.gtf \
  alignment/DEX_shMED1_1/DEX_shMED1_1.sorted.mdup.hardClip.bam
cuffquant.DEX_shMED1_1.14082c65efb9e33dcebebbb735372436.mugqic.done
)
cuffquant_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cuffquant_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cuffquant_3_JOB_ID: cuffquant.DEX_shNIPBL_1
#-------------------------------------------------------------------------------
JOB_NAME=cuffquant.DEX_shNIPBL_1
JOB_DEPENDENCIES=$bam_hard_clip_3_JOB_ID:$cuffmerge_1_JOB_ID
JOB_DONE=job_output/cuffquant/cuffquant.DEX_shNIPBL_1.23326824bc102779a3d26e21a74c6eb1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffquant.DEX_shNIPBL_1.23326824bc102779a3d26e21a74c6eb1.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/DEX_shNIPBL_1 && \
cuffquant -q  \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/DEX_shNIPBL_1 \
  --num-threads 12 \
  cufflinks/AllSamples/merged.gtf \
  alignment/DEX_shNIPBL_1/DEX_shNIPBL_1.sorted.mdup.hardClip.bam
cuffquant.DEX_shNIPBL_1.23326824bc102779a3d26e21a74c6eb1.mugqic.done
)
cuffquant_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cuffquant_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cuffquant_4_JOB_ID: cuffquant.DEX_shSMC1A_3
#-------------------------------------------------------------------------------
JOB_NAME=cuffquant.DEX_shSMC1A_3
JOB_DEPENDENCIES=$bam_hard_clip_4_JOB_ID:$cuffmerge_1_JOB_ID
JOB_DONE=job_output/cuffquant/cuffquant.DEX_shSMC1A_3.e1fb14aa3392cf97f33b4192a1e2de7f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffquant.DEX_shSMC1A_3.e1fb14aa3392cf97f33b4192a1e2de7f.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/DEX_shSMC1A_3 && \
cuffquant -q  \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/DEX_shSMC1A_3 \
  --num-threads 12 \
  cufflinks/AllSamples/merged.gtf \
  alignment/DEX_shSMC1A_3/DEX_shSMC1A_3.sorted.mdup.hardClip.bam
cuffquant.DEX_shSMC1A_3.e1fb14aa3392cf97f33b4192a1e2de7f.mugqic.done
)
cuffquant_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cuffquant_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cuffquant_5_JOB_ID: cuffquant.ETOH_nonMamm_1
#-------------------------------------------------------------------------------
JOB_NAME=cuffquant.ETOH_nonMamm_1
JOB_DEPENDENCIES=$bam_hard_clip_5_JOB_ID:$cuffmerge_1_JOB_ID
JOB_DONE=job_output/cuffquant/cuffquant.ETOH_nonMamm_1.4ca7c3ed0e2d8262539bca90c2b40ea3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffquant.ETOH_nonMamm_1.4ca7c3ed0e2d8262539bca90c2b40ea3.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/ETOH_nonMamm_1 && \
cuffquant -q  \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/ETOH_nonMamm_1 \
  --num-threads 12 \
  cufflinks/AllSamples/merged.gtf \
  alignment/ETOH_nonMamm_1/ETOH_nonMamm_1.sorted.mdup.hardClip.bam
cuffquant.ETOH_nonMamm_1.4ca7c3ed0e2d8262539bca90c2b40ea3.mugqic.done
)
cuffquant_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cuffquant_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cuffquant_6_JOB_ID: cuffquant.ETOH_nonMamm_2
#-------------------------------------------------------------------------------
JOB_NAME=cuffquant.ETOH_nonMamm_2
JOB_DEPENDENCIES=$bam_hard_clip_6_JOB_ID:$cuffmerge_1_JOB_ID
JOB_DONE=job_output/cuffquant/cuffquant.ETOH_nonMamm_2.dcb4d502ad25056f4e908a5bb4ce2b4f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffquant.ETOH_nonMamm_2.dcb4d502ad25056f4e908a5bb4ce2b4f.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/ETOH_nonMamm_2 && \
cuffquant -q  \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/ETOH_nonMamm_2 \
  --num-threads 12 \
  cufflinks/AllSamples/merged.gtf \
  alignment/ETOH_nonMamm_2/ETOH_nonMamm_2.sorted.mdup.hardClip.bam
cuffquant.ETOH_nonMamm_2.dcb4d502ad25056f4e908a5bb4ce2b4f.mugqic.done
)
cuffquant_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cuffquant_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cuffquant_7_JOB_ID: cuffquant.ETOH_shMED1_1
#-------------------------------------------------------------------------------
JOB_NAME=cuffquant.ETOH_shMED1_1
JOB_DEPENDENCIES=$bam_hard_clip_7_JOB_ID:$cuffmerge_1_JOB_ID
JOB_DONE=job_output/cuffquant/cuffquant.ETOH_shMED1_1.ca70f67ec216b5f313a692a6c402c0d5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffquant.ETOH_shMED1_1.ca70f67ec216b5f313a692a6c402c0d5.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/ETOH_shMED1_1 && \
cuffquant -q  \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/ETOH_shMED1_1 \
  --num-threads 12 \
  cufflinks/AllSamples/merged.gtf \
  alignment/ETOH_shMED1_1/ETOH_shMED1_1.sorted.mdup.hardClip.bam
cuffquant.ETOH_shMED1_1.ca70f67ec216b5f313a692a6c402c0d5.mugqic.done
)
cuffquant_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cuffquant_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cuffquant_8_JOB_ID: cuffquant.ETOH_shNIPBL_2
#-------------------------------------------------------------------------------
JOB_NAME=cuffquant.ETOH_shNIPBL_2
JOB_DEPENDENCIES=$bam_hard_clip_8_JOB_ID:$cuffmerge_1_JOB_ID
JOB_DONE=job_output/cuffquant/cuffquant.ETOH_shNIPBL_2.81a44c82dd9a1f200a78925541ec37d1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffquant.ETOH_shNIPBL_2.81a44c82dd9a1f200a78925541ec37d1.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/ETOH_shNIPBL_2 && \
cuffquant -q  \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/ETOH_shNIPBL_2 \
  --num-threads 12 \
  cufflinks/AllSamples/merged.gtf \
  alignment/ETOH_shNIPBL_2/ETOH_shNIPBL_2.sorted.mdup.hardClip.bam
cuffquant.ETOH_shNIPBL_2.81a44c82dd9a1f200a78925541ec37d1.mugqic.done
)
cuffquant_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cuffquant_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cuffquant_9_JOB_ID: cuffquant.ETOH_shNIPBL_1
#-------------------------------------------------------------------------------
JOB_NAME=cuffquant.ETOH_shNIPBL_1
JOB_DEPENDENCIES=$bam_hard_clip_9_JOB_ID:$cuffmerge_1_JOB_ID
JOB_DONE=job_output/cuffquant/cuffquant.ETOH_shNIPBL_1.f3dd4d817bd61b51d37b1e936cd48314.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffquant.ETOH_shNIPBL_1.f3dd4d817bd61b51d37b1e936cd48314.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/ETOH_shNIPBL_1 && \
cuffquant -q  \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/ETOH_shNIPBL_1 \
  --num-threads 12 \
  cufflinks/AllSamples/merged.gtf \
  alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1.sorted.mdup.hardClip.bam
cuffquant.ETOH_shNIPBL_1.f3dd4d817bd61b51d37b1e936cd48314.mugqic.done
)
cuffquant_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cuffquant_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cuffquant_10_JOB_ID: cuffquant.ETOH_shMED1_2
#-------------------------------------------------------------------------------
JOB_NAME=cuffquant.ETOH_shMED1_2
JOB_DEPENDENCIES=$bam_hard_clip_10_JOB_ID:$cuffmerge_1_JOB_ID
JOB_DONE=job_output/cuffquant/cuffquant.ETOH_shMED1_2.94307be1dfc655a4a4a3eb208be03351.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffquant.ETOH_shMED1_2.94307be1dfc655a4a4a3eb208be03351.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/ETOH_shMED1_2 && \
cuffquant -q  \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/ETOH_shMED1_2 \
  --num-threads 12 \
  cufflinks/AllSamples/merged.gtf \
  alignment/ETOH_shMED1_2/ETOH_shMED1_2.sorted.mdup.hardClip.bam
cuffquant.ETOH_shMED1_2.94307be1dfc655a4a4a3eb208be03351.mugqic.done
)
cuffquant_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cuffquant_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cuffquant_11_JOB_ID: cuffquant.ETOH_shSMC1A_3
#-------------------------------------------------------------------------------
JOB_NAME=cuffquant.ETOH_shSMC1A_3
JOB_DEPENDENCIES=$bam_hard_clip_11_JOB_ID:$cuffmerge_1_JOB_ID
JOB_DONE=job_output/cuffquant/cuffquant.ETOH_shSMC1A_3.c02cd2a2cbe1323cf416d76d4ab9739e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffquant.ETOH_shSMC1A_3.c02cd2a2cbe1323cf416d76d4ab9739e.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/ETOH_shSMC1A_3 && \
cuffquant -q  \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/ETOH_shSMC1A_3 \
  --num-threads 12 \
  cufflinks/AllSamples/merged.gtf \
  alignment/ETOH_shSMC1A_3/ETOH_shSMC1A_3.sorted.mdup.hardClip.bam
cuffquant.ETOH_shSMC1A_3.c02cd2a2cbe1323cf416d76d4ab9739e.mugqic.done
)
cuffquant_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cuffquant_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cuffquant_12_JOB_ID: cuffquant.ETOH_shSMC1A_2
#-------------------------------------------------------------------------------
JOB_NAME=cuffquant.ETOH_shSMC1A_2
JOB_DEPENDENCIES=$bam_hard_clip_12_JOB_ID:$cuffmerge_1_JOB_ID
JOB_DONE=job_output/cuffquant/cuffquant.ETOH_shSMC1A_2.53689cbe1e9e5556c2433a0bbb4faf66.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffquant.ETOH_shSMC1A_2.53689cbe1e9e5556c2433a0bbb4faf66.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/ETOH_shSMC1A_2 && \
cuffquant -q  \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/ETOH_shSMC1A_2 \
  --num-threads 12 \
  cufflinks/AllSamples/merged.gtf \
  alignment/ETOH_shSMC1A_2/ETOH_shSMC1A_2.sorted.mdup.hardClip.bam
cuffquant.ETOH_shSMC1A_2.53689cbe1e9e5556c2433a0bbb4faf66.mugqic.done
)
cuffquant_12_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cuffquant_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: cuffdiff
#-------------------------------------------------------------------------------
STEP=cuffdiff
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: cuffdiff_1_JOB_ID: cuffdiff.Dex_vs_EtOH-All
#-------------------------------------------------------------------------------
JOB_NAME=cuffdiff.Dex_vs_EtOH-All
JOB_DEPENDENCIES=$cuffmerge_1_JOB_ID:$cuffquant_1_JOB_ID:$cuffquant_2_JOB_ID:$cuffquant_3_JOB_ID:$cuffquant_4_JOB_ID:$cuffquant_5_JOB_ID:$cuffquant_6_JOB_ID:$cuffquant_7_JOB_ID:$cuffquant_8_JOB_ID:$cuffquant_9_JOB_ID:$cuffquant_10_JOB_ID:$cuffquant_11_JOB_ID:$cuffquant_12_JOB_ID
JOB_DONE=job_output/cuffdiff/cuffdiff.Dex_vs_EtOH-All.5e1fea96a2a8c966270cfa83c1ad2c9f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffdiff.Dex_vs_EtOH-All.5e1fea96a2a8c966270cfa83c1ad2c9f.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cuffdiff/Dex_vs_EtOH-All && \
cuffdiff -u \
  --frag-bias-correct /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --library-type fr-firststrand \
  --output-dir cuffdiff/Dex_vs_EtOH-All \
  --num-threads 12 \
  cufflinks/AllSamples/merged.gtf \
  cufflinks/ETOH_nonMamm_1/abundances.cxb,cufflinks/ETOH_nonMamm_2/abundances.cxb,cufflinks/ETOH_shMED1_1/abundances.cxb,cufflinks/ETOH_shMED1_2/abundances.cxb,cufflinks/ETOH_shNIPBL_1/abundances.cxb,cufflinks/ETOH_shNIPBL_2/abundances.cxb,cufflinks/ETOH_shSMC1A_3/abundances.cxb,cufflinks/ETOH_shSMC1A_2/abundances.cxb \
  cufflinks/DEX_nonMamm_1/abundances.cxb,cufflinks/DEX_shMED1_1/abundances.cxb,cufflinks/DEX_shNIPBL_1/abundances.cxb,cufflinks/DEX_shSMC1A_3/abundances.cxb
cuffdiff.Dex_vs_EtOH-All.5e1fea96a2a8c966270cfa83c1ad2c9f.mugqic.done
)
cuffdiff_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cuffdiff_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cuffdiff_2_JOB_ID: cuffdiff.Dex_vs_EtOH-nonMamm
#-------------------------------------------------------------------------------
JOB_NAME=cuffdiff.Dex_vs_EtOH-nonMamm
JOB_DEPENDENCIES=$cuffmerge_1_JOB_ID:$cuffquant_1_JOB_ID:$cuffquant_5_JOB_ID:$cuffquant_6_JOB_ID
JOB_DONE=job_output/cuffdiff/cuffdiff.Dex_vs_EtOH-nonMamm.b7d095516f72d7a83a330768097a9ee1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffdiff.Dex_vs_EtOH-nonMamm.b7d095516f72d7a83a330768097a9ee1.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cuffdiff/Dex_vs_EtOH-nonMamm && \
cuffdiff -u \
  --frag-bias-correct /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --library-type fr-firststrand \
  --output-dir cuffdiff/Dex_vs_EtOH-nonMamm \
  --num-threads 12 \
  cufflinks/AllSamples/merged.gtf \
  cufflinks/ETOH_nonMamm_1/abundances.cxb,cufflinks/ETOH_nonMamm_2/abundances.cxb \
  cufflinks/DEX_nonMamm_1/abundances.cxb
cuffdiff.Dex_vs_EtOH-nonMamm.b7d095516f72d7a83a330768097a9ee1.mugqic.done
)
cuffdiff_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cuffdiff_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cuffdiff_3_JOB_ID: cuffdiff.Dex_vs_EtOH-shMED1
#-------------------------------------------------------------------------------
JOB_NAME=cuffdiff.Dex_vs_EtOH-shMED1
JOB_DEPENDENCIES=$cuffmerge_1_JOB_ID:$cuffquant_2_JOB_ID:$cuffquant_7_JOB_ID:$cuffquant_10_JOB_ID
JOB_DONE=job_output/cuffdiff/cuffdiff.Dex_vs_EtOH-shMED1.382f6da8688c2d081389fff0f4db2cbf.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffdiff.Dex_vs_EtOH-shMED1.382f6da8688c2d081389fff0f4db2cbf.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cuffdiff/Dex_vs_EtOH-shMED1 && \
cuffdiff -u \
  --frag-bias-correct /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --library-type fr-firststrand \
  --output-dir cuffdiff/Dex_vs_EtOH-shMED1 \
  --num-threads 12 \
  cufflinks/AllSamples/merged.gtf \
  cufflinks/ETOH_shMED1_1/abundances.cxb,cufflinks/ETOH_shMED1_2/abundances.cxb \
  cufflinks/DEX_shMED1_1/abundances.cxb
cuffdiff.Dex_vs_EtOH-shMED1.382f6da8688c2d081389fff0f4db2cbf.mugqic.done
)
cuffdiff_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cuffdiff_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cuffdiff_4_JOB_ID: cuffdiff.Dex_vs_EtOH-shNIPBL
#-------------------------------------------------------------------------------
JOB_NAME=cuffdiff.Dex_vs_EtOH-shNIPBL
JOB_DEPENDENCIES=$cuffmerge_1_JOB_ID:$cuffquant_3_JOB_ID:$cuffquant_8_JOB_ID:$cuffquant_9_JOB_ID
JOB_DONE=job_output/cuffdiff/cuffdiff.Dex_vs_EtOH-shNIPBL.2cfd91152ff64a728431bf8bd657e3c0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffdiff.Dex_vs_EtOH-shNIPBL.2cfd91152ff64a728431bf8bd657e3c0.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cuffdiff/Dex_vs_EtOH-shNIPBL && \
cuffdiff -u \
  --frag-bias-correct /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --library-type fr-firststrand \
  --output-dir cuffdiff/Dex_vs_EtOH-shNIPBL \
  --num-threads 12 \
  cufflinks/AllSamples/merged.gtf \
  cufflinks/ETOH_shNIPBL_1/abundances.cxb,cufflinks/ETOH_shNIPBL_2/abundances.cxb \
  cufflinks/DEX_shNIPBL_1/abundances.cxb
cuffdiff.Dex_vs_EtOH-shNIPBL.2cfd91152ff64a728431bf8bd657e3c0.mugqic.done
)
cuffdiff_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cuffdiff_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cuffdiff_5_JOB_ID: cuffdiff.Dex_vs_EtOH-shSMC1A
#-------------------------------------------------------------------------------
JOB_NAME=cuffdiff.Dex_vs_EtOH-shSMC1A
JOB_DEPENDENCIES=$cuffmerge_1_JOB_ID:$cuffquant_4_JOB_ID:$cuffquant_11_JOB_ID:$cuffquant_12_JOB_ID
JOB_DONE=job_output/cuffdiff/cuffdiff.Dex_vs_EtOH-shSMC1A.2b73f7479258c38e96120e2bf7e8b6a7.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffdiff.Dex_vs_EtOH-shSMC1A.2b73f7479258c38e96120e2bf7e8b6a7.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cuffdiff/Dex_vs_EtOH-shSMC1A && \
cuffdiff -u \
  --frag-bias-correct /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --library-type fr-firststrand \
  --output-dir cuffdiff/Dex_vs_EtOH-shSMC1A \
  --num-threads 12 \
  cufflinks/AllSamples/merged.gtf \
  cufflinks/ETOH_shSMC1A_3/abundances.cxb,cufflinks/ETOH_shSMC1A_2/abundances.cxb \
  cufflinks/DEX_shSMC1A_3/abundances.cxb
cuffdiff.Dex_vs_EtOH-shSMC1A.2b73f7479258c38e96120e2bf7e8b6a7.mugqic.done
)
cuffdiff_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cuffdiff_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cuffdiff_6_JOB_ID: cuffdiff.Dex_vs_EtOH-shSMC1A-3
#-------------------------------------------------------------------------------
JOB_NAME=cuffdiff.Dex_vs_EtOH-shSMC1A-3
JOB_DEPENDENCIES=$cuffmerge_1_JOB_ID:$cuffquant_4_JOB_ID:$cuffquant_11_JOB_ID
JOB_DONE=job_output/cuffdiff/cuffdiff.Dex_vs_EtOH-shSMC1A-3.b013f525f89afee535b6383dfe510194.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffdiff.Dex_vs_EtOH-shSMC1A-3.b013f525f89afee535b6383dfe510194.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cuffdiff/Dex_vs_EtOH-shSMC1A-3 && \
cuffdiff -u \
  --frag-bias-correct /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --library-type fr-firststrand \
  --output-dir cuffdiff/Dex_vs_EtOH-shSMC1A-3 \
  --num-threads 12 \
  cufflinks/AllSamples/merged.gtf \
  cufflinks/ETOH_shSMC1A_3/abundances.cxb \
  cufflinks/DEX_shSMC1A_3/abundances.cxb
cuffdiff.Dex_vs_EtOH-shSMC1A-3.b013f525f89afee535b6383dfe510194.mugqic.done
)
cuffdiff_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cuffdiff_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cuffdiff_7_JOB_ID: cuffdiff.shMED1_vs_shMamm-Dex
#-------------------------------------------------------------------------------
JOB_NAME=cuffdiff.shMED1_vs_shMamm-Dex
JOB_DEPENDENCIES=$cuffmerge_1_JOB_ID:$cuffquant_1_JOB_ID:$cuffquant_2_JOB_ID
JOB_DONE=job_output/cuffdiff/cuffdiff.shMED1_vs_shMamm-Dex.67ba35334ceeab17584c85394f3816dd.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffdiff.shMED1_vs_shMamm-Dex.67ba35334ceeab17584c85394f3816dd.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cuffdiff/shMED1_vs_shMamm-Dex && \
cuffdiff -u \
  --frag-bias-correct /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --library-type fr-firststrand \
  --output-dir cuffdiff/shMED1_vs_shMamm-Dex \
  --num-threads 12 \
  cufflinks/AllSamples/merged.gtf \
  cufflinks/DEX_nonMamm_1/abundances.cxb \
  cufflinks/DEX_shMED1_1/abundances.cxb
cuffdiff.shMED1_vs_shMamm-Dex.67ba35334ceeab17584c85394f3816dd.mugqic.done
)
cuffdiff_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cuffdiff_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cuffdiff_8_JOB_ID: cuffdiff.shMED1_vs_shMamm-EtOH
#-------------------------------------------------------------------------------
JOB_NAME=cuffdiff.shMED1_vs_shMamm-EtOH
JOB_DEPENDENCIES=$cuffmerge_1_JOB_ID:$cuffquant_5_JOB_ID:$cuffquant_6_JOB_ID:$cuffquant_7_JOB_ID:$cuffquant_10_JOB_ID
JOB_DONE=job_output/cuffdiff/cuffdiff.shMED1_vs_shMamm-EtOH.427aef01d832a7b7b1cdb84f6df10c5d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffdiff.shMED1_vs_shMamm-EtOH.427aef01d832a7b7b1cdb84f6df10c5d.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cuffdiff/shMED1_vs_shMamm-EtOH && \
cuffdiff -u \
  --frag-bias-correct /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --library-type fr-firststrand \
  --output-dir cuffdiff/shMED1_vs_shMamm-EtOH \
  --num-threads 12 \
  cufflinks/AllSamples/merged.gtf \
  cufflinks/ETOH_nonMamm_1/abundances.cxb,cufflinks/ETOH_nonMamm_2/abundances.cxb \
  cufflinks/ETOH_shMED1_1/abundances.cxb,cufflinks/ETOH_shMED1_2/abundances.cxb
cuffdiff.shMED1_vs_shMamm-EtOH.427aef01d832a7b7b1cdb84f6df10c5d.mugqic.done
)
cuffdiff_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cuffdiff_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cuffdiff_9_JOB_ID: cuffdiff.shNIPBL_vs_shMamm-Dex
#-------------------------------------------------------------------------------
JOB_NAME=cuffdiff.shNIPBL_vs_shMamm-Dex
JOB_DEPENDENCIES=$cuffmerge_1_JOB_ID:$cuffquant_1_JOB_ID:$cuffquant_3_JOB_ID
JOB_DONE=job_output/cuffdiff/cuffdiff.shNIPBL_vs_shMamm-Dex.b19217fa4d8011c84a5de5ee6a36dafb.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffdiff.shNIPBL_vs_shMamm-Dex.b19217fa4d8011c84a5de5ee6a36dafb.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cuffdiff/shNIPBL_vs_shMamm-Dex && \
cuffdiff -u \
  --frag-bias-correct /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --library-type fr-firststrand \
  --output-dir cuffdiff/shNIPBL_vs_shMamm-Dex \
  --num-threads 12 \
  cufflinks/AllSamples/merged.gtf \
  cufflinks/DEX_nonMamm_1/abundances.cxb \
  cufflinks/DEX_shNIPBL_1/abundances.cxb
cuffdiff.shNIPBL_vs_shMamm-Dex.b19217fa4d8011c84a5de5ee6a36dafb.mugqic.done
)
cuffdiff_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cuffdiff_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cuffdiff_10_JOB_ID: cuffdiff.shNIPBL_vs_shMamm-EtOH
#-------------------------------------------------------------------------------
JOB_NAME=cuffdiff.shNIPBL_vs_shMamm-EtOH
JOB_DEPENDENCIES=$cuffmerge_1_JOB_ID:$cuffquant_5_JOB_ID:$cuffquant_6_JOB_ID:$cuffquant_8_JOB_ID:$cuffquant_9_JOB_ID
JOB_DONE=job_output/cuffdiff/cuffdiff.shNIPBL_vs_shMamm-EtOH.27aaa961efa077e8b23c0093b997eec0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffdiff.shNIPBL_vs_shMamm-EtOH.27aaa961efa077e8b23c0093b997eec0.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cuffdiff/shNIPBL_vs_shMamm-EtOH && \
cuffdiff -u \
  --frag-bias-correct /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --library-type fr-firststrand \
  --output-dir cuffdiff/shNIPBL_vs_shMamm-EtOH \
  --num-threads 12 \
  cufflinks/AllSamples/merged.gtf \
  cufflinks/ETOH_nonMamm_1/abundances.cxb,cufflinks/ETOH_nonMamm_2/abundances.cxb \
  cufflinks/ETOH_shNIPBL_1/abundances.cxb,cufflinks/ETOH_shNIPBL_2/abundances.cxb
cuffdiff.shNIPBL_vs_shMamm-EtOH.27aaa961efa077e8b23c0093b997eec0.mugqic.done
)
cuffdiff_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cuffdiff_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cuffdiff_11_JOB_ID: cuffdiff.shSMC1A-3_vs_shMamm-Dex
#-------------------------------------------------------------------------------
JOB_NAME=cuffdiff.shSMC1A-3_vs_shMamm-Dex
JOB_DEPENDENCIES=$cuffmerge_1_JOB_ID:$cuffquant_1_JOB_ID:$cuffquant_4_JOB_ID
JOB_DONE=job_output/cuffdiff/cuffdiff.shSMC1A-3_vs_shMamm-Dex.4e046fcdc1987a9f88d80e054cba0ea4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffdiff.shSMC1A-3_vs_shMamm-Dex.4e046fcdc1987a9f88d80e054cba0ea4.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cuffdiff/shSMC1A-3_vs_shMamm-Dex && \
cuffdiff -u \
  --frag-bias-correct /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --library-type fr-firststrand \
  --output-dir cuffdiff/shSMC1A-3_vs_shMamm-Dex \
  --num-threads 12 \
  cufflinks/AllSamples/merged.gtf \
  cufflinks/DEX_nonMamm_1/abundances.cxb \
  cufflinks/DEX_shSMC1A_3/abundances.cxb
cuffdiff.shSMC1A-3_vs_shMamm-Dex.4e046fcdc1987a9f88d80e054cba0ea4.mugqic.done
)
cuffdiff_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cuffdiff_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cuffdiff_12_JOB_ID: cuffdiff.shSMC1A_vs_shMamm-EtOH
#-------------------------------------------------------------------------------
JOB_NAME=cuffdiff.shSMC1A_vs_shMamm-EtOH
JOB_DEPENDENCIES=$cuffmerge_1_JOB_ID:$cuffquant_5_JOB_ID:$cuffquant_6_JOB_ID:$cuffquant_11_JOB_ID:$cuffquant_12_JOB_ID
JOB_DONE=job_output/cuffdiff/cuffdiff.shSMC1A_vs_shMamm-EtOH.ac043f622d513b96674d6e30e10b691e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffdiff.shSMC1A_vs_shMamm-EtOH.ac043f622d513b96674d6e30e10b691e.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cuffdiff/shSMC1A_vs_shMamm-EtOH && \
cuffdiff -u \
  --frag-bias-correct /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --library-type fr-firststrand \
  --output-dir cuffdiff/shSMC1A_vs_shMamm-EtOH \
  --num-threads 12 \
  cufflinks/AllSamples/merged.gtf \
  cufflinks/ETOH_nonMamm_1/abundances.cxb,cufflinks/ETOH_nonMamm_2/abundances.cxb \
  cufflinks/ETOH_shSMC1A_3/abundances.cxb,cufflinks/ETOH_shSMC1A_2/abundances.cxb
cuffdiff.shSMC1A_vs_shMamm-EtOH.ac043f622d513b96674d6e30e10b691e.mugqic.done
)
cuffdiff_12_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cuffdiff_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cuffdiff_13_JOB_ID: cuffdiff.shSMC1A-3_vs_shMamm-EtOH
#-------------------------------------------------------------------------------
JOB_NAME=cuffdiff.shSMC1A-3_vs_shMamm-EtOH
JOB_DEPENDENCIES=$cuffmerge_1_JOB_ID:$cuffquant_5_JOB_ID:$cuffquant_6_JOB_ID:$cuffquant_11_JOB_ID
JOB_DONE=job_output/cuffdiff/cuffdiff.shSMC1A-3_vs_shMamm-EtOH.d3ba9b041d80f4faea15625bab14f878.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffdiff.shSMC1A-3_vs_shMamm-EtOH.d3ba9b041d80f4faea15625bab14f878.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cuffdiff/shSMC1A-3_vs_shMamm-EtOH && \
cuffdiff -u \
  --frag-bias-correct /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --library-type fr-firststrand \
  --output-dir cuffdiff/shSMC1A-3_vs_shMamm-EtOH \
  --num-threads 12 \
  cufflinks/AllSamples/merged.gtf \
  cufflinks/ETOH_nonMamm_1/abundances.cxb,cufflinks/ETOH_nonMamm_2/abundances.cxb \
  cufflinks/ETOH_shSMC1A_3/abundances.cxb
cuffdiff.shSMC1A-3_vs_shMamm-EtOH.d3ba9b041d80f4faea15625bab14f878.mugqic.done
)
cuffdiff_13_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cuffdiff_13_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: cuffnorm
#-------------------------------------------------------------------------------
STEP=cuffnorm
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: cuffnorm_1_JOB_ID: cuffnorm
#-------------------------------------------------------------------------------
JOB_NAME=cuffnorm
JOB_DEPENDENCIES=$cuffmerge_1_JOB_ID:$cuffquant_1_JOB_ID:$cuffquant_2_JOB_ID:$cuffquant_3_JOB_ID:$cuffquant_4_JOB_ID:$cuffquant_5_JOB_ID:$cuffquant_6_JOB_ID:$cuffquant_7_JOB_ID:$cuffquant_8_JOB_ID:$cuffquant_9_JOB_ID:$cuffquant_10_JOB_ID:$cuffquant_11_JOB_ID:$cuffquant_12_JOB_ID
JOB_DONE=job_output/cuffnorm/cuffnorm.7b6762ea1a70300d967550081de1336c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffnorm.7b6762ea1a70300d967550081de1336c.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cuffnorm && \
cuffnorm -q  \
  --library-type fr-firststrand \
  --output-dir cuffnorm \
  --num-threads 12 \
  --labels DEX_nonMamm_1,DEX_shMED1_1,DEX_shNIPBL_1,DEX_shSMC1A_3,ETOH_nonMamm_1,ETOH_nonMamm_2,ETOH_shMED1_1,ETOH_shNIPBL_2,ETOH_shNIPBL_1,ETOH_shMED1_2,ETOH_shSMC1A_3,ETOH_shSMC1A_2 \
  cufflinks/AllSamples/merged.gtf \
  cufflinks/DEX_nonMamm_1/abundances.cxb \
  cufflinks/DEX_shMED1_1/abundances.cxb \
  cufflinks/DEX_shNIPBL_1/abundances.cxb \
  cufflinks/DEX_shSMC1A_3/abundances.cxb \
  cufflinks/ETOH_nonMamm_1/abundances.cxb \
  cufflinks/ETOH_nonMamm_2/abundances.cxb \
  cufflinks/ETOH_shMED1_1/abundances.cxb \
  cufflinks/ETOH_shNIPBL_2/abundances.cxb \
  cufflinks/ETOH_shNIPBL_1/abundances.cxb \
  cufflinks/ETOH_shMED1_2/abundances.cxb \
  cufflinks/ETOH_shSMC1A_3/abundances.cxb \
  cufflinks/ETOH_shSMC1A_2/abundances.cxb
cuffnorm.7b6762ea1a70300d967550081de1336c.mugqic.done
)
cuffnorm_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cuffnorm_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: fpkm_correlation_matrix
#-------------------------------------------------------------------------------
STEP=fpkm_correlation_matrix
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: fpkm_correlation_matrix_1_JOB_ID: fpkm_correlation_matrix_transcript
#-------------------------------------------------------------------------------
JOB_NAME=fpkm_correlation_matrix_transcript
JOB_DEPENDENCIES=$cuffnorm_1_JOB_ID
JOB_DONE=job_output/fpkm_correlation_matrix/fpkm_correlation_matrix_transcript.2d56bdad701a899cf3a678871c049eff.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'fpkm_correlation_matrix_transcript.2d56bdad701a899cf3a678871c049eff.mugqic.done'
module load mugqic/R_Bioconductor/3.2.3_3.2 && \
mkdir -p metrics && \
R --no-save --no-restore <<-EOF
dataFile=read.table("cuffnorm/isoforms.fpkm_table",header=T,check.names=F)
fpkm=cbind(dataFile[,2:ncol(dataFile)])
corTable=cor(log2(fpkm+0.1))
corTableOut=rbind(c('Vs.',colnames(corTable)),cbind(rownames(corTable),round(corTable,3)))
write.table(corTableOut,file="metrics/transcripts_fpkm_correlation_matrix.tsv",col.names=F,row.names=F,sep="	",quote=F)
print("done.")

EOF
fpkm_correlation_matrix_transcript.2d56bdad701a899cf3a678871c049eff.mugqic.done
)
fpkm_correlation_matrix_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$fpkm_correlation_matrix_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: fpkm_correlation_matrix_2_JOB_ID: fpkm_correlation_matrix_gene
#-------------------------------------------------------------------------------
JOB_NAME=fpkm_correlation_matrix_gene
JOB_DEPENDENCIES=$cuffnorm_1_JOB_ID
JOB_DONE=job_output/fpkm_correlation_matrix/fpkm_correlation_matrix_gene.ed1ff067871efb85a076c7f045782017.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'fpkm_correlation_matrix_gene.ed1ff067871efb85a076c7f045782017.mugqic.done'
module load mugqic/R_Bioconductor/3.2.3_3.2 && \
R --no-save --no-restore <<-EOF
dataFile=read.table("cuffnorm/genes.fpkm_table",header=T,check.names=F)
fpkm=cbind(dataFile[,2:ncol(dataFile)])
corTable=cor(log2(fpkm+0.1))
corTableOut=rbind(c('Vs.',colnames(corTable)),cbind(rownames(corTable),round(corTable,3)))
write.table(corTableOut,file="metrics/gene_fpkm_correlation_matrix.tsv",col.names=F,row.names=F,sep="	",quote=F)
print("done.")

EOF
fpkm_correlation_matrix_gene.ed1ff067871efb85a076c7f045782017.mugqic.done
)
fpkm_correlation_matrix_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$fpkm_correlation_matrix_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: gq_seq_utils_exploratory_analysis_rnaseq
#-------------------------------------------------------------------------------
STEP=gq_seq_utils_exploratory_analysis_rnaseq
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: gq_seq_utils_exploratory_analysis_rnaseq_1_JOB_ID: gq_seq_utils_exploratory_analysis_rnaseq
#-------------------------------------------------------------------------------
JOB_NAME=gq_seq_utils_exploratory_analysis_rnaseq
JOB_DEPENDENCIES=$raw_counts_metrics_1_JOB_ID:$cuffnorm_1_JOB_ID
JOB_DONE=job_output/gq_seq_utils_exploratory_analysis_rnaseq/gq_seq_utils_exploratory_analysis_rnaseq.12b0410e480afc76cb1b0b732cc2ff1e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'gq_seq_utils_exploratory_analysis_rnaseq.12b0410e480afc76cb1b0b732cc2ff1e.mugqic.done'
module load mugqic/R_Bioconductor/3.2.3_3.2 mugqic/mugqic_R_packages/1.0.4 && \
mkdir -p exploratory && \
R --no-save --no-restore <<-EOF
suppressPackageStartupMessages(library(gqSeqUtils))

exploratoryAnalysisRNAseq(htseq.counts.path="DGE/rawCountMatrix.csv", cuffnorm.fpkms.dir="cuffnorm", genes.path="/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.genes.tsv", output.dir="exploratory")
desc = readRDS(file.path("exploratory","index.RData"))
write.table(desc,file=file.path("exploratory","index.tsv"),sep='	',quote=F,col.names=T,row.names=F)
print("done.")

EOF
gq_seq_utils_exploratory_analysis_rnaseq.12b0410e480afc76cb1b0b732cc2ff1e.mugqic.done
)
gq_seq_utils_exploratory_analysis_rnaseq_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=00:30:0 -q metaq -l nodes=1:ppn=1 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$gq_seq_utils_exploratory_analysis_rnaseq_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: gq_seq_utils_exploratory_analysis_rnaseq_2_JOB_ID: gq_seq_utils_exploratory_analysis_rnaseq_report
#-------------------------------------------------------------------------------
JOB_NAME=gq_seq_utils_exploratory_analysis_rnaseq_report
JOB_DEPENDENCIES=$gq_seq_utils_exploratory_analysis_rnaseq_1_JOB_ID
JOB_DONE=job_output/gq_seq_utils_exploratory_analysis_rnaseq/gq_seq_utils_exploratory_analysis_rnaseq_report.e5d39d481ff0e98e6b12fbacde49f832.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'gq_seq_utils_exploratory_analysis_rnaseq_report.e5d39d481ff0e98e6b12fbacde49f832.mugqic.done'
module load mugqic/R_Bioconductor/3.2.3_3.2 mugqic/pandoc/1.15.2 && \
R --no-save --no-restore <<-'EOF'
report_dir="report";
input_rmarkdown_file = '/cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.3.0/bfx/report/RnaSeq.gq_seq_utils_exploratory_analysis_rnaseq.Rmd'
render_output_dir    = 'report'
rmarkdown_file       = basename(input_rmarkdown_file) # honoring a different WD that location of Rmd file in knitr is problematic
file.copy(from = input_rmarkdown_file, to = rmarkdown_file, overwrite = T)
rmarkdown::render(input = rmarkdown_file, output_format = c("html_document","md_document"), output_dir = render_output_dir  )
file.remove(rmarkdown_file)
EOF
gq_seq_utils_exploratory_analysis_rnaseq_report.e5d39d481ff0e98e6b12fbacde49f832.mugqic.done
)
gq_seq_utils_exploratory_analysis_rnaseq_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$gq_seq_utils_exploratory_analysis_rnaseq_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: gq_seq_utils_exploratory_analysis_rnaseq_3_JOB_ID: cuffnorm_report
#-------------------------------------------------------------------------------
JOB_NAME=cuffnorm_report
JOB_DEPENDENCIES=$cuffmerge_1_JOB_ID
JOB_DONE=job_output/gq_seq_utils_exploratory_analysis_rnaseq/cuffnorm_report.1dc3c217cf66333221bbfbb81befd479.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffnorm_report.1dc3c217cf66333221bbfbb81befd479.mugqic.done'
mkdir -p report && \
zip -r report/cuffAnalysis.zip cufflinks/ cuffdiff/ cuffnorm/ && \
cp \
  /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.3.0/bfx/report/RnaSeq.cuffnorm.md \
  report/RnaSeq.cuffnorm.md
cuffnorm_report.1dc3c217cf66333221bbfbb81befd479.mugqic.done
)
gq_seq_utils_exploratory_analysis_rnaseq_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$gq_seq_utils_exploratory_analysis_rnaseq_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: differential_expression
#-------------------------------------------------------------------------------
STEP=differential_expression
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: differential_expression_1_JOB_ID: differential_expression
#-------------------------------------------------------------------------------
JOB_NAME=differential_expression
JOB_DEPENDENCIES=$raw_counts_metrics_1_JOB_ID
JOB_DONE=job_output/differential_expression/differential_expression.471a56e71dafee5bb143cd0a2091d23c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'differential_expression.471a56e71dafee5bb143cd0a2091d23c.mugqic.done'
module load mugqic/mugqic_tools/2.1.7 mugqic/R_Bioconductor/3.1.2_3.0 && \
mkdir -p DGE && \
Rscript $R_TOOLS/edger.R \
  -d ../../raw/rna-seq/design.txt \
  -c DGE/rawCountMatrix.csv \
  -o DGE && \
Rscript $R_TOOLS/deseq.R \
  -d ../../raw/rna-seq/design.txt \
  -c DGE/rawCountMatrix.csv \
  -o DGE \
  -l
differential_expression.471a56e71dafee5bb143cd0a2091d23c.mugqic.done
)
differential_expression_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=10:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$differential_expression_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# Call home with pipeline statistics
#-------------------------------------------------------------------------------
wget "http://mugqic.hpc.mcgill.ca/cgi-bin/pipeline.cgi?hostname=lg-1r17-n03&ip=10.241.129.13&pipeline=RnaSeq&steps=picard_sam_to_fastq,trimmomatic,merge_trimmomatic_stats,star,picard_merge_sam_files,picard_sort_sam,picard_mark_duplicates,picard_rna_metrics,estimate_ribosomal_rna,bam_hard_clip,rnaseqc,wiggle,raw_counts,raw_counts_metrics,cufflinks,cuffmerge,cuffquant,cuffdiff,cuffnorm,fpkm_correlation_matrix,gq_seq_utils_exploratory_analysis_rnaseq,differential_expression&samples=12&AnonymizedList=dddd03cb6e65c1183657e5e5bf559f48,444c16ed9974bf0f01c9f5de6751c17b,15b56eae2cdf5203323e847cafca2478,9bb6d10c02e361fb068595ab4db4aed6,b55b936bea7aa8fa3185729524325737,cff800af6cafbbb70c28346717710b25,ba0cee9f1a8fbbefae102e42e845ca2b,277b756c8c3e10902132376d11c4df12,ac8379cdfc0119f8bb100d8925e216fd,6aa81e81b0a21275efefd2698d4b38d7,e8f8bd8dfd005d92c279b40cd9c7e216,da8d2067113079941b7b9b4ded748341" --quiet --output-document=/dev/null

