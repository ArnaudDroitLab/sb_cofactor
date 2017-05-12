#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

#-------------------------------------------------------------------------------
# RnaSeq PBSScheduler Job Submission Bash script
# Version: 2.2.1-beta
# Created on: 2017-04-21T20:09:05
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
#   cuffdiff: 11 jobs
#   cuffnorm: 1 job
#   fpkm_correlation_matrix: 2 jobs
#   gq_seq_utils_exploratory_analysis_rnaseq: 3 jobs
#   differential_expression: 1 job
#   TOTAL: 268 jobs
#-------------------------------------------------------------------------------

OUTPUT_DIR=/gs/scratch/efournier/CofactorHR/A549/output/pipeline
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
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.DEX_nonMamm_1_1.a6ef7e66bfd0ca3c2c9d20a30c40c43a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.DEX_nonMamm_1_1.a6ef7e66bfd0ca3c2c9d20a30c40c43a.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.002.Index_3.SB_A549-DEX-nonMamm.bam \
 FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.002.Index_3.SB_A549-DEX-nonMamm.pair1.fastq.gz \
  SECOND_END_FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.002.Index_3.SB_A549-DEX-nonMamm.pair2.fastq.gz
picard_sam_to_fastq.DEX_nonMamm_1_1.a6ef7e66bfd0ca3c2c9d20a30c40c43a.mugqic.done
)
picard_sam_to_fastq_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m | grep "[0-9]")
echo "$picard_sam_to_fastq_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_2_JOB_ID: picard_sam_to_fastq.DEX_nonMamm_1_2
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.DEX_nonMamm_1_2
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.DEX_nonMamm_1_2.8f68807bf7bdb4b77da4111f144e8e01.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.DEX_nonMamm_1_2.8f68807bf7bdb4b77da4111f144e8e01.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.003.Index_3.SB_A549-DEX-nonMamm.bam \
 FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.003.Index_3.SB_A549-DEX-nonMamm.pair1.fastq.gz \
  SECOND_END_FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.003.Index_3.SB_A549-DEX-nonMamm.pair2.fastq.gz
picard_sam_to_fastq.DEX_nonMamm_1_2.8f68807bf7bdb4b77da4111f144e8e01.mugqic.done
)
picard_sam_to_fastq_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m | grep "[0-9]")
echo "$picard_sam_to_fastq_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_3_JOB_ID: picard_sam_to_fastq.DEX_shMED1_1_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.DEX_shMED1_1_1
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.DEX_shMED1_1_1.2f4951deb94389bfcb13dde9103078c6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.DEX_shMED1_1_1.2f4951deb94389bfcb13dde9103078c6.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.002.Index_22.SB_A549-DEX-shMED1-2.bam \
 FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.002.Index_22.SB_A549-DEX-shMED1-2.pair1.fastq.gz \
  SECOND_END_FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.002.Index_22.SB_A549-DEX-shMED1-2.pair2.fastq.gz
picard_sam_to_fastq.DEX_shMED1_1_1.2f4951deb94389bfcb13dde9103078c6.mugqic.done
)
picard_sam_to_fastq_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m | grep "[0-9]")
echo "$picard_sam_to_fastq_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_4_JOB_ID: picard_sam_to_fastq.DEX_shMED1_1_2
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.DEX_shMED1_1_2
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.DEX_shMED1_1_2.58adc6c83555c4883d2fa674afd7171f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.DEX_shMED1_1_2.58adc6c83555c4883d2fa674afd7171f.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.003.Index_22.SB_A549-DEX-shMED1-2.bam \
 FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.003.Index_22.SB_A549-DEX-shMED1-2.pair1.fastq.gz \
  SECOND_END_FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.003.Index_22.SB_A549-DEX-shMED1-2.pair2.fastq.gz
picard_sam_to_fastq.DEX_shMED1_1_2.58adc6c83555c4883d2fa674afd7171f.mugqic.done
)
picard_sam_to_fastq_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m | grep "[0-9]")
echo "$picard_sam_to_fastq_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_5_JOB_ID: picard_sam_to_fastq.DEX_shNIPBL_1_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.DEX_shNIPBL_1_1
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.DEX_shNIPBL_1_1.6cd754f9d17b17437ae51f27aa86aa28.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.DEX_shNIPBL_1_1.6cd754f9d17b17437ae51f27aa86aa28.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.002.Index_9.SB_A549-DEX-shNIPBL-3.bam \
 FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.002.Index_9.SB_A549-DEX-shNIPBL-3.pair1.fastq.gz \
  SECOND_END_FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.002.Index_9.SB_A549-DEX-shNIPBL-3.pair2.fastq.gz
picard_sam_to_fastq.DEX_shNIPBL_1_1.6cd754f9d17b17437ae51f27aa86aa28.mugqic.done
)
picard_sam_to_fastq_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m | grep "[0-9]")
echo "$picard_sam_to_fastq_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_6_JOB_ID: picard_sam_to_fastq.DEX_shNIPBL_1_2
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.DEX_shNIPBL_1_2
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.DEX_shNIPBL_1_2.1af94fadcfe8fb0597f4d4b3db7fa803.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.DEX_shNIPBL_1_2.1af94fadcfe8fb0597f4d4b3db7fa803.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.003.Index_9.SB_A549-DEX-shNIPBL-3.bam \
 FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.003.Index_9.SB_A549-DEX-shNIPBL-3.pair1.fastq.gz \
  SECOND_END_FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.003.Index_9.SB_A549-DEX-shNIPBL-3.pair2.fastq.gz
picard_sam_to_fastq.DEX_shNIPBL_1_2.1af94fadcfe8fb0597f4d4b3db7fa803.mugqic.done
)
picard_sam_to_fastq_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m | grep "[0-9]")
echo "$picard_sam_to_fastq_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_7_JOB_ID: picard_sam_to_fastq.DEX_shSMC1A_1_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.DEX_shSMC1A_1_1
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.DEX_shSMC1A_1_1.765c66546c9cfb917284b073eba42c84.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.DEX_shSMC1A_1_1.765c66546c9cfb917284b073eba42c84.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.002.Index_27.SB_A549-DEX-shSMC1A-3.bam \
 FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.002.Index_27.SB_A549-DEX-shSMC1A-3.pair1.fastq.gz \
  SECOND_END_FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.002.Index_27.SB_A549-DEX-shSMC1A-3.pair2.fastq.gz
picard_sam_to_fastq.DEX_shSMC1A_1_1.765c66546c9cfb917284b073eba42c84.mugqic.done
)
picard_sam_to_fastq_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m | grep "[0-9]")
echo "$picard_sam_to_fastq_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_8_JOB_ID: picard_sam_to_fastq.DEX_shSMC1A_1_2
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.DEX_shSMC1A_1_2
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.DEX_shSMC1A_1_2.cf5f1f39c1a130a47030d6cc79c749e1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.DEX_shSMC1A_1_2.cf5f1f39c1a130a47030d6cc79c749e1.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.003.Index_27.SB_A549-DEX-shSMC1A-3.bam \
 FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.003.Index_27.SB_A549-DEX-shSMC1A-3.pair1.fastq.gz \
  SECOND_END_FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.003.Index_27.SB_A549-DEX-shSMC1A-3.pair2.fastq.gz
picard_sam_to_fastq.DEX_shSMC1A_1_2.cf5f1f39c1a130a47030d6cc79c749e1.mugqic.done
)
picard_sam_to_fastq_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m | grep "[0-9]")
echo "$picard_sam_to_fastq_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_9_JOB_ID: picard_sam_to_fastq.ETOH_nonMamm_1_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.ETOH_nonMamm_1_1
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.ETOH_nonMamm_1_1.3e7e4d3c5a222c5f8370b2626417ec64.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.ETOH_nonMamm_1_1.3e7e4d3c5a222c5f8370b2626417ec64.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.002.Index_1.SB_A549-nonMamm.bam \
 FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.002.Index_1.SB_A549-nonMamm.pair1.fastq.gz \
  SECOND_END_FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.002.Index_1.SB_A549-nonMamm.pair2.fastq.gz
picard_sam_to_fastq.ETOH_nonMamm_1_1.3e7e4d3c5a222c5f8370b2626417ec64.mugqic.done
)
picard_sam_to_fastq_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m | grep "[0-9]")
echo "$picard_sam_to_fastq_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_10_JOB_ID: picard_sam_to_fastq.ETOH_nonMamm_1_2
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.ETOH_nonMamm_1_2
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.ETOH_nonMamm_1_2.90e9da94a1d5820f1d6f762308e6dcf3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.ETOH_nonMamm_1_2.90e9da94a1d5820f1d6f762308e6dcf3.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.003.Index_1.SB_A549-nonMamm.bam \
 FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.003.Index_1.SB_A549-nonMamm.pair1.fastq.gz \
  SECOND_END_FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.003.Index_1.SB_A549-nonMamm.pair2.fastq.gz
picard_sam_to_fastq.ETOH_nonMamm_1_2.90e9da94a1d5820f1d6f762308e6dcf3.mugqic.done
)
picard_sam_to_fastq_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m | grep "[0-9]")
echo "$picard_sam_to_fastq_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_11_JOB_ID: picard_sam_to_fastq.ETOH_nonMamm_2_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.ETOH_nonMamm_2_1
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.ETOH_nonMamm_2_1.582b478d7e369ccaa8fa3371b955e8e3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.ETOH_nonMamm_2_1.582b478d7e369ccaa8fa3371b955e8e3.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2430.006.Index_10.A549_ETOH_nonMAMM.bam \
 FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2430.006.Index_10.A549_ETOH_nonMAMM.pair1.fastq.gz \
  SECOND_END_FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2430.006.Index_10.A549_ETOH_nonMAMM.pair2.fastq.gz
picard_sam_to_fastq.ETOH_nonMamm_2_1.582b478d7e369ccaa8fa3371b955e8e3.mugqic.done
)
picard_sam_to_fastq_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m | grep "[0-9]")
echo "$picard_sam_to_fastq_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_12_JOB_ID: picard_sam_to_fastq.ETOH_shMED1_1_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.ETOH_shMED1_1_1
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.ETOH_shMED1_1_1.140859b1b2c6cc06034dafeb97a9ff41.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.ETOH_shMED1_1_1.140859b1b2c6cc06034dafeb97a9ff41.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.002.Index_10.SB_A549-shMED1-2.bam \
 FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.002.Index_10.SB_A549-shMED1-2.pair1.fastq.gz \
  SECOND_END_FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.002.Index_10.SB_A549-shMED1-2.pair2.fastq.gz
picard_sam_to_fastq.ETOH_shMED1_1_1.140859b1b2c6cc06034dafeb97a9ff41.mugqic.done
)
picard_sam_to_fastq_12_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m | grep "[0-9]")
echo "$picard_sam_to_fastq_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_13_JOB_ID: picard_sam_to_fastq.ETOH_shMED1_1_2
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.ETOH_shMED1_1_2
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.ETOH_shMED1_1_2.88ef8fe3ffa8f572a96837ecfb461a39.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.ETOH_shMED1_1_2.88ef8fe3ffa8f572a96837ecfb461a39.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.003.Index_10.SB_A549-shMED1-2.bam \
 FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.003.Index_10.SB_A549-shMED1-2.pair1.fastq.gz \
  SECOND_END_FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.003.Index_10.SB_A549-shMED1-2.pair2.fastq.gz
picard_sam_to_fastq.ETOH_shMED1_1_2.88ef8fe3ffa8f572a96837ecfb461a39.mugqic.done
)
picard_sam_to_fastq_13_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m | grep "[0-9]")
echo "$picard_sam_to_fastq_13_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_14_JOB_ID: picard_sam_to_fastq.ETOH_shNIPBL_2_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.ETOH_shNIPBL_2_1
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.ETOH_shNIPBL_2_1.b489d8c64a491352c641924d551609b4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.ETOH_shNIPBL_2_1.b489d8c64a491352c641924d551609b4.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2430.006.Index_11.A549_ETOH_shMED1-2.bam \
 FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2430.006.Index_11.A549_ETOH_shMED1-2.pair1.fastq.gz \
  SECOND_END_FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2430.006.Index_11.A549_ETOH_shMED1-2.pair2.fastq.gz
picard_sam_to_fastq.ETOH_shNIPBL_2_1.b489d8c64a491352c641924d551609b4.mugqic.done
)
picard_sam_to_fastq_14_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m | grep "[0-9]")
echo "$picard_sam_to_fastq_14_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_15_JOB_ID: picard_sam_to_fastq.ETOH_shNIPBL_1_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.ETOH_shNIPBL_1_1
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.ETOH_shNIPBL_1_1.068201d911d6b420ed30405e43d527a9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.ETOH_shNIPBL_1_1.068201d911d6b420ed30405e43d527a9.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.002.Index_8.SB_A549-shNIPBL-3.bam \
 FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.002.Index_8.SB_A549-shNIPBL-3.pair1.fastq.gz \
  SECOND_END_FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.002.Index_8.SB_A549-shNIPBL-3.pair2.fastq.gz
picard_sam_to_fastq.ETOH_shNIPBL_1_1.068201d911d6b420ed30405e43d527a9.mugqic.done
)
picard_sam_to_fastq_15_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m | grep "[0-9]")
echo "$picard_sam_to_fastq_15_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_16_JOB_ID: picard_sam_to_fastq.ETOH_shNIPBL_1_2
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.ETOH_shNIPBL_1_2
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.ETOH_shNIPBL_1_2.ffc86a8b2d53b7b01f53b66b24e96a41.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.ETOH_shNIPBL_1_2.ffc86a8b2d53b7b01f53b66b24e96a41.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.003.Index_8.SB_A549-shNIPBL-3.bam \
 FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.003.Index_8.SB_A549-shNIPBL-3.pair1.fastq.gz \
  SECOND_END_FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.003.Index_8.SB_A549-shNIPBL-3.pair2.fastq.gz
picard_sam_to_fastq.ETOH_shNIPBL_1_2.ffc86a8b2d53b7b01f53b66b24e96a41.mugqic.done
)
picard_sam_to_fastq_16_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m | grep "[0-9]")
echo "$picard_sam_to_fastq_16_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_17_JOB_ID: picard_sam_to_fastq.ETOH_shMED1_2_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.ETOH_shMED1_2_1
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.ETOH_shMED1_2_1.6f49f284ffd57729a570f23c2c4a9e20.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.ETOH_shMED1_2_1.6f49f284ffd57729a570f23c2c4a9e20.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2430.006.Index_20.A549_ETOH_shNIPBL-3.bam \
 FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2430.006.Index_20.A549_ETOH_shNIPBL-3.pair1.fastq.gz \
  SECOND_END_FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2430.006.Index_20.A549_ETOH_shNIPBL-3.pair2.fastq.gz
picard_sam_to_fastq.ETOH_shMED1_2_1.6f49f284ffd57729a570f23c2c4a9e20.mugqic.done
)
picard_sam_to_fastq_17_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m | grep "[0-9]")
echo "$picard_sam_to_fastq_17_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_18_JOB_ID: picard_sam_to_fastq.ETOH_shSMC1A_1_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.ETOH_shSMC1A_1_1
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.ETOH_shSMC1A_1_1.729ede24c09243c00b038d0876d735d2.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.ETOH_shSMC1A_1_1.729ede24c09243c00b038d0876d735d2.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.002.Index_11.SB_A549-shSMC1A-3.bam \
 FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.002.Index_11.SB_A549-shSMC1A-3.pair1.fastq.gz \
  SECOND_END_FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.002.Index_11.SB_A549-shSMC1A-3.pair2.fastq.gz
picard_sam_to_fastq.ETOH_shSMC1A_1_1.729ede24c09243c00b038d0876d735d2.mugqic.done
)
picard_sam_to_fastq_18_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m | grep "[0-9]")
echo "$picard_sam_to_fastq_18_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_19_JOB_ID: picard_sam_to_fastq.ETOH_shSMC1A_1_2
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.ETOH_shSMC1A_1_2
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.ETOH_shSMC1A_1_2.a0019708b61801bf1509d58ea7614afe.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.ETOH_shSMC1A_1_2.a0019708b61801bf1509d58ea7614afe.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.003.Index_11.SB_A549-shSMC1A-3.bam \
 FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.003.Index_11.SB_A549-shSMC1A-3.pair1.fastq.gz \
  SECOND_END_FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.003.Index_11.SB_A549-shSMC1A-3.pair2.fastq.gz
picard_sam_to_fastq.ETOH_shSMC1A_1_2.a0019708b61801bf1509d58ea7614afe.mugqic.done
)
picard_sam_to_fastq_19_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m | grep "[0-9]")
echo "$picard_sam_to_fastq_19_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_20_JOB_ID: picard_sam_to_fastq.ETOH_shSMC1A_2_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.ETOH_shSMC1A_2_1
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.ETOH_shSMC1A_2_1.5636d2d54b9e3e2d7b8a7cbb41e9846f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.ETOH_shSMC1A_2_1.5636d2d54b9e3e2d7b8a7cbb41e9846f.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/picard.jar SamToFastq \
 VALIDATION_STRINGENCY=LENIENT \
 INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2430.006.Index_22.A549_ETOH_shSMC1A-2.bam \
 FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2430.006.Index_22.A549_ETOH_shSMC1A-2.pair1.fastq.gz \
  SECOND_END_FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/HI.2430.006.Index_22.A549_ETOH_shSMC1A-2.pair2.fastq.gz
picard_sam_to_fastq.ETOH_shSMC1A_2_1.5636d2d54b9e3e2d7b8a7cbb41e9846f.mugqic.done
)
picard_sam_to_fastq_20_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m | grep "[0-9]")
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
JOB_DONE=job_output/trimmomatic/trimmomatic.DEX_nonMamm_1_1.6e56ee3d795ce43fd98e41ba6f60a350.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.DEX_nonMamm_1_1.6e56ee3d795ce43fd98e41ba6f60a350.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/DEX_nonMamm_1 && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR PE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.002.Index_3.SB_A549-DEX-nonMamm.pair1.fastq.gz \
  /gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.002.Index_3.SB_A549-DEX-nonMamm.pair2.fastq.gz \
  trim/DEX_nonMamm_1/DEX_nonMamm_1_1.trim.pair1.fastq.gz \
  trim/DEX_nonMamm_1/DEX_nonMamm_1_1.trim.single1.fastq.gz \
  trim/DEX_nonMamm_1/DEX_nonMamm_1_1.trim.pair2.fastq.gz \
  trim/DEX_nonMamm_1/DEX_nonMamm_1_1.trim.single2.fastq.gz \
  ILLUMINACLIP:/cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.0/bfx/adapters-truseq.fa:2:30:15:8:true \
  TRAILING:30 \
  MINLEN:32 \
  2> trim/DEX_nonMamm_1/DEX_nonMamm_1_1.trim.log
trimmomatic.DEX_nonMamm_1_1.6e56ee3d795ce43fd98e41ba6f60a350.mugqic.done
)
trimmomatic_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_2_JOB_ID: trimmomatic.DEX_nonMamm_1_2
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.DEX_nonMamm_1_2
JOB_DEPENDENCIES=$picard_sam_to_fastq_2_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.DEX_nonMamm_1_2.8f6e842a8e0373133b1e1df0d1d0ff17.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.DEX_nonMamm_1_2.8f6e842a8e0373133b1e1df0d1d0ff17.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/DEX_nonMamm_1 && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR PE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.003.Index_3.SB_A549-DEX-nonMamm.pair1.fastq.gz \
  /gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.003.Index_3.SB_A549-DEX-nonMamm.pair2.fastq.gz \
  trim/DEX_nonMamm_1/DEX_nonMamm_1_2.trim.pair1.fastq.gz \
  trim/DEX_nonMamm_1/DEX_nonMamm_1_2.trim.single1.fastq.gz \
  trim/DEX_nonMamm_1/DEX_nonMamm_1_2.trim.pair2.fastq.gz \
  trim/DEX_nonMamm_1/DEX_nonMamm_1_2.trim.single2.fastq.gz \
  ILLUMINACLIP:/cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.0/bfx/adapters-truseq.fa:2:30:15:8:true \
  TRAILING:30 \
  MINLEN:32 \
  2> trim/DEX_nonMamm_1/DEX_nonMamm_1_2.trim.log
trimmomatic.DEX_nonMamm_1_2.8f6e842a8e0373133b1e1df0d1d0ff17.mugqic.done
)
trimmomatic_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_3_JOB_ID: trimmomatic.DEX_shMED1_1_1
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.DEX_shMED1_1_1
JOB_DEPENDENCIES=$picard_sam_to_fastq_3_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.DEX_shMED1_1_1.686374d6156437711be0f31290a3584c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.DEX_shMED1_1_1.686374d6156437711be0f31290a3584c.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/DEX_shMED1_1 && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR PE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.002.Index_22.SB_A549-DEX-shMED1-2.pair1.fastq.gz \
  /gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.002.Index_22.SB_A549-DEX-shMED1-2.pair2.fastq.gz \
  trim/DEX_shMED1_1/DEX_shMED1_1_1.trim.pair1.fastq.gz \
  trim/DEX_shMED1_1/DEX_shMED1_1_1.trim.single1.fastq.gz \
  trim/DEX_shMED1_1/DEX_shMED1_1_1.trim.pair2.fastq.gz \
  trim/DEX_shMED1_1/DEX_shMED1_1_1.trim.single2.fastq.gz \
  ILLUMINACLIP:/cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.0/bfx/adapters-truseq.fa:2:30:15:8:true \
  TRAILING:30 \
  MINLEN:32 \
  2> trim/DEX_shMED1_1/DEX_shMED1_1_1.trim.log
trimmomatic.DEX_shMED1_1_1.686374d6156437711be0f31290a3584c.mugqic.done
)
trimmomatic_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_4_JOB_ID: trimmomatic.DEX_shMED1_1_2
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.DEX_shMED1_1_2
JOB_DEPENDENCIES=$picard_sam_to_fastq_4_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.DEX_shMED1_1_2.919b03b6fb7303530207e857c43d3d2b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.DEX_shMED1_1_2.919b03b6fb7303530207e857c43d3d2b.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/DEX_shMED1_1 && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR PE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.003.Index_22.SB_A549-DEX-shMED1-2.pair1.fastq.gz \
  /gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.003.Index_22.SB_A549-DEX-shMED1-2.pair2.fastq.gz \
  trim/DEX_shMED1_1/DEX_shMED1_1_2.trim.pair1.fastq.gz \
  trim/DEX_shMED1_1/DEX_shMED1_1_2.trim.single1.fastq.gz \
  trim/DEX_shMED1_1/DEX_shMED1_1_2.trim.pair2.fastq.gz \
  trim/DEX_shMED1_1/DEX_shMED1_1_2.trim.single2.fastq.gz \
  ILLUMINACLIP:/cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.0/bfx/adapters-truseq.fa:2:30:15:8:true \
  TRAILING:30 \
  MINLEN:32 \
  2> trim/DEX_shMED1_1/DEX_shMED1_1_2.trim.log
trimmomatic.DEX_shMED1_1_2.919b03b6fb7303530207e857c43d3d2b.mugqic.done
)
trimmomatic_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_5_JOB_ID: trimmomatic.DEX_shNIPBL_1_1
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.DEX_shNIPBL_1_1
JOB_DEPENDENCIES=$picard_sam_to_fastq_5_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.DEX_shNIPBL_1_1.effd510105035e5c5c69eb6ed987e447.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.DEX_shNIPBL_1_1.effd510105035e5c5c69eb6ed987e447.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/DEX_shNIPBL_1 && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR PE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.002.Index_9.SB_A549-DEX-shNIPBL-3.pair1.fastq.gz \
  /gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.002.Index_9.SB_A549-DEX-shNIPBL-3.pair2.fastq.gz \
  trim/DEX_shNIPBL_1/DEX_shNIPBL_1_1.trim.pair1.fastq.gz \
  trim/DEX_shNIPBL_1/DEX_shNIPBL_1_1.trim.single1.fastq.gz \
  trim/DEX_shNIPBL_1/DEX_shNIPBL_1_1.trim.pair2.fastq.gz \
  trim/DEX_shNIPBL_1/DEX_shNIPBL_1_1.trim.single2.fastq.gz \
  ILLUMINACLIP:/cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.0/bfx/adapters-truseq.fa:2:30:15:8:true \
  TRAILING:30 \
  MINLEN:32 \
  2> trim/DEX_shNIPBL_1/DEX_shNIPBL_1_1.trim.log
trimmomatic.DEX_shNIPBL_1_1.effd510105035e5c5c69eb6ed987e447.mugqic.done
)
trimmomatic_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_6_JOB_ID: trimmomatic.DEX_shNIPBL_1_2
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.DEX_shNIPBL_1_2
JOB_DEPENDENCIES=$picard_sam_to_fastq_6_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.DEX_shNIPBL_1_2.5f3e0e7f23c85a90cea01ce01fd33da7.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.DEX_shNIPBL_1_2.5f3e0e7f23c85a90cea01ce01fd33da7.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/DEX_shNIPBL_1 && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR PE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.003.Index_9.SB_A549-DEX-shNIPBL-3.pair1.fastq.gz \
  /gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.003.Index_9.SB_A549-DEX-shNIPBL-3.pair2.fastq.gz \
  trim/DEX_shNIPBL_1/DEX_shNIPBL_1_2.trim.pair1.fastq.gz \
  trim/DEX_shNIPBL_1/DEX_shNIPBL_1_2.trim.single1.fastq.gz \
  trim/DEX_shNIPBL_1/DEX_shNIPBL_1_2.trim.pair2.fastq.gz \
  trim/DEX_shNIPBL_1/DEX_shNIPBL_1_2.trim.single2.fastq.gz \
  ILLUMINACLIP:/cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.0/bfx/adapters-truseq.fa:2:30:15:8:true \
  TRAILING:30 \
  MINLEN:32 \
  2> trim/DEX_shNIPBL_1/DEX_shNIPBL_1_2.trim.log
trimmomatic.DEX_shNIPBL_1_2.5f3e0e7f23c85a90cea01ce01fd33da7.mugqic.done
)
trimmomatic_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_7_JOB_ID: trimmomatic.DEX_shSMC1A_1_1
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.DEX_shSMC1A_1_1
JOB_DEPENDENCIES=$picard_sam_to_fastq_7_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.DEX_shSMC1A_1_1.26442b84152a5e1369427561ca77bdfc.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.DEX_shSMC1A_1_1.26442b84152a5e1369427561ca77bdfc.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/DEX_shSMC1A_1 && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR PE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.002.Index_27.SB_A549-DEX-shSMC1A-3.pair1.fastq.gz \
  /gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.002.Index_27.SB_A549-DEX-shSMC1A-3.pair2.fastq.gz \
  trim/DEX_shSMC1A_1/DEX_shSMC1A_1_1.trim.pair1.fastq.gz \
  trim/DEX_shSMC1A_1/DEX_shSMC1A_1_1.trim.single1.fastq.gz \
  trim/DEX_shSMC1A_1/DEX_shSMC1A_1_1.trim.pair2.fastq.gz \
  trim/DEX_shSMC1A_1/DEX_shSMC1A_1_1.trim.single2.fastq.gz \
  ILLUMINACLIP:/cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.0/bfx/adapters-truseq.fa:2:30:15:8:true \
  TRAILING:30 \
  MINLEN:32 \
  2> trim/DEX_shSMC1A_1/DEX_shSMC1A_1_1.trim.log
trimmomatic.DEX_shSMC1A_1_1.26442b84152a5e1369427561ca77bdfc.mugqic.done
)
trimmomatic_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_8_JOB_ID: trimmomatic.DEX_shSMC1A_1_2
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.DEX_shSMC1A_1_2
JOB_DEPENDENCIES=$picard_sam_to_fastq_8_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.DEX_shSMC1A_1_2.d20636cd5f7462fb39093dc9bf81bb0a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.DEX_shSMC1A_1_2.d20636cd5f7462fb39093dc9bf81bb0a.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/DEX_shSMC1A_1 && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR PE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.003.Index_27.SB_A549-DEX-shSMC1A-3.pair1.fastq.gz \
  /gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.003.Index_27.SB_A549-DEX-shSMC1A-3.pair2.fastq.gz \
  trim/DEX_shSMC1A_1/DEX_shSMC1A_1_2.trim.pair1.fastq.gz \
  trim/DEX_shSMC1A_1/DEX_shSMC1A_1_2.trim.single1.fastq.gz \
  trim/DEX_shSMC1A_1/DEX_shSMC1A_1_2.trim.pair2.fastq.gz \
  trim/DEX_shSMC1A_1/DEX_shSMC1A_1_2.trim.single2.fastq.gz \
  ILLUMINACLIP:/cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.0/bfx/adapters-truseq.fa:2:30:15:8:true \
  TRAILING:30 \
  MINLEN:32 \
  2> trim/DEX_shSMC1A_1/DEX_shSMC1A_1_2.trim.log
trimmomatic.DEX_shSMC1A_1_2.d20636cd5f7462fb39093dc9bf81bb0a.mugqic.done
)
trimmomatic_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_9_JOB_ID: trimmomatic.ETOH_nonMamm_1_1
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.ETOH_nonMamm_1_1
JOB_DEPENDENCIES=$picard_sam_to_fastq_9_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.ETOH_nonMamm_1_1.42171e8f7ef2b2ef60a546495b9c914e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.ETOH_nonMamm_1_1.42171e8f7ef2b2ef60a546495b9c914e.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/ETOH_nonMamm_1 && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR PE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.002.Index_1.SB_A549-nonMamm.pair1.fastq.gz \
  /gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.002.Index_1.SB_A549-nonMamm.pair2.fastq.gz \
  trim/ETOH_nonMamm_1/ETOH_nonMamm_1_1.trim.pair1.fastq.gz \
  trim/ETOH_nonMamm_1/ETOH_nonMamm_1_1.trim.single1.fastq.gz \
  trim/ETOH_nonMamm_1/ETOH_nonMamm_1_1.trim.pair2.fastq.gz \
  trim/ETOH_nonMamm_1/ETOH_nonMamm_1_1.trim.single2.fastq.gz \
  ILLUMINACLIP:/cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.0/bfx/adapters-truseq.fa:2:30:15:8:true \
  TRAILING:30 \
  MINLEN:32 \
  2> trim/ETOH_nonMamm_1/ETOH_nonMamm_1_1.trim.log
trimmomatic.ETOH_nonMamm_1_1.42171e8f7ef2b2ef60a546495b9c914e.mugqic.done
)
trimmomatic_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_10_JOB_ID: trimmomatic.ETOH_nonMamm_1_2
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.ETOH_nonMamm_1_2
JOB_DEPENDENCIES=$picard_sam_to_fastq_10_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.ETOH_nonMamm_1_2.1f0e452d17f5d596f2db0cff88a1d823.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.ETOH_nonMamm_1_2.1f0e452d17f5d596f2db0cff88a1d823.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/ETOH_nonMamm_1 && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR PE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.003.Index_1.SB_A549-nonMamm.pair1.fastq.gz \
  /gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.003.Index_1.SB_A549-nonMamm.pair2.fastq.gz \
  trim/ETOH_nonMamm_1/ETOH_nonMamm_1_2.trim.pair1.fastq.gz \
  trim/ETOH_nonMamm_1/ETOH_nonMamm_1_2.trim.single1.fastq.gz \
  trim/ETOH_nonMamm_1/ETOH_nonMamm_1_2.trim.pair2.fastq.gz \
  trim/ETOH_nonMamm_1/ETOH_nonMamm_1_2.trim.single2.fastq.gz \
  ILLUMINACLIP:/cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.0/bfx/adapters-truseq.fa:2:30:15:8:true \
  TRAILING:30 \
  MINLEN:32 \
  2> trim/ETOH_nonMamm_1/ETOH_nonMamm_1_2.trim.log
trimmomatic.ETOH_nonMamm_1_2.1f0e452d17f5d596f2db0cff88a1d823.mugqic.done
)
trimmomatic_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_11_JOB_ID: trimmomatic.ETOH_nonMamm_2_1
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.ETOH_nonMamm_2_1
JOB_DEPENDENCIES=$picard_sam_to_fastq_11_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.ETOH_nonMamm_2_1.ff1d9e14a971b5b7c53052345d3a9139.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.ETOH_nonMamm_2_1.ff1d9e14a971b5b7c53052345d3a9139.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/ETOH_nonMamm_2 && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR PE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/HI.2430.006.Index_10.A549_ETOH_nonMAMM.pair1.fastq.gz \
  /gs/scratch/efournier/CofactorHR/A549/raw/HI.2430.006.Index_10.A549_ETOH_nonMAMM.pair2.fastq.gz \
  trim/ETOH_nonMamm_2/ETOH_nonMamm_2_1.trim.pair1.fastq.gz \
  trim/ETOH_nonMamm_2/ETOH_nonMamm_2_1.trim.single1.fastq.gz \
  trim/ETOH_nonMamm_2/ETOH_nonMamm_2_1.trim.pair2.fastq.gz \
  trim/ETOH_nonMamm_2/ETOH_nonMamm_2_1.trim.single2.fastq.gz \
  ILLUMINACLIP:/cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.0/bfx/adapters-truseq.fa:2:30:15:8:true \
  TRAILING:30 \
  MINLEN:32 \
  2> trim/ETOH_nonMamm_2/ETOH_nonMamm_2_1.trim.log
trimmomatic.ETOH_nonMamm_2_1.ff1d9e14a971b5b7c53052345d3a9139.mugqic.done
)
trimmomatic_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_12_JOB_ID: trimmomatic.ETOH_shMED1_1_1
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.ETOH_shMED1_1_1
JOB_DEPENDENCIES=$picard_sam_to_fastq_12_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.ETOH_shMED1_1_1.0a2eef97f909509d1a288ed38946e302.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.ETOH_shMED1_1_1.0a2eef97f909509d1a288ed38946e302.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/ETOH_shMED1_1 && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR PE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.002.Index_10.SB_A549-shMED1-2.pair1.fastq.gz \
  /gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.002.Index_10.SB_A549-shMED1-2.pair2.fastq.gz \
  trim/ETOH_shMED1_1/ETOH_shMED1_1_1.trim.pair1.fastq.gz \
  trim/ETOH_shMED1_1/ETOH_shMED1_1_1.trim.single1.fastq.gz \
  trim/ETOH_shMED1_1/ETOH_shMED1_1_1.trim.pair2.fastq.gz \
  trim/ETOH_shMED1_1/ETOH_shMED1_1_1.trim.single2.fastq.gz \
  ILLUMINACLIP:/cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.0/bfx/adapters-truseq.fa:2:30:15:8:true \
  TRAILING:30 \
  MINLEN:32 \
  2> trim/ETOH_shMED1_1/ETOH_shMED1_1_1.trim.log
trimmomatic.ETOH_shMED1_1_1.0a2eef97f909509d1a288ed38946e302.mugqic.done
)
trimmomatic_12_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_13_JOB_ID: trimmomatic.ETOH_shMED1_1_2
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.ETOH_shMED1_1_2
JOB_DEPENDENCIES=$picard_sam_to_fastq_13_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.ETOH_shMED1_1_2.3e5487e7a47f656b9cae0cbb94198b6f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.ETOH_shMED1_1_2.3e5487e7a47f656b9cae0cbb94198b6f.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/ETOH_shMED1_1 && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR PE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.003.Index_10.SB_A549-shMED1-2.pair1.fastq.gz \
  /gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.003.Index_10.SB_A549-shMED1-2.pair2.fastq.gz \
  trim/ETOH_shMED1_1/ETOH_shMED1_1_2.trim.pair1.fastq.gz \
  trim/ETOH_shMED1_1/ETOH_shMED1_1_2.trim.single1.fastq.gz \
  trim/ETOH_shMED1_1/ETOH_shMED1_1_2.trim.pair2.fastq.gz \
  trim/ETOH_shMED1_1/ETOH_shMED1_1_2.trim.single2.fastq.gz \
  ILLUMINACLIP:/cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.0/bfx/adapters-truseq.fa:2:30:15:8:true \
  TRAILING:30 \
  MINLEN:32 \
  2> trim/ETOH_shMED1_1/ETOH_shMED1_1_2.trim.log
trimmomatic.ETOH_shMED1_1_2.3e5487e7a47f656b9cae0cbb94198b6f.mugqic.done
)
trimmomatic_13_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_13_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_14_JOB_ID: trimmomatic.ETOH_shNIPBL_2_1
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.ETOH_shNIPBL_2_1
JOB_DEPENDENCIES=$picard_sam_to_fastq_14_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.ETOH_shNIPBL_2_1.c655a0ae5be4f8b497d10ed6032a38d5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.ETOH_shNIPBL_2_1.c655a0ae5be4f8b497d10ed6032a38d5.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/ETOH_shNIPBL_2 && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR PE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/HI.2430.006.Index_11.A549_ETOH_shMED1-2.pair1.fastq.gz \
  /gs/scratch/efournier/CofactorHR/A549/raw/HI.2430.006.Index_11.A549_ETOH_shMED1-2.pair2.fastq.gz \
  trim/ETOH_shNIPBL_2/ETOH_shNIPBL_2_1.trim.pair1.fastq.gz \
  trim/ETOH_shNIPBL_2/ETOH_shNIPBL_2_1.trim.single1.fastq.gz \
  trim/ETOH_shNIPBL_2/ETOH_shNIPBL_2_1.trim.pair2.fastq.gz \
  trim/ETOH_shNIPBL_2/ETOH_shNIPBL_2_1.trim.single2.fastq.gz \
  ILLUMINACLIP:/cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.0/bfx/adapters-truseq.fa:2:30:15:8:true \
  TRAILING:30 \
  MINLEN:32 \
  2> trim/ETOH_shNIPBL_2/ETOH_shNIPBL_2_1.trim.log
trimmomatic.ETOH_shNIPBL_2_1.c655a0ae5be4f8b497d10ed6032a38d5.mugqic.done
)
trimmomatic_14_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_14_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_15_JOB_ID: trimmomatic.ETOH_shNIPBL_1_1
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.ETOH_shNIPBL_1_1
JOB_DEPENDENCIES=$picard_sam_to_fastq_15_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.ETOH_shNIPBL_1_1.604f5cf3467eef887c9681ab75066770.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.ETOH_shNIPBL_1_1.604f5cf3467eef887c9681ab75066770.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/ETOH_shNIPBL_1 && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR PE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.002.Index_8.SB_A549-shNIPBL-3.pair1.fastq.gz \
  /gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.002.Index_8.SB_A549-shNIPBL-3.pair2.fastq.gz \
  trim/ETOH_shNIPBL_1/ETOH_shNIPBL_1_1.trim.pair1.fastq.gz \
  trim/ETOH_shNIPBL_1/ETOH_shNIPBL_1_1.trim.single1.fastq.gz \
  trim/ETOH_shNIPBL_1/ETOH_shNIPBL_1_1.trim.pair2.fastq.gz \
  trim/ETOH_shNIPBL_1/ETOH_shNIPBL_1_1.trim.single2.fastq.gz \
  ILLUMINACLIP:/cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.0/bfx/adapters-truseq.fa:2:30:15:8:true \
  TRAILING:30 \
  MINLEN:32 \
  2> trim/ETOH_shNIPBL_1/ETOH_shNIPBL_1_1.trim.log
trimmomatic.ETOH_shNIPBL_1_1.604f5cf3467eef887c9681ab75066770.mugqic.done
)
trimmomatic_15_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_15_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_16_JOB_ID: trimmomatic.ETOH_shNIPBL_1_2
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.ETOH_shNIPBL_1_2
JOB_DEPENDENCIES=$picard_sam_to_fastq_16_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.ETOH_shNIPBL_1_2.ad6aeb3d07e4c3e9f922aebf27420656.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.ETOH_shNIPBL_1_2.ad6aeb3d07e4c3e9f922aebf27420656.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/ETOH_shNIPBL_1 && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR PE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.003.Index_8.SB_A549-shNIPBL-3.pair1.fastq.gz \
  /gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.003.Index_8.SB_A549-shNIPBL-3.pair2.fastq.gz \
  trim/ETOH_shNIPBL_1/ETOH_shNIPBL_1_2.trim.pair1.fastq.gz \
  trim/ETOH_shNIPBL_1/ETOH_shNIPBL_1_2.trim.single1.fastq.gz \
  trim/ETOH_shNIPBL_1/ETOH_shNIPBL_1_2.trim.pair2.fastq.gz \
  trim/ETOH_shNIPBL_1/ETOH_shNIPBL_1_2.trim.single2.fastq.gz \
  ILLUMINACLIP:/cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.0/bfx/adapters-truseq.fa:2:30:15:8:true \
  TRAILING:30 \
  MINLEN:32 \
  2> trim/ETOH_shNIPBL_1/ETOH_shNIPBL_1_2.trim.log
trimmomatic.ETOH_shNIPBL_1_2.ad6aeb3d07e4c3e9f922aebf27420656.mugqic.done
)
trimmomatic_16_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_16_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_17_JOB_ID: trimmomatic.ETOH_shMED1_2_1
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.ETOH_shMED1_2_1
JOB_DEPENDENCIES=$picard_sam_to_fastq_17_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.ETOH_shMED1_2_1.b3c996fe2e075a50a7743edae83d671f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.ETOH_shMED1_2_1.b3c996fe2e075a50a7743edae83d671f.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/ETOH_shMED1_2 && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR PE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/HI.2430.006.Index_20.A549_ETOH_shNIPBL-3.pair1.fastq.gz \
  /gs/scratch/efournier/CofactorHR/A549/raw/HI.2430.006.Index_20.A549_ETOH_shNIPBL-3.pair2.fastq.gz \
  trim/ETOH_shMED1_2/ETOH_shMED1_2_1.trim.pair1.fastq.gz \
  trim/ETOH_shMED1_2/ETOH_shMED1_2_1.trim.single1.fastq.gz \
  trim/ETOH_shMED1_2/ETOH_shMED1_2_1.trim.pair2.fastq.gz \
  trim/ETOH_shMED1_2/ETOH_shMED1_2_1.trim.single2.fastq.gz \
  ILLUMINACLIP:/cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.0/bfx/adapters-truseq.fa:2:30:15:8:true \
  TRAILING:30 \
  MINLEN:32 \
  2> trim/ETOH_shMED1_2/ETOH_shMED1_2_1.trim.log
trimmomatic.ETOH_shMED1_2_1.b3c996fe2e075a50a7743edae83d671f.mugqic.done
)
trimmomatic_17_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_17_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_18_JOB_ID: trimmomatic.ETOH_shSMC1A_1_1
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.ETOH_shSMC1A_1_1
JOB_DEPENDENCIES=$picard_sam_to_fastq_18_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.ETOH_shSMC1A_1_1.5024423d42ea13ae746f3769b0e62014.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.ETOH_shSMC1A_1_1.5024423d42ea13ae746f3769b0e62014.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/ETOH_shSMC1A_1 && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR PE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.002.Index_11.SB_A549-shSMC1A-3.pair1.fastq.gz \
  /gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.002.Index_11.SB_A549-shSMC1A-3.pair2.fastq.gz \
  trim/ETOH_shSMC1A_1/ETOH_shSMC1A_1_1.trim.pair1.fastq.gz \
  trim/ETOH_shSMC1A_1/ETOH_shSMC1A_1_1.trim.single1.fastq.gz \
  trim/ETOH_shSMC1A_1/ETOH_shSMC1A_1_1.trim.pair2.fastq.gz \
  trim/ETOH_shSMC1A_1/ETOH_shSMC1A_1_1.trim.single2.fastq.gz \
  ILLUMINACLIP:/cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.0/bfx/adapters-truseq.fa:2:30:15:8:true \
  TRAILING:30 \
  MINLEN:32 \
  2> trim/ETOH_shSMC1A_1/ETOH_shSMC1A_1_1.trim.log
trimmomatic.ETOH_shSMC1A_1_1.5024423d42ea13ae746f3769b0e62014.mugqic.done
)
trimmomatic_18_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_18_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_19_JOB_ID: trimmomatic.ETOH_shSMC1A_1_2
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.ETOH_shSMC1A_1_2
JOB_DEPENDENCIES=$picard_sam_to_fastq_19_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.ETOH_shSMC1A_1_2.6d76daf6306f4e286f65d99a432a3256.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.ETOH_shSMC1A_1_2.6d76daf6306f4e286f65d99a432a3256.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/ETOH_shSMC1A_1 && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR PE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.003.Index_11.SB_A549-shSMC1A-3.pair1.fastq.gz \
  /gs/scratch/efournier/CofactorHR/A549/raw/HI.2031.003.Index_11.SB_A549-shSMC1A-3.pair2.fastq.gz \
  trim/ETOH_shSMC1A_1/ETOH_shSMC1A_1_2.trim.pair1.fastq.gz \
  trim/ETOH_shSMC1A_1/ETOH_shSMC1A_1_2.trim.single1.fastq.gz \
  trim/ETOH_shSMC1A_1/ETOH_shSMC1A_1_2.trim.pair2.fastq.gz \
  trim/ETOH_shSMC1A_1/ETOH_shSMC1A_1_2.trim.single2.fastq.gz \
  ILLUMINACLIP:/cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.0/bfx/adapters-truseq.fa:2:30:15:8:true \
  TRAILING:30 \
  MINLEN:32 \
  2> trim/ETOH_shSMC1A_1/ETOH_shSMC1A_1_2.trim.log
trimmomatic.ETOH_shSMC1A_1_2.6d76daf6306f4e286f65d99a432a3256.mugqic.done
)
trimmomatic_19_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_19_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_20_JOB_ID: trimmomatic.ETOH_shSMC1A_2_1
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.ETOH_shSMC1A_2_1
JOB_DEPENDENCIES=$picard_sam_to_fastq_20_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.ETOH_shSMC1A_2_1.03c725f43c6b63ed29bdd5098c41cf96.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.ETOH_shSMC1A_2_1.03c725f43c6b63ed29bdd5098c41cf96.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.36 && \
mkdir -p trim/ETOH_shSMC1A_2 && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR PE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/HI.2430.006.Index_22.A549_ETOH_shSMC1A-2.pair1.fastq.gz \
  /gs/scratch/efournier/CofactorHR/A549/raw/HI.2430.006.Index_22.A549_ETOH_shSMC1A-2.pair2.fastq.gz \
  trim/ETOH_shSMC1A_2/ETOH_shSMC1A_2_1.trim.pair1.fastq.gz \
  trim/ETOH_shSMC1A_2/ETOH_shSMC1A_2_1.trim.single1.fastq.gz \
  trim/ETOH_shSMC1A_2/ETOH_shSMC1A_2_1.trim.pair2.fastq.gz \
  trim/ETOH_shSMC1A_2/ETOH_shSMC1A_2_1.trim.single2.fastq.gz \
  ILLUMINACLIP:/cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.1.0/bfx/adapters-truseq.fa:2:30:15:8:true \
  TRAILING:30 \
  MINLEN:32 \
  2> trim/ETOH_shSMC1A_2/ETOH_shSMC1A_2_1.trim.log
trimmomatic.ETOH_shSMC1A_2_1.03c725f43c6b63ed29bdd5098c41cf96.mugqic.done
)
trimmomatic_20_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
JOB_DONE=job_output/merge_trimmomatic_stats/merge_trimmomatic_stats.c86c0871420cc36358cd01ab8fab5a99.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'merge_trimmomatic_stats.c86c0871420cc36358cd01ab8fab5a99.mugqic.done'
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
grep ^Input trim/DEX_shSMC1A_1/DEX_shSMC1A_1_1.trim.log | \
perl -pe 's/^Input Read Pairs: (\d+).*Both Surviving: (\d+).*Forward Only Surviving: (\d+).*$/DEX_shSMC1A_1	DEX_shSMC1A_1_1	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/DEX_shSMC1A_1/DEX_shSMC1A_1_2.trim.log | \
perl -pe 's/^Input Read Pairs: (\d+).*Both Surviving: (\d+).*Forward Only Surviving: (\d+).*$/DEX_shSMC1A_1	DEX_shSMC1A_1_2	\1	\2/' | \
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
grep ^Input trim/ETOH_shSMC1A_1/ETOH_shSMC1A_1_1.trim.log | \
perl -pe 's/^Input Read Pairs: (\d+).*Both Surviving: (\d+).*Forward Only Surviving: (\d+).*$/ETOH_shSMC1A_1	ETOH_shSMC1A_1_1	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/ETOH_shSMC1A_1/ETOH_shSMC1A_1_2.trim.log | \
perl -pe 's/^Input Read Pairs: (\d+).*Both Surviving: (\d+).*Forward Only Surviving: (\d+).*$/ETOH_shSMC1A_1	ETOH_shSMC1A_1_2	\1	\2/' | \
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
merge_trimmomatic_stats.c86c0871420cc36358cd01ab8fab5a99.mugqic.done
)
merge_trimmomatic_stats_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
JOB_DONE=job_output/star/star_align.1.DEX_nonMamm_1_1.d22d684e9e3235bc35c0985dd5a73c3f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.1.DEX_nonMamm_1_1.d22d684e9e3235bc35c0985dd5a73c3f.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment_1stPass/DEX_nonMamm_1/DEX_nonMamm_1_1 && \
STAR --runMode alignReads \
  --genomeDir /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/star_index/UCSC2009-03-08.sjdbOverhang99 \
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
star_align.1.DEX_nonMamm_1_1.d22d684e9e3235bc35c0985dd5a73c3f.mugqic.done
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
JOB_DONE=job_output/star/star_align.1.DEX_nonMamm_1_2.3efdc45f11d655546ed64f4506e3af92.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.1.DEX_nonMamm_1_2.3efdc45f11d655546ed64f4506e3af92.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment_1stPass/DEX_nonMamm_1/DEX_nonMamm_1_2 && \
STAR --runMode alignReads \
  --genomeDir /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/star_index/UCSC2009-03-08.sjdbOverhang99 \
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
star_align.1.DEX_nonMamm_1_2.3efdc45f11d655546ed64f4506e3af92.mugqic.done
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
JOB_DONE=job_output/star/star_align.1.DEX_shMED1_1_1.1d9e06a300facc1a0c6367e840e068cb.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.1.DEX_shMED1_1_1.1d9e06a300facc1a0c6367e840e068cb.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment_1stPass/DEX_shMED1_1/DEX_shMED1_1_1 && \
STAR --runMode alignReads \
  --genomeDir /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/star_index/UCSC2009-03-08.sjdbOverhang99 \
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
star_align.1.DEX_shMED1_1_1.1d9e06a300facc1a0c6367e840e068cb.mugqic.done
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
JOB_DONE=job_output/star/star_align.1.DEX_shMED1_1_2.04be0c0399a20243fa565b6148b67555.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.1.DEX_shMED1_1_2.04be0c0399a20243fa565b6148b67555.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment_1stPass/DEX_shMED1_1/DEX_shMED1_1_2 && \
STAR --runMode alignReads \
  --genomeDir /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/star_index/UCSC2009-03-08.sjdbOverhang99 \
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
star_align.1.DEX_shMED1_1_2.04be0c0399a20243fa565b6148b67555.mugqic.done
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
JOB_DONE=job_output/star/star_align.1.DEX_shNIPBL_1_1.c74325805efa72ebce0d6939e8f70ac1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.1.DEX_shNIPBL_1_1.c74325805efa72ebce0d6939e8f70ac1.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment_1stPass/DEX_shNIPBL_1/DEX_shNIPBL_1_1 && \
STAR --runMode alignReads \
  --genomeDir /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/star_index/UCSC2009-03-08.sjdbOverhang99 \
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
star_align.1.DEX_shNIPBL_1_1.c74325805efa72ebce0d6939e8f70ac1.mugqic.done
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
JOB_DONE=job_output/star/star_align.1.DEX_shNIPBL_1_2.3da7bd06ee7f3fd4e58fabb703a4cd23.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.1.DEX_shNIPBL_1_2.3da7bd06ee7f3fd4e58fabb703a4cd23.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment_1stPass/DEX_shNIPBL_1/DEX_shNIPBL_1_2 && \
STAR --runMode alignReads \
  --genomeDir /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/star_index/UCSC2009-03-08.sjdbOverhang99 \
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
star_align.1.DEX_shNIPBL_1_2.3da7bd06ee7f3fd4e58fabb703a4cd23.mugqic.done
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
JOB_DONE=job_output/star/star_align.1.DEX_shSMC1A_1_1.e7608c9c9946926a558e182596f4487a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.1.DEX_shSMC1A_1_1.e7608c9c9946926a558e182596f4487a.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment_1stPass/DEX_shSMC1A_1/DEX_shSMC1A_1_1 && \
STAR --runMode alignReads \
  --genomeDir /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/star_index/UCSC2009-03-08.sjdbOverhang99 \
  --readFilesIn \
    trim/DEX_shSMC1A_1/DEX_shSMC1A_1_1.trim.pair1.fastq.gz \
    trim/DEX_shSMC1A_1/DEX_shSMC1A_1_1.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix alignment_1stPass/DEX_shSMC1A_1/DEX_shSMC1A_1_1/ \
  --outSAMattrRGline ID:"DEX_shSMC1A_1_1" 	PL:"ILLUMINA" 			SM:"DEX_shSMC1A_1" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 70000000000 \
  --limitIObufferSize 1000000000
star_align.1.DEX_shSMC1A_1_1.e7608c9c9946926a558e182596f4487a.mugqic.done
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
JOB_DONE=job_output/star/star_align.1.DEX_shSMC1A_1_2.561501d76a768d75116a553429f1bd96.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.1.DEX_shSMC1A_1_2.561501d76a768d75116a553429f1bd96.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment_1stPass/DEX_shSMC1A_1/DEX_shSMC1A_1_2 && \
STAR --runMode alignReads \
  --genomeDir /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/star_index/UCSC2009-03-08.sjdbOverhang99 \
  --readFilesIn \
    trim/DEX_shSMC1A_1/DEX_shSMC1A_1_2.trim.pair1.fastq.gz \
    trim/DEX_shSMC1A_1/DEX_shSMC1A_1_2.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix alignment_1stPass/DEX_shSMC1A_1/DEX_shSMC1A_1_2/ \
  --outSAMattrRGline ID:"DEX_shSMC1A_1_2" 	PL:"ILLUMINA" 			SM:"DEX_shSMC1A_1" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 70000000000 \
  --limitIObufferSize 1000000000
star_align.1.DEX_shSMC1A_1_2.561501d76a768d75116a553429f1bd96.mugqic.done
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
JOB_DONE=job_output/star/star_align.1.ETOH_nonMamm_1_1.c499371783b1e469ceb7f96401568c21.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.1.ETOH_nonMamm_1_1.c499371783b1e469ceb7f96401568c21.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment_1stPass/ETOH_nonMamm_1/ETOH_nonMamm_1_1 && \
STAR --runMode alignReads \
  --genomeDir /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/star_index/UCSC2009-03-08.sjdbOverhang99 \
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
star_align.1.ETOH_nonMamm_1_1.c499371783b1e469ceb7f96401568c21.mugqic.done
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
JOB_DONE=job_output/star/star_align.1.ETOH_nonMamm_1_2.0a8d115398b24f4f500e7dbb8c0cc099.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.1.ETOH_nonMamm_1_2.0a8d115398b24f4f500e7dbb8c0cc099.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment_1stPass/ETOH_nonMamm_1/ETOH_nonMamm_1_2 && \
STAR --runMode alignReads \
  --genomeDir /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/star_index/UCSC2009-03-08.sjdbOverhang99 \
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
star_align.1.ETOH_nonMamm_1_2.0a8d115398b24f4f500e7dbb8c0cc099.mugqic.done
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
JOB_DONE=job_output/star/star_align.1.ETOH_nonMamm_2_1.39d99292db9df9cb39507e0d7a619975.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.1.ETOH_nonMamm_2_1.39d99292db9df9cb39507e0d7a619975.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment_1stPass/ETOH_nonMamm_2/ETOH_nonMamm_2_1 && \
STAR --runMode alignReads \
  --genomeDir /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/star_index/UCSC2009-03-08.sjdbOverhang99 \
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
star_align.1.ETOH_nonMamm_2_1.39d99292db9df9cb39507e0d7a619975.mugqic.done
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
JOB_DONE=job_output/star/star_align.1.ETOH_shMED1_1_1.a26c57a97a1b1ed481eb4bb6dc276cf1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.1.ETOH_shMED1_1_1.a26c57a97a1b1ed481eb4bb6dc276cf1.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment_1stPass/ETOH_shMED1_1/ETOH_shMED1_1_1 && \
STAR --runMode alignReads \
  --genomeDir /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/star_index/UCSC2009-03-08.sjdbOverhang99 \
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
star_align.1.ETOH_shMED1_1_1.a26c57a97a1b1ed481eb4bb6dc276cf1.mugqic.done
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
JOB_DONE=job_output/star/star_align.1.ETOH_shMED1_1_2.39888ac3b587369abb730ca38d45a5e6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.1.ETOH_shMED1_1_2.39888ac3b587369abb730ca38d45a5e6.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment_1stPass/ETOH_shMED1_1/ETOH_shMED1_1_2 && \
STAR --runMode alignReads \
  --genomeDir /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/star_index/UCSC2009-03-08.sjdbOverhang99 \
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
star_align.1.ETOH_shMED1_1_2.39888ac3b587369abb730ca38d45a5e6.mugqic.done
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
JOB_DONE=job_output/star/star_align.1.ETOH_shNIPBL_2_1.8fec688cd08c9996ed4fc33d764689b6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.1.ETOH_shNIPBL_2_1.8fec688cd08c9996ed4fc33d764689b6.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment_1stPass/ETOH_shNIPBL_2/ETOH_shNIPBL_2_1 && \
STAR --runMode alignReads \
  --genomeDir /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/star_index/UCSC2009-03-08.sjdbOverhang99 \
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
star_align.1.ETOH_shNIPBL_2_1.8fec688cd08c9996ed4fc33d764689b6.mugqic.done
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
JOB_DONE=job_output/star/star_align.1.ETOH_shNIPBL_1_1.3ca02f7e038f898d53db1f4e4ec96794.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.1.ETOH_shNIPBL_1_1.3ca02f7e038f898d53db1f4e4ec96794.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment_1stPass/ETOH_shNIPBL_1/ETOH_shNIPBL_1_1 && \
STAR --runMode alignReads \
  --genomeDir /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/star_index/UCSC2009-03-08.sjdbOverhang99 \
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
star_align.1.ETOH_shNIPBL_1_1.3ca02f7e038f898d53db1f4e4ec96794.mugqic.done
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
JOB_DONE=job_output/star/star_align.1.ETOH_shNIPBL_1_2.aa4266547969084816c19d941d4cfe68.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.1.ETOH_shNIPBL_1_2.aa4266547969084816c19d941d4cfe68.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment_1stPass/ETOH_shNIPBL_1/ETOH_shNIPBL_1_2 && \
STAR --runMode alignReads \
  --genomeDir /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/star_index/UCSC2009-03-08.sjdbOverhang99 \
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
star_align.1.ETOH_shNIPBL_1_2.aa4266547969084816c19d941d4cfe68.mugqic.done
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
JOB_DONE=job_output/star/star_align.1.ETOH_shMED1_2_1.dc0404a969c178a9d3aa93a1627a540f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.1.ETOH_shMED1_2_1.dc0404a969c178a9d3aa93a1627a540f.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment_1stPass/ETOH_shMED1_2/ETOH_shMED1_2_1 && \
STAR --runMode alignReads \
  --genomeDir /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/star_index/UCSC2009-03-08.sjdbOverhang99 \
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
star_align.1.ETOH_shMED1_2_1.dc0404a969c178a9d3aa93a1627a540f.mugqic.done
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
JOB_DONE=job_output/star/star_align.1.ETOH_shSMC1A_1_1.7d438a063998d200485c7ae4a2c47a18.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.1.ETOH_shSMC1A_1_1.7d438a063998d200485c7ae4a2c47a18.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment_1stPass/ETOH_shSMC1A_1/ETOH_shSMC1A_1_1 && \
STAR --runMode alignReads \
  --genomeDir /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/star_index/UCSC2009-03-08.sjdbOverhang99 \
  --readFilesIn \
    trim/ETOH_shSMC1A_1/ETOH_shSMC1A_1_1.trim.pair1.fastq.gz \
    trim/ETOH_shSMC1A_1/ETOH_shSMC1A_1_1.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix alignment_1stPass/ETOH_shSMC1A_1/ETOH_shSMC1A_1_1/ \
  --outSAMattrRGline ID:"ETOH_shSMC1A_1_1" 	PL:"ILLUMINA" 			SM:"ETOH_shSMC1A_1" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 70000000000 \
  --limitIObufferSize 1000000000
star_align.1.ETOH_shSMC1A_1_1.7d438a063998d200485c7ae4a2c47a18.mugqic.done
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
JOB_DONE=job_output/star/star_align.1.ETOH_shSMC1A_1_2.0b66c0d4800af13707986a90c37f1eb0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.1.ETOH_shSMC1A_1_2.0b66c0d4800af13707986a90c37f1eb0.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment_1stPass/ETOH_shSMC1A_1/ETOH_shSMC1A_1_2 && \
STAR --runMode alignReads \
  --genomeDir /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/star_index/UCSC2009-03-08.sjdbOverhang99 \
  --readFilesIn \
    trim/ETOH_shSMC1A_1/ETOH_shSMC1A_1_2.trim.pair1.fastq.gz \
    trim/ETOH_shSMC1A_1/ETOH_shSMC1A_1_2.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix alignment_1stPass/ETOH_shSMC1A_1/ETOH_shSMC1A_1_2/ \
  --outSAMattrRGline ID:"ETOH_shSMC1A_1_2" 	PL:"ILLUMINA" 			SM:"ETOH_shSMC1A_1" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 70000000000 \
  --limitIObufferSize 1000000000
star_align.1.ETOH_shSMC1A_1_2.0b66c0d4800af13707986a90c37f1eb0.mugqic.done
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
JOB_DONE=job_output/star/star_align.1.ETOH_shSMC1A_2_1.38698d80682f49e8c9a7a73984974981.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.1.ETOH_shSMC1A_2_1.38698d80682f49e8c9a7a73984974981.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment_1stPass/ETOH_shSMC1A_2/ETOH_shSMC1A_2_1 && \
STAR --runMode alignReads \
  --genomeDir /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/star_index/UCSC2009-03-08.sjdbOverhang99 \
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
star_align.1.ETOH_shSMC1A_2_1.38698d80682f49e8c9a7a73984974981.mugqic.done
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
JOB_DONE=job_output/star/star_index.AllSamples.7044bb24fbc5acd6d3de84474301f5a1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_index.AllSamples.7044bb24fbc5acd6d3de84474301f5a1.mugqic.done'
module load mugqic/star/2.5.1b && \
cat \
  alignment_1stPass/DEX_nonMamm_1/DEX_nonMamm_1_1/SJ.out.tab \
  alignment_1stPass/DEX_nonMamm_1/DEX_nonMamm_1_2/SJ.out.tab \
  alignment_1stPass/DEX_shMED1_1/DEX_shMED1_1_1/SJ.out.tab \
  alignment_1stPass/DEX_shMED1_1/DEX_shMED1_1_2/SJ.out.tab \
  alignment_1stPass/DEX_shNIPBL_1/DEX_shNIPBL_1_1/SJ.out.tab \
  alignment_1stPass/DEX_shNIPBL_1/DEX_shNIPBL_1_2/SJ.out.tab \
  alignment_1stPass/DEX_shSMC1A_1/DEX_shSMC1A_1_1/SJ.out.tab \
  alignment_1stPass/DEX_shSMC1A_1/DEX_shSMC1A_1_2/SJ.out.tab \
  alignment_1stPass/ETOH_nonMamm_1/ETOH_nonMamm_1_1/SJ.out.tab \
  alignment_1stPass/ETOH_nonMamm_1/ETOH_nonMamm_1_2/SJ.out.tab \
  alignment_1stPass/ETOH_nonMamm_2/ETOH_nonMamm_2_1/SJ.out.tab \
  alignment_1stPass/ETOH_shMED1_1/ETOH_shMED1_1_1/SJ.out.tab \
  alignment_1stPass/ETOH_shMED1_1/ETOH_shMED1_1_2/SJ.out.tab \
  alignment_1stPass/ETOH_shNIPBL_2/ETOH_shNIPBL_2_1/SJ.out.tab \
  alignment_1stPass/ETOH_shNIPBL_1/ETOH_shNIPBL_1_1/SJ.out.tab \
  alignment_1stPass/ETOH_shNIPBL_1/ETOH_shNIPBL_1_2/SJ.out.tab \
  alignment_1stPass/ETOH_shMED1_2/ETOH_shMED1_2_1/SJ.out.tab \
  alignment_1stPass/ETOH_shSMC1A_1/ETOH_shSMC1A_1_1/SJ.out.tab \
  alignment_1stPass/ETOH_shSMC1A_1/ETOH_shSMC1A_1_2/SJ.out.tab \
  alignment_1stPass/ETOH_shSMC1A_2/ETOH_shSMC1A_2_1/SJ.out.tab | \
awk 'BEGIN {OFS="	"; strChar[0]="."; strChar[1]="+"; strChar[2]="-"} {if($5>0){print $1,$2,$3,strChar[$4]}}' | sort -k1,1h -k2,2n > alignment_1stPass/AllSamples.SJ.out.tab && \
mkdir -p reference.Merged && \
STAR --runMode genomeGenerate \
  --genomeDir reference.Merged \
  --genomeFastaFiles /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa \
  --runThreadN 12 \
  --limitGenomeGenerateRAM 50000000000 \
  --sjdbFileChrStartEnd alignment_1stPass/AllSamples.SJ.out.tab \
  --limitIObufferSize 1000000000 \
  --sjdbOverhang 99
star_index.AllSamples.7044bb24fbc5acd6d3de84474301f5a1.mugqic.done
)
star_21_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q lm -l nodes=1:ppn=12 -l pmem=5700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
JOB_DONE=job_output/star/star_align.2.DEX_shSMC1A_1_1.0ec6d5680146a7cb90da2f087367f168.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.2.DEX_shSMC1A_1_1.0ec6d5680146a7cb90da2f087367f168.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment/DEX_shSMC1A_1/DEX_shSMC1A_1_1 && \
STAR --runMode alignReads \
  --genomeDir reference.Merged \
  --readFilesIn \
    trim/DEX_shSMC1A_1/DEX_shSMC1A_1_1.trim.pair1.fastq.gz \
    trim/DEX_shSMC1A_1/DEX_shSMC1A_1_1.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix alignment/DEX_shSMC1A_1/DEX_shSMC1A_1_1/ \
  --outSAMattrRGline ID:"DEX_shSMC1A_1_1" 	PL:"ILLUMINA" 			SM:"DEX_shSMC1A_1" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 70000000000 \
  --limitBAMsortRAM 70000000000 \
  --limitIObufferSize 1000000000 \
  --outWigType wiggle read1_5p --outWigStrand Stranded --outWigReferencesPrefix chr \
  --chimSegmentMin 21
star_align.2.DEX_shSMC1A_1_1.0ec6d5680146a7cb90da2f087367f168.mugqic.done
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
JOB_DONE=job_output/star/star_align.2.DEX_shSMC1A_1_2.7b7c40227ae9eabac454d5b2520aeac8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.2.DEX_shSMC1A_1_2.7b7c40227ae9eabac454d5b2520aeac8.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment/DEX_shSMC1A_1/DEX_shSMC1A_1_2 && \
STAR --runMode alignReads \
  --genomeDir reference.Merged \
  --readFilesIn \
    trim/DEX_shSMC1A_1/DEX_shSMC1A_1_2.trim.pair1.fastq.gz \
    trim/DEX_shSMC1A_1/DEX_shSMC1A_1_2.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix alignment/DEX_shSMC1A_1/DEX_shSMC1A_1_2/ \
  --outSAMattrRGline ID:"DEX_shSMC1A_1_2" 	PL:"ILLUMINA" 			SM:"DEX_shSMC1A_1" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 70000000000 \
  --limitBAMsortRAM 70000000000 \
  --limitIObufferSize 1000000000 \
  --outWigType wiggle read1_5p --outWigStrand Stranded --outWigReferencesPrefix chr \
  --chimSegmentMin 21
star_align.2.DEX_shSMC1A_1_2.7b7c40227ae9eabac454d5b2520aeac8.mugqic.done
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
JOB_DONE=job_output/star/star_align.2.ETOH_shSMC1A_1_1.6fa69bf38ae740c75a2bce9bda7f7d71.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.2.ETOH_shSMC1A_1_1.6fa69bf38ae740c75a2bce9bda7f7d71.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment/ETOH_shSMC1A_1/ETOH_shSMC1A_1_1 && \
STAR --runMode alignReads \
  --genomeDir reference.Merged \
  --readFilesIn \
    trim/ETOH_shSMC1A_1/ETOH_shSMC1A_1_1.trim.pair1.fastq.gz \
    trim/ETOH_shSMC1A_1/ETOH_shSMC1A_1_1.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix alignment/ETOH_shSMC1A_1/ETOH_shSMC1A_1_1/ \
  --outSAMattrRGline ID:"ETOH_shSMC1A_1_1" 	PL:"ILLUMINA" 			SM:"ETOH_shSMC1A_1" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 70000000000 \
  --limitBAMsortRAM 70000000000 \
  --limitIObufferSize 1000000000 \
  --outWigType wiggle read1_5p --outWigStrand Stranded --outWigReferencesPrefix chr \
  --chimSegmentMin 21
star_align.2.ETOH_shSMC1A_1_1.6fa69bf38ae740c75a2bce9bda7f7d71.mugqic.done
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
JOB_DONE=job_output/star/star_align.2.ETOH_shSMC1A_1_2.041445f16b03ccd06ecc3beabb0c4330.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.2.ETOH_shSMC1A_1_2.041445f16b03ccd06ecc3beabb0c4330.mugqic.done'
module load mugqic/star/2.5.1b && \
mkdir -p alignment/ETOH_shSMC1A_1/ETOH_shSMC1A_1_2 && \
STAR --runMode alignReads \
  --genomeDir reference.Merged \
  --readFilesIn \
    trim/ETOH_shSMC1A_1/ETOH_shSMC1A_1_2.trim.pair1.fastq.gz \
    trim/ETOH_shSMC1A_1/ETOH_shSMC1A_1_2.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix alignment/ETOH_shSMC1A_1/ETOH_shSMC1A_1_2/ \
  --outSAMattrRGline ID:"ETOH_shSMC1A_1_2" 	PL:"ILLUMINA" 			SM:"ETOH_shSMC1A_1" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 70000000000 \
  --limitBAMsortRAM 70000000000 \
  --limitIObufferSize 1000000000 \
  --outWigType wiggle read1_5p --outWigStrand Stranded --outWigReferencesPrefix chr \
  --chimSegmentMin 21
star_align.2.ETOH_shSMC1A_1_2.041445f16b03ccd06ecc3beabb0c4330.mugqic.done
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
JOB_DONE=job_output/star/star_report.91fb8ed6ea061fef22348cd2adfdba9e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_report.91fb8ed6ea061fef22348cd2adfdba9e.mugqic.done'
module load mugqic/pandoc/1.15.2 && \
mkdir -p report && \
pandoc --to=markdown \
  --template /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.3.0/bfx/report/RnaSeq.star.md \
  --variable scientific_name="Homo_sapiens" \
  --variable assembly="hg19" \
  /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.3.0/bfx/report/RnaSeq.star.md \
  > report/RnaSeq.star.md
star_report.91fb8ed6ea061fef22348cd2adfdba9e.mugqic.done
)
star_42_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_4_JOB_ID: picard_merge_sam_files.DEX_shSMC1A_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_merge_sam_files.DEX_shSMC1A_1
JOB_DEPENDENCIES=$star_28_JOB_ID:$star_29_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/picard_merge_sam_files.DEX_shSMC1A_1.4552066ba1c15d18a136681b23143643.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_merge_sam_files.DEX_shSMC1A_1.4552066ba1c15d18a136681b23143643.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/picard.jar MergeSamFiles \
 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/DEX_shSMC1A_1/DEX_shSMC1A_1_1/Aligned.sortedByCoord.out.bam \
  INPUT=alignment/DEX_shSMC1A_1/DEX_shSMC1A_1_2/Aligned.sortedByCoord.out.bam \
 OUTPUT=alignment/DEX_shSMC1A_1/DEX_shSMC1A_1.sorted.bam \
 MAX_RECORDS_IN_RAM=5750000
picard_merge_sam_files.DEX_shSMC1A_1.4552066ba1c15d18a136681b23143643.mugqic.done
)
picard_merge_sam_files_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_8_JOB_ID: picard_merge_sam_files.ETOH_shSMC1A_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_merge_sam_files.ETOH_shSMC1A_1
JOB_DEPENDENCIES=$star_39_JOB_ID:$star_40_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/picard_merge_sam_files.ETOH_shSMC1A_1.3be9305fa6d61e24b4a63fa3d0dc66d9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_merge_sam_files.ETOH_shSMC1A_1.3be9305fa6d61e24b4a63fa3d0dc66d9.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/picard.jar MergeSamFiles \
 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/ETOH_shSMC1A_1/ETOH_shSMC1A_1_1/Aligned.sortedByCoord.out.bam \
  INPUT=alignment/ETOH_shSMC1A_1/ETOH_shSMC1A_1_2/Aligned.sortedByCoord.out.bam \
 OUTPUT=alignment/ETOH_shSMC1A_1/ETOH_shSMC1A_1.sorted.bam \
 MAX_RECORDS_IN_RAM=5750000
picard_merge_sam_files.ETOH_shSMC1A_1.3be9305fa6d61e24b4a63fa3d0dc66d9.mugqic.done
)
picard_merge_sam_files_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_sort_sam_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sort_sam_4_JOB_ID: picard_sort_sam.DEX_shSMC1A_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_sort_sam.DEX_shSMC1A_1
JOB_DEPENDENCIES=$picard_merge_sam_files_4_JOB_ID
JOB_DONE=job_output/picard_sort_sam/picard_sort_sam.DEX_shSMC1A_1.620c7d34f6a1f32844d8762038cf225b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sort_sam.DEX_shSMC1A_1.620c7d34f6a1f32844d8762038cf225b.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
 java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/DEX_shSMC1A_1/DEX_shSMC1A_1.sorted.bam \
 OUTPUT=alignment/DEX_shSMC1A_1/DEX_shSMC1A_1.QueryNameSorted.bam \
 SORT_ORDER=queryname \
 MAX_RECORDS_IN_RAM=5750000
picard_sort_sam.DEX_shSMC1A_1.620c7d34f6a1f32844d8762038cf225b.mugqic.done
)
picard_sort_sam_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_sort_sam_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sort_sam_11_JOB_ID: picard_sort_sam.ETOH_shSMC1A_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_sort_sam.ETOH_shSMC1A_1
JOB_DEPENDENCIES=$picard_merge_sam_files_8_JOB_ID
JOB_DONE=job_output/picard_sort_sam/picard_sort_sam.ETOH_shSMC1A_1.ea56cd1ef106ddeb6005fa878d965b0c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sort_sam.ETOH_shSMC1A_1.ea56cd1ef106ddeb6005fa878d965b0c.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
 java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/ETOH_shSMC1A_1/ETOH_shSMC1A_1.sorted.bam \
 OUTPUT=alignment/ETOH_shSMC1A_1/ETOH_shSMC1A_1.QueryNameSorted.bam \
 SORT_ORDER=queryname \
 MAX_RECORDS_IN_RAM=5750000
picard_sort_sam.ETOH_shSMC1A_1.ea56cd1ef106ddeb6005fa878d965b0c.mugqic.done
)
picard_sort_sam_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_4_JOB_ID: picard_mark_duplicates.DEX_shSMC1A_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.DEX_shSMC1A_1
JOB_DEPENDENCIES=$picard_merge_sam_files_4_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.DEX_shSMC1A_1.7ca283ad4ac7964342845cd906651ae1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.DEX_shSMC1A_1.7ca283ad4ac7964342845cd906651ae1.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx14G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/DEX_shSMC1A_1/DEX_shSMC1A_1.sorted.bam \
 OUTPUT=alignment/DEX_shSMC1A_1/DEX_shSMC1A_1.sorted.mdup.bam \
 METRICS_FILE=alignment/DEX_shSMC1A_1/DEX_shSMC1A_1.sorted.mdup.metrics \
 MAX_RECORDS_IN_RAM=3500000
picard_mark_duplicates.DEX_shSMC1A_1.7ca283ad4ac7964342845cd906651ae1.mugqic.done
)
picard_mark_duplicates_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_11_JOB_ID: picard_mark_duplicates.ETOH_shSMC1A_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.ETOH_shSMC1A_1
JOB_DEPENDENCIES=$picard_merge_sam_files_8_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.ETOH_shSMC1A_1.f832b06c8e103b733e272f6d2ea8887f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.ETOH_shSMC1A_1.f832b06c8e103b733e272f6d2ea8887f.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx14G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/ETOH_shSMC1A_1/ETOH_shSMC1A_1.sorted.bam \
 OUTPUT=alignment/ETOH_shSMC1A_1/ETOH_shSMC1A_1.sorted.mdup.bam \
 METRICS_FILE=alignment/ETOH_shSMC1A_1/ETOH_shSMC1A_1.sorted.mdup.metrics \
 MAX_RECORDS_IN_RAM=3500000
picard_mark_duplicates.ETOH_shSMC1A_1.f832b06c8e103b733e272f6d2ea8887f.mugqic.done
)
picard_mark_duplicates_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=6 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
JOB_DONE=job_output/picard_rna_metrics/picard_rna_metrics.DEX_nonMamm_1.a63bd72892b731f9325d68e4763e7259.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_rna_metrics.DEX_nonMamm_1.a63bd72892b731f9325d68e4763e7259.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 mugqic/R_Bioconductor/3.2.3_3.2 && \
mkdir -p metrics/DEX_nonMamm_1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=4 -Xmx27G -jar $PICARD_HOME/picard.jar CollectMultipleMetrics \
 PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT \
 TMP_DIR=/gs/scratch/$USER \
 REFERENCE_SEQUENCE=/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa \
 INPUT=alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.bam \
 OUTPUT=metrics/DEX_nonMamm_1/DEX_nonMamm_1 \
 MAX_RECORDS_IN_RAM=5750000 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=4 -Xmx27G -jar $PICARD_HOME/picard.jar CollectRnaSeqMetrics \
 VALIDATION_STRINGENCY=SILENT  \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.bam \
 OUTPUT=metrics/DEX_nonMamm_1/DEX_nonMamm_1.picard_rna_metrics \
 REF_FLAT=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.ref_flat.tsv \
 STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
 MINIMUM_LENGTH=200 \
 REFERENCE_SEQUENCE=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa \
 MAX_RECORDS_IN_RAM=5750000
picard_rna_metrics.DEX_nonMamm_1.a63bd72892b731f9325d68e4763e7259.mugqic.done
)
picard_rna_metrics_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_rna_metrics_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_rna_metrics_2_JOB_ID: picard_rna_metrics.DEX_shMED1_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_rna_metrics.DEX_shMED1_1
JOB_DEPENDENCIES=$picard_mark_duplicates_2_JOB_ID
JOB_DONE=job_output/picard_rna_metrics/picard_rna_metrics.DEX_shMED1_1.8c884b2a04f9b665661e22b3b1a256bc.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_rna_metrics.DEX_shMED1_1.8c884b2a04f9b665661e22b3b1a256bc.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 mugqic/R_Bioconductor/3.2.3_3.2 && \
mkdir -p metrics/DEX_shMED1_1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=4 -Xmx27G -jar $PICARD_HOME/picard.jar CollectMultipleMetrics \
 PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT \
 TMP_DIR=/gs/scratch/$USER \
 REFERENCE_SEQUENCE=/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa \
 INPUT=alignment/DEX_shMED1_1/DEX_shMED1_1.sorted.mdup.bam \
 OUTPUT=metrics/DEX_shMED1_1/DEX_shMED1_1 \
 MAX_RECORDS_IN_RAM=5750000 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=4 -Xmx27G -jar $PICARD_HOME/picard.jar CollectRnaSeqMetrics \
 VALIDATION_STRINGENCY=SILENT  \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/DEX_shMED1_1/DEX_shMED1_1.sorted.mdup.bam \
 OUTPUT=metrics/DEX_shMED1_1/DEX_shMED1_1.picard_rna_metrics \
 REF_FLAT=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.ref_flat.tsv \
 STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
 MINIMUM_LENGTH=200 \
 REFERENCE_SEQUENCE=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa \
 MAX_RECORDS_IN_RAM=5750000
picard_rna_metrics.DEX_shMED1_1.8c884b2a04f9b665661e22b3b1a256bc.mugqic.done
)
picard_rna_metrics_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_rna_metrics_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_rna_metrics_3_JOB_ID: picard_rna_metrics.DEX_shNIPBL_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_rna_metrics.DEX_shNIPBL_1
JOB_DEPENDENCIES=$picard_mark_duplicates_3_JOB_ID
JOB_DONE=job_output/picard_rna_metrics/picard_rna_metrics.DEX_shNIPBL_1.1b36b9c3e8c7e97f59826bd65496c085.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_rna_metrics.DEX_shNIPBL_1.1b36b9c3e8c7e97f59826bd65496c085.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 mugqic/R_Bioconductor/3.2.3_3.2 && \
mkdir -p metrics/DEX_shNIPBL_1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=4 -Xmx27G -jar $PICARD_HOME/picard.jar CollectMultipleMetrics \
 PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT \
 TMP_DIR=/gs/scratch/$USER \
 REFERENCE_SEQUENCE=/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa \
 INPUT=alignment/DEX_shNIPBL_1/DEX_shNIPBL_1.sorted.mdup.bam \
 OUTPUT=metrics/DEX_shNIPBL_1/DEX_shNIPBL_1 \
 MAX_RECORDS_IN_RAM=5750000 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=4 -Xmx27G -jar $PICARD_HOME/picard.jar CollectRnaSeqMetrics \
 VALIDATION_STRINGENCY=SILENT  \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/DEX_shNIPBL_1/DEX_shNIPBL_1.sorted.mdup.bam \
 OUTPUT=metrics/DEX_shNIPBL_1/DEX_shNIPBL_1.picard_rna_metrics \
 REF_FLAT=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.ref_flat.tsv \
 STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
 MINIMUM_LENGTH=200 \
 REFERENCE_SEQUENCE=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa \
 MAX_RECORDS_IN_RAM=5750000
picard_rna_metrics.DEX_shNIPBL_1.1b36b9c3e8c7e97f59826bd65496c085.mugqic.done
)
picard_rna_metrics_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_rna_metrics_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_rna_metrics_4_JOB_ID: picard_rna_metrics.DEX_shSMC1A_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_rna_metrics.DEX_shSMC1A_1
JOB_DEPENDENCIES=$picard_mark_duplicates_4_JOB_ID
JOB_DONE=job_output/picard_rna_metrics/picard_rna_metrics.DEX_shSMC1A_1.2a74113af923af6d8d4d558c42c6a48d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_rna_metrics.DEX_shSMC1A_1.2a74113af923af6d8d4d558c42c6a48d.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 mugqic/R_Bioconductor/3.2.3_3.2 && \
mkdir -p metrics/DEX_shSMC1A_1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=4 -Xmx27G -jar $PICARD_HOME/picard.jar CollectMultipleMetrics \
 PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT \
 TMP_DIR=/gs/scratch/$USER \
 REFERENCE_SEQUENCE=/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa \
 INPUT=alignment/DEX_shSMC1A_1/DEX_shSMC1A_1.sorted.mdup.bam \
 OUTPUT=metrics/DEX_shSMC1A_1/DEX_shSMC1A_1 \
 MAX_RECORDS_IN_RAM=5750000 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=4 -Xmx27G -jar $PICARD_HOME/picard.jar CollectRnaSeqMetrics \
 VALIDATION_STRINGENCY=SILENT  \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/DEX_shSMC1A_1/DEX_shSMC1A_1.sorted.mdup.bam \
 OUTPUT=metrics/DEX_shSMC1A_1/DEX_shSMC1A_1.picard_rna_metrics \
 REF_FLAT=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.ref_flat.tsv \
 STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
 MINIMUM_LENGTH=200 \
 REFERENCE_SEQUENCE=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa \
 MAX_RECORDS_IN_RAM=5750000
picard_rna_metrics.DEX_shSMC1A_1.2a74113af923af6d8d4d558c42c6a48d.mugqic.done
)
picard_rna_metrics_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_rna_metrics_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_rna_metrics_5_JOB_ID: picard_rna_metrics.ETOH_nonMamm_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_rna_metrics.ETOH_nonMamm_1
JOB_DEPENDENCIES=$picard_mark_duplicates_5_JOB_ID
JOB_DONE=job_output/picard_rna_metrics/picard_rna_metrics.ETOH_nonMamm_1.801b77fc97df695139fa5f4921127308.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_rna_metrics.ETOH_nonMamm_1.801b77fc97df695139fa5f4921127308.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 mugqic/R_Bioconductor/3.2.3_3.2 && \
mkdir -p metrics/ETOH_nonMamm_1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=4 -Xmx27G -jar $PICARD_HOME/picard.jar CollectMultipleMetrics \
 PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT \
 TMP_DIR=/gs/scratch/$USER \
 REFERENCE_SEQUENCE=/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa \
 INPUT=alignment/ETOH_nonMamm_1/ETOH_nonMamm_1.sorted.mdup.bam \
 OUTPUT=metrics/ETOH_nonMamm_1/ETOH_nonMamm_1 \
 MAX_RECORDS_IN_RAM=5750000 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=4 -Xmx27G -jar $PICARD_HOME/picard.jar CollectRnaSeqMetrics \
 VALIDATION_STRINGENCY=SILENT  \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/ETOH_nonMamm_1/ETOH_nonMamm_1.sorted.mdup.bam \
 OUTPUT=metrics/ETOH_nonMamm_1/ETOH_nonMamm_1.picard_rna_metrics \
 REF_FLAT=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.ref_flat.tsv \
 STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
 MINIMUM_LENGTH=200 \
 REFERENCE_SEQUENCE=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa \
 MAX_RECORDS_IN_RAM=5750000
picard_rna_metrics.ETOH_nonMamm_1.801b77fc97df695139fa5f4921127308.mugqic.done
)
picard_rna_metrics_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_rna_metrics_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_rna_metrics_6_JOB_ID: picard_rna_metrics.ETOH_nonMamm_2
#-------------------------------------------------------------------------------
JOB_NAME=picard_rna_metrics.ETOH_nonMamm_2
JOB_DEPENDENCIES=$picard_mark_duplicates_6_JOB_ID
JOB_DONE=job_output/picard_rna_metrics/picard_rna_metrics.ETOH_nonMamm_2.716c94c9994c52f3a04f16c2aac19b0a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_rna_metrics.ETOH_nonMamm_2.716c94c9994c52f3a04f16c2aac19b0a.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 mugqic/R_Bioconductor/3.2.3_3.2 && \
mkdir -p metrics/ETOH_nonMamm_2 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=4 -Xmx27G -jar $PICARD_HOME/picard.jar CollectMultipleMetrics \
 PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT \
 TMP_DIR=/gs/scratch/$USER \
 REFERENCE_SEQUENCE=/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa \
 INPUT=alignment/ETOH_nonMamm_2/ETOH_nonMamm_2.sorted.mdup.bam \
 OUTPUT=metrics/ETOH_nonMamm_2/ETOH_nonMamm_2 \
 MAX_RECORDS_IN_RAM=5750000 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=4 -Xmx27G -jar $PICARD_HOME/picard.jar CollectRnaSeqMetrics \
 VALIDATION_STRINGENCY=SILENT  \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/ETOH_nonMamm_2/ETOH_nonMamm_2.sorted.mdup.bam \
 OUTPUT=metrics/ETOH_nonMamm_2/ETOH_nonMamm_2.picard_rna_metrics \
 REF_FLAT=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.ref_flat.tsv \
 STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
 MINIMUM_LENGTH=200 \
 REFERENCE_SEQUENCE=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa \
 MAX_RECORDS_IN_RAM=5750000
picard_rna_metrics.ETOH_nonMamm_2.716c94c9994c52f3a04f16c2aac19b0a.mugqic.done
)
picard_rna_metrics_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_rna_metrics_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_rna_metrics_7_JOB_ID: picard_rna_metrics.ETOH_shMED1_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_rna_metrics.ETOH_shMED1_1
JOB_DEPENDENCIES=$picard_mark_duplicates_7_JOB_ID
JOB_DONE=job_output/picard_rna_metrics/picard_rna_metrics.ETOH_shMED1_1.e0fe59bb1eea394cdd89a1414927893e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_rna_metrics.ETOH_shMED1_1.e0fe59bb1eea394cdd89a1414927893e.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 mugqic/R_Bioconductor/3.2.3_3.2 && \
mkdir -p metrics/ETOH_shMED1_1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=4 -Xmx27G -jar $PICARD_HOME/picard.jar CollectMultipleMetrics \
 PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT \
 TMP_DIR=/gs/scratch/$USER \
 REFERENCE_SEQUENCE=/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa \
 INPUT=alignment/ETOH_shMED1_1/ETOH_shMED1_1.sorted.mdup.bam \
 OUTPUT=metrics/ETOH_shMED1_1/ETOH_shMED1_1 \
 MAX_RECORDS_IN_RAM=5750000 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=4 -Xmx27G -jar $PICARD_HOME/picard.jar CollectRnaSeqMetrics \
 VALIDATION_STRINGENCY=SILENT  \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/ETOH_shMED1_1/ETOH_shMED1_1.sorted.mdup.bam \
 OUTPUT=metrics/ETOH_shMED1_1/ETOH_shMED1_1.picard_rna_metrics \
 REF_FLAT=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.ref_flat.tsv \
 STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
 MINIMUM_LENGTH=200 \
 REFERENCE_SEQUENCE=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa \
 MAX_RECORDS_IN_RAM=5750000
picard_rna_metrics.ETOH_shMED1_1.e0fe59bb1eea394cdd89a1414927893e.mugqic.done
)
picard_rna_metrics_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_rna_metrics_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_rna_metrics_8_JOB_ID: picard_rna_metrics.ETOH_shNIPBL_2
#-------------------------------------------------------------------------------
JOB_NAME=picard_rna_metrics.ETOH_shNIPBL_2
JOB_DEPENDENCIES=$picard_mark_duplicates_8_JOB_ID
JOB_DONE=job_output/picard_rna_metrics/picard_rna_metrics.ETOH_shNIPBL_2.31f88d45ab6a1bc7e42ff6d366837b6e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_rna_metrics.ETOH_shNIPBL_2.31f88d45ab6a1bc7e42ff6d366837b6e.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 mugqic/R_Bioconductor/3.2.3_3.2 && \
mkdir -p metrics/ETOH_shNIPBL_2 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=4 -Xmx27G -jar $PICARD_HOME/picard.jar CollectMultipleMetrics \
 PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT \
 TMP_DIR=/gs/scratch/$USER \
 REFERENCE_SEQUENCE=/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa \
 INPUT=alignment/ETOH_shNIPBL_2/ETOH_shNIPBL_2.sorted.mdup.bam \
 OUTPUT=metrics/ETOH_shNIPBL_2/ETOH_shNIPBL_2 \
 MAX_RECORDS_IN_RAM=5750000 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=4 -Xmx27G -jar $PICARD_HOME/picard.jar CollectRnaSeqMetrics \
 VALIDATION_STRINGENCY=SILENT  \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/ETOH_shNIPBL_2/ETOH_shNIPBL_2.sorted.mdup.bam \
 OUTPUT=metrics/ETOH_shNIPBL_2/ETOH_shNIPBL_2.picard_rna_metrics \
 REF_FLAT=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.ref_flat.tsv \
 STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
 MINIMUM_LENGTH=200 \
 REFERENCE_SEQUENCE=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa \
 MAX_RECORDS_IN_RAM=5750000
picard_rna_metrics.ETOH_shNIPBL_2.31f88d45ab6a1bc7e42ff6d366837b6e.mugqic.done
)
picard_rna_metrics_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_rna_metrics_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_rna_metrics_9_JOB_ID: picard_rna_metrics.ETOH_shNIPBL_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_rna_metrics.ETOH_shNIPBL_1
JOB_DEPENDENCIES=$picard_mark_duplicates_9_JOB_ID
JOB_DONE=job_output/picard_rna_metrics/picard_rna_metrics.ETOH_shNIPBL_1.38cb6ea10580725a411fb8efd05f3778.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_rna_metrics.ETOH_shNIPBL_1.38cb6ea10580725a411fb8efd05f3778.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 mugqic/R_Bioconductor/3.2.3_3.2 && \
mkdir -p metrics/ETOH_shNIPBL_1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=4 -Xmx27G -jar $PICARD_HOME/picard.jar CollectMultipleMetrics \
 PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT \
 TMP_DIR=/gs/scratch/$USER \
 REFERENCE_SEQUENCE=/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa \
 INPUT=alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1.sorted.mdup.bam \
 OUTPUT=metrics/ETOH_shNIPBL_1/ETOH_shNIPBL_1 \
 MAX_RECORDS_IN_RAM=5750000 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=4 -Xmx27G -jar $PICARD_HOME/picard.jar CollectRnaSeqMetrics \
 VALIDATION_STRINGENCY=SILENT  \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1.sorted.mdup.bam \
 OUTPUT=metrics/ETOH_shNIPBL_1/ETOH_shNIPBL_1.picard_rna_metrics \
 REF_FLAT=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.ref_flat.tsv \
 STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
 MINIMUM_LENGTH=200 \
 REFERENCE_SEQUENCE=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa \
 MAX_RECORDS_IN_RAM=5750000
picard_rna_metrics.ETOH_shNIPBL_1.38cb6ea10580725a411fb8efd05f3778.mugqic.done
)
picard_rna_metrics_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_rna_metrics_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_rna_metrics_10_JOB_ID: picard_rna_metrics.ETOH_shMED1_2
#-------------------------------------------------------------------------------
JOB_NAME=picard_rna_metrics.ETOH_shMED1_2
JOB_DEPENDENCIES=$picard_mark_duplicates_10_JOB_ID
JOB_DONE=job_output/picard_rna_metrics/picard_rna_metrics.ETOH_shMED1_2.be67e57a9ac693b20340fb6c6fdddbf3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_rna_metrics.ETOH_shMED1_2.be67e57a9ac693b20340fb6c6fdddbf3.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 mugqic/R_Bioconductor/3.2.3_3.2 && \
mkdir -p metrics/ETOH_shMED1_2 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=4 -Xmx27G -jar $PICARD_HOME/picard.jar CollectMultipleMetrics \
 PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT \
 TMP_DIR=/gs/scratch/$USER \
 REFERENCE_SEQUENCE=/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa \
 INPUT=alignment/ETOH_shMED1_2/ETOH_shMED1_2.sorted.mdup.bam \
 OUTPUT=metrics/ETOH_shMED1_2/ETOH_shMED1_2 \
 MAX_RECORDS_IN_RAM=5750000 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=4 -Xmx27G -jar $PICARD_HOME/picard.jar CollectRnaSeqMetrics \
 VALIDATION_STRINGENCY=SILENT  \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/ETOH_shMED1_2/ETOH_shMED1_2.sorted.mdup.bam \
 OUTPUT=metrics/ETOH_shMED1_2/ETOH_shMED1_2.picard_rna_metrics \
 REF_FLAT=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.ref_flat.tsv \
 STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
 MINIMUM_LENGTH=200 \
 REFERENCE_SEQUENCE=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa \
 MAX_RECORDS_IN_RAM=5750000
picard_rna_metrics.ETOH_shMED1_2.be67e57a9ac693b20340fb6c6fdddbf3.mugqic.done
)
picard_rna_metrics_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_rna_metrics_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_rna_metrics_11_JOB_ID: picard_rna_metrics.ETOH_shSMC1A_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_rna_metrics.ETOH_shSMC1A_1
JOB_DEPENDENCIES=$picard_mark_duplicates_11_JOB_ID
JOB_DONE=job_output/picard_rna_metrics/picard_rna_metrics.ETOH_shSMC1A_1.6ec31c191c420a687e81e5c92b3b1788.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_rna_metrics.ETOH_shSMC1A_1.6ec31c191c420a687e81e5c92b3b1788.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 mugqic/R_Bioconductor/3.2.3_3.2 && \
mkdir -p metrics/ETOH_shSMC1A_1 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=4 -Xmx27G -jar $PICARD_HOME/picard.jar CollectMultipleMetrics \
 PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT \
 TMP_DIR=/gs/scratch/$USER \
 REFERENCE_SEQUENCE=/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa \
 INPUT=alignment/ETOH_shSMC1A_1/ETOH_shSMC1A_1.sorted.mdup.bam \
 OUTPUT=metrics/ETOH_shSMC1A_1/ETOH_shSMC1A_1 \
 MAX_RECORDS_IN_RAM=5750000 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=4 -Xmx27G -jar $PICARD_HOME/picard.jar CollectRnaSeqMetrics \
 VALIDATION_STRINGENCY=SILENT  \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/ETOH_shSMC1A_1/ETOH_shSMC1A_1.sorted.mdup.bam \
 OUTPUT=metrics/ETOH_shSMC1A_1/ETOH_shSMC1A_1.picard_rna_metrics \
 REF_FLAT=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.ref_flat.tsv \
 STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
 MINIMUM_LENGTH=200 \
 REFERENCE_SEQUENCE=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa \
 MAX_RECORDS_IN_RAM=5750000
picard_rna_metrics.ETOH_shSMC1A_1.6ec31c191c420a687e81e5c92b3b1788.mugqic.done
)
picard_rna_metrics_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_rna_metrics_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_rna_metrics_12_JOB_ID: picard_rna_metrics.ETOH_shSMC1A_2
#-------------------------------------------------------------------------------
JOB_NAME=picard_rna_metrics.ETOH_shSMC1A_2
JOB_DEPENDENCIES=$picard_mark_duplicates_12_JOB_ID
JOB_DONE=job_output/picard_rna_metrics/picard_rna_metrics.ETOH_shSMC1A_2.45ee3a6794c9743caddf26df4b8801f0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_rna_metrics.ETOH_shSMC1A_2.45ee3a6794c9743caddf26df4b8801f0.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 mugqic/R_Bioconductor/3.2.3_3.2 && \
mkdir -p metrics/ETOH_shSMC1A_2 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=4 -Xmx27G -jar $PICARD_HOME/picard.jar CollectMultipleMetrics \
 PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT \
 TMP_DIR=/gs/scratch/$USER \
 REFERENCE_SEQUENCE=/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa \
 INPUT=alignment/ETOH_shSMC1A_2/ETOH_shSMC1A_2.sorted.mdup.bam \
 OUTPUT=metrics/ETOH_shSMC1A_2/ETOH_shSMC1A_2 \
 MAX_RECORDS_IN_RAM=5750000 && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=4 -Xmx27G -jar $PICARD_HOME/picard.jar CollectRnaSeqMetrics \
 VALIDATION_STRINGENCY=SILENT  \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/ETOH_shSMC1A_2/ETOH_shSMC1A_2.sorted.mdup.bam \
 OUTPUT=metrics/ETOH_shSMC1A_2/ETOH_shSMC1A_2.picard_rna_metrics \
 REF_FLAT=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.ref_flat.tsv \
 STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
 MINIMUM_LENGTH=200 \
 REFERENCE_SEQUENCE=$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa \
 MAX_RECORDS_IN_RAM=5750000
picard_rna_metrics.ETOH_shSMC1A_2.45ee3a6794c9743caddf26df4b8801f0.mugqic.done
)
picard_rna_metrics_12_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
JOB_DONE=job_output/estimate_ribosomal_rna/bwa_mem_rRNA.DEX_nonMamm_1_1.0ab847dd3960ef33d2855fd88662561f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_rRNA.DEX_nonMamm_1_1.0ab847dd3960ef33d2855fd88662561f.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/bvatools/1.6 mugqic/bwa/0.7.12 mugqic/picard/2.0.1 mugqic/mugqic_tools/2.1.7 mugqic/python/2.7.13 && \
mkdir -p alignment/DEX_nonMamm_1/DEX_nonMamm_1_1 metrics/DEX_nonMamm_1/DEX_nonMamm_1_1 && \
java -XX:ParallelGCThreads=1 -Xmx10G -jar $BVATOOLS_JAR \
  bam2fq --mapped ONLY \
  --bam alignment/DEX_nonMamm_1/DEX_nonMamm_1_1/Aligned.sortedByCoord.out.bam    | \
bwa mem  \
  -M -t 10 \
  -R '@RG	ID:DEX_nonMamm_1_1	SM:DEX_nonMamm_1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/annotations/rrna_bwa_index/Homo_sapiens.hg19.UCSC2009-03-08.rrna.fa \
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
  -g $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.gtf \
  -o metrics/DEX_nonMamm_1/DEX_nonMamm_1_1/DEX_nonMamm_1_1rRNA.stats.tsv \
  -t transcript
bwa_mem_rRNA.DEX_nonMamm_1_1.0ab847dd3960ef33d2855fd88662561f.mugqic.done
)
estimate_ribosomal_rna_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$estimate_ribosomal_rna_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: estimate_ribosomal_rna_2_JOB_ID: bwa_mem_rRNA.DEX_nonMamm_1_2
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_rRNA.DEX_nonMamm_1_2
JOB_DEPENDENCIES=$star_23_JOB_ID
JOB_DONE=job_output/estimate_ribosomal_rna/bwa_mem_rRNA.DEX_nonMamm_1_2.0b005259235f11672d0338c3d05c953e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_rRNA.DEX_nonMamm_1_2.0b005259235f11672d0338c3d05c953e.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/bvatools/1.6 mugqic/bwa/0.7.12 mugqic/picard/2.0.1 mugqic/mugqic_tools/2.1.7 mugqic/python/2.7.13 && \
mkdir -p alignment/DEX_nonMamm_1/DEX_nonMamm_1_2 metrics/DEX_nonMamm_1/DEX_nonMamm_1_2 && \
java -XX:ParallelGCThreads=1 -Xmx10G -jar $BVATOOLS_JAR \
  bam2fq --mapped ONLY \
  --bam alignment/DEX_nonMamm_1/DEX_nonMamm_1_2/Aligned.sortedByCoord.out.bam    | \
bwa mem  \
  -M -t 10 \
  -R '@RG	ID:DEX_nonMamm_1_2	SM:DEX_nonMamm_1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/annotations/rrna_bwa_index/Homo_sapiens.hg19.UCSC2009-03-08.rrna.fa \
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
  -g $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.gtf \
  -o metrics/DEX_nonMamm_1/DEX_nonMamm_1_2/DEX_nonMamm_1_2rRNA.stats.tsv \
  -t transcript
bwa_mem_rRNA.DEX_nonMamm_1_2.0b005259235f11672d0338c3d05c953e.mugqic.done
)
estimate_ribosomal_rna_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$estimate_ribosomal_rna_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: estimate_ribosomal_rna_3_JOB_ID: bwa_mem_rRNA.DEX_shMED1_1_1
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_rRNA.DEX_shMED1_1_1
JOB_DEPENDENCIES=$star_24_JOB_ID
JOB_DONE=job_output/estimate_ribosomal_rna/bwa_mem_rRNA.DEX_shMED1_1_1.4addca27bdeceb5a6c67123d5b12525f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_rRNA.DEX_shMED1_1_1.4addca27bdeceb5a6c67123d5b12525f.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/bvatools/1.6 mugqic/bwa/0.7.12 mugqic/picard/2.0.1 mugqic/mugqic_tools/2.1.7 mugqic/python/2.7.13 && \
mkdir -p alignment/DEX_shMED1_1/DEX_shMED1_1_1 metrics/DEX_shMED1_1/DEX_shMED1_1_1 && \
java -XX:ParallelGCThreads=1 -Xmx10G -jar $BVATOOLS_JAR \
  bam2fq --mapped ONLY \
  --bam alignment/DEX_shMED1_1/DEX_shMED1_1_1/Aligned.sortedByCoord.out.bam    | \
bwa mem  \
  -M -t 10 \
  -R '@RG	ID:DEX_shMED1_1_1	SM:DEX_shMED1_1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/annotations/rrna_bwa_index/Homo_sapiens.hg19.UCSC2009-03-08.rrna.fa \
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
  -g $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.gtf \
  -o metrics/DEX_shMED1_1/DEX_shMED1_1_1/DEX_shMED1_1_1rRNA.stats.tsv \
  -t transcript
bwa_mem_rRNA.DEX_shMED1_1_1.4addca27bdeceb5a6c67123d5b12525f.mugqic.done
)
estimate_ribosomal_rna_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$estimate_ribosomal_rna_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: estimate_ribosomal_rna_4_JOB_ID: bwa_mem_rRNA.DEX_shMED1_1_2
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_rRNA.DEX_shMED1_1_2
JOB_DEPENDENCIES=$star_25_JOB_ID
JOB_DONE=job_output/estimate_ribosomal_rna/bwa_mem_rRNA.DEX_shMED1_1_2.a63710e3ee0fb14f2d18f738c65588b4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_rRNA.DEX_shMED1_1_2.a63710e3ee0fb14f2d18f738c65588b4.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/bvatools/1.6 mugqic/bwa/0.7.12 mugqic/picard/2.0.1 mugqic/mugqic_tools/2.1.7 mugqic/python/2.7.13 && \
mkdir -p alignment/DEX_shMED1_1/DEX_shMED1_1_2 metrics/DEX_shMED1_1/DEX_shMED1_1_2 && \
java -XX:ParallelGCThreads=1 -Xmx10G -jar $BVATOOLS_JAR \
  bam2fq --mapped ONLY \
  --bam alignment/DEX_shMED1_1/DEX_shMED1_1_2/Aligned.sortedByCoord.out.bam    | \
bwa mem  \
  -M -t 10 \
  -R '@RG	ID:DEX_shMED1_1_2	SM:DEX_shMED1_1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/annotations/rrna_bwa_index/Homo_sapiens.hg19.UCSC2009-03-08.rrna.fa \
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
  -g $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.gtf \
  -o metrics/DEX_shMED1_1/DEX_shMED1_1_2/DEX_shMED1_1_2rRNA.stats.tsv \
  -t transcript
bwa_mem_rRNA.DEX_shMED1_1_2.a63710e3ee0fb14f2d18f738c65588b4.mugqic.done
)
estimate_ribosomal_rna_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$estimate_ribosomal_rna_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: estimate_ribosomal_rna_5_JOB_ID: bwa_mem_rRNA.DEX_shNIPBL_1_1
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_rRNA.DEX_shNIPBL_1_1
JOB_DEPENDENCIES=$star_26_JOB_ID
JOB_DONE=job_output/estimate_ribosomal_rna/bwa_mem_rRNA.DEX_shNIPBL_1_1.0695eaaec1c7386903543f6e8fbeb39f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_rRNA.DEX_shNIPBL_1_1.0695eaaec1c7386903543f6e8fbeb39f.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/bvatools/1.6 mugqic/bwa/0.7.12 mugqic/picard/2.0.1 mugqic/mugqic_tools/2.1.7 mugqic/python/2.7.13 && \
mkdir -p alignment/DEX_shNIPBL_1/DEX_shNIPBL_1_1 metrics/DEX_shNIPBL_1/DEX_shNIPBL_1_1 && \
java -XX:ParallelGCThreads=1 -Xmx10G -jar $BVATOOLS_JAR \
  bam2fq --mapped ONLY \
  --bam alignment/DEX_shNIPBL_1/DEX_shNIPBL_1_1/Aligned.sortedByCoord.out.bam    | \
bwa mem  \
  -M -t 10 \
  -R '@RG	ID:DEX_shNIPBL_1_1	SM:DEX_shNIPBL_1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/annotations/rrna_bwa_index/Homo_sapiens.hg19.UCSC2009-03-08.rrna.fa \
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
  -g $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.gtf \
  -o metrics/DEX_shNIPBL_1/DEX_shNIPBL_1_1/DEX_shNIPBL_1_1rRNA.stats.tsv \
  -t transcript
bwa_mem_rRNA.DEX_shNIPBL_1_1.0695eaaec1c7386903543f6e8fbeb39f.mugqic.done
)
estimate_ribosomal_rna_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$estimate_ribosomal_rna_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: estimate_ribosomal_rna_6_JOB_ID: bwa_mem_rRNA.DEX_shNIPBL_1_2
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_rRNA.DEX_shNIPBL_1_2
JOB_DEPENDENCIES=$star_27_JOB_ID
JOB_DONE=job_output/estimate_ribosomal_rna/bwa_mem_rRNA.DEX_shNIPBL_1_2.cd65de41daf8390ef22c8f0d1addc0b3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_rRNA.DEX_shNIPBL_1_2.cd65de41daf8390ef22c8f0d1addc0b3.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/bvatools/1.6 mugqic/bwa/0.7.12 mugqic/picard/2.0.1 mugqic/mugqic_tools/2.1.7 mugqic/python/2.7.13 && \
mkdir -p alignment/DEX_shNIPBL_1/DEX_shNIPBL_1_2 metrics/DEX_shNIPBL_1/DEX_shNIPBL_1_2 && \
java -XX:ParallelGCThreads=1 -Xmx10G -jar $BVATOOLS_JAR \
  bam2fq --mapped ONLY \
  --bam alignment/DEX_shNIPBL_1/DEX_shNIPBL_1_2/Aligned.sortedByCoord.out.bam    | \
bwa mem  \
  -M -t 10 \
  -R '@RG	ID:DEX_shNIPBL_1_2	SM:DEX_shNIPBL_1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/annotations/rrna_bwa_index/Homo_sapiens.hg19.UCSC2009-03-08.rrna.fa \
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
  -g $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.gtf \
  -o metrics/DEX_shNIPBL_1/DEX_shNIPBL_1_2/DEX_shNIPBL_1_2rRNA.stats.tsv \
  -t transcript
bwa_mem_rRNA.DEX_shNIPBL_1_2.cd65de41daf8390ef22c8f0d1addc0b3.mugqic.done
)
estimate_ribosomal_rna_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$estimate_ribosomal_rna_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: estimate_ribosomal_rna_7_JOB_ID: bwa_mem_rRNA.DEX_shSMC1A_1_1
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_rRNA.DEX_shSMC1A_1_1
JOB_DEPENDENCIES=$star_28_JOB_ID
JOB_DONE=job_output/estimate_ribosomal_rna/bwa_mem_rRNA.DEX_shSMC1A_1_1.07e2895047478291b653bcd2d813ba1e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_rRNA.DEX_shSMC1A_1_1.07e2895047478291b653bcd2d813ba1e.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/bvatools/1.6 mugqic/bwa/0.7.12 mugqic/picard/2.0.1 mugqic/mugqic_tools/2.1.7 mugqic/python/2.7.13 && \
mkdir -p alignment/DEX_shSMC1A_1/DEX_shSMC1A_1_1 metrics/DEX_shSMC1A_1/DEX_shSMC1A_1_1 && \
java -XX:ParallelGCThreads=1 -Xmx10G -jar $BVATOOLS_JAR \
  bam2fq --mapped ONLY \
  --bam alignment/DEX_shSMC1A_1/DEX_shSMC1A_1_1/Aligned.sortedByCoord.out.bam    | \
bwa mem  \
  -M -t 10 \
  -R '@RG	ID:DEX_shSMC1A_1_1	SM:DEX_shSMC1A_1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/annotations/rrna_bwa_index/Homo_sapiens.hg19.UCSC2009-03-08.rrna.fa \
  /dev/stdin | \
 java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx7G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=/dev/stdin \
 OUTPUT=metrics/DEX_shSMC1A_1/DEX_shSMC1A_1_1/DEX_shSMC1A_1_1rRNA.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=1750000 && \
python $PYTHON_TOOLS/rrnaBAMcounter.py \
  -i metrics/DEX_shSMC1A_1/DEX_shSMC1A_1_1/DEX_shSMC1A_1_1rRNA.bam \
  -g $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.gtf \
  -o metrics/DEX_shSMC1A_1/DEX_shSMC1A_1_1/DEX_shSMC1A_1_1rRNA.stats.tsv \
  -t transcript
bwa_mem_rRNA.DEX_shSMC1A_1_1.07e2895047478291b653bcd2d813ba1e.mugqic.done
)
estimate_ribosomal_rna_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$estimate_ribosomal_rna_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: estimate_ribosomal_rna_8_JOB_ID: bwa_mem_rRNA.DEX_shSMC1A_1_2
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_rRNA.DEX_shSMC1A_1_2
JOB_DEPENDENCIES=$star_29_JOB_ID
JOB_DONE=job_output/estimate_ribosomal_rna/bwa_mem_rRNA.DEX_shSMC1A_1_2.b551610311338eb084d782918b1e6f80.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_rRNA.DEX_shSMC1A_1_2.b551610311338eb084d782918b1e6f80.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/bvatools/1.6 mugqic/bwa/0.7.12 mugqic/picard/2.0.1 mugqic/mugqic_tools/2.1.7 mugqic/python/2.7.13 && \
mkdir -p alignment/DEX_shSMC1A_1/DEX_shSMC1A_1_2 metrics/DEX_shSMC1A_1/DEX_shSMC1A_1_2 && \
java -XX:ParallelGCThreads=1 -Xmx10G -jar $BVATOOLS_JAR \
  bam2fq --mapped ONLY \
  --bam alignment/DEX_shSMC1A_1/DEX_shSMC1A_1_2/Aligned.sortedByCoord.out.bam    | \
bwa mem  \
  -M -t 10 \
  -R '@RG	ID:DEX_shSMC1A_1_2	SM:DEX_shSMC1A_1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/annotations/rrna_bwa_index/Homo_sapiens.hg19.UCSC2009-03-08.rrna.fa \
  /dev/stdin | \
 java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx7G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=/dev/stdin \
 OUTPUT=metrics/DEX_shSMC1A_1/DEX_shSMC1A_1_2/DEX_shSMC1A_1_2rRNA.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=1750000 && \
python $PYTHON_TOOLS/rrnaBAMcounter.py \
  -i metrics/DEX_shSMC1A_1/DEX_shSMC1A_1_2/DEX_shSMC1A_1_2rRNA.bam \
  -g $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.gtf \
  -o metrics/DEX_shSMC1A_1/DEX_shSMC1A_1_2/DEX_shSMC1A_1_2rRNA.stats.tsv \
  -t transcript
bwa_mem_rRNA.DEX_shSMC1A_1_2.b551610311338eb084d782918b1e6f80.mugqic.done
)
estimate_ribosomal_rna_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$estimate_ribosomal_rna_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: estimate_ribosomal_rna_9_JOB_ID: bwa_mem_rRNA.ETOH_nonMamm_1_1
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_rRNA.ETOH_nonMamm_1_1
JOB_DEPENDENCIES=$star_30_JOB_ID
JOB_DONE=job_output/estimate_ribosomal_rna/bwa_mem_rRNA.ETOH_nonMamm_1_1.f848c88cf3bd7ed1f89a2ef2af256ea6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_rRNA.ETOH_nonMamm_1_1.f848c88cf3bd7ed1f89a2ef2af256ea6.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/bvatools/1.6 mugqic/bwa/0.7.12 mugqic/picard/2.0.1 mugqic/mugqic_tools/2.1.7 mugqic/python/2.7.13 && \
mkdir -p alignment/ETOH_nonMamm_1/ETOH_nonMamm_1_1 metrics/ETOH_nonMamm_1/ETOH_nonMamm_1_1 && \
java -XX:ParallelGCThreads=1 -Xmx10G -jar $BVATOOLS_JAR \
  bam2fq --mapped ONLY \
  --bam alignment/ETOH_nonMamm_1/ETOH_nonMamm_1_1/Aligned.sortedByCoord.out.bam    | \
bwa mem  \
  -M -t 10 \
  -R '@RG	ID:ETOH_nonMamm_1_1	SM:ETOH_nonMamm_1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/annotations/rrna_bwa_index/Homo_sapiens.hg19.UCSC2009-03-08.rrna.fa \
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
  -g $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.gtf \
  -o metrics/ETOH_nonMamm_1/ETOH_nonMamm_1_1/ETOH_nonMamm_1_1rRNA.stats.tsv \
  -t transcript
bwa_mem_rRNA.ETOH_nonMamm_1_1.f848c88cf3bd7ed1f89a2ef2af256ea6.mugqic.done
)
estimate_ribosomal_rna_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$estimate_ribosomal_rna_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: estimate_ribosomal_rna_10_JOB_ID: bwa_mem_rRNA.ETOH_nonMamm_1_2
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_rRNA.ETOH_nonMamm_1_2
JOB_DEPENDENCIES=$star_31_JOB_ID
JOB_DONE=job_output/estimate_ribosomal_rna/bwa_mem_rRNA.ETOH_nonMamm_1_2.511059832aa47a88574114e9eb4a56b8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_rRNA.ETOH_nonMamm_1_2.511059832aa47a88574114e9eb4a56b8.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/bvatools/1.6 mugqic/bwa/0.7.12 mugqic/picard/2.0.1 mugqic/mugqic_tools/2.1.7 mugqic/python/2.7.13 && \
mkdir -p alignment/ETOH_nonMamm_1/ETOH_nonMamm_1_2 metrics/ETOH_nonMamm_1/ETOH_nonMamm_1_2 && \
java -XX:ParallelGCThreads=1 -Xmx10G -jar $BVATOOLS_JAR \
  bam2fq --mapped ONLY \
  --bam alignment/ETOH_nonMamm_1/ETOH_nonMamm_1_2/Aligned.sortedByCoord.out.bam    | \
bwa mem  \
  -M -t 10 \
  -R '@RG	ID:ETOH_nonMamm_1_2	SM:ETOH_nonMamm_1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/annotations/rrna_bwa_index/Homo_sapiens.hg19.UCSC2009-03-08.rrna.fa \
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
  -g $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.gtf \
  -o metrics/ETOH_nonMamm_1/ETOH_nonMamm_1_2/ETOH_nonMamm_1_2rRNA.stats.tsv \
  -t transcript
bwa_mem_rRNA.ETOH_nonMamm_1_2.511059832aa47a88574114e9eb4a56b8.mugqic.done
)
estimate_ribosomal_rna_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$estimate_ribosomal_rna_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: estimate_ribosomal_rna_11_JOB_ID: bwa_mem_rRNA.ETOH_nonMamm_2_1
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_rRNA.ETOH_nonMamm_2_1
JOB_DEPENDENCIES=$star_32_JOB_ID
JOB_DONE=job_output/estimate_ribosomal_rna/bwa_mem_rRNA.ETOH_nonMamm_2_1.cc1bc2c6e90d9925ca37de13075d512e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_rRNA.ETOH_nonMamm_2_1.cc1bc2c6e90d9925ca37de13075d512e.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/bvatools/1.6 mugqic/bwa/0.7.12 mugqic/picard/2.0.1 mugqic/mugqic_tools/2.1.7 mugqic/python/2.7.13 && \
mkdir -p alignment/ETOH_nonMamm_2/ETOH_nonMamm_2_1 metrics/ETOH_nonMamm_2/ETOH_nonMamm_2_1 && \
java -XX:ParallelGCThreads=1 -Xmx10G -jar $BVATOOLS_JAR \
  bam2fq --mapped ONLY \
  --bam alignment/ETOH_nonMamm_2/ETOH_nonMamm_2_1/Aligned.sortedByCoord.out.bam    | \
bwa mem  \
  -M -t 10 \
  -R '@RG	ID:ETOH_nonMamm_2_1	SM:ETOH_nonMamm_2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/annotations/rrna_bwa_index/Homo_sapiens.hg19.UCSC2009-03-08.rrna.fa \
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
  -g $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.gtf \
  -o metrics/ETOH_nonMamm_2/ETOH_nonMamm_2_1/ETOH_nonMamm_2_1rRNA.stats.tsv \
  -t transcript
bwa_mem_rRNA.ETOH_nonMamm_2_1.cc1bc2c6e90d9925ca37de13075d512e.mugqic.done
)
estimate_ribosomal_rna_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$estimate_ribosomal_rna_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: estimate_ribosomal_rna_12_JOB_ID: bwa_mem_rRNA.ETOH_shMED1_1_1
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_rRNA.ETOH_shMED1_1_1
JOB_DEPENDENCIES=$star_33_JOB_ID
JOB_DONE=job_output/estimate_ribosomal_rna/bwa_mem_rRNA.ETOH_shMED1_1_1.16b50cf4fa3d3d8b2c2b8555fb13b3d8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_rRNA.ETOH_shMED1_1_1.16b50cf4fa3d3d8b2c2b8555fb13b3d8.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/bvatools/1.6 mugqic/bwa/0.7.12 mugqic/picard/2.0.1 mugqic/mugqic_tools/2.1.7 mugqic/python/2.7.13 && \
mkdir -p alignment/ETOH_shMED1_1/ETOH_shMED1_1_1 metrics/ETOH_shMED1_1/ETOH_shMED1_1_1 && \
java -XX:ParallelGCThreads=1 -Xmx10G -jar $BVATOOLS_JAR \
  bam2fq --mapped ONLY \
  --bam alignment/ETOH_shMED1_1/ETOH_shMED1_1_1/Aligned.sortedByCoord.out.bam    | \
bwa mem  \
  -M -t 10 \
  -R '@RG	ID:ETOH_shMED1_1_1	SM:ETOH_shMED1_1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/annotations/rrna_bwa_index/Homo_sapiens.hg19.UCSC2009-03-08.rrna.fa \
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
  -g $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.gtf \
  -o metrics/ETOH_shMED1_1/ETOH_shMED1_1_1/ETOH_shMED1_1_1rRNA.stats.tsv \
  -t transcript
bwa_mem_rRNA.ETOH_shMED1_1_1.16b50cf4fa3d3d8b2c2b8555fb13b3d8.mugqic.done
)
estimate_ribosomal_rna_12_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$estimate_ribosomal_rna_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: estimate_ribosomal_rna_13_JOB_ID: bwa_mem_rRNA.ETOH_shMED1_1_2
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_rRNA.ETOH_shMED1_1_2
JOB_DEPENDENCIES=$star_34_JOB_ID
JOB_DONE=job_output/estimate_ribosomal_rna/bwa_mem_rRNA.ETOH_shMED1_1_2.34a2851182ce70839435d05c0bfdfb3c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_rRNA.ETOH_shMED1_1_2.34a2851182ce70839435d05c0bfdfb3c.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/bvatools/1.6 mugqic/bwa/0.7.12 mugqic/picard/2.0.1 mugqic/mugqic_tools/2.1.7 mugqic/python/2.7.13 && \
mkdir -p alignment/ETOH_shMED1_1/ETOH_shMED1_1_2 metrics/ETOH_shMED1_1/ETOH_shMED1_1_2 && \
java -XX:ParallelGCThreads=1 -Xmx10G -jar $BVATOOLS_JAR \
  bam2fq --mapped ONLY \
  --bam alignment/ETOH_shMED1_1/ETOH_shMED1_1_2/Aligned.sortedByCoord.out.bam    | \
bwa mem  \
  -M -t 10 \
  -R '@RG	ID:ETOH_shMED1_1_2	SM:ETOH_shMED1_1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/annotations/rrna_bwa_index/Homo_sapiens.hg19.UCSC2009-03-08.rrna.fa \
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
  -g $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.gtf \
  -o metrics/ETOH_shMED1_1/ETOH_shMED1_1_2/ETOH_shMED1_1_2rRNA.stats.tsv \
  -t transcript
bwa_mem_rRNA.ETOH_shMED1_1_2.34a2851182ce70839435d05c0bfdfb3c.mugqic.done
)
estimate_ribosomal_rna_13_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$estimate_ribosomal_rna_13_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: estimate_ribosomal_rna_14_JOB_ID: bwa_mem_rRNA.ETOH_shNIPBL_2_1
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_rRNA.ETOH_shNIPBL_2_1
JOB_DEPENDENCIES=$star_35_JOB_ID
JOB_DONE=job_output/estimate_ribosomal_rna/bwa_mem_rRNA.ETOH_shNIPBL_2_1.ce7b61672c7b38cdf8d39a4b721e9fd1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_rRNA.ETOH_shNIPBL_2_1.ce7b61672c7b38cdf8d39a4b721e9fd1.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/bvatools/1.6 mugqic/bwa/0.7.12 mugqic/picard/2.0.1 mugqic/mugqic_tools/2.1.7 mugqic/python/2.7.13 && \
mkdir -p alignment/ETOH_shNIPBL_2/ETOH_shNIPBL_2_1 metrics/ETOH_shNIPBL_2/ETOH_shNIPBL_2_1 && \
java -XX:ParallelGCThreads=1 -Xmx10G -jar $BVATOOLS_JAR \
  bam2fq --mapped ONLY \
  --bam alignment/ETOH_shNIPBL_2/ETOH_shNIPBL_2_1/Aligned.sortedByCoord.out.bam    | \
bwa mem  \
  -M -t 10 \
  -R '@RG	ID:ETOH_shNIPBL_2_1	SM:ETOH_shNIPBL_2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/annotations/rrna_bwa_index/Homo_sapiens.hg19.UCSC2009-03-08.rrna.fa \
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
  -g $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.gtf \
  -o metrics/ETOH_shNIPBL_2/ETOH_shNIPBL_2_1/ETOH_shNIPBL_2_1rRNA.stats.tsv \
  -t transcript
bwa_mem_rRNA.ETOH_shNIPBL_2_1.ce7b61672c7b38cdf8d39a4b721e9fd1.mugqic.done
)
estimate_ribosomal_rna_14_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$estimate_ribosomal_rna_14_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: estimate_ribosomal_rna_15_JOB_ID: bwa_mem_rRNA.ETOH_shNIPBL_1_1
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_rRNA.ETOH_shNIPBL_1_1
JOB_DEPENDENCIES=$star_36_JOB_ID
JOB_DONE=job_output/estimate_ribosomal_rna/bwa_mem_rRNA.ETOH_shNIPBL_1_1.b65d6695556f18e74907ae85438efaf0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_rRNA.ETOH_shNIPBL_1_1.b65d6695556f18e74907ae85438efaf0.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/bvatools/1.6 mugqic/bwa/0.7.12 mugqic/picard/2.0.1 mugqic/mugqic_tools/2.1.7 mugqic/python/2.7.13 && \
mkdir -p alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1_1 metrics/ETOH_shNIPBL_1/ETOH_shNIPBL_1_1 && \
java -XX:ParallelGCThreads=1 -Xmx10G -jar $BVATOOLS_JAR \
  bam2fq --mapped ONLY \
  --bam alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1_1/Aligned.sortedByCoord.out.bam    | \
bwa mem  \
  -M -t 10 \
  -R '@RG	ID:ETOH_shNIPBL_1_1	SM:ETOH_shNIPBL_1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/annotations/rrna_bwa_index/Homo_sapiens.hg19.UCSC2009-03-08.rrna.fa \
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
  -g $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.gtf \
  -o metrics/ETOH_shNIPBL_1/ETOH_shNIPBL_1_1/ETOH_shNIPBL_1_1rRNA.stats.tsv \
  -t transcript
bwa_mem_rRNA.ETOH_shNIPBL_1_1.b65d6695556f18e74907ae85438efaf0.mugqic.done
)
estimate_ribosomal_rna_15_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$estimate_ribosomal_rna_15_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: estimate_ribosomal_rna_16_JOB_ID: bwa_mem_rRNA.ETOH_shNIPBL_1_2
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_rRNA.ETOH_shNIPBL_1_2
JOB_DEPENDENCIES=$star_37_JOB_ID
JOB_DONE=job_output/estimate_ribosomal_rna/bwa_mem_rRNA.ETOH_shNIPBL_1_2.6887ca561b3df00f9bebcd4e7c9c08bc.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_rRNA.ETOH_shNIPBL_1_2.6887ca561b3df00f9bebcd4e7c9c08bc.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/bvatools/1.6 mugqic/bwa/0.7.12 mugqic/picard/2.0.1 mugqic/mugqic_tools/2.1.7 mugqic/python/2.7.13 && \
mkdir -p alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1_2 metrics/ETOH_shNIPBL_1/ETOH_shNIPBL_1_2 && \
java -XX:ParallelGCThreads=1 -Xmx10G -jar $BVATOOLS_JAR \
  bam2fq --mapped ONLY \
  --bam alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1_2/Aligned.sortedByCoord.out.bam    | \
bwa mem  \
  -M -t 10 \
  -R '@RG	ID:ETOH_shNIPBL_1_2	SM:ETOH_shNIPBL_1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/annotations/rrna_bwa_index/Homo_sapiens.hg19.UCSC2009-03-08.rrna.fa \
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
  -g $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.gtf \
  -o metrics/ETOH_shNIPBL_1/ETOH_shNIPBL_1_2/ETOH_shNIPBL_1_2rRNA.stats.tsv \
  -t transcript
bwa_mem_rRNA.ETOH_shNIPBL_1_2.6887ca561b3df00f9bebcd4e7c9c08bc.mugqic.done
)
estimate_ribosomal_rna_16_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$estimate_ribosomal_rna_16_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: estimate_ribosomal_rna_17_JOB_ID: bwa_mem_rRNA.ETOH_shMED1_2_1
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_rRNA.ETOH_shMED1_2_1
JOB_DEPENDENCIES=$star_38_JOB_ID
JOB_DONE=job_output/estimate_ribosomal_rna/bwa_mem_rRNA.ETOH_shMED1_2_1.ad1f283928a3ca2ee7dca6ae67f5d56a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_rRNA.ETOH_shMED1_2_1.ad1f283928a3ca2ee7dca6ae67f5d56a.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/bvatools/1.6 mugqic/bwa/0.7.12 mugqic/picard/2.0.1 mugqic/mugqic_tools/2.1.7 mugqic/python/2.7.13 && \
mkdir -p alignment/ETOH_shMED1_2/ETOH_shMED1_2_1 metrics/ETOH_shMED1_2/ETOH_shMED1_2_1 && \
java -XX:ParallelGCThreads=1 -Xmx10G -jar $BVATOOLS_JAR \
  bam2fq --mapped ONLY \
  --bam alignment/ETOH_shMED1_2/ETOH_shMED1_2_1/Aligned.sortedByCoord.out.bam    | \
bwa mem  \
  -M -t 10 \
  -R '@RG	ID:ETOH_shMED1_2_1	SM:ETOH_shMED1_2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/annotations/rrna_bwa_index/Homo_sapiens.hg19.UCSC2009-03-08.rrna.fa \
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
  -g $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.gtf \
  -o metrics/ETOH_shMED1_2/ETOH_shMED1_2_1/ETOH_shMED1_2_1rRNA.stats.tsv \
  -t transcript
bwa_mem_rRNA.ETOH_shMED1_2_1.ad1f283928a3ca2ee7dca6ae67f5d56a.mugqic.done
)
estimate_ribosomal_rna_17_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$estimate_ribosomal_rna_17_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: estimate_ribosomal_rna_18_JOB_ID: bwa_mem_rRNA.ETOH_shSMC1A_1_1
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_rRNA.ETOH_shSMC1A_1_1
JOB_DEPENDENCIES=$star_39_JOB_ID
JOB_DONE=job_output/estimate_ribosomal_rna/bwa_mem_rRNA.ETOH_shSMC1A_1_1.766b30f43746e82efc78bdd0bf6ec738.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_rRNA.ETOH_shSMC1A_1_1.766b30f43746e82efc78bdd0bf6ec738.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/bvatools/1.6 mugqic/bwa/0.7.12 mugqic/picard/2.0.1 mugqic/mugqic_tools/2.1.7 mugqic/python/2.7.13 && \
mkdir -p alignment/ETOH_shSMC1A_1/ETOH_shSMC1A_1_1 metrics/ETOH_shSMC1A_1/ETOH_shSMC1A_1_1 && \
java -XX:ParallelGCThreads=1 -Xmx10G -jar $BVATOOLS_JAR \
  bam2fq --mapped ONLY \
  --bam alignment/ETOH_shSMC1A_1/ETOH_shSMC1A_1_1/Aligned.sortedByCoord.out.bam    | \
bwa mem  \
  -M -t 10 \
  -R '@RG	ID:ETOH_shSMC1A_1_1	SM:ETOH_shSMC1A_1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/annotations/rrna_bwa_index/Homo_sapiens.hg19.UCSC2009-03-08.rrna.fa \
  /dev/stdin | \
 java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx7G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=/dev/stdin \
 OUTPUT=metrics/ETOH_shSMC1A_1/ETOH_shSMC1A_1_1/ETOH_shSMC1A_1_1rRNA.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=1750000 && \
python $PYTHON_TOOLS/rrnaBAMcounter.py \
  -i metrics/ETOH_shSMC1A_1/ETOH_shSMC1A_1_1/ETOH_shSMC1A_1_1rRNA.bam \
  -g $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.gtf \
  -o metrics/ETOH_shSMC1A_1/ETOH_shSMC1A_1_1/ETOH_shSMC1A_1_1rRNA.stats.tsv \
  -t transcript
bwa_mem_rRNA.ETOH_shSMC1A_1_1.766b30f43746e82efc78bdd0bf6ec738.mugqic.done
)
estimate_ribosomal_rna_18_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$estimate_ribosomal_rna_18_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: estimate_ribosomal_rna_19_JOB_ID: bwa_mem_rRNA.ETOH_shSMC1A_1_2
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_rRNA.ETOH_shSMC1A_1_2
JOB_DEPENDENCIES=$star_40_JOB_ID
JOB_DONE=job_output/estimate_ribosomal_rna/bwa_mem_rRNA.ETOH_shSMC1A_1_2.895198658907358b33babed988d7dbde.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_rRNA.ETOH_shSMC1A_1_2.895198658907358b33babed988d7dbde.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/bvatools/1.6 mugqic/bwa/0.7.12 mugqic/picard/2.0.1 mugqic/mugqic_tools/2.1.7 mugqic/python/2.7.13 && \
mkdir -p alignment/ETOH_shSMC1A_1/ETOH_shSMC1A_1_2 metrics/ETOH_shSMC1A_1/ETOH_shSMC1A_1_2 && \
java -XX:ParallelGCThreads=1 -Xmx10G -jar $BVATOOLS_JAR \
  bam2fq --mapped ONLY \
  --bam alignment/ETOH_shSMC1A_1/ETOH_shSMC1A_1_2/Aligned.sortedByCoord.out.bam    | \
bwa mem  \
  -M -t 10 \
  -R '@RG	ID:ETOH_shSMC1A_1_2	SM:ETOH_shSMC1A_1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/annotations/rrna_bwa_index/Homo_sapiens.hg19.UCSC2009-03-08.rrna.fa \
  /dev/stdin | \
 java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx7G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=/dev/stdin \
 OUTPUT=metrics/ETOH_shSMC1A_1/ETOH_shSMC1A_1_2/ETOH_shSMC1A_1_2rRNA.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=1750000 && \
python $PYTHON_TOOLS/rrnaBAMcounter.py \
  -i metrics/ETOH_shSMC1A_1/ETOH_shSMC1A_1_2/ETOH_shSMC1A_1_2rRNA.bam \
  -g $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.gtf \
  -o metrics/ETOH_shSMC1A_1/ETOH_shSMC1A_1_2/ETOH_shSMC1A_1_2rRNA.stats.tsv \
  -t transcript
bwa_mem_rRNA.ETOH_shSMC1A_1_2.895198658907358b33babed988d7dbde.mugqic.done
)
estimate_ribosomal_rna_19_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$estimate_ribosomal_rna_19_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: estimate_ribosomal_rna_20_JOB_ID: bwa_mem_rRNA.ETOH_shSMC1A_2_1
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_rRNA.ETOH_shSMC1A_2_1
JOB_DEPENDENCIES=$star_41_JOB_ID
JOB_DONE=job_output/estimate_ribosomal_rna/bwa_mem_rRNA.ETOH_shSMC1A_2_1.d36e5063521e517942b0707a91baf5d0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_rRNA.ETOH_shSMC1A_2_1.d36e5063521e517942b0707a91baf5d0.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/bvatools/1.6 mugqic/bwa/0.7.12 mugqic/picard/2.0.1 mugqic/mugqic_tools/2.1.7 mugqic/python/2.7.13 && \
mkdir -p alignment/ETOH_shSMC1A_2/ETOH_shSMC1A_2_1 metrics/ETOH_shSMC1A_2/ETOH_shSMC1A_2_1 && \
java -XX:ParallelGCThreads=1 -Xmx10G -jar $BVATOOLS_JAR \
  bam2fq --mapped ONLY \
  --bam alignment/ETOH_shSMC1A_2/ETOH_shSMC1A_2_1/Aligned.sortedByCoord.out.bam    | \
bwa mem  \
  -M -t 10 \
  -R '@RG	ID:ETOH_shSMC1A_2_1	SM:ETOH_shSMC1A_2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/annotations/rrna_bwa_index/Homo_sapiens.hg19.UCSC2009-03-08.rrna.fa \
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
  -g $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.gtf \
  -o metrics/ETOH_shSMC1A_2/ETOH_shSMC1A_2_1/ETOH_shSMC1A_2_1rRNA.stats.tsv \
  -t transcript
bwa_mem_rRNA.ETOH_shSMC1A_2_1.d36e5063521e517942b0707a91baf5d0.mugqic.done
)
estimate_ribosomal_rna_20_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bam_hard_clip_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bam_hard_clip_4_JOB_ID: tuxedo_hard_clip.DEX_shSMC1A_1
#-------------------------------------------------------------------------------
JOB_NAME=tuxedo_hard_clip.DEX_shSMC1A_1
JOB_DEPENDENCIES=$picard_mark_duplicates_4_JOB_ID
JOB_DONE=job_output/bam_hard_clip/tuxedo_hard_clip.DEX_shSMC1A_1.574bf4ff47a6a82460112955f51443fd.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'tuxedo_hard_clip.DEX_shSMC1A_1.574bf4ff47a6a82460112955f51443fd.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -h \
  alignment/DEX_shSMC1A_1/DEX_shSMC1A_1.sorted.mdup.bam | \
awk 'BEGIN {OFS="\t"} {if (substr($1,1,1)=="@") {print;next}; split($6,C,/[0-9]*/); split($6,L,/[SMDIN]/); if (C[2]=="S") {$10=substr($10,L[1]+1); $11=substr($11,L[1]+1)}; if (C[length(C)]=="S") {L1=length($10)-L[length(L)-1]; $10=substr($10,1,L1); $11=substr($11,1,L1); }; gsub(/[0-9]*S/,"",$6); print}'  | \
samtools view -hbS \
  - \
  > alignment/DEX_shSMC1A_1/DEX_shSMC1A_1.sorted.mdup.hardClip.bam
tuxedo_hard_clip.DEX_shSMC1A_1.574bf4ff47a6a82460112955f51443fd.mugqic.done
)
bam_hard_clip_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bam_hard_clip_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bam_hard_clip_11_JOB_ID: tuxedo_hard_clip.ETOH_shSMC1A_1
#-------------------------------------------------------------------------------
JOB_NAME=tuxedo_hard_clip.ETOH_shSMC1A_1
JOB_DEPENDENCIES=$picard_mark_duplicates_11_JOB_ID
JOB_DONE=job_output/bam_hard_clip/tuxedo_hard_clip.ETOH_shSMC1A_1.29cb8140950b16436e94ea8afc7f62a9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'tuxedo_hard_clip.ETOH_shSMC1A_1.29cb8140950b16436e94ea8afc7f62a9.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -h \
  alignment/ETOH_shSMC1A_1/ETOH_shSMC1A_1.sorted.mdup.bam | \
awk 'BEGIN {OFS="\t"} {if (substr($1,1,1)=="@") {print;next}; split($6,C,/[0-9]*/); split($6,L,/[SMDIN]/); if (C[2]=="S") {$10=substr($10,L[1]+1); $11=substr($11,L[1]+1)}; if (C[length(C)]=="S") {L1=length($10)-L[length(L)-1]; $10=substr($10,1,L1); $11=substr($11,1,L1); }; gsub(/[0-9]*S/,"",$6); print}'  | \
samtools view -hbS \
  - \
  > alignment/ETOH_shSMC1A_1/ETOH_shSMC1A_1.sorted.mdup.hardClip.bam
tuxedo_hard_clip.ETOH_shSMC1A_1.29cb8140950b16436e94ea8afc7f62a9.mugqic.done
)
bam_hard_clip_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
JOB_DONE=job_output/rnaseqc/rnaseqc.0314459992eca561cfa0678d6ac5145d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'rnaseqc.0314459992eca561cfa0678d6ac5145d.mugqic.done'
module load mugqic/java/openjdk-jdk1.7.0_60 mugqic/bwa/0.7.12 mugqic/rnaseqc/1.1.8 && \
mkdir -p metrics/rnaseqRep && \
echo "Sample	BamFile	Note
DEX_nonMamm_1	alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.bam	RNAseq
DEX_shMED1_1	alignment/DEX_shMED1_1/DEX_shMED1_1.sorted.mdup.bam	RNAseq
DEX_shNIPBL_1	alignment/DEX_shNIPBL_1/DEX_shNIPBL_1.sorted.mdup.bam	RNAseq
DEX_shSMC1A_1	alignment/DEX_shSMC1A_1/DEX_shSMC1A_1.sorted.mdup.bam	RNAseq
ETOH_nonMamm_1	alignment/ETOH_nonMamm_1/ETOH_nonMamm_1.sorted.mdup.bam	RNAseq
ETOH_nonMamm_2	alignment/ETOH_nonMamm_2/ETOH_nonMamm_2.sorted.mdup.bam	RNAseq
ETOH_shMED1_1	alignment/ETOH_shMED1_1/ETOH_shMED1_1.sorted.mdup.bam	RNAseq
ETOH_shNIPBL_2	alignment/ETOH_shNIPBL_2/ETOH_shNIPBL_2.sorted.mdup.bam	RNAseq
ETOH_shNIPBL_1	alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1.sorted.mdup.bam	RNAseq
ETOH_shMED1_2	alignment/ETOH_shMED1_2/ETOH_shMED1_2.sorted.mdup.bam	RNAseq
ETOH_shSMC1A_1	alignment/ETOH_shSMC1A_1/ETOH_shSMC1A_1.sorted.mdup.bam	RNAseq
ETOH_shSMC1A_2	alignment/ETOH_shSMC1A_2/ETOH_shSMC1A_2.sorted.mdup.bam	RNAseq" \
  > alignment/rnaseqc.samples.txt && \
touch dummy_rRNA.fa && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=4 -Xmx27G -jar $RNASEQC_JAR \
  -n 1000 \
  -o metrics/rnaseqRep \
  -r /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa \
  -s alignment/rnaseqc.samples.txt \
  -t /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.transcript_id.gtf \
  -ttype 2\
  -BWArRNA dummy_rRNA.fa && \
zip -r metrics/rnaseqRep.zip metrics/rnaseqRep
rnaseqc.0314459992eca561cfa0678d6ac5145d.mugqic.done
)
rnaseqc_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:00 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
JOB_DONE=job_output/wiggle/wiggle.DEX_nonMamm_1.forward.63840fd1f61e1eacfb96971af7cc6177.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.DEX_nonMamm_1.forward.63840fd1f61e1eacfb96971af7cc6177.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/bedtools/2.25.0 mugqic/ucsc/v326 && \
mkdir -p tracks/DEX_nonMamm_1 tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 81  alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.forward.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.forward.bam \
  -g /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa.fai \
  > tracks/DEX_nonMamm_1/DEX_nonMamm_1.forward.bedGraph && \
sort -k1,1 -k2,2n tracks/DEX_nonMamm_1/DEX_nonMamm_1.forward.bedGraph > tracks/DEX_nonMamm_1/DEX_nonMamm_1.forward.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/DEX_nonMamm_1/DEX_nonMamm_1.forward.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa.fai \
  tracks/bigWig/DEX_nonMamm_1.forward.bw
wiggle.DEX_nonMamm_1.forward.63840fd1f61e1eacfb96971af7cc6177.mugqic.done
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
JOB_DONE=job_output/wiggle/wiggle.DEX_nonMamm_1.reverse.e559b6a49cd836631ec43c8bf9c53ab2.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.DEX_nonMamm_1.reverse.e559b6a49cd836631ec43c8bf9c53ab2.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/bedtools/2.25.0 mugqic/ucsc/v326 && \
mkdir -p tracks/DEX_nonMamm_1 tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 97  alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.reverse.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.reverse.bam \
  -g /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa.fai \
  > tracks/DEX_nonMamm_1/DEX_nonMamm_1.reverse.bedGraph && \
sort -k1,1 -k2,2n tracks/DEX_nonMamm_1/DEX_nonMamm_1.reverse.bedGraph > tracks/DEX_nonMamm_1/DEX_nonMamm_1.reverse.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/DEX_nonMamm_1/DEX_nonMamm_1.reverse.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa.fai \
  tracks/bigWig/DEX_nonMamm_1.reverse.bw
wiggle.DEX_nonMamm_1.reverse.e559b6a49cd836631ec43c8bf9c53ab2.mugqic.done
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
JOB_DONE=job_output/wiggle/wiggle.DEX_shMED1_1.forward.9bbea9dd0f211cee045782d4c2f44fe0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.DEX_shMED1_1.forward.9bbea9dd0f211cee045782d4c2f44fe0.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/bedtools/2.25.0 mugqic/ucsc/v326 && \
mkdir -p tracks/DEX_shMED1_1 tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 81  alignment/DEX_shMED1_1/DEX_shMED1_1.sorted.mdup.forward.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/DEX_shMED1_1/DEX_shMED1_1.sorted.mdup.forward.bam \
  -g /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa.fai \
  > tracks/DEX_shMED1_1/DEX_shMED1_1.forward.bedGraph && \
sort -k1,1 -k2,2n tracks/DEX_shMED1_1/DEX_shMED1_1.forward.bedGraph > tracks/DEX_shMED1_1/DEX_shMED1_1.forward.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/DEX_shMED1_1/DEX_shMED1_1.forward.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa.fai \
  tracks/bigWig/DEX_shMED1_1.forward.bw
wiggle.DEX_shMED1_1.forward.9bbea9dd0f211cee045782d4c2f44fe0.mugqic.done
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
JOB_DONE=job_output/wiggle/wiggle.DEX_shMED1_1.reverse.af28b06832ae349402703410fcbf098a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.DEX_shMED1_1.reverse.af28b06832ae349402703410fcbf098a.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/bedtools/2.25.0 mugqic/ucsc/v326 && \
mkdir -p tracks/DEX_shMED1_1 tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 97  alignment/DEX_shMED1_1/DEX_shMED1_1.sorted.mdup.reverse.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/DEX_shMED1_1/DEX_shMED1_1.sorted.mdup.reverse.bam \
  -g /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa.fai \
  > tracks/DEX_shMED1_1/DEX_shMED1_1.reverse.bedGraph && \
sort -k1,1 -k2,2n tracks/DEX_shMED1_1/DEX_shMED1_1.reverse.bedGraph > tracks/DEX_shMED1_1/DEX_shMED1_1.reverse.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/DEX_shMED1_1/DEX_shMED1_1.reverse.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa.fai \
  tracks/bigWig/DEX_shMED1_1.reverse.bw
wiggle.DEX_shMED1_1.reverse.af28b06832ae349402703410fcbf098a.mugqic.done
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
JOB_DONE=job_output/wiggle/wiggle.DEX_shNIPBL_1.forward.ac0c6303600353804eb44deb026b9c31.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.DEX_shNIPBL_1.forward.ac0c6303600353804eb44deb026b9c31.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/bedtools/2.25.0 mugqic/ucsc/v326 && \
mkdir -p tracks/DEX_shNIPBL_1 tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 81  alignment/DEX_shNIPBL_1/DEX_shNIPBL_1.sorted.mdup.forward.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/DEX_shNIPBL_1/DEX_shNIPBL_1.sorted.mdup.forward.bam \
  -g /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa.fai \
  > tracks/DEX_shNIPBL_1/DEX_shNIPBL_1.forward.bedGraph && \
sort -k1,1 -k2,2n tracks/DEX_shNIPBL_1/DEX_shNIPBL_1.forward.bedGraph > tracks/DEX_shNIPBL_1/DEX_shNIPBL_1.forward.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/DEX_shNIPBL_1/DEX_shNIPBL_1.forward.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa.fai \
  tracks/bigWig/DEX_shNIPBL_1.forward.bw
wiggle.DEX_shNIPBL_1.forward.ac0c6303600353804eb44deb026b9c31.mugqic.done
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
JOB_DONE=job_output/wiggle/wiggle.DEX_shNIPBL_1.reverse.c50015e045a5ae1b78636f14aa691e06.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.DEX_shNIPBL_1.reverse.c50015e045a5ae1b78636f14aa691e06.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/bedtools/2.25.0 mugqic/ucsc/v326 && \
mkdir -p tracks/DEX_shNIPBL_1 tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 97  alignment/DEX_shNIPBL_1/DEX_shNIPBL_1.sorted.mdup.reverse.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/DEX_shNIPBL_1/DEX_shNIPBL_1.sorted.mdup.reverse.bam \
  -g /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa.fai \
  > tracks/DEX_shNIPBL_1/DEX_shNIPBL_1.reverse.bedGraph && \
sort -k1,1 -k2,2n tracks/DEX_shNIPBL_1/DEX_shNIPBL_1.reverse.bedGraph > tracks/DEX_shNIPBL_1/DEX_shNIPBL_1.reverse.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/DEX_shNIPBL_1/DEX_shNIPBL_1.reverse.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa.fai \
  tracks/bigWig/DEX_shNIPBL_1.reverse.bw
wiggle.DEX_shNIPBL_1.reverse.c50015e045a5ae1b78636f14aa691e06.mugqic.done
)
wiggle_12_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_13_JOB_ID: wiggle.DEX_shSMC1A_1.forward_strandspec
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.DEX_shSMC1A_1.forward_strandspec
JOB_DEPENDENCIES=$picard_mark_duplicates_4_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.DEX_shSMC1A_1.forward_strandspec.8076f37a19bb1db9b670b49cbf007ee5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.DEX_shSMC1A_1.forward_strandspec.8076f37a19bb1db9b670b49cbf007ee5.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
samtools view -bh -F 256 -f 81 \
  alignment/DEX_shSMC1A_1/DEX_shSMC1A_1.sorted.mdup.bam \
  > alignment/DEX_shSMC1A_1/DEX_shSMC1A_1.sorted.mdup.tmp1.forward.bam && \
samtools view -bh -F 256 -f 161 \
  alignment/DEX_shSMC1A_1/DEX_shSMC1A_1.sorted.mdup.bam \
  > alignment/DEX_shSMC1A_1/DEX_shSMC1A_1.sorted.mdup.tmp2.forward.bam && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/picard.jar MergeSamFiles \
 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/DEX_shSMC1A_1/DEX_shSMC1A_1.sorted.mdup.tmp1.forward.bam \
  INPUT=alignment/DEX_shSMC1A_1/DEX_shSMC1A_1.sorted.mdup.tmp2.forward.bam \
 OUTPUT=alignment/DEX_shSMC1A_1/DEX_shSMC1A_1.sorted.mdup.forward.bam \
 MAX_RECORDS_IN_RAM=5750000 && \
rm alignment/DEX_shSMC1A_1/DEX_shSMC1A_1.sorted.mdup.tmp1.forward.bam alignment/DEX_shSMC1A_1/DEX_shSMC1A_1.sorted.mdup.tmp2.forward.bam
wiggle.DEX_shSMC1A_1.forward_strandspec.8076f37a19bb1db9b670b49cbf007ee5.mugqic.done
)
wiggle_13_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_13_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_14_JOB_ID: wiggle.DEX_shSMC1A_1.reverse_strandspec
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.DEX_shSMC1A_1.reverse_strandspec
JOB_DEPENDENCIES=$picard_mark_duplicates_4_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.DEX_shSMC1A_1.reverse_strandspec.36fa8c7b202948904194f24eb3f38268.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.DEX_shSMC1A_1.reverse_strandspec.36fa8c7b202948904194f24eb3f38268.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p tracks/DEX_shSMC1A_1 tracks/bigWig && \
samtools view -bh -F 256 -f 97 \
  alignment/DEX_shSMC1A_1/DEX_shSMC1A_1.sorted.mdup.bam \
  > alignment/DEX_shSMC1A_1/DEX_shSMC1A_1.sorted.mdup.tmp1.reverse.bam && \
samtools view -bh -F 256 -f 145 \
  alignment/DEX_shSMC1A_1/DEX_shSMC1A_1.sorted.mdup.bam \
  > alignment/DEX_shSMC1A_1/DEX_shSMC1A_1.sorted.mdup.tmp2.reverse.bam && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/picard.jar MergeSamFiles \
 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/DEX_shSMC1A_1/DEX_shSMC1A_1.sorted.mdup.tmp1.reverse.bam \
  INPUT=alignment/DEX_shSMC1A_1/DEX_shSMC1A_1.sorted.mdup.tmp2.reverse.bam \
 OUTPUT=alignment/DEX_shSMC1A_1/DEX_shSMC1A_1.sorted.mdup.reverse.bam \
 MAX_RECORDS_IN_RAM=5750000 && \
rm alignment/DEX_shSMC1A_1/DEX_shSMC1A_1.sorted.mdup.tmp1.reverse.bam alignment/DEX_shSMC1A_1/DEX_shSMC1A_1.sorted.mdup.tmp2.reverse.bam
wiggle.DEX_shSMC1A_1.reverse_strandspec.36fa8c7b202948904194f24eb3f38268.mugqic.done
)
wiggle_14_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_14_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_15_JOB_ID: wiggle.DEX_shSMC1A_1.forward
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.DEX_shSMC1A_1.forward
JOB_DEPENDENCIES=$wiggle_13_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.DEX_shSMC1A_1.forward.59bfbd313cbd238bc98a0358396e82ee.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.DEX_shSMC1A_1.forward.59bfbd313cbd238bc98a0358396e82ee.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/bedtools/2.25.0 mugqic/ucsc/v326 && \
mkdir -p tracks/DEX_shSMC1A_1 tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 81  alignment/DEX_shSMC1A_1/DEX_shSMC1A_1.sorted.mdup.forward.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/DEX_shSMC1A_1/DEX_shSMC1A_1.sorted.mdup.forward.bam \
  -g /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa.fai \
  > tracks/DEX_shSMC1A_1/DEX_shSMC1A_1.forward.bedGraph && \
sort -k1,1 -k2,2n tracks/DEX_shSMC1A_1/DEX_shSMC1A_1.forward.bedGraph > tracks/DEX_shSMC1A_1/DEX_shSMC1A_1.forward.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/DEX_shSMC1A_1/DEX_shSMC1A_1.forward.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa.fai \
  tracks/bigWig/DEX_shSMC1A_1.forward.bw
wiggle.DEX_shSMC1A_1.forward.59bfbd313cbd238bc98a0358396e82ee.mugqic.done
)
wiggle_15_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_15_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_16_JOB_ID: wiggle.DEX_shSMC1A_1.reverse
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.DEX_shSMC1A_1.reverse
JOB_DEPENDENCIES=$wiggle_14_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.DEX_shSMC1A_1.reverse.9ff78611e519cfdddf46ff715340666b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.DEX_shSMC1A_1.reverse.9ff78611e519cfdddf46ff715340666b.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/bedtools/2.25.0 mugqic/ucsc/v326 && \
mkdir -p tracks/DEX_shSMC1A_1 tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 97  alignment/DEX_shSMC1A_1/DEX_shSMC1A_1.sorted.mdup.reverse.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/DEX_shSMC1A_1/DEX_shSMC1A_1.sorted.mdup.reverse.bam \
  -g /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa.fai \
  > tracks/DEX_shSMC1A_1/DEX_shSMC1A_1.reverse.bedGraph && \
sort -k1,1 -k2,2n tracks/DEX_shSMC1A_1/DEX_shSMC1A_1.reverse.bedGraph > tracks/DEX_shSMC1A_1/DEX_shSMC1A_1.reverse.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/DEX_shSMC1A_1/DEX_shSMC1A_1.reverse.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa.fai \
  tracks/bigWig/DEX_shSMC1A_1.reverse.bw
wiggle.DEX_shSMC1A_1.reverse.9ff78611e519cfdddf46ff715340666b.mugqic.done
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
JOB_DONE=job_output/wiggle/wiggle.ETOH_nonMamm_1.forward.016658d29b0f9783affa16a8b1b6a04e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.ETOH_nonMamm_1.forward.016658d29b0f9783affa16a8b1b6a04e.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/bedtools/2.25.0 mugqic/ucsc/v326 && \
mkdir -p tracks/ETOH_nonMamm_1 tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 81  alignment/ETOH_nonMamm_1/ETOH_nonMamm_1.sorted.mdup.forward.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/ETOH_nonMamm_1/ETOH_nonMamm_1.sorted.mdup.forward.bam \
  -g /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa.fai \
  > tracks/ETOH_nonMamm_1/ETOH_nonMamm_1.forward.bedGraph && \
sort -k1,1 -k2,2n tracks/ETOH_nonMamm_1/ETOH_nonMamm_1.forward.bedGraph > tracks/ETOH_nonMamm_1/ETOH_nonMamm_1.forward.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/ETOH_nonMamm_1/ETOH_nonMamm_1.forward.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa.fai \
  tracks/bigWig/ETOH_nonMamm_1.forward.bw
wiggle.ETOH_nonMamm_1.forward.016658d29b0f9783affa16a8b1b6a04e.mugqic.done
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
JOB_DONE=job_output/wiggle/wiggle.ETOH_nonMamm_1.reverse.eb45c5345891148b7ae990b9dd3c950d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.ETOH_nonMamm_1.reverse.eb45c5345891148b7ae990b9dd3c950d.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/bedtools/2.25.0 mugqic/ucsc/v326 && \
mkdir -p tracks/ETOH_nonMamm_1 tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 97  alignment/ETOH_nonMamm_1/ETOH_nonMamm_1.sorted.mdup.reverse.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/ETOH_nonMamm_1/ETOH_nonMamm_1.sorted.mdup.reverse.bam \
  -g /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa.fai \
  > tracks/ETOH_nonMamm_1/ETOH_nonMamm_1.reverse.bedGraph && \
sort -k1,1 -k2,2n tracks/ETOH_nonMamm_1/ETOH_nonMamm_1.reverse.bedGraph > tracks/ETOH_nonMamm_1/ETOH_nonMamm_1.reverse.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/ETOH_nonMamm_1/ETOH_nonMamm_1.reverse.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa.fai \
  tracks/bigWig/ETOH_nonMamm_1.reverse.bw
wiggle.ETOH_nonMamm_1.reverse.eb45c5345891148b7ae990b9dd3c950d.mugqic.done
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
JOB_DONE=job_output/wiggle/wiggle.ETOH_nonMamm_2.forward.3c1332ea38a755a16cfec3b918af562d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.ETOH_nonMamm_2.forward.3c1332ea38a755a16cfec3b918af562d.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/bedtools/2.25.0 mugqic/ucsc/v326 && \
mkdir -p tracks/ETOH_nonMamm_2 tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 81  alignment/ETOH_nonMamm_2/ETOH_nonMamm_2.sorted.mdup.forward.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/ETOH_nonMamm_2/ETOH_nonMamm_2.sorted.mdup.forward.bam \
  -g /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa.fai \
  > tracks/ETOH_nonMamm_2/ETOH_nonMamm_2.forward.bedGraph && \
sort -k1,1 -k2,2n tracks/ETOH_nonMamm_2/ETOH_nonMamm_2.forward.bedGraph > tracks/ETOH_nonMamm_2/ETOH_nonMamm_2.forward.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/ETOH_nonMamm_2/ETOH_nonMamm_2.forward.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa.fai \
  tracks/bigWig/ETOH_nonMamm_2.forward.bw
wiggle.ETOH_nonMamm_2.forward.3c1332ea38a755a16cfec3b918af562d.mugqic.done
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
JOB_DONE=job_output/wiggle/wiggle.ETOH_nonMamm_2.reverse.dd8cd64a12a5fa8a408dd462bc3b6889.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.ETOH_nonMamm_2.reverse.dd8cd64a12a5fa8a408dd462bc3b6889.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/bedtools/2.25.0 mugqic/ucsc/v326 && \
mkdir -p tracks/ETOH_nonMamm_2 tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 97  alignment/ETOH_nonMamm_2/ETOH_nonMamm_2.sorted.mdup.reverse.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/ETOH_nonMamm_2/ETOH_nonMamm_2.sorted.mdup.reverse.bam \
  -g /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa.fai \
  > tracks/ETOH_nonMamm_2/ETOH_nonMamm_2.reverse.bedGraph && \
sort -k1,1 -k2,2n tracks/ETOH_nonMamm_2/ETOH_nonMamm_2.reverse.bedGraph > tracks/ETOH_nonMamm_2/ETOH_nonMamm_2.reverse.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/ETOH_nonMamm_2/ETOH_nonMamm_2.reverse.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa.fai \
  tracks/bigWig/ETOH_nonMamm_2.reverse.bw
wiggle.ETOH_nonMamm_2.reverse.dd8cd64a12a5fa8a408dd462bc3b6889.mugqic.done
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
JOB_DONE=job_output/wiggle/wiggle.ETOH_shMED1_1.forward.b6b0f2a4b1a7fbf0a250472f0d976dc6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.ETOH_shMED1_1.forward.b6b0f2a4b1a7fbf0a250472f0d976dc6.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/bedtools/2.25.0 mugqic/ucsc/v326 && \
mkdir -p tracks/ETOH_shMED1_1 tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 81  alignment/ETOH_shMED1_1/ETOH_shMED1_1.sorted.mdup.forward.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/ETOH_shMED1_1/ETOH_shMED1_1.sorted.mdup.forward.bam \
  -g /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa.fai \
  > tracks/ETOH_shMED1_1/ETOH_shMED1_1.forward.bedGraph && \
sort -k1,1 -k2,2n tracks/ETOH_shMED1_1/ETOH_shMED1_1.forward.bedGraph > tracks/ETOH_shMED1_1/ETOH_shMED1_1.forward.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/ETOH_shMED1_1/ETOH_shMED1_1.forward.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa.fai \
  tracks/bigWig/ETOH_shMED1_1.forward.bw
wiggle.ETOH_shMED1_1.forward.b6b0f2a4b1a7fbf0a250472f0d976dc6.mugqic.done
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
JOB_DONE=job_output/wiggle/wiggle.ETOH_shMED1_1.reverse.e7b41a7a9c48849b0accde54b124c5d6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.ETOH_shMED1_1.reverse.e7b41a7a9c48849b0accde54b124c5d6.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/bedtools/2.25.0 mugqic/ucsc/v326 && \
mkdir -p tracks/ETOH_shMED1_1 tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 97  alignment/ETOH_shMED1_1/ETOH_shMED1_1.sorted.mdup.reverse.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/ETOH_shMED1_1/ETOH_shMED1_1.sorted.mdup.reverse.bam \
  -g /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa.fai \
  > tracks/ETOH_shMED1_1/ETOH_shMED1_1.reverse.bedGraph && \
sort -k1,1 -k2,2n tracks/ETOH_shMED1_1/ETOH_shMED1_1.reverse.bedGraph > tracks/ETOH_shMED1_1/ETOH_shMED1_1.reverse.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/ETOH_shMED1_1/ETOH_shMED1_1.reverse.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa.fai \
  tracks/bigWig/ETOH_shMED1_1.reverse.bw
wiggle.ETOH_shMED1_1.reverse.e7b41a7a9c48849b0accde54b124c5d6.mugqic.done
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
JOB_DONE=job_output/wiggle/wiggle.ETOH_shNIPBL_2.forward.4648ed8549223378099c48bdffa98b3e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.ETOH_shNIPBL_2.forward.4648ed8549223378099c48bdffa98b3e.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/bedtools/2.25.0 mugqic/ucsc/v326 && \
mkdir -p tracks/ETOH_shNIPBL_2 tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 81  alignment/ETOH_shNIPBL_2/ETOH_shNIPBL_2.sorted.mdup.forward.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/ETOH_shNIPBL_2/ETOH_shNIPBL_2.sorted.mdup.forward.bam \
  -g /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa.fai \
  > tracks/ETOH_shNIPBL_2/ETOH_shNIPBL_2.forward.bedGraph && \
sort -k1,1 -k2,2n tracks/ETOH_shNIPBL_2/ETOH_shNIPBL_2.forward.bedGraph > tracks/ETOH_shNIPBL_2/ETOH_shNIPBL_2.forward.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/ETOH_shNIPBL_2/ETOH_shNIPBL_2.forward.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa.fai \
  tracks/bigWig/ETOH_shNIPBL_2.forward.bw
wiggle.ETOH_shNIPBL_2.forward.4648ed8549223378099c48bdffa98b3e.mugqic.done
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
JOB_DONE=job_output/wiggle/wiggle.ETOH_shNIPBL_2.reverse.8a529c61d1f7c537cb72b1ed5b192fb0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.ETOH_shNIPBL_2.reverse.8a529c61d1f7c537cb72b1ed5b192fb0.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/bedtools/2.25.0 mugqic/ucsc/v326 && \
mkdir -p tracks/ETOH_shNIPBL_2 tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 97  alignment/ETOH_shNIPBL_2/ETOH_shNIPBL_2.sorted.mdup.reverse.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/ETOH_shNIPBL_2/ETOH_shNIPBL_2.sorted.mdup.reverse.bam \
  -g /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa.fai \
  > tracks/ETOH_shNIPBL_2/ETOH_shNIPBL_2.reverse.bedGraph && \
sort -k1,1 -k2,2n tracks/ETOH_shNIPBL_2/ETOH_shNIPBL_2.reverse.bedGraph > tracks/ETOH_shNIPBL_2/ETOH_shNIPBL_2.reverse.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/ETOH_shNIPBL_2/ETOH_shNIPBL_2.reverse.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa.fai \
  tracks/bigWig/ETOH_shNIPBL_2.reverse.bw
wiggle.ETOH_shNIPBL_2.reverse.8a529c61d1f7c537cb72b1ed5b192fb0.mugqic.done
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
JOB_DONE=job_output/wiggle/wiggle.ETOH_shNIPBL_1.forward.a0adf405454644d3ee04ba087976b62a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.ETOH_shNIPBL_1.forward.a0adf405454644d3ee04ba087976b62a.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/bedtools/2.25.0 mugqic/ucsc/v326 && \
mkdir -p tracks/ETOH_shNIPBL_1 tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 81  alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1.sorted.mdup.forward.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1.sorted.mdup.forward.bam \
  -g /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa.fai \
  > tracks/ETOH_shNIPBL_1/ETOH_shNIPBL_1.forward.bedGraph && \
sort -k1,1 -k2,2n tracks/ETOH_shNIPBL_1/ETOH_shNIPBL_1.forward.bedGraph > tracks/ETOH_shNIPBL_1/ETOH_shNIPBL_1.forward.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/ETOH_shNIPBL_1/ETOH_shNIPBL_1.forward.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa.fai \
  tracks/bigWig/ETOH_shNIPBL_1.forward.bw
wiggle.ETOH_shNIPBL_1.forward.a0adf405454644d3ee04ba087976b62a.mugqic.done
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
JOB_DONE=job_output/wiggle/wiggle.ETOH_shNIPBL_1.reverse.777ef3767fa23999d51aa54f2943b85f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.ETOH_shNIPBL_1.reverse.777ef3767fa23999d51aa54f2943b85f.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/bedtools/2.25.0 mugqic/ucsc/v326 && \
mkdir -p tracks/ETOH_shNIPBL_1 tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 97  alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1.sorted.mdup.reverse.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1.sorted.mdup.reverse.bam \
  -g /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa.fai \
  > tracks/ETOH_shNIPBL_1/ETOH_shNIPBL_1.reverse.bedGraph && \
sort -k1,1 -k2,2n tracks/ETOH_shNIPBL_1/ETOH_shNIPBL_1.reverse.bedGraph > tracks/ETOH_shNIPBL_1/ETOH_shNIPBL_1.reverse.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/ETOH_shNIPBL_1/ETOH_shNIPBL_1.reverse.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa.fai \
  tracks/bigWig/ETOH_shNIPBL_1.reverse.bw
wiggle.ETOH_shNIPBL_1.reverse.777ef3767fa23999d51aa54f2943b85f.mugqic.done
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
JOB_DONE=job_output/wiggle/wiggle.ETOH_shMED1_2.forward.b85c19e2a717ba6d61dba65c7218a507.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.ETOH_shMED1_2.forward.b85c19e2a717ba6d61dba65c7218a507.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/bedtools/2.25.0 mugqic/ucsc/v326 && \
mkdir -p tracks/ETOH_shMED1_2 tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 81  alignment/ETOH_shMED1_2/ETOH_shMED1_2.sorted.mdup.forward.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/ETOH_shMED1_2/ETOH_shMED1_2.sorted.mdup.forward.bam \
  -g /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa.fai \
  > tracks/ETOH_shMED1_2/ETOH_shMED1_2.forward.bedGraph && \
sort -k1,1 -k2,2n tracks/ETOH_shMED1_2/ETOH_shMED1_2.forward.bedGraph > tracks/ETOH_shMED1_2/ETOH_shMED1_2.forward.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/ETOH_shMED1_2/ETOH_shMED1_2.forward.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa.fai \
  tracks/bigWig/ETOH_shMED1_2.forward.bw
wiggle.ETOH_shMED1_2.forward.b85c19e2a717ba6d61dba65c7218a507.mugqic.done
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
JOB_DONE=job_output/wiggle/wiggle.ETOH_shMED1_2.reverse.1e180227d9ca666433c1d1f9ea090111.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.ETOH_shMED1_2.reverse.1e180227d9ca666433c1d1f9ea090111.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/bedtools/2.25.0 mugqic/ucsc/v326 && \
mkdir -p tracks/ETOH_shMED1_2 tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 97  alignment/ETOH_shMED1_2/ETOH_shMED1_2.sorted.mdup.reverse.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/ETOH_shMED1_2/ETOH_shMED1_2.sorted.mdup.reverse.bam \
  -g /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa.fai \
  > tracks/ETOH_shMED1_2/ETOH_shMED1_2.reverse.bedGraph && \
sort -k1,1 -k2,2n tracks/ETOH_shMED1_2/ETOH_shMED1_2.reverse.bedGraph > tracks/ETOH_shMED1_2/ETOH_shMED1_2.reverse.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/ETOH_shMED1_2/ETOH_shMED1_2.reverse.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa.fai \
  tracks/bigWig/ETOH_shMED1_2.reverse.bw
wiggle.ETOH_shMED1_2.reverse.1e180227d9ca666433c1d1f9ea090111.mugqic.done
)
wiggle_40_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_40_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_41_JOB_ID: wiggle.ETOH_shSMC1A_1.forward_strandspec
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.ETOH_shSMC1A_1.forward_strandspec
JOB_DEPENDENCIES=$picard_mark_duplicates_11_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.ETOH_shSMC1A_1.forward_strandspec.86b34dc981bd6148bce7df2c5ae58de6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.ETOH_shSMC1A_1.forward_strandspec.86b34dc981bd6148bce7df2c5ae58de6.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
samtools view -bh -F 256 -f 81 \
  alignment/ETOH_shSMC1A_1/ETOH_shSMC1A_1.sorted.mdup.bam \
  > alignment/ETOH_shSMC1A_1/ETOH_shSMC1A_1.sorted.mdup.tmp1.forward.bam && \
samtools view -bh -F 256 -f 161 \
  alignment/ETOH_shSMC1A_1/ETOH_shSMC1A_1.sorted.mdup.bam \
  > alignment/ETOH_shSMC1A_1/ETOH_shSMC1A_1.sorted.mdup.tmp2.forward.bam && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/picard.jar MergeSamFiles \
 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/ETOH_shSMC1A_1/ETOH_shSMC1A_1.sorted.mdup.tmp1.forward.bam \
  INPUT=alignment/ETOH_shSMC1A_1/ETOH_shSMC1A_1.sorted.mdup.tmp2.forward.bam \
 OUTPUT=alignment/ETOH_shSMC1A_1/ETOH_shSMC1A_1.sorted.mdup.forward.bam \
 MAX_RECORDS_IN_RAM=5750000 && \
rm alignment/ETOH_shSMC1A_1/ETOH_shSMC1A_1.sorted.mdup.tmp1.forward.bam alignment/ETOH_shSMC1A_1/ETOH_shSMC1A_1.sorted.mdup.tmp2.forward.bam
wiggle.ETOH_shSMC1A_1.forward_strandspec.86b34dc981bd6148bce7df2c5ae58de6.mugqic.done
)
wiggle_41_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_41_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_42_JOB_ID: wiggle.ETOH_shSMC1A_1.reverse_strandspec
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.ETOH_shSMC1A_1.reverse_strandspec
JOB_DEPENDENCIES=$picard_mark_duplicates_11_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.ETOH_shSMC1A_1.reverse_strandspec.f14e4607f81b16e619c87e7a61e989b2.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.ETOH_shSMC1A_1.reverse_strandspec.f14e4607f81b16e619c87e7a61e989b2.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/2.0.1 && \
mkdir -p tracks/ETOH_shSMC1A_1 tracks/bigWig && \
samtools view -bh -F 256 -f 97 \
  alignment/ETOH_shSMC1A_1/ETOH_shSMC1A_1.sorted.mdup.bam \
  > alignment/ETOH_shSMC1A_1/ETOH_shSMC1A_1.sorted.mdup.tmp1.reverse.bam && \
samtools view -bh -F 256 -f 145 \
  alignment/ETOH_shSMC1A_1/ETOH_shSMC1A_1.sorted.mdup.bam \
  > alignment/ETOH_shSMC1A_1/ETOH_shSMC1A_1.sorted.mdup.tmp2.reverse.bam && \
java -Djava.io.tmpdir=/gs/scratch/$USER -XX:ParallelGCThreads=1 -Xmx27G -jar $PICARD_HOME/picard.jar MergeSamFiles \
 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
 TMP_DIR=/gs/scratch/$USER \
 INPUT=alignment/ETOH_shSMC1A_1/ETOH_shSMC1A_1.sorted.mdup.tmp1.reverse.bam \
  INPUT=alignment/ETOH_shSMC1A_1/ETOH_shSMC1A_1.sorted.mdup.tmp2.reverse.bam \
 OUTPUT=alignment/ETOH_shSMC1A_1/ETOH_shSMC1A_1.sorted.mdup.reverse.bam \
 MAX_RECORDS_IN_RAM=5750000 && \
rm alignment/ETOH_shSMC1A_1/ETOH_shSMC1A_1.sorted.mdup.tmp1.reverse.bam alignment/ETOH_shSMC1A_1/ETOH_shSMC1A_1.sorted.mdup.tmp2.reverse.bam
wiggle.ETOH_shSMC1A_1.reverse_strandspec.f14e4607f81b16e619c87e7a61e989b2.mugqic.done
)
wiggle_42_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_42_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_43_JOB_ID: wiggle.ETOH_shSMC1A_1.forward
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.ETOH_shSMC1A_1.forward
JOB_DEPENDENCIES=$wiggle_41_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.ETOH_shSMC1A_1.forward.6d3ee9039c7c4c2d6b4d5cbbe2262d4b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.ETOH_shSMC1A_1.forward.6d3ee9039c7c4c2d6b4d5cbbe2262d4b.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/bedtools/2.25.0 mugqic/ucsc/v326 && \
mkdir -p tracks/ETOH_shSMC1A_1 tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 81  alignment/ETOH_shSMC1A_1/ETOH_shSMC1A_1.sorted.mdup.forward.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/ETOH_shSMC1A_1/ETOH_shSMC1A_1.sorted.mdup.forward.bam \
  -g /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa.fai \
  > tracks/ETOH_shSMC1A_1/ETOH_shSMC1A_1.forward.bedGraph && \
sort -k1,1 -k2,2n tracks/ETOH_shSMC1A_1/ETOH_shSMC1A_1.forward.bedGraph > tracks/ETOH_shSMC1A_1/ETOH_shSMC1A_1.forward.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/ETOH_shSMC1A_1/ETOH_shSMC1A_1.forward.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa.fai \
  tracks/bigWig/ETOH_shSMC1A_1.forward.bw
wiggle.ETOH_shSMC1A_1.forward.6d3ee9039c7c4c2d6b4d5cbbe2262d4b.mugqic.done
)
wiggle_43_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$wiggle_43_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: wiggle_44_JOB_ID: wiggle.ETOH_shSMC1A_1.reverse
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.ETOH_shSMC1A_1.reverse
JOB_DEPENDENCIES=$wiggle_42_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.ETOH_shSMC1A_1.reverse.c9e005bc7e6c3b6417fcda7af7159424.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.ETOH_shSMC1A_1.reverse.c9e005bc7e6c3b6417fcda7af7159424.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/bedtools/2.25.0 mugqic/ucsc/v326 && \
mkdir -p tracks/ETOH_shSMC1A_1 tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 97  alignment/ETOH_shSMC1A_1/ETOH_shSMC1A_1.sorted.mdup.reverse.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/ETOH_shSMC1A_1/ETOH_shSMC1A_1.sorted.mdup.reverse.bam \
  -g /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa.fai \
  > tracks/ETOH_shSMC1A_1/ETOH_shSMC1A_1.reverse.bedGraph && \
sort -k1,1 -k2,2n tracks/ETOH_shSMC1A_1/ETOH_shSMC1A_1.reverse.bedGraph > tracks/ETOH_shSMC1A_1/ETOH_shSMC1A_1.reverse.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/ETOH_shSMC1A_1/ETOH_shSMC1A_1.reverse.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa.fai \
  tracks/bigWig/ETOH_shSMC1A_1.reverse.bw
wiggle.ETOH_shSMC1A_1.reverse.c9e005bc7e6c3b6417fcda7af7159424.mugqic.done
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
JOB_DONE=job_output/wiggle/wiggle.ETOH_shSMC1A_2.forward.1d23dad2de92b6db72e3c6d3ebb633b8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.ETOH_shSMC1A_2.forward.1d23dad2de92b6db72e3c6d3ebb633b8.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/bedtools/2.25.0 mugqic/ucsc/v326 && \
mkdir -p tracks/ETOH_shSMC1A_2 tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 81  alignment/ETOH_shSMC1A_2/ETOH_shSMC1A_2.sorted.mdup.forward.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/ETOH_shSMC1A_2/ETOH_shSMC1A_2.sorted.mdup.forward.bam \
  -g /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa.fai \
  > tracks/ETOH_shSMC1A_2/ETOH_shSMC1A_2.forward.bedGraph && \
sort -k1,1 -k2,2n tracks/ETOH_shSMC1A_2/ETOH_shSMC1A_2.forward.bedGraph > tracks/ETOH_shSMC1A_2/ETOH_shSMC1A_2.forward.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/ETOH_shSMC1A_2/ETOH_shSMC1A_2.forward.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa.fai \
  tracks/bigWig/ETOH_shSMC1A_2.forward.bw
wiggle.ETOH_shSMC1A_2.forward.1d23dad2de92b6db72e3c6d3ebb633b8.mugqic.done
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
JOB_DONE=job_output/wiggle/wiggle.ETOH_shSMC1A_2.reverse.c37850b5994b6cc87c59fe92b9587322.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.ETOH_shSMC1A_2.reverse.c37850b5994b6cc87c59fe92b9587322.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/bedtools/2.25.0 mugqic/ucsc/v326 && \
mkdir -p tracks/ETOH_shSMC1A_2 tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 97  alignment/ETOH_shSMC1A_2/ETOH_shSMC1A_2.sorted.mdup.reverse.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/ETOH_shSMC1A_2/ETOH_shSMC1A_2.sorted.mdup.reverse.bam \
  -g /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa.fai \
  > tracks/ETOH_shSMC1A_2/ETOH_shSMC1A_2.reverse.bedGraph && \
sort -k1,1 -k2,2n tracks/ETOH_shSMC1A_2/ETOH_shSMC1A_2.reverse.bedGraph > tracks/ETOH_shSMC1A_2/ETOH_shSMC1A_2.reverse.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/ETOH_shSMC1A_2/ETOH_shSMC1A_2.reverse.bedGraph.sorted \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa.fai \
  tracks/bigWig/ETOH_shSMC1A_2.reverse.bw
wiggle.ETOH_shSMC1A_2.reverse.c37850b5994b6cc87c59fe92b9587322.mugqic.done
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
JOB_DONE=job_output/raw_counts/htseq_count.DEX_nonMamm_1.cda09854fa5c0036c14132a63939f746.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'htseq_count.DEX_nonMamm_1.cda09854fa5c0036c14132a63939f746.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/python/2.7.13 && \
mkdir -p raw_counts && \
samtools view -F 4 \
  alignment/DEX_nonMamm_1/DEX_nonMamm_1.QueryNameSorted.bam | \
htseq-count - \
  -m intersection-nonempty \
  --stranded=reverse \
  --format=sam \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.gtf \
  > raw_counts/DEX_nonMamm_1.readcounts.csv
htseq_count.DEX_nonMamm_1.cda09854fa5c0036c14132a63939f746.mugqic.done
)
raw_counts_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$raw_counts_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: raw_counts_2_JOB_ID: htseq_count.DEX_shMED1_1
#-------------------------------------------------------------------------------
JOB_NAME=htseq_count.DEX_shMED1_1
JOB_DEPENDENCIES=$picard_sort_sam_2_JOB_ID
JOB_DONE=job_output/raw_counts/htseq_count.DEX_shMED1_1.20649cac2eaf3c487e3806d360495b6c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'htseq_count.DEX_shMED1_1.20649cac2eaf3c487e3806d360495b6c.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/python/2.7.13 && \
mkdir -p raw_counts && \
samtools view -F 4 \
  alignment/DEX_shMED1_1/DEX_shMED1_1.QueryNameSorted.bam | \
htseq-count - \
  -m intersection-nonempty \
  --stranded=reverse \
  --format=sam \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.gtf \
  > raw_counts/DEX_shMED1_1.readcounts.csv
htseq_count.DEX_shMED1_1.20649cac2eaf3c487e3806d360495b6c.mugqic.done
)
raw_counts_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$raw_counts_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: raw_counts_3_JOB_ID: htseq_count.DEX_shNIPBL_1
#-------------------------------------------------------------------------------
JOB_NAME=htseq_count.DEX_shNIPBL_1
JOB_DEPENDENCIES=$picard_sort_sam_3_JOB_ID
JOB_DONE=job_output/raw_counts/htseq_count.DEX_shNIPBL_1.b4c69e6af7975293fb13169f1f651f5d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'htseq_count.DEX_shNIPBL_1.b4c69e6af7975293fb13169f1f651f5d.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/python/2.7.13 && \
mkdir -p raw_counts && \
samtools view -F 4 \
  alignment/DEX_shNIPBL_1/DEX_shNIPBL_1.QueryNameSorted.bam | \
htseq-count - \
  -m intersection-nonempty \
  --stranded=reverse \
  --format=sam \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.gtf \
  > raw_counts/DEX_shNIPBL_1.readcounts.csv
htseq_count.DEX_shNIPBL_1.b4c69e6af7975293fb13169f1f651f5d.mugqic.done
)
raw_counts_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$raw_counts_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: raw_counts_4_JOB_ID: htseq_count.DEX_shSMC1A_1
#-------------------------------------------------------------------------------
JOB_NAME=htseq_count.DEX_shSMC1A_1
JOB_DEPENDENCIES=$picard_sort_sam_4_JOB_ID
JOB_DONE=job_output/raw_counts/htseq_count.DEX_shSMC1A_1.249a222d0812a0290850fb498e0f5a6c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'htseq_count.DEX_shSMC1A_1.249a222d0812a0290850fb498e0f5a6c.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/python/2.7.13 && \
mkdir -p raw_counts && \
samtools view -F 4 \
  alignment/DEX_shSMC1A_1/DEX_shSMC1A_1.QueryNameSorted.bam | \
htseq-count - \
  -m intersection-nonempty \
  --stranded=reverse \
  --format=sam \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.gtf \
  > raw_counts/DEX_shSMC1A_1.readcounts.csv
htseq_count.DEX_shSMC1A_1.249a222d0812a0290850fb498e0f5a6c.mugqic.done
)
raw_counts_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$raw_counts_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: raw_counts_5_JOB_ID: htseq_count.ETOH_nonMamm_1
#-------------------------------------------------------------------------------
JOB_NAME=htseq_count.ETOH_nonMamm_1
JOB_DEPENDENCIES=$picard_sort_sam_5_JOB_ID
JOB_DONE=job_output/raw_counts/htseq_count.ETOH_nonMamm_1.c002606d033d98f723d9720fcb62c27d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'htseq_count.ETOH_nonMamm_1.c002606d033d98f723d9720fcb62c27d.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/python/2.7.13 && \
mkdir -p raw_counts && \
samtools view -F 4 \
  alignment/ETOH_nonMamm_1/ETOH_nonMamm_1.QueryNameSorted.bam | \
htseq-count - \
  -m intersection-nonempty \
  --stranded=reverse \
  --format=sam \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.gtf \
  > raw_counts/ETOH_nonMamm_1.readcounts.csv
htseq_count.ETOH_nonMamm_1.c002606d033d98f723d9720fcb62c27d.mugqic.done
)
raw_counts_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$raw_counts_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: raw_counts_6_JOB_ID: htseq_count.ETOH_nonMamm_2
#-------------------------------------------------------------------------------
JOB_NAME=htseq_count.ETOH_nonMamm_2
JOB_DEPENDENCIES=$picard_sort_sam_6_JOB_ID
JOB_DONE=job_output/raw_counts/htseq_count.ETOH_nonMamm_2.2f4b039f189ead7f2c504bc5b3c82267.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'htseq_count.ETOH_nonMamm_2.2f4b039f189ead7f2c504bc5b3c82267.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/python/2.7.13 && \
mkdir -p raw_counts && \
samtools view -F 4 \
  alignment/ETOH_nonMamm_2/ETOH_nonMamm_2.QueryNameSorted.bam | \
htseq-count - \
  -m intersection-nonempty \
  --stranded=reverse \
  --format=sam \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.gtf \
  > raw_counts/ETOH_nonMamm_2.readcounts.csv
htseq_count.ETOH_nonMamm_2.2f4b039f189ead7f2c504bc5b3c82267.mugqic.done
)
raw_counts_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$raw_counts_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: raw_counts_7_JOB_ID: htseq_count.ETOH_shMED1_1
#-------------------------------------------------------------------------------
JOB_NAME=htseq_count.ETOH_shMED1_1
JOB_DEPENDENCIES=$picard_sort_sam_7_JOB_ID
JOB_DONE=job_output/raw_counts/htseq_count.ETOH_shMED1_1.1805706205cc6eefcf397f33fc39af06.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'htseq_count.ETOH_shMED1_1.1805706205cc6eefcf397f33fc39af06.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/python/2.7.13 && \
mkdir -p raw_counts && \
samtools view -F 4 \
  alignment/ETOH_shMED1_1/ETOH_shMED1_1.QueryNameSorted.bam | \
htseq-count - \
  -m intersection-nonempty \
  --stranded=reverse \
  --format=sam \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.gtf \
  > raw_counts/ETOH_shMED1_1.readcounts.csv
htseq_count.ETOH_shMED1_1.1805706205cc6eefcf397f33fc39af06.mugqic.done
)
raw_counts_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$raw_counts_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: raw_counts_8_JOB_ID: htseq_count.ETOH_shNIPBL_2
#-------------------------------------------------------------------------------
JOB_NAME=htseq_count.ETOH_shNIPBL_2
JOB_DEPENDENCIES=$picard_sort_sam_8_JOB_ID
JOB_DONE=job_output/raw_counts/htseq_count.ETOH_shNIPBL_2.64d8757ddc896dac291d14f05c36fb4f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'htseq_count.ETOH_shNIPBL_2.64d8757ddc896dac291d14f05c36fb4f.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/python/2.7.13 && \
mkdir -p raw_counts && \
samtools view -F 4 \
  alignment/ETOH_shNIPBL_2/ETOH_shNIPBL_2.QueryNameSorted.bam | \
htseq-count - \
  -m intersection-nonempty \
  --stranded=reverse \
  --format=sam \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.gtf \
  > raw_counts/ETOH_shNIPBL_2.readcounts.csv
htseq_count.ETOH_shNIPBL_2.64d8757ddc896dac291d14f05c36fb4f.mugqic.done
)
raw_counts_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$raw_counts_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: raw_counts_9_JOB_ID: htseq_count.ETOH_shNIPBL_1
#-------------------------------------------------------------------------------
JOB_NAME=htseq_count.ETOH_shNIPBL_1
JOB_DEPENDENCIES=$picard_sort_sam_9_JOB_ID
JOB_DONE=job_output/raw_counts/htseq_count.ETOH_shNIPBL_1.ec5717e8bfb63aab5a8b0a3885573b27.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'htseq_count.ETOH_shNIPBL_1.ec5717e8bfb63aab5a8b0a3885573b27.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/python/2.7.13 && \
mkdir -p raw_counts && \
samtools view -F 4 \
  alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1.QueryNameSorted.bam | \
htseq-count - \
  -m intersection-nonempty \
  --stranded=reverse \
  --format=sam \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.gtf \
  > raw_counts/ETOH_shNIPBL_1.readcounts.csv
htseq_count.ETOH_shNIPBL_1.ec5717e8bfb63aab5a8b0a3885573b27.mugqic.done
)
raw_counts_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$raw_counts_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: raw_counts_10_JOB_ID: htseq_count.ETOH_shMED1_2
#-------------------------------------------------------------------------------
JOB_NAME=htseq_count.ETOH_shMED1_2
JOB_DEPENDENCIES=$picard_sort_sam_10_JOB_ID
JOB_DONE=job_output/raw_counts/htseq_count.ETOH_shMED1_2.728a1efbfc6b094f685e20206001f2f8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'htseq_count.ETOH_shMED1_2.728a1efbfc6b094f685e20206001f2f8.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/python/2.7.13 && \
mkdir -p raw_counts && \
samtools view -F 4 \
  alignment/ETOH_shMED1_2/ETOH_shMED1_2.QueryNameSorted.bam | \
htseq-count - \
  -m intersection-nonempty \
  --stranded=reverse \
  --format=sam \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.gtf \
  > raw_counts/ETOH_shMED1_2.readcounts.csv
htseq_count.ETOH_shMED1_2.728a1efbfc6b094f685e20206001f2f8.mugqic.done
)
raw_counts_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$raw_counts_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: raw_counts_11_JOB_ID: htseq_count.ETOH_shSMC1A_1
#-------------------------------------------------------------------------------
JOB_NAME=htseq_count.ETOH_shSMC1A_1
JOB_DEPENDENCIES=$picard_sort_sam_11_JOB_ID
JOB_DONE=job_output/raw_counts/htseq_count.ETOH_shSMC1A_1.b01c9a24249f2f5489b7685970beb0e3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'htseq_count.ETOH_shSMC1A_1.b01c9a24249f2f5489b7685970beb0e3.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/python/2.7.13 && \
mkdir -p raw_counts && \
samtools view -F 4 \
  alignment/ETOH_shSMC1A_1/ETOH_shSMC1A_1.QueryNameSorted.bam | \
htseq-count - \
  -m intersection-nonempty \
  --stranded=reverse \
  --format=sam \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.gtf \
  > raw_counts/ETOH_shSMC1A_1.readcounts.csv
htseq_count.ETOH_shSMC1A_1.b01c9a24249f2f5489b7685970beb0e3.mugqic.done
)
raw_counts_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$raw_counts_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: raw_counts_12_JOB_ID: htseq_count.ETOH_shSMC1A_2
#-------------------------------------------------------------------------------
JOB_NAME=htseq_count.ETOH_shSMC1A_2
JOB_DEPENDENCIES=$picard_sort_sam_12_JOB_ID
JOB_DONE=job_output/raw_counts/htseq_count.ETOH_shSMC1A_2.bd31be1e8d841acb792d73fc68d648e3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'htseq_count.ETOH_shSMC1A_2.bd31be1e8d841acb792d73fc68d648e3.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/python/2.7.13 && \
mkdir -p raw_counts && \
samtools view -F 4 \
  alignment/ETOH_shSMC1A_2/ETOH_shSMC1A_2.QueryNameSorted.bam | \
htseq-count - \
  -m intersection-nonempty \
  --stranded=reverse \
  --format=sam \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.gtf \
  > raw_counts/ETOH_shSMC1A_2.readcounts.csv
htseq_count.ETOH_shSMC1A_2.bd31be1e8d841acb792d73fc68d648e3.mugqic.done
)
raw_counts_12_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
JOB_DONE=job_output/raw_counts_metrics/metrics.matrix.917566f24290226944d5860151cc320b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'metrics.matrix.917566f24290226944d5860151cc320b.mugqic.done'
module load mugqic/mugqic_tools/2.1.7 && \
mkdir -p DGE && \
gtf2tmpMatrix.awk \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.gtf \
  DGE/tmpMatrix.txt && \
HEAD='Gene\tSymbol' && \
for read_count_file in \
  raw_counts/DEX_nonMamm_1.readcounts.csv \
  raw_counts/DEX_shMED1_1.readcounts.csv \
  raw_counts/DEX_shNIPBL_1.readcounts.csv \
  raw_counts/DEX_shSMC1A_1.readcounts.csv \
  raw_counts/ETOH_nonMamm_1.readcounts.csv \
  raw_counts/ETOH_nonMamm_2.readcounts.csv \
  raw_counts/ETOH_shMED1_1.readcounts.csv \
  raw_counts/ETOH_shNIPBL_2.readcounts.csv \
  raw_counts/ETOH_shNIPBL_1.readcounts.csv \
  raw_counts/ETOH_shMED1_2.readcounts.csv \
  raw_counts/ETOH_shSMC1A_1.readcounts.csv \
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
metrics.matrix.917566f24290226944d5860151cc320b.mugqic.done
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
JOB_DONE=job_output/raw_counts_metrics/rpkm_saturation.c0eb91b36e1bda91fd4b23dd50be7f18.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'rpkm_saturation.c0eb91b36e1bda91fd4b23dd50be7f18.mugqic.done'
module load mugqic/R_Bioconductor/3.2.3_3.2 mugqic/mugqic_tools/2.1.7 && \
mkdir -p metrics/saturation && \
Rscript $R_TOOLS/rpkmSaturation.R \
  DGE/rawCountMatrix.csv \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.genes.length.tsv \
  raw_counts \
  metrics/saturation \
  11 \
  1 && \
zip -r metrics/saturation.zip metrics/saturation
rpkm_saturation.c0eb91b36e1bda91fd4b23dd50be7f18.mugqic.done
)
raw_counts_metrics_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q lm -l nodes=1:ppn=12 -l pmem=5700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
JOB_DONE=job_output/cufflinks/cufflinks.DEX_nonMamm_1.ca8b3bd243d0d7ec1227d05138860091.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cufflinks.DEX_nonMamm_1.ca8b3bd243d0d7ec1227d05138860091.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/DEX_nonMamm_1 && \
cufflinks -q  \
  --GTF-guide /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.gtf \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/DEX_nonMamm_1 \
  --num-threads 12 \
  alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.hardClip.bam
cufflinks.DEX_nonMamm_1.ca8b3bd243d0d7ec1227d05138860091.mugqic.done
)
cufflinks_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cufflinks_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cufflinks_2_JOB_ID: cufflinks.DEX_shMED1_1
#-------------------------------------------------------------------------------
JOB_NAME=cufflinks.DEX_shMED1_1
JOB_DEPENDENCIES=$bam_hard_clip_2_JOB_ID
JOB_DONE=job_output/cufflinks/cufflinks.DEX_shMED1_1.e5ef7b9d7e4a5015f98cdb7225773dbf.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cufflinks.DEX_shMED1_1.e5ef7b9d7e4a5015f98cdb7225773dbf.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/DEX_shMED1_1 && \
cufflinks -q  \
  --GTF-guide /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.gtf \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/DEX_shMED1_1 \
  --num-threads 12 \
  alignment/DEX_shMED1_1/DEX_shMED1_1.sorted.mdup.hardClip.bam
cufflinks.DEX_shMED1_1.e5ef7b9d7e4a5015f98cdb7225773dbf.mugqic.done
)
cufflinks_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cufflinks_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cufflinks_3_JOB_ID: cufflinks.DEX_shNIPBL_1
#-------------------------------------------------------------------------------
JOB_NAME=cufflinks.DEX_shNIPBL_1
JOB_DEPENDENCIES=$bam_hard_clip_3_JOB_ID
JOB_DONE=job_output/cufflinks/cufflinks.DEX_shNIPBL_1.bde4374866671374e3d0b15f5b8bf4db.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cufflinks.DEX_shNIPBL_1.bde4374866671374e3d0b15f5b8bf4db.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/DEX_shNIPBL_1 && \
cufflinks -q  \
  --GTF-guide /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.gtf \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/DEX_shNIPBL_1 \
  --num-threads 12 \
  alignment/DEX_shNIPBL_1/DEX_shNIPBL_1.sorted.mdup.hardClip.bam
cufflinks.DEX_shNIPBL_1.bde4374866671374e3d0b15f5b8bf4db.mugqic.done
)
cufflinks_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cufflinks_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cufflinks_4_JOB_ID: cufflinks.DEX_shSMC1A_1
#-------------------------------------------------------------------------------
JOB_NAME=cufflinks.DEX_shSMC1A_1
JOB_DEPENDENCIES=$bam_hard_clip_4_JOB_ID
JOB_DONE=job_output/cufflinks/cufflinks.DEX_shSMC1A_1.c0971089d83286ac7fa1a153e998509d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cufflinks.DEX_shSMC1A_1.c0971089d83286ac7fa1a153e998509d.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/DEX_shSMC1A_1 && \
cufflinks -q  \
  --GTF-guide /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.gtf \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/DEX_shSMC1A_1 \
  --num-threads 12 \
  alignment/DEX_shSMC1A_1/DEX_shSMC1A_1.sorted.mdup.hardClip.bam
cufflinks.DEX_shSMC1A_1.c0971089d83286ac7fa1a153e998509d.mugqic.done
)
cufflinks_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cufflinks_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cufflinks_5_JOB_ID: cufflinks.ETOH_nonMamm_1
#-------------------------------------------------------------------------------
JOB_NAME=cufflinks.ETOH_nonMamm_1
JOB_DEPENDENCIES=$bam_hard_clip_5_JOB_ID
JOB_DONE=job_output/cufflinks/cufflinks.ETOH_nonMamm_1.d7d7a7f9be7e2a4eba8e7ca5f37614a7.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cufflinks.ETOH_nonMamm_1.d7d7a7f9be7e2a4eba8e7ca5f37614a7.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/ETOH_nonMamm_1 && \
cufflinks -q  \
  --GTF-guide /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.gtf \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/ETOH_nonMamm_1 \
  --num-threads 12 \
  alignment/ETOH_nonMamm_1/ETOH_nonMamm_1.sorted.mdup.hardClip.bam
cufflinks.ETOH_nonMamm_1.d7d7a7f9be7e2a4eba8e7ca5f37614a7.mugqic.done
)
cufflinks_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cufflinks_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cufflinks_6_JOB_ID: cufflinks.ETOH_nonMamm_2
#-------------------------------------------------------------------------------
JOB_NAME=cufflinks.ETOH_nonMamm_2
JOB_DEPENDENCIES=$bam_hard_clip_6_JOB_ID
JOB_DONE=job_output/cufflinks/cufflinks.ETOH_nonMamm_2.e0b9c328f565f4fac39eb2336172ca44.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cufflinks.ETOH_nonMamm_2.e0b9c328f565f4fac39eb2336172ca44.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/ETOH_nonMamm_2 && \
cufflinks -q  \
  --GTF-guide /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.gtf \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/ETOH_nonMamm_2 \
  --num-threads 12 \
  alignment/ETOH_nonMamm_2/ETOH_nonMamm_2.sorted.mdup.hardClip.bam
cufflinks.ETOH_nonMamm_2.e0b9c328f565f4fac39eb2336172ca44.mugqic.done
)
cufflinks_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cufflinks_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cufflinks_7_JOB_ID: cufflinks.ETOH_shMED1_1
#-------------------------------------------------------------------------------
JOB_NAME=cufflinks.ETOH_shMED1_1
JOB_DEPENDENCIES=$bam_hard_clip_7_JOB_ID
JOB_DONE=job_output/cufflinks/cufflinks.ETOH_shMED1_1.327a5c9c118d386856234ca2b4e4b87b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cufflinks.ETOH_shMED1_1.327a5c9c118d386856234ca2b4e4b87b.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/ETOH_shMED1_1 && \
cufflinks -q  \
  --GTF-guide /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.gtf \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/ETOH_shMED1_1 \
  --num-threads 12 \
  alignment/ETOH_shMED1_1/ETOH_shMED1_1.sorted.mdup.hardClip.bam
cufflinks.ETOH_shMED1_1.327a5c9c118d386856234ca2b4e4b87b.mugqic.done
)
cufflinks_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cufflinks_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cufflinks_8_JOB_ID: cufflinks.ETOH_shNIPBL_2
#-------------------------------------------------------------------------------
JOB_NAME=cufflinks.ETOH_shNIPBL_2
JOB_DEPENDENCIES=$bam_hard_clip_8_JOB_ID
JOB_DONE=job_output/cufflinks/cufflinks.ETOH_shNIPBL_2.8e8f785bb88181d9b645c0686bbe48f9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cufflinks.ETOH_shNIPBL_2.8e8f785bb88181d9b645c0686bbe48f9.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/ETOH_shNIPBL_2 && \
cufflinks -q  \
  --GTF-guide /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.gtf \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/ETOH_shNIPBL_2 \
  --num-threads 12 \
  alignment/ETOH_shNIPBL_2/ETOH_shNIPBL_2.sorted.mdup.hardClip.bam
cufflinks.ETOH_shNIPBL_2.8e8f785bb88181d9b645c0686bbe48f9.mugqic.done
)
cufflinks_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cufflinks_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cufflinks_9_JOB_ID: cufflinks.ETOH_shNIPBL_1
#-------------------------------------------------------------------------------
JOB_NAME=cufflinks.ETOH_shNIPBL_1
JOB_DEPENDENCIES=$bam_hard_clip_9_JOB_ID
JOB_DONE=job_output/cufflinks/cufflinks.ETOH_shNIPBL_1.7c6486119a7d1bf771a8ec2b633ddc07.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cufflinks.ETOH_shNIPBL_1.7c6486119a7d1bf771a8ec2b633ddc07.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/ETOH_shNIPBL_1 && \
cufflinks -q  \
  --GTF-guide /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.gtf \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/ETOH_shNIPBL_1 \
  --num-threads 12 \
  alignment/ETOH_shNIPBL_1/ETOH_shNIPBL_1.sorted.mdup.hardClip.bam
cufflinks.ETOH_shNIPBL_1.7c6486119a7d1bf771a8ec2b633ddc07.mugqic.done
)
cufflinks_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cufflinks_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cufflinks_10_JOB_ID: cufflinks.ETOH_shMED1_2
#-------------------------------------------------------------------------------
JOB_NAME=cufflinks.ETOH_shMED1_2
JOB_DEPENDENCIES=$bam_hard_clip_10_JOB_ID
JOB_DONE=job_output/cufflinks/cufflinks.ETOH_shMED1_2.b3d7445fb71a8f34c4a35229a6d9b90d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cufflinks.ETOH_shMED1_2.b3d7445fb71a8f34c4a35229a6d9b90d.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/ETOH_shMED1_2 && \
cufflinks -q  \
  --GTF-guide /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.gtf \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/ETOH_shMED1_2 \
  --num-threads 12 \
  alignment/ETOH_shMED1_2/ETOH_shMED1_2.sorted.mdup.hardClip.bam
cufflinks.ETOH_shMED1_2.b3d7445fb71a8f34c4a35229a6d9b90d.mugqic.done
)
cufflinks_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cufflinks_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cufflinks_11_JOB_ID: cufflinks.ETOH_shSMC1A_1
#-------------------------------------------------------------------------------
JOB_NAME=cufflinks.ETOH_shSMC1A_1
JOB_DEPENDENCIES=$bam_hard_clip_11_JOB_ID
JOB_DONE=job_output/cufflinks/cufflinks.ETOH_shSMC1A_1.49975403a84fbcf007c4502417aef96d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cufflinks.ETOH_shSMC1A_1.49975403a84fbcf007c4502417aef96d.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/ETOH_shSMC1A_1 && \
cufflinks -q  \
  --GTF-guide /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.gtf \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/ETOH_shSMC1A_1 \
  --num-threads 12 \
  alignment/ETOH_shSMC1A_1/ETOH_shSMC1A_1.sorted.mdup.hardClip.bam
cufflinks.ETOH_shSMC1A_1.49975403a84fbcf007c4502417aef96d.mugqic.done
)
cufflinks_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cufflinks_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cufflinks_12_JOB_ID: cufflinks.ETOH_shSMC1A_2
#-------------------------------------------------------------------------------
JOB_NAME=cufflinks.ETOH_shSMC1A_2
JOB_DEPENDENCIES=$bam_hard_clip_12_JOB_ID
JOB_DONE=job_output/cufflinks/cufflinks.ETOH_shSMC1A_2.73b6ce54bffaf4dacaee117ba9352698.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cufflinks.ETOH_shSMC1A_2.73b6ce54bffaf4dacaee117ba9352698.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/ETOH_shSMC1A_2 && \
cufflinks -q  \
  --GTF-guide /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.gtf \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/ETOH_shSMC1A_2 \
  --num-threads 12 \
  alignment/ETOH_shSMC1A_2/ETOH_shSMC1A_2.sorted.mdup.hardClip.bam
cufflinks.ETOH_shSMC1A_2.73b6ce54bffaf4dacaee117ba9352698.mugqic.done
)
cufflinks_12_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
JOB_DONE=job_output/cuffmerge/cuffmerge.f360a3fcd3da46574d0f234151a8a932.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffmerge.f360a3fcd3da46574d0f234151a8a932.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/AllSamples && \
`cat > cufflinks/cuffmerge.samples.txt << END
cufflinks/DEX_nonMamm_1/transcripts.gtf
cufflinks/DEX_shMED1_1/transcripts.gtf
cufflinks/DEX_shNIPBL_1/transcripts.gtf
cufflinks/DEX_shSMC1A_1/transcripts.gtf
cufflinks/ETOH_nonMamm_1/transcripts.gtf
cufflinks/ETOH_nonMamm_2/transcripts.gtf
cufflinks/ETOH_shMED1_1/transcripts.gtf
cufflinks/ETOH_shNIPBL_2/transcripts.gtf
cufflinks/ETOH_shNIPBL_1/transcripts.gtf
cufflinks/ETOH_shMED1_2/transcripts.gtf
cufflinks/ETOH_shSMC1A_1/transcripts.gtf
cufflinks/ETOH_shSMC1A_2/transcripts.gtf
END

` && \
cuffmerge  \
  --ref-gtf /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.gtf \
  --ref-sequence /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa \
  -o cufflinks/AllSamples \
  --num-threads 12 \
  cufflinks/cuffmerge.samples.txt
cuffmerge.f360a3fcd3da46574d0f234151a8a932.mugqic.done
)
cuffmerge_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cuffquant_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cuffquant_4_JOB_ID: cuffquant.DEX_shSMC1A_1
#-------------------------------------------------------------------------------
JOB_NAME=cuffquant.DEX_shSMC1A_1
JOB_DEPENDENCIES=$bam_hard_clip_4_JOB_ID:$cuffmerge_1_JOB_ID
JOB_DONE=job_output/cuffquant/cuffquant.DEX_shSMC1A_1.b2566942ac972d3b85ba1a8082ef9555.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffquant.DEX_shSMC1A_1.b2566942ac972d3b85ba1a8082ef9555.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/DEX_shSMC1A_1 && \
cuffquant -q  \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/DEX_shSMC1A_1 \
  --num-threads 12 \
  cufflinks/AllSamples/merged.gtf \
  alignment/DEX_shSMC1A_1/DEX_shSMC1A_1.sorted.mdup.hardClip.bam
cuffquant.DEX_shSMC1A_1.b2566942ac972d3b85ba1a8082ef9555.mugqic.done
)
cuffquant_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cuffquant_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cuffquant_11_JOB_ID: cuffquant.ETOH_shSMC1A_1
#-------------------------------------------------------------------------------
JOB_NAME=cuffquant.ETOH_shSMC1A_1
JOB_DEPENDENCIES=$bam_hard_clip_11_JOB_ID:$cuffmerge_1_JOB_ID
JOB_DONE=job_output/cuffquant/cuffquant.ETOH_shSMC1A_1.ec4b0d5ff6a15301be9f89834584b983.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffquant.ETOH_shSMC1A_1.ec4b0d5ff6a15301be9f89834584b983.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/ETOH_shSMC1A_1 && \
cuffquant -q  \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/ETOH_shSMC1A_1 \
  --num-threads 12 \
  cufflinks/AllSamples/merged.gtf \
  alignment/ETOH_shSMC1A_1/ETOH_shSMC1A_1.sorted.mdup.hardClip.bam
cuffquant.ETOH_shSMC1A_1.ec4b0d5ff6a15301be9f89834584b983.mugqic.done
)
cuffquant_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
JOB_DONE=job_output/cuffdiff/cuffdiff.Dex_vs_EtOH-All.967a4c401487d24abbc36467d59c7fbd.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffdiff.Dex_vs_EtOH-All.967a4c401487d24abbc36467d59c7fbd.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cuffdiff/Dex_vs_EtOH-All && \
cuffdiff -u \
  --frag-bias-correct /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa \
  --library-type fr-firststrand \
  --output-dir cuffdiff/Dex_vs_EtOH-All \
  --num-threads 12 \
  cufflinks/AllSamples/merged.gtf \
  cufflinks/ETOH_nonMamm_1/abundances.cxb,cufflinks/ETOH_nonMamm_2/abundances.cxb,cufflinks/ETOH_shMED1_1/abundances.cxb,cufflinks/ETOH_shMED1_2/abundances.cxb,cufflinks/ETOH_shNIPBL_1/abundances.cxb,cufflinks/ETOH_shNIPBL_2/abundances.cxb,cufflinks/ETOH_shSMC1A_1/abundances.cxb,cufflinks/ETOH_shSMC1A_2/abundances.cxb \
  cufflinks/DEX_nonMamm_1/abundances.cxb,cufflinks/DEX_shMED1_1/abundances.cxb,cufflinks/DEX_shNIPBL_1/abundances.cxb,cufflinks/DEX_shSMC1A_1/abundances.cxb
cuffdiff.Dex_vs_EtOH-All.967a4c401487d24abbc36467d59c7fbd.mugqic.done
)
cuffdiff_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cuffdiff_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cuffdiff_2_JOB_ID: cuffdiff.Dex_vs_EtOH-nonMamm
#-------------------------------------------------------------------------------
JOB_NAME=cuffdiff.Dex_vs_EtOH-nonMamm
JOB_DEPENDENCIES=$cuffmerge_1_JOB_ID:$cuffquant_1_JOB_ID:$cuffquant_5_JOB_ID:$cuffquant_6_JOB_ID
JOB_DONE=job_output/cuffdiff/cuffdiff.Dex_vs_EtOH-nonMamm.e80031d1d2a9368f1395b83365898f49.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffdiff.Dex_vs_EtOH-nonMamm.e80031d1d2a9368f1395b83365898f49.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cuffdiff/Dex_vs_EtOH-nonMamm && \
cuffdiff -u \
  --frag-bias-correct /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa \
  --library-type fr-firststrand \
  --output-dir cuffdiff/Dex_vs_EtOH-nonMamm \
  --num-threads 12 \
  cufflinks/AllSamples/merged.gtf \
  cufflinks/ETOH_nonMamm_1/abundances.cxb,cufflinks/ETOH_nonMamm_2/abundances.cxb \
  cufflinks/DEX_nonMamm_1/abundances.cxb
cuffdiff.Dex_vs_EtOH-nonMamm.e80031d1d2a9368f1395b83365898f49.mugqic.done
)
cuffdiff_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cuffdiff_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cuffdiff_3_JOB_ID: cuffdiff.Dex_vs_EtOH-shMED1
#-------------------------------------------------------------------------------
JOB_NAME=cuffdiff.Dex_vs_EtOH-shMED1
JOB_DEPENDENCIES=$cuffmerge_1_JOB_ID:$cuffquant_2_JOB_ID:$cuffquant_7_JOB_ID:$cuffquant_10_JOB_ID
JOB_DONE=job_output/cuffdiff/cuffdiff.Dex_vs_EtOH-shMED1.7ff9489bcf9db0028653ed30bb1106e8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffdiff.Dex_vs_EtOH-shMED1.7ff9489bcf9db0028653ed30bb1106e8.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cuffdiff/Dex_vs_EtOH-shMED1 && \
cuffdiff -u \
  --frag-bias-correct /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa \
  --library-type fr-firststrand \
  --output-dir cuffdiff/Dex_vs_EtOH-shMED1 \
  --num-threads 12 \
  cufflinks/AllSamples/merged.gtf \
  cufflinks/ETOH_shMED1_1/abundances.cxb,cufflinks/ETOH_shMED1_2/abundances.cxb \
  cufflinks/DEX_shMED1_1/abundances.cxb
cuffdiff.Dex_vs_EtOH-shMED1.7ff9489bcf9db0028653ed30bb1106e8.mugqic.done
)
cuffdiff_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cuffdiff_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cuffdiff_4_JOB_ID: cuffdiff.Dex_vs_EtOH-shNIPBL
#-------------------------------------------------------------------------------
JOB_NAME=cuffdiff.Dex_vs_EtOH-shNIPBL
JOB_DEPENDENCIES=$cuffmerge_1_JOB_ID:$cuffquant_3_JOB_ID:$cuffquant_8_JOB_ID:$cuffquant_9_JOB_ID
JOB_DONE=job_output/cuffdiff/cuffdiff.Dex_vs_EtOH-shNIPBL.373fa5432e07574aeabdbf0fa0bc71e4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffdiff.Dex_vs_EtOH-shNIPBL.373fa5432e07574aeabdbf0fa0bc71e4.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cuffdiff/Dex_vs_EtOH-shNIPBL && \
cuffdiff -u \
  --frag-bias-correct /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa \
  --library-type fr-firststrand \
  --output-dir cuffdiff/Dex_vs_EtOH-shNIPBL \
  --num-threads 12 \
  cufflinks/AllSamples/merged.gtf \
  cufflinks/ETOH_shNIPBL_1/abundances.cxb,cufflinks/ETOH_shNIPBL_2/abundances.cxb \
  cufflinks/DEX_shNIPBL_1/abundances.cxb
cuffdiff.Dex_vs_EtOH-shNIPBL.373fa5432e07574aeabdbf0fa0bc71e4.mugqic.done
)
cuffdiff_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cuffdiff_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cuffdiff_5_JOB_ID: cuffdiff.Dex_vs_EtOH-shSMC1A
#-------------------------------------------------------------------------------
JOB_NAME=cuffdiff.Dex_vs_EtOH-shSMC1A
JOB_DEPENDENCIES=$cuffmerge_1_JOB_ID:$cuffquant_4_JOB_ID:$cuffquant_11_JOB_ID:$cuffquant_12_JOB_ID
JOB_DONE=job_output/cuffdiff/cuffdiff.Dex_vs_EtOH-shSMC1A.6186837f92401626756a3b63396efc35.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffdiff.Dex_vs_EtOH-shSMC1A.6186837f92401626756a3b63396efc35.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cuffdiff/Dex_vs_EtOH-shSMC1A && \
cuffdiff -u \
  --frag-bias-correct /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa \
  --library-type fr-firststrand \
  --output-dir cuffdiff/Dex_vs_EtOH-shSMC1A \
  --num-threads 12 \
  cufflinks/AllSamples/merged.gtf \
  cufflinks/ETOH_shSMC1A_1/abundances.cxb,cufflinks/ETOH_shSMC1A_2/abundances.cxb \
  cufflinks/DEX_shSMC1A_1/abundances.cxb
cuffdiff.Dex_vs_EtOH-shSMC1A.6186837f92401626756a3b63396efc35.mugqic.done
)
cuffdiff_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cuffdiff_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cuffdiff_6_JOB_ID: cuffdiff.shMED1_vs_shMamm-Dex
#-------------------------------------------------------------------------------
JOB_NAME=cuffdiff.shMED1_vs_shMamm-Dex
JOB_DEPENDENCIES=$cuffmerge_1_JOB_ID:$cuffquant_1_JOB_ID:$cuffquant_2_JOB_ID
JOB_DONE=job_output/cuffdiff/cuffdiff.shMED1_vs_shMamm-Dex.d716efd20d507e8a7039704970cf64a1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffdiff.shMED1_vs_shMamm-Dex.d716efd20d507e8a7039704970cf64a1.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cuffdiff/shMED1_vs_shMamm-Dex && \
cuffdiff -u \
  --frag-bias-correct /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa \
  --library-type fr-firststrand \
  --output-dir cuffdiff/shMED1_vs_shMamm-Dex \
  --num-threads 12 \
  cufflinks/AllSamples/merged.gtf \
  cufflinks/DEX_nonMamm_1/abundances.cxb \
  cufflinks/DEX_shMED1_1/abundances.cxb
cuffdiff.shMED1_vs_shMamm-Dex.d716efd20d507e8a7039704970cf64a1.mugqic.done
)
cuffdiff_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cuffdiff_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cuffdiff_7_JOB_ID: cuffdiff.shMED1_vs_shMamm-EtOH
#-------------------------------------------------------------------------------
JOB_NAME=cuffdiff.shMED1_vs_shMamm-EtOH
JOB_DEPENDENCIES=$cuffmerge_1_JOB_ID:$cuffquant_1_JOB_ID:$cuffquant_3_JOB_ID
JOB_DONE=job_output/cuffdiff/cuffdiff.shMED1_vs_shMamm-EtOH.582c53ac90f805833a2dbf31fafe8213.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffdiff.shMED1_vs_shMamm-EtOH.582c53ac90f805833a2dbf31fafe8213.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cuffdiff/shMED1_vs_shMamm-EtOH && \
cuffdiff -u \
  --frag-bias-correct /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa \
  --library-type fr-firststrand \
  --output-dir cuffdiff/shMED1_vs_shMamm-EtOH \
  --num-threads 12 \
  cufflinks/AllSamples/merged.gtf \
  cufflinks/DEX_nonMamm_1/abundances.cxb \
  cufflinks/DEX_shNIPBL_1/abundances.cxb
cuffdiff.shMED1_vs_shMamm-EtOH.582c53ac90f805833a2dbf31fafe8213.mugqic.done
)
cuffdiff_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cuffdiff_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cuffdiff_8_JOB_ID: cuffdiff.shNIPBL_vs_shMamm-Dex
#-------------------------------------------------------------------------------
JOB_NAME=cuffdiff.shNIPBL_vs_shMamm-Dex
JOB_DEPENDENCIES=$cuffmerge_1_JOB_ID:$cuffquant_1_JOB_ID:$cuffquant_4_JOB_ID
JOB_DONE=job_output/cuffdiff/cuffdiff.shNIPBL_vs_shMamm-Dex.e6c247accaa370eff941e68ad8958530.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffdiff.shNIPBL_vs_shMamm-Dex.e6c247accaa370eff941e68ad8958530.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cuffdiff/shNIPBL_vs_shMamm-Dex && \
cuffdiff -u \
  --frag-bias-correct /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa \
  --library-type fr-firststrand \
  --output-dir cuffdiff/shNIPBL_vs_shMamm-Dex \
  --num-threads 12 \
  cufflinks/AllSamples/merged.gtf \
  cufflinks/DEX_nonMamm_1/abundances.cxb \
  cufflinks/DEX_shSMC1A_1/abundances.cxb
cuffdiff.shNIPBL_vs_shMamm-Dex.e6c247accaa370eff941e68ad8958530.mugqic.done
)
cuffdiff_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cuffdiff_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cuffdiff_9_JOB_ID: cuffdiff.shNIPBL_vs_shMamm-EtOH
#-------------------------------------------------------------------------------
JOB_NAME=cuffdiff.shNIPBL_vs_shMamm-EtOH
JOB_DEPENDENCIES=$cuffmerge_1_JOB_ID:$cuffquant_5_JOB_ID:$cuffquant_6_JOB_ID:$cuffquant_7_JOB_ID:$cuffquant_10_JOB_ID
JOB_DONE=job_output/cuffdiff/cuffdiff.shNIPBL_vs_shMamm-EtOH.4fe3bac07b624cc00e477810a6d0d94d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffdiff.shNIPBL_vs_shMamm-EtOH.4fe3bac07b624cc00e477810a6d0d94d.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cuffdiff/shNIPBL_vs_shMamm-EtOH && \
cuffdiff -u \
  --frag-bias-correct /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa \
  --library-type fr-firststrand \
  --output-dir cuffdiff/shNIPBL_vs_shMamm-EtOH \
  --num-threads 12 \
  cufflinks/AllSamples/merged.gtf \
  cufflinks/ETOH_nonMamm_1/abundances.cxb,cufflinks/ETOH_nonMamm_2/abundances.cxb \
  cufflinks/ETOH_shMED1_1/abundances.cxb,cufflinks/ETOH_shMED1_2/abundances.cxb
cuffdiff.shNIPBL_vs_shMamm-EtOH.4fe3bac07b624cc00e477810a6d0d94d.mugqic.done
)
cuffdiff_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cuffdiff_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cuffdiff_10_JOB_ID: cuffdiff.shSMC1A_vs_shMamm-Dex
#-------------------------------------------------------------------------------
JOB_NAME=cuffdiff.shSMC1A_vs_shMamm-Dex
JOB_DEPENDENCIES=$cuffmerge_1_JOB_ID:$cuffquant_5_JOB_ID:$cuffquant_6_JOB_ID:$cuffquant_8_JOB_ID:$cuffquant_9_JOB_ID
JOB_DONE=job_output/cuffdiff/cuffdiff.shSMC1A_vs_shMamm-Dex.445ae9e298b9be44c70f0029adbbf7f9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffdiff.shSMC1A_vs_shMamm-Dex.445ae9e298b9be44c70f0029adbbf7f9.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cuffdiff/shSMC1A_vs_shMamm-Dex && \
cuffdiff -u \
  --frag-bias-correct /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa \
  --library-type fr-firststrand \
  --output-dir cuffdiff/shSMC1A_vs_shMamm-Dex \
  --num-threads 12 \
  cufflinks/AllSamples/merged.gtf \
  cufflinks/ETOH_nonMamm_1/abundances.cxb,cufflinks/ETOH_nonMamm_2/abundances.cxb \
  cufflinks/ETOH_shNIPBL_1/abundances.cxb,cufflinks/ETOH_shNIPBL_2/abundances.cxb
cuffdiff.shSMC1A_vs_shMamm-Dex.445ae9e298b9be44c70f0029adbbf7f9.mugqic.done
)
cuffdiff_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cuffdiff_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: cuffdiff_11_JOB_ID: cuffdiff.shSMC1A_vs_shMamm-EtOH
#-------------------------------------------------------------------------------
JOB_NAME=cuffdiff.shSMC1A_vs_shMamm-EtOH
JOB_DEPENDENCIES=$cuffmerge_1_JOB_ID:$cuffquant_5_JOB_ID:$cuffquant_6_JOB_ID:$cuffquant_11_JOB_ID:$cuffquant_12_JOB_ID
JOB_DONE=job_output/cuffdiff/cuffdiff.shSMC1A_vs_shMamm-EtOH.ab3749ac0d58e693561633e4592ff885.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffdiff.shSMC1A_vs_shMamm-EtOH.ab3749ac0d58e693561633e4592ff885.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cuffdiff/shSMC1A_vs_shMamm-EtOH && \
cuffdiff -u \
  --frag-bias-correct /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/Homo_sapiens.hg19.fa \
  --library-type fr-firststrand \
  --output-dir cuffdiff/shSMC1A_vs_shMamm-EtOH \
  --num-threads 12 \
  cufflinks/AllSamples/merged.gtf \
  cufflinks/ETOH_nonMamm_1/abundances.cxb,cufflinks/ETOH_nonMamm_2/abundances.cxb \
  cufflinks/ETOH_shSMC1A_1/abundances.cxb,cufflinks/ETOH_shSMC1A_2/abundances.cxb
cuffdiff.shSMC1A_vs_shMamm-EtOH.ab3749ac0d58e693561633e4592ff885.mugqic.done
)
cuffdiff_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$cuffdiff_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


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
JOB_DONE=job_output/cuffnorm/cuffnorm.994bba553af9832278a2ef26b4e6692c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffnorm.994bba553af9832278a2ef26b4e6692c.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cuffnorm && \
cuffnorm -q  \
  --library-type fr-firststrand \
  --output-dir cuffnorm \
  --num-threads 12 \
  --labels DEX_nonMamm_1,DEX_shMED1_1,DEX_shNIPBL_1,DEX_shSMC1A_1,ETOH_nonMamm_1,ETOH_nonMamm_2,ETOH_shMED1_1,ETOH_shNIPBL_2,ETOH_shNIPBL_1,ETOH_shMED1_2,ETOH_shSMC1A_1,ETOH_shSMC1A_2 \
  cufflinks/AllSamples/merged.gtf \
  cufflinks/DEX_nonMamm_1/abundances.cxb \
  cufflinks/DEX_shMED1_1/abundances.cxb \
  cufflinks/DEX_shNIPBL_1/abundances.cxb \
  cufflinks/DEX_shSMC1A_1/abundances.cxb \
  cufflinks/ETOH_nonMamm_1/abundances.cxb \
  cufflinks/ETOH_nonMamm_2/abundances.cxb \
  cufflinks/ETOH_shMED1_1/abundances.cxb \
  cufflinks/ETOH_shNIPBL_2/abundances.cxb \
  cufflinks/ETOH_shNIPBL_1/abundances.cxb \
  cufflinks/ETOH_shMED1_2/abundances.cxb \
  cufflinks/ETOH_shSMC1A_1/abundances.cxb \
  cufflinks/ETOH_shSMC1A_2/abundances.cxb
cuffnorm.994bba553af9832278a2ef26b4e6692c.mugqic.done
)
cuffnorm_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=12 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
JOB_DONE=job_output/gq_seq_utils_exploratory_analysis_rnaseq/gq_seq_utils_exploratory_analysis_rnaseq.1d827a2be63be1d458ce17ba7cd9303e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'gq_seq_utils_exploratory_analysis_rnaseq.1d827a2be63be1d458ce17ba7cd9303e.mugqic.done'
module load mugqic/R_Bioconductor/3.2.3_3.2 mugqic/mugqic_R_packages/1.0.4 && \
mkdir -p exploratory && \
R --no-save --no-restore <<-EOF
suppressPackageStartupMessages(library(gqSeqUtils))

exploratoryAnalysisRNAseq(htseq.counts.path="DGE/rawCountMatrix.csv", cuffnorm.fpkms.dir="cuffnorm", genes.path="/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/annotations/Homo_sapiens.hg19.UCSC2009-03-08.genes.tsv", output.dir="exploratory")
desc = readRDS(file.path("exploratory","index.RData"))
write.table(desc,file=file.path("exploratory","index.tsv"),sep='	',quote=F,col.names=T,row.names=F)
print("done.")

EOF
gq_seq_utils_exploratory_analysis_rnaseq.1d827a2be63be1d458ce17ba7cd9303e.mugqic.done
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=23:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
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
JOB_DONE=job_output/differential_expression/differential_expression.099e216402a405450f3bc26982cb8aee.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'differential_expression.099e216402a405450f3bc26982cb8aee.mugqic.done'
module load mugqic/mugqic_tools/2.1.7 mugqic/R_Bioconductor/3.1.2_3.0 && \
mkdir -p DGE && \
Rscript $R_TOOLS/edger.R \
  -d ../../raw/design.txt \
  -c DGE/rawCountMatrix.csv \
  -o DGE && \
Rscript $R_TOOLS/deseq.R \
  -d ../../raw/design.txt \
  -c DGE/rawCountMatrix.csv \
  -o DGE \
  -l
differential_expression.099e216402a405450f3bc26982cb8aee.mugqic.done
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
wget "http://mugqic.hpc.mcgill.ca/cgi-bin/pipeline.cgi?hostname=lg-1r17-n04&ip=10.241.129.14&pipeline=RnaSeq&steps=picard_sam_to_fastq,trimmomatic,merge_trimmomatic_stats,star,picard_merge_sam_files,picard_sort_sam,picard_mark_duplicates,picard_rna_metrics,estimate_ribosomal_rna,bam_hard_clip,rnaseqc,wiggle,raw_counts,raw_counts_metrics,cufflinks,cuffmerge,cuffquant,cuffdiff,cuffnorm,fpkm_correlation_matrix,gq_seq_utils_exploratory_analysis_rnaseq,differential_expression&samples=12&AnonymizedList=dddd03cb6e65c1183657e5e5bf559f48,444c16ed9974bf0f01c9f5de6751c17b,15b56eae2cdf5203323e847cafca2478,9bb6d10c02e361fb068595ab4db4aed6,b55b936bea7aa8fa3185729524325737,3da4051f49d673a80b1ef99172a9627b,ba0cee9f1a8fbbefae102e42e845ca2b,277b756c8c3e10902132376d11c4df12,ac8379cdfc0119f8bb100d8925e216fd,6aa81e81b0a21275efefd2698d4b38d7,e8f8bd8dfd005d92c279b40cd9c7e216,717ac814ab6a2eb5d20d9f23b67ca76a" --quiet --output-document=/dev/null

