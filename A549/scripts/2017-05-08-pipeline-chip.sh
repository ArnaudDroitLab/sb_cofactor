#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

#-------------------------------------------------------------------------------
# ChipSeq PBSScheduler Job Submission Bash script
# Version: 2.2.0
# Created on: 2017-05-08T19:51:29
# Steps:
#   picard_sam_to_fastq: 6 jobs
#   trimmomatic: 6 jobs
#   merge_trimmomatic_stats: 1 job
#   bwa_mem_picard_sort_sam: 7 jobs
#   samtools_view_filter: 7 jobs
#   picard_merge_sam_files: 6 jobs
#   picard_mark_duplicates: 7 jobs
#   metrics: 2 jobs
#   homer_make_tag_directory: 6 jobs
#   qc_metrics: 1 job
#   homer_make_ucsc_file: 7 jobs
#   macs2_callpeak: 6 jobs
#   homer_annotate_peaks: 6 jobs
#   homer_find_motifs_genome: 6 jobs
#   annotation_graphs: 1 job
#   TOTAL: 75 jobs
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
# JOB: picard_sam_to_fastq_1_JOB_ID: picard_sam_to_fastq.A549_CHIP_DEX_WCE_MF1_rep1_RS
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.A549_CHIP_DEX_WCE_MF1_rep1_RS
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.A549_CHIP_DEX_WCE_MF1_rep1_RS.4974ade2f734eeb1cadfec5d4944a52d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.A549_CHIP_DEX_WCE_MF1_rep1_RS.4974ade2f734eeb1cadfec5d4944a52d.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/SamToFastq.jar \
  VALIDATION_STRINGENCY=LENIENT \
  INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.1613.003.Index_4.A549_Dex_WCE_MF_ChIP1_8_rep1.bam \
  FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.1613.003.Index_4.A549_Dex_WCE_MF_ChIP1_8_rep1.single.fastq.gz
picard_sam_to_fastq.A549_CHIP_DEX_WCE_MF1_rep1_RS.4974ade2f734eeb1cadfec5d4944a52d.mugqic.done
)
picard_sam_to_fastq_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=3 | grep "[0-9]")
echo "$picard_sam_to_fastq_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_2_JOB_ID: picard_sam_to_fastq.A549_CHIP_DEX_NIPBL_MF1_rep1_RS
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.A549_CHIP_DEX_NIPBL_MF1_rep1_RS
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.A549_CHIP_DEX_NIPBL_MF1_rep1_RS.b1acebc16a106634a250edb8382255b0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.A549_CHIP_DEX_NIPBL_MF1_rep1_RS.b1acebc16a106634a250edb8382255b0.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/SamToFastq.jar \
  VALIDATION_STRINGENCY=LENIENT \
  INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.1613.003.Index_5.A549_Dex_NIPBL_MF_ChIP1_1_rep1.bam \
  FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.1613.003.Index_5.A549_Dex_NIPBL_MF_ChIP1_1_rep1.single.fastq.gz
picard_sam_to_fastq.A549_CHIP_DEX_NIPBL_MF1_rep1_RS.b1acebc16a106634a250edb8382255b0.mugqic.done
)
picard_sam_to_fastq_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=3 | grep "[0-9]")
echo "$picard_sam_to_fastq_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_3_JOB_ID: picard_sam_to_fastq.A549_CHIP_DEX_MED1_MF1_rep1_RS
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.A549_CHIP_DEX_MED1_MF1_rep1_RS
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.A549_CHIP_DEX_MED1_MF1_rep1_RS.aeafad905dccabd75c553773f5044a44.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.A549_CHIP_DEX_MED1_MF1_rep1_RS.aeafad905dccabd75c553773f5044a44.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/SamToFastq.jar \
  VALIDATION_STRINGENCY=LENIENT \
  INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.1613.003.Index_12.A549_Dex_MED1_MF_ChIP1_2_rep1.bam \
  FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.1613.003.Index_12.A549_Dex_MED1_MF_ChIP1_2_rep1.single.fastq.gz
picard_sam_to_fastq.A549_CHIP_DEX_MED1_MF1_rep1_RS.aeafad905dccabd75c553773f5044a44.mugqic.done
)
picard_sam_to_fastq_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=3 | grep "[0-9]")
echo "$picard_sam_to_fastq_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_4_JOB_ID: picard_sam_to_fastq.A549_CHIP_DEX_SMC1A_MF1_rep1_RS
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.A549_CHIP_DEX_SMC1A_MF1_rep1_RS
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.A549_CHIP_DEX_SMC1A_MF1_rep1_RS.5670d5edf63e14aaa77e68fde4a5eacc.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.A549_CHIP_DEX_SMC1A_MF1_rep1_RS.5670d5edf63e14aaa77e68fde4a5eacc.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/SamToFastq.jar \
  VALIDATION_STRINGENCY=LENIENT \
  INPUT=/gs/project/fhq-091-aa/Raw_Data/ChIP_2014-06-10/HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam \
  FASTQ=/gs/project/fhq-091-aa/Raw_Data/ChIP_2014-06-10/HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.single.fastq.gz
picard_sam_to_fastq.A549_CHIP_DEX_SMC1A_MF1_rep1_RS.5670d5edf63e14aaa77e68fde4a5eacc.mugqic.done
)
picard_sam_to_fastq_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=3 | grep "[0-9]")
echo "$picard_sam_to_fastq_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_5_JOB_ID: picard_sam_to_fastq.A549_CHIP_DEX_CDK9_MF1_rep1_RS
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.A549_CHIP_DEX_CDK9_MF1_rep1_RS
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.A549_CHIP_DEX_CDK9_MF1_rep1_RS.e3dda43eb2d16190c596aa563f9e298e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.A549_CHIP_DEX_CDK9_MF1_rep1_RS.e3dda43eb2d16190c596aa563f9e298e.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/SamToFastq.jar \
  VALIDATION_STRINGENCY=LENIENT \
  INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.1613.003.Index_1.A549_Dex_CDK9_MF_ChIP1_4_rep1.bam \
  FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.1613.003.Index_1.A549_Dex_CDK9_MF_ChIP1_4_rep1.single.fastq.gz
picard_sam_to_fastq.A549_CHIP_DEX_CDK9_MF1_rep1_RS.e3dda43eb2d16190c596aa563f9e298e.mugqic.done
)
picard_sam_to_fastq_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=3 | grep "[0-9]")
echo "$picard_sam_to_fastq_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_sam_to_fastq_6_JOB_ID: picard_sam_to_fastq.A549_CHIP_DEX_BRD4_MF1_rep1_RS
#-------------------------------------------------------------------------------
JOB_NAME=picard_sam_to_fastq.A549_CHIP_DEX_BRD4_MF1_rep1_RS
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_sam_to_fastq/picard_sam_to_fastq.A549_CHIP_DEX_BRD4_MF1_rep1_RS.3c948273b033bfdb7fc008e8b6c4371f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sam_to_fastq.A549_CHIP_DEX_BRD4_MF1_rep1_RS.3c948273b033bfdb7fc008e8b6c4371f.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx10G -jar $PICARD_HOME/SamToFastq.jar \
  VALIDATION_STRINGENCY=LENIENT \
  INPUT=/gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.1613.003.Index_19.A549_Dex_SMC1_MF_ChIP1_3_rep1.bam \
  FASTQ=/gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.1613.003.Index_19.A549_Dex_SMC1_MF_ChIP1_3_rep1.single.fastq.gz
picard_sam_to_fastq.A549_CHIP_DEX_BRD4_MF1_rep1_RS.3c948273b033bfdb7fc008e8b6c4371f.mugqic.done
)
picard_sam_to_fastq_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=3 | grep "[0-9]")
echo "$picard_sam_to_fastq_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: trimmomatic
#-------------------------------------------------------------------------------
STEP=trimmomatic
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: trimmomatic_1_JOB_ID: trimmomatic.A549_CHIP_DEX_WCE_MF1_rep1_RS
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.A549_CHIP_DEX_WCE_MF1_rep1_RS
JOB_DEPENDENCIES=$picard_sam_to_fastq_1_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.A549_CHIP_DEX_WCE_MF1_rep1_RS.795829b80ec9119a37e520f6185fae79.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.A549_CHIP_DEX_WCE_MF1_rep1_RS.795829b80ec9119a37e520f6185fae79.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/A549_CHIP_DEX_WCE_MF1_rep1 && \
`cat > trim/A549_CHIP_DEX_WCE_MF1_rep1/A549_CHIP_DEX_WCE_MF1_rep1_RS.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.1613.003.Index_4.A549_Dex_WCE_MF_ChIP1_8_rep1.single.fastq.gz \
  trim/A549_CHIP_DEX_WCE_MF1_rep1/A549_CHIP_DEX_WCE_MF1_rep1_RS.trim.single.fastq.gz \
  ILLUMINACLIP:trim/A549_CHIP_DEX_WCE_MF1_rep1/A549_CHIP_DEX_WCE_MF1_rep1_RS.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/A549_CHIP_DEX_WCE_MF1_rep1/A549_CHIP_DEX_WCE_MF1_rep1_RS.trim.log
trimmomatic.A549_CHIP_DEX_WCE_MF1_rep1_RS.795829b80ec9119a37e520f6185fae79.mugqic.done
)
trimmomatic_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_2_JOB_ID: trimmomatic.A549_CHIP_DEX_NIPBL_MF1_rep1_RS
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.A549_CHIP_DEX_NIPBL_MF1_rep1_RS
JOB_DEPENDENCIES=$picard_sam_to_fastq_2_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.A549_CHIP_DEX_NIPBL_MF1_rep1_RS.6d05fb41dd30aec50660ffd2364bc864.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.A549_CHIP_DEX_NIPBL_MF1_rep1_RS.6d05fb41dd30aec50660ffd2364bc864.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/A549_CHIP_DEX_NIPBL_MF1_rep1 && \
`cat > trim/A549_CHIP_DEX_NIPBL_MF1_rep1/A549_CHIP_DEX_NIPBL_MF1_rep1_RS.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.1613.003.Index_5.A549_Dex_NIPBL_MF_ChIP1_1_rep1.single.fastq.gz \
  trim/A549_CHIP_DEX_NIPBL_MF1_rep1/A549_CHIP_DEX_NIPBL_MF1_rep1_RS.trim.single.fastq.gz \
  ILLUMINACLIP:trim/A549_CHIP_DEX_NIPBL_MF1_rep1/A549_CHIP_DEX_NIPBL_MF1_rep1_RS.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/A549_CHIP_DEX_NIPBL_MF1_rep1/A549_CHIP_DEX_NIPBL_MF1_rep1_RS.trim.log
trimmomatic.A549_CHIP_DEX_NIPBL_MF1_rep1_RS.6d05fb41dd30aec50660ffd2364bc864.mugqic.done
)
trimmomatic_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_3_JOB_ID: trimmomatic.A549_CHIP_DEX_MED1_MF1_rep1_RS
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.A549_CHIP_DEX_MED1_MF1_rep1_RS
JOB_DEPENDENCIES=$picard_sam_to_fastq_3_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.A549_CHIP_DEX_MED1_MF1_rep1_RS.a8496602b06e40a5f643fd068215ba27.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.A549_CHIP_DEX_MED1_MF1_rep1_RS.a8496602b06e40a5f643fd068215ba27.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/A549_CHIP_DEX_MED1_MF1_rep1 && \
`cat > trim/A549_CHIP_DEX_MED1_MF1_rep1/A549_CHIP_DEX_MED1_MF1_rep1_RS.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.1613.003.Index_12.A549_Dex_MED1_MF_ChIP1_2_rep1.single.fastq.gz \
  trim/A549_CHIP_DEX_MED1_MF1_rep1/A549_CHIP_DEX_MED1_MF1_rep1_RS.trim.single.fastq.gz \
  ILLUMINACLIP:trim/A549_CHIP_DEX_MED1_MF1_rep1/A549_CHIP_DEX_MED1_MF1_rep1_RS.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/A549_CHIP_DEX_MED1_MF1_rep1/A549_CHIP_DEX_MED1_MF1_rep1_RS.trim.log
trimmomatic.A549_CHIP_DEX_MED1_MF1_rep1_RS.a8496602b06e40a5f643fd068215ba27.mugqic.done
)
trimmomatic_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_4_JOB_ID: trimmomatic.A549_CHIP_DEX_SMC1A_MF1_rep1_RS
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.A549_CHIP_DEX_SMC1A_MF1_rep1_RS
JOB_DEPENDENCIES=$picard_sam_to_fastq_4_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.A549_CHIP_DEX_SMC1A_MF1_rep1_RS.017661099c032ac8e73b0fa1cbc46d5a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.A549_CHIP_DEX_SMC1A_MF1_rep1_RS.017661099c032ac8e73b0fa1cbc46d5a.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/A549_CHIP_DEX_SMC1A_MF1_rep1 && \
`cat > trim/A549_CHIP_DEX_SMC1A_MF1_rep1/A549_CHIP_DEX_SMC1A_MF1_rep1_RS.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/project/fhq-091-aa/Raw_Data/ChIP_2014-06-10/HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.single.fastq.gz \
  trim/A549_CHIP_DEX_SMC1A_MF1_rep1/A549_CHIP_DEX_SMC1A_MF1_rep1_RS.trim.single.fastq.gz \
  ILLUMINACLIP:trim/A549_CHIP_DEX_SMC1A_MF1_rep1/A549_CHIP_DEX_SMC1A_MF1_rep1_RS.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/A549_CHIP_DEX_SMC1A_MF1_rep1/A549_CHIP_DEX_SMC1A_MF1_rep1_RS.trim.log
trimmomatic.A549_CHIP_DEX_SMC1A_MF1_rep1_RS.017661099c032ac8e73b0fa1cbc46d5a.mugqic.done
)
trimmomatic_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_5_JOB_ID: trimmomatic.A549_CHIP_DEX_CDK9_MF1_rep1_RS
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.A549_CHIP_DEX_CDK9_MF1_rep1_RS
JOB_DEPENDENCIES=$picard_sam_to_fastq_5_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.A549_CHIP_DEX_CDK9_MF1_rep1_RS.45fea0e0ec0c0b3bb58dc1ae3b77750f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.A549_CHIP_DEX_CDK9_MF1_rep1_RS.45fea0e0ec0c0b3bb58dc1ae3b77750f.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/A549_CHIP_DEX_CDK9_MF1_rep1 && \
`cat > trim/A549_CHIP_DEX_CDK9_MF1_rep1/A549_CHIP_DEX_CDK9_MF1_rep1_RS.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.1613.003.Index_1.A549_Dex_CDK9_MF_ChIP1_4_rep1.single.fastq.gz \
  trim/A549_CHIP_DEX_CDK9_MF1_rep1/A549_CHIP_DEX_CDK9_MF1_rep1_RS.trim.single.fastq.gz \
  ILLUMINACLIP:trim/A549_CHIP_DEX_CDK9_MF1_rep1/A549_CHIP_DEX_CDK9_MF1_rep1_RS.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/A549_CHIP_DEX_CDK9_MF1_rep1/A549_CHIP_DEX_CDK9_MF1_rep1_RS.trim.log
trimmomatic.A549_CHIP_DEX_CDK9_MF1_rep1_RS.45fea0e0ec0c0b3bb58dc1ae3b77750f.mugqic.done
)
trimmomatic_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_6_JOB_ID: trimmomatic.A549_CHIP_DEX_BRD4_MF1_rep1_RS
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.A549_CHIP_DEX_BRD4_MF1_rep1_RS
JOB_DEPENDENCIES=$picard_sam_to_fastq_6_JOB_ID
JOB_DONE=job_output/trimmomatic/trimmomatic.A549_CHIP_DEX_BRD4_MF1_rep1_RS.83a1a7ca9c924c5fd1bfddbfac0f3b00.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.A549_CHIP_DEX_BRD4_MF1_rep1_RS.83a1a7ca9c924c5fd1bfddbfac0f3b00.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/A549_CHIP_DEX_BRD4_MF1_rep1 && \
`cat > trim/A549_CHIP_DEX_BRD4_MF1_rep1/A549_CHIP_DEX_BRD4_MF1_rep1_RS.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.1613.003.Index_19.A549_Dex_SMC1_MF_ChIP1_3_rep1.single.fastq.gz \
  trim/A549_CHIP_DEX_BRD4_MF1_rep1/A549_CHIP_DEX_BRD4_MF1_rep1_RS.trim.single.fastq.gz \
  ILLUMINACLIP:trim/A549_CHIP_DEX_BRD4_MF1_rep1/A549_CHIP_DEX_BRD4_MF1_rep1_RS.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/A549_CHIP_DEX_BRD4_MF1_rep1/A549_CHIP_DEX_BRD4_MF1_rep1_RS.trim.log
trimmomatic.A549_CHIP_DEX_BRD4_MF1_rep1_RS.83a1a7ca9c924c5fd1bfddbfac0f3b00.mugqic.done
)
trimmomatic_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$trimmomatic_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: merge_trimmomatic_stats
#-------------------------------------------------------------------------------
STEP=merge_trimmomatic_stats
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: merge_trimmomatic_stats_1_JOB_ID: merge_trimmomatic_stats
#-------------------------------------------------------------------------------
JOB_NAME=merge_trimmomatic_stats
JOB_DEPENDENCIES=$trimmomatic_1_JOB_ID:$trimmomatic_2_JOB_ID:$trimmomatic_3_JOB_ID:$trimmomatic_4_JOB_ID:$trimmomatic_5_JOB_ID:$trimmomatic_6_JOB_ID
JOB_DONE=job_output/merge_trimmomatic_stats/merge_trimmomatic_stats.df05b927a25262779b4eedf9327b7a91.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'merge_trimmomatic_stats.df05b927a25262779b4eedf9327b7a91.mugqic.done'
module load mugqic/pandoc/1.15.2 && \
mkdir -p metrics && \
echo 'Sample	Readset	Raw Single Reads #	Surviving Single Reads #	Surviving Single Reads %' > metrics/trimReadsetTable.tsv && \
grep ^Input trim/A549_CHIP_DEX_WCE_MF1_rep1/A549_CHIP_DEX_WCE_MF1_rep1_RS.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/A549_CHIP_DEX_WCE_MF1_rep1	A549_CHIP_DEX_WCE_MF1_rep1_RS	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/A549_CHIP_DEX_NIPBL_MF1_rep1/A549_CHIP_DEX_NIPBL_MF1_rep1_RS.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/A549_CHIP_DEX_NIPBL_MF1_rep1	A549_CHIP_DEX_NIPBL_MF1_rep1_RS	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/A549_CHIP_DEX_MED1_MF1_rep1/A549_CHIP_DEX_MED1_MF1_rep1_RS.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/A549_CHIP_DEX_MED1_MF1_rep1	A549_CHIP_DEX_MED1_MF1_rep1_RS	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/A549_CHIP_DEX_SMC1A_MF1_rep1/A549_CHIP_DEX_SMC1A_MF1_rep1_RS.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/A549_CHIP_DEX_SMC1A_MF1_rep1	A549_CHIP_DEX_SMC1A_MF1_rep1_RS	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/A549_CHIP_DEX_CDK9_MF1_rep1/A549_CHIP_DEX_CDK9_MF1_rep1_RS.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/A549_CHIP_DEX_CDK9_MF1_rep1	A549_CHIP_DEX_CDK9_MF1_rep1_RS	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/A549_CHIP_DEX_BRD4_MF1_rep1/A549_CHIP_DEX_BRD4_MF1_rep1_RS.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/A549_CHIP_DEX_BRD4_MF1_rep1	A549_CHIP_DEX_BRD4_MF1_rep1_RS	\1	\2/' | \
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
merge_trimmomatic_stats.df05b927a25262779b4eedf9327b7a91.mugqic.done
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
# JOB: bwa_mem_picard_sort_sam_1_JOB_ID: bwa_mem_picard_sort_sam.A549_CHIP_DEX_WCE_MF1_rep1_RS
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.A549_CHIP_DEX_WCE_MF1_rep1_RS
JOB_DEPENDENCIES=$trimmomatic_1_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.A549_CHIP_DEX_WCE_MF1_rep1_RS.1add0ae8a6aca5a230b02cea753ff11c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.A549_CHIP_DEX_WCE_MF1_rep1_RS.1add0ae8a6aca5a230b02cea753ff11c.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/A549_CHIP_DEX_WCE_MF1_rep1/A549_CHIP_DEX_WCE_MF1_rep1_RS && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:A549_CHIP_DEX_WCE_MF1_rep1_RS	SM:A549_CHIP_DEX_WCE_MF1_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Mus_musculus.mm10/genome/bwa_index/Mus_musculus.mm10.fa \
  trim/A549_CHIP_DEX_WCE_MF1_rep1/A549_CHIP_DEX_WCE_MF1_rep1_RS.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/A549_CHIP_DEX_WCE_MF1_rep1/A549_CHIP_DEX_WCE_MF1_rep1_RS/A549_CHIP_DEX_WCE_MF1_rep1_RS.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.A549_CHIP_DEX_WCE_MF1_rep1_RS.1add0ae8a6aca5a230b02cea753ff11c.mugqic.done
)
bwa_mem_picard_sort_sam_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_2_JOB_ID: bwa_mem_picard_sort_sam.A549_CHIP_DEX_NIPBL_MF1_rep1_RS
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.A549_CHIP_DEX_NIPBL_MF1_rep1_RS
JOB_DEPENDENCIES=$trimmomatic_2_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.A549_CHIP_DEX_NIPBL_MF1_rep1_RS.78b0d4e21f1ff05433296c0aac9e29c1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.A549_CHIP_DEX_NIPBL_MF1_rep1_RS.78b0d4e21f1ff05433296c0aac9e29c1.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/A549_CHIP_DEX_NIPBL_MF1_rep1/A549_CHIP_DEX_NIPBL_MF1_rep1_RS && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:A549_CHIP_DEX_NIPBL_MF1_rep1_RS	SM:A549_CHIP_DEX_NIPBL_MF1_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Mus_musculus.mm10/genome/bwa_index/Mus_musculus.mm10.fa \
  trim/A549_CHIP_DEX_NIPBL_MF1_rep1/A549_CHIP_DEX_NIPBL_MF1_rep1_RS.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/A549_CHIP_DEX_NIPBL_MF1_rep1/A549_CHIP_DEX_NIPBL_MF1_rep1_RS/A549_CHIP_DEX_NIPBL_MF1_rep1_RS.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.A549_CHIP_DEX_NIPBL_MF1_rep1_RS.78b0d4e21f1ff05433296c0aac9e29c1.mugqic.done
)
bwa_mem_picard_sort_sam_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_3_JOB_ID: bwa_mem_picard_sort_sam.A549_CHIP_DEX_MED1_MF1_rep1_RS
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.A549_CHIP_DEX_MED1_MF1_rep1_RS
JOB_DEPENDENCIES=$trimmomatic_3_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.A549_CHIP_DEX_MED1_MF1_rep1_RS.63045c0f8c99407b0040234321df06c1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.A549_CHIP_DEX_MED1_MF1_rep1_RS.63045c0f8c99407b0040234321df06c1.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/A549_CHIP_DEX_MED1_MF1_rep1/A549_CHIP_DEX_MED1_MF1_rep1_RS && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:A549_CHIP_DEX_MED1_MF1_rep1_RS	SM:A549_CHIP_DEX_MED1_MF1_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Mus_musculus.mm10/genome/bwa_index/Mus_musculus.mm10.fa \
  trim/A549_CHIP_DEX_MED1_MF1_rep1/A549_CHIP_DEX_MED1_MF1_rep1_RS.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/A549_CHIP_DEX_MED1_MF1_rep1/A549_CHIP_DEX_MED1_MF1_rep1_RS/A549_CHIP_DEX_MED1_MF1_rep1_RS.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.A549_CHIP_DEX_MED1_MF1_rep1_RS.63045c0f8c99407b0040234321df06c1.mugqic.done
)
bwa_mem_picard_sort_sam_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_4_JOB_ID: bwa_mem_picard_sort_sam.A549_CHIP_DEX_SMC1A_MF1_rep1_RS
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.A549_CHIP_DEX_SMC1A_MF1_rep1_RS
JOB_DEPENDENCIES=$trimmomatic_4_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.A549_CHIP_DEX_SMC1A_MF1_rep1_RS.11d911cb5eb1d97188c436c0a3ebbc84.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.A549_CHIP_DEX_SMC1A_MF1_rep1_RS.11d911cb5eb1d97188c436c0a3ebbc84.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/A549_CHIP_DEX_SMC1A_MF1_rep1/A549_CHIP_DEX_SMC1A_MF1_rep1_RS && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:A549_CHIP_DEX_SMC1A_MF1_rep1_RS	SM:A549_CHIP_DEX_SMC1A_MF1_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Mus_musculus.mm10/genome/bwa_index/Mus_musculus.mm10.fa \
  trim/A549_CHIP_DEX_SMC1A_MF1_rep1/A549_CHIP_DEX_SMC1A_MF1_rep1_RS.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/A549_CHIP_DEX_SMC1A_MF1_rep1/A549_CHIP_DEX_SMC1A_MF1_rep1_RS/A549_CHIP_DEX_SMC1A_MF1_rep1_RS.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.A549_CHIP_DEX_SMC1A_MF1_rep1_RS.11d911cb5eb1d97188c436c0a3ebbc84.mugqic.done
)
bwa_mem_picard_sort_sam_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_5_JOB_ID: bwa_mem_picard_sort_sam.A549_CHIP_DEX_CDK9_MF1_rep1_RS
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.A549_CHIP_DEX_CDK9_MF1_rep1_RS
JOB_DEPENDENCIES=$trimmomatic_5_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.A549_CHIP_DEX_CDK9_MF1_rep1_RS.68b1c3de94a345bd24d5db277ee007e2.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.A549_CHIP_DEX_CDK9_MF1_rep1_RS.68b1c3de94a345bd24d5db277ee007e2.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/A549_CHIP_DEX_CDK9_MF1_rep1/A549_CHIP_DEX_CDK9_MF1_rep1_RS && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:A549_CHIP_DEX_CDK9_MF1_rep1_RS	SM:A549_CHIP_DEX_CDK9_MF1_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Mus_musculus.mm10/genome/bwa_index/Mus_musculus.mm10.fa \
  trim/A549_CHIP_DEX_CDK9_MF1_rep1/A549_CHIP_DEX_CDK9_MF1_rep1_RS.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/A549_CHIP_DEX_CDK9_MF1_rep1/A549_CHIP_DEX_CDK9_MF1_rep1_RS/A549_CHIP_DEX_CDK9_MF1_rep1_RS.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.A549_CHIP_DEX_CDK9_MF1_rep1_RS.68b1c3de94a345bd24d5db277ee007e2.mugqic.done
)
bwa_mem_picard_sort_sam_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_6_JOB_ID: bwa_mem_picard_sort_sam.A549_CHIP_DEX_BRD4_MF1_rep1_RS
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.A549_CHIP_DEX_BRD4_MF1_rep1_RS
JOB_DEPENDENCIES=$trimmomatic_6_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.A549_CHIP_DEX_BRD4_MF1_rep1_RS.f7aaffe809db66c6e50b28c757d329f7.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.A549_CHIP_DEX_BRD4_MF1_rep1_RS.f7aaffe809db66c6e50b28c757d329f7.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/A549_CHIP_DEX_BRD4_MF1_rep1/A549_CHIP_DEX_BRD4_MF1_rep1_RS && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:A549_CHIP_DEX_BRD4_MF1_rep1_RS	SM:A549_CHIP_DEX_BRD4_MF1_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Mus_musculus.mm10/genome/bwa_index/Mus_musculus.mm10.fa \
  trim/A549_CHIP_DEX_BRD4_MF1_rep1/A549_CHIP_DEX_BRD4_MF1_rep1_RS.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/A549_CHIP_DEX_BRD4_MF1_rep1/A549_CHIP_DEX_BRD4_MF1_rep1_RS/A549_CHIP_DEX_BRD4_MF1_rep1_RS.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.A549_CHIP_DEX_BRD4_MF1_rep1_RS.f7aaffe809db66c6e50b28c757d329f7.mugqic.done
)
bwa_mem_picard_sort_sam_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_7_JOB_ID: bwa_mem_picard_sort_sam_report
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam_report
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_1_JOB_ID:$bwa_mem_picard_sort_sam_2_JOB_ID:$bwa_mem_picard_sort_sam_3_JOB_ID:$bwa_mem_picard_sort_sam_4_JOB_ID:$bwa_mem_picard_sort_sam_5_JOB_ID:$bwa_mem_picard_sort_sam_6_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam_report.1353c26fc9a055a24e53e09763a8940d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam_report.1353c26fc9a055a24e53e09763a8940d.mugqic.done'
module load mugqic/pandoc/1.15.2 && \
mkdir -p report && \
pandoc --to=markdown \
  --template /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/DnaSeq.bwa_mem_picard_sort_sam.md \
  --variable scientific_name="Mus_musculus" \
  --variable assembly="mm10" \
  /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/DnaSeq.bwa_mem_picard_sort_sam.md \
  > report/DnaSeq.bwa_mem_picard_sort_sam.md
bwa_mem_picard_sort_sam_report.1353c26fc9a055a24e53e09763a8940d.mugqic.done
)
bwa_mem_picard_sort_sam_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: samtools_view_filter
#-------------------------------------------------------------------------------
STEP=samtools_view_filter
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_1_JOB_ID: samtools_view_filter.A549_CHIP_DEX_WCE_MF1_rep1_RS
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.A549_CHIP_DEX_WCE_MF1_rep1_RS
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_1_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.A549_CHIP_DEX_WCE_MF1_rep1_RS.19e30d66c882dbfe5919ed10234bfeeb.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.A549_CHIP_DEX_WCE_MF1_rep1_RS.19e30d66c882dbfe5919ed10234bfeeb.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/A549_CHIP_DEX_WCE_MF1_rep1/A549_CHIP_DEX_WCE_MF1_rep1_RS/A549_CHIP_DEX_WCE_MF1_rep1_RS.sorted.bam \
  > alignment/A549_CHIP_DEX_WCE_MF1_rep1/A549_CHIP_DEX_WCE_MF1_rep1_RS/A549_CHIP_DEX_WCE_MF1_rep1_RS.sorted.filtered.bam
samtools_view_filter.A549_CHIP_DEX_WCE_MF1_rep1_RS.19e30d66c882dbfe5919ed10234bfeeb.mugqic.done
)
samtools_view_filter_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_2_JOB_ID: samtools_view_filter.A549_CHIP_DEX_NIPBL_MF1_rep1_RS
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.A549_CHIP_DEX_NIPBL_MF1_rep1_RS
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_2_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.A549_CHIP_DEX_NIPBL_MF1_rep1_RS.b9fd67b7a207c50b11dcd812507861e9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.A549_CHIP_DEX_NIPBL_MF1_rep1_RS.b9fd67b7a207c50b11dcd812507861e9.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/A549_CHIP_DEX_NIPBL_MF1_rep1/A549_CHIP_DEX_NIPBL_MF1_rep1_RS/A549_CHIP_DEX_NIPBL_MF1_rep1_RS.sorted.bam \
  > alignment/A549_CHIP_DEX_NIPBL_MF1_rep1/A549_CHIP_DEX_NIPBL_MF1_rep1_RS/A549_CHIP_DEX_NIPBL_MF1_rep1_RS.sorted.filtered.bam
samtools_view_filter.A549_CHIP_DEX_NIPBL_MF1_rep1_RS.b9fd67b7a207c50b11dcd812507861e9.mugqic.done
)
samtools_view_filter_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_3_JOB_ID: samtools_view_filter.A549_CHIP_DEX_MED1_MF1_rep1_RS
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.A549_CHIP_DEX_MED1_MF1_rep1_RS
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_3_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.A549_CHIP_DEX_MED1_MF1_rep1_RS.b75f7742807e1c04ee2b2eef895fd3ed.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.A549_CHIP_DEX_MED1_MF1_rep1_RS.b75f7742807e1c04ee2b2eef895fd3ed.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/A549_CHIP_DEX_MED1_MF1_rep1/A549_CHIP_DEX_MED1_MF1_rep1_RS/A549_CHIP_DEX_MED1_MF1_rep1_RS.sorted.bam \
  > alignment/A549_CHIP_DEX_MED1_MF1_rep1/A549_CHIP_DEX_MED1_MF1_rep1_RS/A549_CHIP_DEX_MED1_MF1_rep1_RS.sorted.filtered.bam
samtools_view_filter.A549_CHIP_DEX_MED1_MF1_rep1_RS.b75f7742807e1c04ee2b2eef895fd3ed.mugqic.done
)
samtools_view_filter_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_4_JOB_ID: samtools_view_filter.A549_CHIP_DEX_SMC1A_MF1_rep1_RS
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.A549_CHIP_DEX_SMC1A_MF1_rep1_RS
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_4_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.A549_CHIP_DEX_SMC1A_MF1_rep1_RS.0c5ed11f7e69c1bc67e0026c84eca57f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.A549_CHIP_DEX_SMC1A_MF1_rep1_RS.0c5ed11f7e69c1bc67e0026c84eca57f.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/A549_CHIP_DEX_SMC1A_MF1_rep1/A549_CHIP_DEX_SMC1A_MF1_rep1_RS/A549_CHIP_DEX_SMC1A_MF1_rep1_RS.sorted.bam \
  > alignment/A549_CHIP_DEX_SMC1A_MF1_rep1/A549_CHIP_DEX_SMC1A_MF1_rep1_RS/A549_CHIP_DEX_SMC1A_MF1_rep1_RS.sorted.filtered.bam
samtools_view_filter.A549_CHIP_DEX_SMC1A_MF1_rep1_RS.0c5ed11f7e69c1bc67e0026c84eca57f.mugqic.done
)
samtools_view_filter_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_5_JOB_ID: samtools_view_filter.A549_CHIP_DEX_CDK9_MF1_rep1_RS
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.A549_CHIP_DEX_CDK9_MF1_rep1_RS
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_5_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.A549_CHIP_DEX_CDK9_MF1_rep1_RS.699e51d92a0d095136b335f5532e25b3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.A549_CHIP_DEX_CDK9_MF1_rep1_RS.699e51d92a0d095136b335f5532e25b3.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/A549_CHIP_DEX_CDK9_MF1_rep1/A549_CHIP_DEX_CDK9_MF1_rep1_RS/A549_CHIP_DEX_CDK9_MF1_rep1_RS.sorted.bam \
  > alignment/A549_CHIP_DEX_CDK9_MF1_rep1/A549_CHIP_DEX_CDK9_MF1_rep1_RS/A549_CHIP_DEX_CDK9_MF1_rep1_RS.sorted.filtered.bam
samtools_view_filter.A549_CHIP_DEX_CDK9_MF1_rep1_RS.699e51d92a0d095136b335f5532e25b3.mugqic.done
)
samtools_view_filter_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_6_JOB_ID: samtools_view_filter.A549_CHIP_DEX_BRD4_MF1_rep1_RS
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.A549_CHIP_DEX_BRD4_MF1_rep1_RS
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_6_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.A549_CHIP_DEX_BRD4_MF1_rep1_RS.a498a7628fa7a4993f5ca6de2988e020.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.A549_CHIP_DEX_BRD4_MF1_rep1_RS.a498a7628fa7a4993f5ca6de2988e020.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/A549_CHIP_DEX_BRD4_MF1_rep1/A549_CHIP_DEX_BRD4_MF1_rep1_RS/A549_CHIP_DEX_BRD4_MF1_rep1_RS.sorted.bam \
  > alignment/A549_CHIP_DEX_BRD4_MF1_rep1/A549_CHIP_DEX_BRD4_MF1_rep1_RS/A549_CHIP_DEX_BRD4_MF1_rep1_RS.sorted.filtered.bam
samtools_view_filter.A549_CHIP_DEX_BRD4_MF1_rep1_RS.a498a7628fa7a4993f5ca6de2988e020.mugqic.done
)
samtools_view_filter_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_7_JOB_ID: samtools_view_filter_report
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter_report
JOB_DEPENDENCIES=$samtools_view_filter_1_JOB_ID:$samtools_view_filter_2_JOB_ID:$samtools_view_filter_3_JOB_ID:$samtools_view_filter_4_JOB_ID:$samtools_view_filter_5_JOB_ID:$samtools_view_filter_6_JOB_ID
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
samtools_view_filter_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: picard_merge_sam_files
#-------------------------------------------------------------------------------
STEP=picard_merge_sam_files
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_1_JOB_ID: symlink_readset_sample_bam.A549_CHIP_DEX_WCE_MF1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.A549_CHIP_DEX_WCE_MF1_rep1
JOB_DEPENDENCIES=$samtools_view_filter_1_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.A549_CHIP_DEX_WCE_MF1_rep1.760c8b5221c1a3c43dd193944e8f4f5f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.A549_CHIP_DEX_WCE_MF1_rep1.760c8b5221c1a3c43dd193944e8f4f5f.mugqic.done'
mkdir -p alignment/A549_CHIP_DEX_WCE_MF1_rep1 && \
ln -s -f A549_CHIP_DEX_WCE_MF1_rep1_RS/A549_CHIP_DEX_WCE_MF1_rep1_RS.sorted.filtered.bam alignment/A549_CHIP_DEX_WCE_MF1_rep1/A549_CHIP_DEX_WCE_MF1_rep1.merged.bam
symlink_readset_sample_bam.A549_CHIP_DEX_WCE_MF1_rep1.760c8b5221c1a3c43dd193944e8f4f5f.mugqic.done
)
picard_merge_sam_files_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_2_JOB_ID: symlink_readset_sample_bam.A549_CHIP_DEX_NIPBL_MF1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.A549_CHIP_DEX_NIPBL_MF1_rep1
JOB_DEPENDENCIES=$samtools_view_filter_2_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.A549_CHIP_DEX_NIPBL_MF1_rep1.7557413415fe02eac81d4bfb16844afb.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.A549_CHIP_DEX_NIPBL_MF1_rep1.7557413415fe02eac81d4bfb16844afb.mugqic.done'
mkdir -p alignment/A549_CHIP_DEX_NIPBL_MF1_rep1 && \
ln -s -f A549_CHIP_DEX_NIPBL_MF1_rep1_RS/A549_CHIP_DEX_NIPBL_MF1_rep1_RS.sorted.filtered.bam alignment/A549_CHIP_DEX_NIPBL_MF1_rep1/A549_CHIP_DEX_NIPBL_MF1_rep1.merged.bam
symlink_readset_sample_bam.A549_CHIP_DEX_NIPBL_MF1_rep1.7557413415fe02eac81d4bfb16844afb.mugqic.done
)
picard_merge_sam_files_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_3_JOB_ID: symlink_readset_sample_bam.A549_CHIP_DEX_MED1_MF1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.A549_CHIP_DEX_MED1_MF1_rep1
JOB_DEPENDENCIES=$samtools_view_filter_3_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.A549_CHIP_DEX_MED1_MF1_rep1.8ae64304e8466456222bb5179e7aa153.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.A549_CHIP_DEX_MED1_MF1_rep1.8ae64304e8466456222bb5179e7aa153.mugqic.done'
mkdir -p alignment/A549_CHIP_DEX_MED1_MF1_rep1 && \
ln -s -f A549_CHIP_DEX_MED1_MF1_rep1_RS/A549_CHIP_DEX_MED1_MF1_rep1_RS.sorted.filtered.bam alignment/A549_CHIP_DEX_MED1_MF1_rep1/A549_CHIP_DEX_MED1_MF1_rep1.merged.bam
symlink_readset_sample_bam.A549_CHIP_DEX_MED1_MF1_rep1.8ae64304e8466456222bb5179e7aa153.mugqic.done
)
picard_merge_sam_files_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_4_JOB_ID: symlink_readset_sample_bam.A549_CHIP_DEX_SMC1A_MF1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.A549_CHIP_DEX_SMC1A_MF1_rep1
JOB_DEPENDENCIES=$samtools_view_filter_4_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.A549_CHIP_DEX_SMC1A_MF1_rep1.d795174e85f8e6e37fedfde0b5300405.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.A549_CHIP_DEX_SMC1A_MF1_rep1.d795174e85f8e6e37fedfde0b5300405.mugqic.done'
mkdir -p alignment/A549_CHIP_DEX_SMC1A_MF1_rep1 && \
ln -s -f A549_CHIP_DEX_SMC1A_MF1_rep1_RS/A549_CHIP_DEX_SMC1A_MF1_rep1_RS.sorted.filtered.bam alignment/A549_CHIP_DEX_SMC1A_MF1_rep1/A549_CHIP_DEX_SMC1A_MF1_rep1.merged.bam
symlink_readset_sample_bam.A549_CHIP_DEX_SMC1A_MF1_rep1.d795174e85f8e6e37fedfde0b5300405.mugqic.done
)
picard_merge_sam_files_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_5_JOB_ID: symlink_readset_sample_bam.A549_CHIP_DEX_CDK9_MF1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.A549_CHIP_DEX_CDK9_MF1_rep1
JOB_DEPENDENCIES=$samtools_view_filter_5_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.A549_CHIP_DEX_CDK9_MF1_rep1.7a3af75bd307c3ad14ed78229a4c9890.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.A549_CHIP_DEX_CDK9_MF1_rep1.7a3af75bd307c3ad14ed78229a4c9890.mugqic.done'
mkdir -p alignment/A549_CHIP_DEX_CDK9_MF1_rep1 && \
ln -s -f A549_CHIP_DEX_CDK9_MF1_rep1_RS/A549_CHIP_DEX_CDK9_MF1_rep1_RS.sorted.filtered.bam alignment/A549_CHIP_DEX_CDK9_MF1_rep1/A549_CHIP_DEX_CDK9_MF1_rep1.merged.bam
symlink_readset_sample_bam.A549_CHIP_DEX_CDK9_MF1_rep1.7a3af75bd307c3ad14ed78229a4c9890.mugqic.done
)
picard_merge_sam_files_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_6_JOB_ID: symlink_readset_sample_bam.A549_CHIP_DEX_BRD4_MF1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.A549_CHIP_DEX_BRD4_MF1_rep1
JOB_DEPENDENCIES=$samtools_view_filter_6_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.A549_CHIP_DEX_BRD4_MF1_rep1.868a27240bc063b68d3169ed92011cbd.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.A549_CHIP_DEX_BRD4_MF1_rep1.868a27240bc063b68d3169ed92011cbd.mugqic.done'
mkdir -p alignment/A549_CHIP_DEX_BRD4_MF1_rep1 && \
ln -s -f A549_CHIP_DEX_BRD4_MF1_rep1_RS/A549_CHIP_DEX_BRD4_MF1_rep1_RS.sorted.filtered.bam alignment/A549_CHIP_DEX_BRD4_MF1_rep1/A549_CHIP_DEX_BRD4_MF1_rep1.merged.bam
symlink_readset_sample_bam.A549_CHIP_DEX_BRD4_MF1_rep1.868a27240bc063b68d3169ed92011cbd.mugqic.done
)
picard_merge_sam_files_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: picard_mark_duplicates
#-------------------------------------------------------------------------------
STEP=picard_mark_duplicates
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_1_JOB_ID: picard_mark_duplicates.A549_CHIP_DEX_WCE_MF1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.A549_CHIP_DEX_WCE_MF1_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_1_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.A549_CHIP_DEX_WCE_MF1_rep1.aaacbabef5efa57931b708e852d94618.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.A549_CHIP_DEX_WCE_MF1_rep1.aaacbabef5efa57931b708e852d94618.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/A549_CHIP_DEX_WCE_MF1_rep1/A549_CHIP_DEX_WCE_MF1_rep1.merged.bam \
  OUTPUT=alignment/A549_CHIP_DEX_WCE_MF1_rep1/A549_CHIP_DEX_WCE_MF1_rep1.sorted.dup.bam \
  METRICS_FILE=alignment/A549_CHIP_DEX_WCE_MF1_rep1/A549_CHIP_DEX_WCE_MF1_rep1.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.A549_CHIP_DEX_WCE_MF1_rep1.aaacbabef5efa57931b708e852d94618.mugqic.done
)
picard_mark_duplicates_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_2_JOB_ID: picard_mark_duplicates.A549_CHIP_DEX_NIPBL_MF1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.A549_CHIP_DEX_NIPBL_MF1_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_2_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.A549_CHIP_DEX_NIPBL_MF1_rep1.1d598203f0030461ecab0f5d0ecffec7.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.A549_CHIP_DEX_NIPBL_MF1_rep1.1d598203f0030461ecab0f5d0ecffec7.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/A549_CHIP_DEX_NIPBL_MF1_rep1/A549_CHIP_DEX_NIPBL_MF1_rep1.merged.bam \
  OUTPUT=alignment/A549_CHIP_DEX_NIPBL_MF1_rep1/A549_CHIP_DEX_NIPBL_MF1_rep1.sorted.dup.bam \
  METRICS_FILE=alignment/A549_CHIP_DEX_NIPBL_MF1_rep1/A549_CHIP_DEX_NIPBL_MF1_rep1.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.A549_CHIP_DEX_NIPBL_MF1_rep1.1d598203f0030461ecab0f5d0ecffec7.mugqic.done
)
picard_mark_duplicates_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_3_JOB_ID: picard_mark_duplicates.A549_CHIP_DEX_MED1_MF1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.A549_CHIP_DEX_MED1_MF1_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_3_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.A549_CHIP_DEX_MED1_MF1_rep1.b4ccd9d6504d647424f2efb945b4fdca.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.A549_CHIP_DEX_MED1_MF1_rep1.b4ccd9d6504d647424f2efb945b4fdca.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/A549_CHIP_DEX_MED1_MF1_rep1/A549_CHIP_DEX_MED1_MF1_rep1.merged.bam \
  OUTPUT=alignment/A549_CHIP_DEX_MED1_MF1_rep1/A549_CHIP_DEX_MED1_MF1_rep1.sorted.dup.bam \
  METRICS_FILE=alignment/A549_CHIP_DEX_MED1_MF1_rep1/A549_CHIP_DEX_MED1_MF1_rep1.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.A549_CHIP_DEX_MED1_MF1_rep1.b4ccd9d6504d647424f2efb945b4fdca.mugqic.done
)
picard_mark_duplicates_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_4_JOB_ID: picard_mark_duplicates.A549_CHIP_DEX_SMC1A_MF1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.A549_CHIP_DEX_SMC1A_MF1_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_4_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.A549_CHIP_DEX_SMC1A_MF1_rep1.5564205ef192528fde9a43412cc01702.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.A549_CHIP_DEX_SMC1A_MF1_rep1.5564205ef192528fde9a43412cc01702.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/A549_CHIP_DEX_SMC1A_MF1_rep1/A549_CHIP_DEX_SMC1A_MF1_rep1.merged.bam \
  OUTPUT=alignment/A549_CHIP_DEX_SMC1A_MF1_rep1/A549_CHIP_DEX_SMC1A_MF1_rep1.sorted.dup.bam \
  METRICS_FILE=alignment/A549_CHIP_DEX_SMC1A_MF1_rep1/A549_CHIP_DEX_SMC1A_MF1_rep1.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.A549_CHIP_DEX_SMC1A_MF1_rep1.5564205ef192528fde9a43412cc01702.mugqic.done
)
picard_mark_duplicates_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_5_JOB_ID: picard_mark_duplicates.A549_CHIP_DEX_CDK9_MF1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.A549_CHIP_DEX_CDK9_MF1_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_5_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.A549_CHIP_DEX_CDK9_MF1_rep1.d06bc32a2a1968583e07ab3053b8ad3e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.A549_CHIP_DEX_CDK9_MF1_rep1.d06bc32a2a1968583e07ab3053b8ad3e.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/A549_CHIP_DEX_CDK9_MF1_rep1/A549_CHIP_DEX_CDK9_MF1_rep1.merged.bam \
  OUTPUT=alignment/A549_CHIP_DEX_CDK9_MF1_rep1/A549_CHIP_DEX_CDK9_MF1_rep1.sorted.dup.bam \
  METRICS_FILE=alignment/A549_CHIP_DEX_CDK9_MF1_rep1/A549_CHIP_DEX_CDK9_MF1_rep1.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.A549_CHIP_DEX_CDK9_MF1_rep1.d06bc32a2a1968583e07ab3053b8ad3e.mugqic.done
)
picard_mark_duplicates_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_6_JOB_ID: picard_mark_duplicates.A549_CHIP_DEX_BRD4_MF1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.A549_CHIP_DEX_BRD4_MF1_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_6_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.A549_CHIP_DEX_BRD4_MF1_rep1.aad8fa05dc641593779a044aacdc5125.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.A549_CHIP_DEX_BRD4_MF1_rep1.aad8fa05dc641593779a044aacdc5125.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/A549_CHIP_DEX_BRD4_MF1_rep1/A549_CHIP_DEX_BRD4_MF1_rep1.merged.bam \
  OUTPUT=alignment/A549_CHIP_DEX_BRD4_MF1_rep1/A549_CHIP_DEX_BRD4_MF1_rep1.sorted.dup.bam \
  METRICS_FILE=alignment/A549_CHIP_DEX_BRD4_MF1_rep1/A549_CHIP_DEX_BRD4_MF1_rep1.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.A549_CHIP_DEX_BRD4_MF1_rep1.aad8fa05dc641593779a044aacdc5125.mugqic.done
)
picard_mark_duplicates_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_7_JOB_ID: picard_mark_duplicates_report
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates_report
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_2_JOB_ID:$picard_mark_duplicates_3_JOB_ID:$picard_mark_duplicates_4_JOB_ID:$picard_mark_duplicates_5_JOB_ID:$picard_mark_duplicates_6_JOB_ID
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
picard_mark_duplicates_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: metrics
#-------------------------------------------------------------------------------
STEP=metrics
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: metrics_1_JOB_ID: metrics.flagstat
#-------------------------------------------------------------------------------
JOB_NAME=metrics.flagstat
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_2_JOB_ID:$picard_mark_duplicates_3_JOB_ID:$picard_mark_duplicates_4_JOB_ID:$picard_mark_duplicates_5_JOB_ID:$picard_mark_duplicates_6_JOB_ID
JOB_DONE=job_output/metrics/metrics.flagstat.87387a3ff30e6afbb6fce3478b80ae87.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'metrics.flagstat.87387a3ff30e6afbb6fce3478b80ae87.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools flagstat \
  alignment/A549_CHIP_DEX_WCE_MF1_rep1/A549_CHIP_DEX_WCE_MF1_rep1.sorted.dup.bam \
  > alignment/A549_CHIP_DEX_WCE_MF1_rep1/A549_CHIP_DEX_WCE_MF1_rep1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/A549_CHIP_DEX_NIPBL_MF1_rep1/A549_CHIP_DEX_NIPBL_MF1_rep1.sorted.dup.bam \
  > alignment/A549_CHIP_DEX_NIPBL_MF1_rep1/A549_CHIP_DEX_NIPBL_MF1_rep1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/A549_CHIP_DEX_MED1_MF1_rep1/A549_CHIP_DEX_MED1_MF1_rep1.sorted.dup.bam \
  > alignment/A549_CHIP_DEX_MED1_MF1_rep1/A549_CHIP_DEX_MED1_MF1_rep1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/A549_CHIP_DEX_SMC1A_MF1_rep1/A549_CHIP_DEX_SMC1A_MF1_rep1.sorted.dup.bam \
  > alignment/A549_CHIP_DEX_SMC1A_MF1_rep1/A549_CHIP_DEX_SMC1A_MF1_rep1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/A549_CHIP_DEX_CDK9_MF1_rep1/A549_CHIP_DEX_CDK9_MF1_rep1.sorted.dup.bam \
  > alignment/A549_CHIP_DEX_CDK9_MF1_rep1/A549_CHIP_DEX_CDK9_MF1_rep1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/A549_CHIP_DEX_BRD4_MF1_rep1/A549_CHIP_DEX_BRD4_MF1_rep1.sorted.dup.bam \
  > alignment/A549_CHIP_DEX_BRD4_MF1_rep1/A549_CHIP_DEX_BRD4_MF1_rep1.sorted.dup.bam.flagstat
metrics.flagstat.87387a3ff30e6afbb6fce3478b80ae87.mugqic.done
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
JOB_DONE=job_output/metrics/metrics_report.bd34eba763c95a3d92fac9ef8a0c7779.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'metrics_report.bd34eba763c95a3d92fac9ef8a0c7779.mugqic.done'
module load mugqic/pandoc/1.15.2 && \
for sample in A549_CHIP_DEX_WCE_MF1_rep1 A549_CHIP_DEX_NIPBL_MF1_rep1 A549_CHIP_DEX_MED1_MF1_rep1 A549_CHIP_DEX_SMC1A_MF1_rep1 A549_CHIP_DEX_CDK9_MF1_rep1 A549_CHIP_DEX_BRD4_MF1_rep1
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

metrics_report.bd34eba763c95a3d92fac9ef8a0c7779.mugqic.done
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
# JOB: homer_make_tag_directory_1_JOB_ID: homer_make_tag_directory.A549_CHIP_DEX_WCE_MF1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.A549_CHIP_DEX_WCE_MF1_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.A549_CHIP_DEX_WCE_MF1_rep1.c9140a2625b3b262e4de88f0555f915e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.A549_CHIP_DEX_WCE_MF1_rep1.c9140a2625b3b262e4de88f0555f915e.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/A549_CHIP_DEX_WCE_MF1_rep1 \
  alignment/A549_CHIP_DEX_WCE_MF1_rep1/A549_CHIP_DEX_WCE_MF1_rep1.sorted.dup.bam \
  -checkGC -genome mm10
homer_make_tag_directory.A549_CHIP_DEX_WCE_MF1_rep1.c9140a2625b3b262e4de88f0555f915e.mugqic.done
)
homer_make_tag_directory_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_2_JOB_ID: homer_make_tag_directory.A549_CHIP_DEX_NIPBL_MF1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.A549_CHIP_DEX_NIPBL_MF1_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_2_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.A549_CHIP_DEX_NIPBL_MF1_rep1.9a700cd40c43e5c97b9dea58387679c8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.A549_CHIP_DEX_NIPBL_MF1_rep1.9a700cd40c43e5c97b9dea58387679c8.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/A549_CHIP_DEX_NIPBL_MF1_rep1 \
  alignment/A549_CHIP_DEX_NIPBL_MF1_rep1/A549_CHIP_DEX_NIPBL_MF1_rep1.sorted.dup.bam \
  -checkGC -genome mm10
homer_make_tag_directory.A549_CHIP_DEX_NIPBL_MF1_rep1.9a700cd40c43e5c97b9dea58387679c8.mugqic.done
)
homer_make_tag_directory_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_3_JOB_ID: homer_make_tag_directory.A549_CHIP_DEX_MED1_MF1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.A549_CHIP_DEX_MED1_MF1_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_3_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.A549_CHIP_DEX_MED1_MF1_rep1.4aca34c8733d02d1fa162fbecae9ec0c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.A549_CHIP_DEX_MED1_MF1_rep1.4aca34c8733d02d1fa162fbecae9ec0c.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/A549_CHIP_DEX_MED1_MF1_rep1 \
  alignment/A549_CHIP_DEX_MED1_MF1_rep1/A549_CHIP_DEX_MED1_MF1_rep1.sorted.dup.bam \
  -checkGC -genome mm10
homer_make_tag_directory.A549_CHIP_DEX_MED1_MF1_rep1.4aca34c8733d02d1fa162fbecae9ec0c.mugqic.done
)
homer_make_tag_directory_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_4_JOB_ID: homer_make_tag_directory.A549_CHIP_DEX_SMC1A_MF1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.A549_CHIP_DEX_SMC1A_MF1_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_4_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.A549_CHIP_DEX_SMC1A_MF1_rep1.50e9cf27677ddad4808f5573b72370c1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.A549_CHIP_DEX_SMC1A_MF1_rep1.50e9cf27677ddad4808f5573b72370c1.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/A549_CHIP_DEX_SMC1A_MF1_rep1 \
  alignment/A549_CHIP_DEX_SMC1A_MF1_rep1/A549_CHIP_DEX_SMC1A_MF1_rep1.sorted.dup.bam \
  -checkGC -genome mm10
homer_make_tag_directory.A549_CHIP_DEX_SMC1A_MF1_rep1.50e9cf27677ddad4808f5573b72370c1.mugqic.done
)
homer_make_tag_directory_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_5_JOB_ID: homer_make_tag_directory.A549_CHIP_DEX_CDK9_MF1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.A549_CHIP_DEX_CDK9_MF1_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_5_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.A549_CHIP_DEX_CDK9_MF1_rep1.b5d53f7de1e64d98c6c08ba0d17537bc.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.A549_CHIP_DEX_CDK9_MF1_rep1.b5d53f7de1e64d98c6c08ba0d17537bc.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/A549_CHIP_DEX_CDK9_MF1_rep1 \
  alignment/A549_CHIP_DEX_CDK9_MF1_rep1/A549_CHIP_DEX_CDK9_MF1_rep1.sorted.dup.bam \
  -checkGC -genome mm10
homer_make_tag_directory.A549_CHIP_DEX_CDK9_MF1_rep1.b5d53f7de1e64d98c6c08ba0d17537bc.mugqic.done
)
homer_make_tag_directory_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_6_JOB_ID: homer_make_tag_directory.A549_CHIP_DEX_BRD4_MF1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.A549_CHIP_DEX_BRD4_MF1_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_6_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.A549_CHIP_DEX_BRD4_MF1_rep1.236654b77dcfe06fca8cd3f130bf6a6b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.A549_CHIP_DEX_BRD4_MF1_rep1.236654b77dcfe06fca8cd3f130bf6a6b.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/A549_CHIP_DEX_BRD4_MF1_rep1 \
  alignment/A549_CHIP_DEX_BRD4_MF1_rep1/A549_CHIP_DEX_BRD4_MF1_rep1.sorted.dup.bam \
  -checkGC -genome mm10
homer_make_tag_directory.A549_CHIP_DEX_BRD4_MF1_rep1.236654b77dcfe06fca8cd3f130bf6a6b.mugqic.done
)
homer_make_tag_directory_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: qc_metrics
#-------------------------------------------------------------------------------
STEP=qc_metrics
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: qc_metrics_1_JOB_ID: qc_plots_R
#-------------------------------------------------------------------------------
JOB_NAME=qc_plots_R
JOB_DEPENDENCIES=$homer_make_tag_directory_1_JOB_ID:$homer_make_tag_directory_2_JOB_ID:$homer_make_tag_directory_3_JOB_ID:$homer_make_tag_directory_4_JOB_ID:$homer_make_tag_directory_5_JOB_ID:$homer_make_tag_directory_6_JOB_ID
JOB_DONE=job_output/qc_metrics/qc_plots_R.8ef22ad92c71d1a312e691295adb7b38.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'qc_plots_R.8ef22ad92c71d1a312e691295adb7b38.mugqic.done'
module load mugqic/mugqic_tools/2.1.5 mugqic/R_Bioconductor/3.2.3_3.2 && \
mkdir -p graphs && \
Rscript $R_TOOLS/chipSeqGenerateQCMetrics.R \
  ../../raw/chip-seq/design.txt \
  /gs/scratch/efournier/CofactorHR/A549/output/chip-pipeline && \
cp /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/ChipSeq.qc_metrics.md report/ChipSeq.qc_metrics.md && \
for sample in A549_CHIP_DEX_WCE_MF1_rep1 A549_CHIP_DEX_NIPBL_MF1_rep1 A549_CHIP_DEX_MED1_MF1_rep1 A549_CHIP_DEX_SMC1A_MF1_rep1 A549_CHIP_DEX_CDK9_MF1_rep1 A549_CHIP_DEX_BRD4_MF1_rep1
do
  cp --parents graphs/${sample}_QC_Metrics.ps report/
  convert -rotate 90 graphs/${sample}_QC_Metrics.ps report/graphs/${sample}_QC_Metrics.png
  echo -e "----

![QC Metrics for Sample $sample ([download high-res image](graphs/${sample}_QC_Metrics.ps))](graphs/${sample}_QC_Metrics.png)
" \
  >> report/ChipSeq.qc_metrics.md
done
qc_plots_R.8ef22ad92c71d1a312e691295adb7b38.mugqic.done
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
# JOB: homer_make_ucsc_file_1_JOB_ID: homer_make_ucsc_file.A549_CHIP_DEX_WCE_MF1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.A549_CHIP_DEX_WCE_MF1_rep1
JOB_DEPENDENCIES=$homer_make_tag_directory_1_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.A549_CHIP_DEX_WCE_MF1_rep1.eb31a0d4c3954a647f76beeb9e3fbae5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.A549_CHIP_DEX_WCE_MF1_rep1.eb31a0d4c3954a647f76beeb9e3fbae5.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/A549_CHIP_DEX_WCE_MF1_rep1 && \
makeUCSCfile \
  tags/A549_CHIP_DEX_WCE_MF1_rep1 | \
gzip -1 -c > tracks/A549_CHIP_DEX_WCE_MF1_rep1/A549_CHIP_DEX_WCE_MF1_rep1.ucsc.bedGraph.gz
homer_make_ucsc_file.A549_CHIP_DEX_WCE_MF1_rep1.eb31a0d4c3954a647f76beeb9e3fbae5.mugqic.done
)
homer_make_ucsc_file_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_2_JOB_ID: homer_make_ucsc_file.A549_CHIP_DEX_NIPBL_MF1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.A549_CHIP_DEX_NIPBL_MF1_rep1
JOB_DEPENDENCIES=$homer_make_tag_directory_2_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.A549_CHIP_DEX_NIPBL_MF1_rep1.b5fdf5496d80a17dc6bf8541e93f444a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.A549_CHIP_DEX_NIPBL_MF1_rep1.b5fdf5496d80a17dc6bf8541e93f444a.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/A549_CHIP_DEX_NIPBL_MF1_rep1 && \
makeUCSCfile \
  tags/A549_CHIP_DEX_NIPBL_MF1_rep1 | \
gzip -1 -c > tracks/A549_CHIP_DEX_NIPBL_MF1_rep1/A549_CHIP_DEX_NIPBL_MF1_rep1.ucsc.bedGraph.gz
homer_make_ucsc_file.A549_CHIP_DEX_NIPBL_MF1_rep1.b5fdf5496d80a17dc6bf8541e93f444a.mugqic.done
)
homer_make_ucsc_file_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_3_JOB_ID: homer_make_ucsc_file.A549_CHIP_DEX_MED1_MF1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.A549_CHIP_DEX_MED1_MF1_rep1
JOB_DEPENDENCIES=$homer_make_tag_directory_3_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.A549_CHIP_DEX_MED1_MF1_rep1.37536a184c4f014609134e00dde1aa66.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.A549_CHIP_DEX_MED1_MF1_rep1.37536a184c4f014609134e00dde1aa66.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/A549_CHIP_DEX_MED1_MF1_rep1 && \
makeUCSCfile \
  tags/A549_CHIP_DEX_MED1_MF1_rep1 | \
gzip -1 -c > tracks/A549_CHIP_DEX_MED1_MF1_rep1/A549_CHIP_DEX_MED1_MF1_rep1.ucsc.bedGraph.gz
homer_make_ucsc_file.A549_CHIP_DEX_MED1_MF1_rep1.37536a184c4f014609134e00dde1aa66.mugqic.done
)
homer_make_ucsc_file_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_4_JOB_ID: homer_make_ucsc_file.A549_CHIP_DEX_SMC1A_MF1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.A549_CHIP_DEX_SMC1A_MF1_rep1
JOB_DEPENDENCIES=$homer_make_tag_directory_4_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.A549_CHIP_DEX_SMC1A_MF1_rep1.f55b84dccb73ac5873a57739fa6814d2.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.A549_CHIP_DEX_SMC1A_MF1_rep1.f55b84dccb73ac5873a57739fa6814d2.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/A549_CHIP_DEX_SMC1A_MF1_rep1 && \
makeUCSCfile \
  tags/A549_CHIP_DEX_SMC1A_MF1_rep1 | \
gzip -1 -c > tracks/A549_CHIP_DEX_SMC1A_MF1_rep1/A549_CHIP_DEX_SMC1A_MF1_rep1.ucsc.bedGraph.gz
homer_make_ucsc_file.A549_CHIP_DEX_SMC1A_MF1_rep1.f55b84dccb73ac5873a57739fa6814d2.mugqic.done
)
homer_make_ucsc_file_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_5_JOB_ID: homer_make_ucsc_file.A549_CHIP_DEX_CDK9_MF1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.A549_CHIP_DEX_CDK9_MF1_rep1
JOB_DEPENDENCIES=$homer_make_tag_directory_5_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.A549_CHIP_DEX_CDK9_MF1_rep1.040796f4f3b22384020fad7d259f3306.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.A549_CHIP_DEX_CDK9_MF1_rep1.040796f4f3b22384020fad7d259f3306.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/A549_CHIP_DEX_CDK9_MF1_rep1 && \
makeUCSCfile \
  tags/A549_CHIP_DEX_CDK9_MF1_rep1 | \
gzip -1 -c > tracks/A549_CHIP_DEX_CDK9_MF1_rep1/A549_CHIP_DEX_CDK9_MF1_rep1.ucsc.bedGraph.gz
homer_make_ucsc_file.A549_CHIP_DEX_CDK9_MF1_rep1.040796f4f3b22384020fad7d259f3306.mugqic.done
)
homer_make_ucsc_file_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_6_JOB_ID: homer_make_ucsc_file.A549_CHIP_DEX_BRD4_MF1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.A549_CHIP_DEX_BRD4_MF1_rep1
JOB_DEPENDENCIES=$homer_make_tag_directory_6_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.A549_CHIP_DEX_BRD4_MF1_rep1.4fbc510b3bc2d47b04c00a7b30e6917a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.A549_CHIP_DEX_BRD4_MF1_rep1.4fbc510b3bc2d47b04c00a7b30e6917a.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/A549_CHIP_DEX_BRD4_MF1_rep1 && \
makeUCSCfile \
  tags/A549_CHIP_DEX_BRD4_MF1_rep1 | \
gzip -1 -c > tracks/A549_CHIP_DEX_BRD4_MF1_rep1/A549_CHIP_DEX_BRD4_MF1_rep1.ucsc.bedGraph.gz
homer_make_ucsc_file.A549_CHIP_DEX_BRD4_MF1_rep1.4fbc510b3bc2d47b04c00a7b30e6917a.mugqic.done
)
homer_make_ucsc_file_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_7_JOB_ID: homer_make_ucsc_file_report
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file_report
JOB_DEPENDENCIES=$homer_make_ucsc_file_6_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file_report.c9cac04b3bc71f34cbcde042cf5b419d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file_report.c9cac04b3bc71f34cbcde042cf5b419d.mugqic.done'
mkdir -p report && \
zip -r report/tracks.zip tracks/*/*.ucsc.bedGraph.gz && \
cp /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/ChipSeq.homer_make_ucsc_file.md report/
homer_make_ucsc_file_report.c9cac04b3bc71f34cbcde042cf5b419d.mugqic.done
)
homer_make_ucsc_file_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: macs2_callpeak
#-------------------------------------------------------------------------------
STEP=macs2_callpeak
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_1_JOB_ID: macs2_callpeak.A549-DEX_NIPBL_MF1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549-DEX_NIPBL_MF1
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_2_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549-DEX_NIPBL_MF1.ea6d0a732f1dab83db3a1e05c7e9fb98.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549-DEX_NIPBL_MF1.ea6d0a732f1dab83db3a1e05c7e9fb98.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/A549-DEX_NIPBL_MF1 && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2184697419.2 \
  --treatment \
  alignment/A549_CHIP_DEX_NIPBL_MF1_rep1/A549_CHIP_DEX_NIPBL_MF1_rep1.sorted.dup.bam \
  --control \
  alignment/A549_CHIP_DEX_WCE_MF1_rep1/A549_CHIP_DEX_WCE_MF1_rep1.sorted.dup.bam \
  --name peak_call/A549-DEX_NIPBL_MF1/A549-DEX_NIPBL_MF1 \
  >& peak_call/A549-DEX_NIPBL_MF1/A549-DEX_NIPBL_MF1.diag.macs.out
macs2_callpeak.A549-DEX_NIPBL_MF1.ea6d0a732f1dab83db3a1e05c7e9fb98.mugqic.done
)
macs2_callpeak_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_2_JOB_ID: macs2_callpeak.A549-DEX_MED1_MF1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549-DEX_MED1_MF1
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_3_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549-DEX_MED1_MF1.e18d9ad15622958f115f48bfa1ee00ac.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549-DEX_MED1_MF1.e18d9ad15622958f115f48bfa1ee00ac.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/A549-DEX_MED1_MF1 && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2184697419.2 \
  --treatment \
  alignment/A549_CHIP_DEX_MED1_MF1_rep1/A549_CHIP_DEX_MED1_MF1_rep1.sorted.dup.bam \
  --control \
  alignment/A549_CHIP_DEX_WCE_MF1_rep1/A549_CHIP_DEX_WCE_MF1_rep1.sorted.dup.bam \
  --name peak_call/A549-DEX_MED1_MF1/A549-DEX_MED1_MF1 \
  >& peak_call/A549-DEX_MED1_MF1/A549-DEX_MED1_MF1.diag.macs.out
macs2_callpeak.A549-DEX_MED1_MF1.e18d9ad15622958f115f48bfa1ee00ac.mugqic.done
)
macs2_callpeak_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_3_JOB_ID: macs2_callpeak.A549-DEX_SMC1A_MF1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549-DEX_SMC1A_MF1
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_4_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549-DEX_SMC1A_MF1.a6f06ddcaabb849d9e35049b5b0c4ade.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549-DEX_SMC1A_MF1.a6f06ddcaabb849d9e35049b5b0c4ade.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/A549-DEX_SMC1A_MF1 && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2184697419.2 \
  --treatment \
  alignment/A549_CHIP_DEX_SMC1A_MF1_rep1/A549_CHIP_DEX_SMC1A_MF1_rep1.sorted.dup.bam \
  --control \
  alignment/A549_CHIP_DEX_WCE_MF1_rep1/A549_CHIP_DEX_WCE_MF1_rep1.sorted.dup.bam \
  --name peak_call/A549-DEX_SMC1A_MF1/A549-DEX_SMC1A_MF1 \
  >& peak_call/A549-DEX_SMC1A_MF1/A549-DEX_SMC1A_MF1.diag.macs.out
macs2_callpeak.A549-DEX_SMC1A_MF1.a6f06ddcaabb849d9e35049b5b0c4ade.mugqic.done
)
macs2_callpeak_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_4_JOB_ID: macs2_callpeak.A549-DEX_CDK9_MF1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549-DEX_CDK9_MF1
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_5_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549-DEX_CDK9_MF1.d9ab2443dfc25b7a729e8127260ff0f6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549-DEX_CDK9_MF1.d9ab2443dfc25b7a729e8127260ff0f6.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/A549-DEX_CDK9_MF1 && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2184697419.2 \
  --treatment \
  alignment/A549_CHIP_DEX_CDK9_MF1_rep1/A549_CHIP_DEX_CDK9_MF1_rep1.sorted.dup.bam \
  --control \
  alignment/A549_CHIP_DEX_WCE_MF1_rep1/A549_CHIP_DEX_WCE_MF1_rep1.sorted.dup.bam \
  --name peak_call/A549-DEX_CDK9_MF1/A549-DEX_CDK9_MF1 \
  >& peak_call/A549-DEX_CDK9_MF1/A549-DEX_CDK9_MF1.diag.macs.out
macs2_callpeak.A549-DEX_CDK9_MF1.d9ab2443dfc25b7a729e8127260ff0f6.mugqic.done
)
macs2_callpeak_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_5_JOB_ID: macs2_callpeak.A549-DEX_BRD4_MF1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549-DEX_BRD4_MF1
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_6_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549-DEX_BRD4_MF1.f5b63b3fe5ca9e05dcf917e8862aa1ee.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549-DEX_BRD4_MF1.f5b63b3fe5ca9e05dcf917e8862aa1ee.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/A549-DEX_BRD4_MF1 && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2184697419.2 \
  --treatment \
  alignment/A549_CHIP_DEX_BRD4_MF1_rep1/A549_CHIP_DEX_BRD4_MF1_rep1.sorted.dup.bam \
  --control \
  alignment/A549_CHIP_DEX_WCE_MF1_rep1/A549_CHIP_DEX_WCE_MF1_rep1.sorted.dup.bam \
  --name peak_call/A549-DEX_BRD4_MF1/A549-DEX_BRD4_MF1 \
  >& peak_call/A549-DEX_BRD4_MF1/A549-DEX_BRD4_MF1.diag.macs.out
macs2_callpeak.A549-DEX_BRD4_MF1.f5b63b3fe5ca9e05dcf917e8862aa1ee.mugqic.done
)
macs2_callpeak_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_6_JOB_ID: macs2_callpeak_report
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_report
JOB_DEPENDENCIES=$macs2_callpeak_1_JOB_ID:$macs2_callpeak_2_JOB_ID:$macs2_callpeak_3_JOB_ID:$macs2_callpeak_4_JOB_ID:$macs2_callpeak_5_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_report.a0e7dca0269f0310058013e66c984ed3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_report.a0e7dca0269f0310058013e66c984ed3.mugqic.done'
mkdir -p report && \
cp /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/ChipSeq.macs2_callpeak.md report/ && \
for contrast in A549-DEX_NIPBL_MF1 A549-DEX_MED1_MF1 A549-DEX_SMC1A_MF1 A549-DEX_CDK9_MF1 A549-DEX_BRD4_MF1
do
  cp -a --parents peak_call/$contrast/ report/ && \
  echo -e "* [Peak Calls File for Design $contrast](peak_call/$contrast/${contrast}_peaks.xls)" \
  >> report/ChipSeq.macs2_callpeak.md
done
macs2_callpeak_report.a0e7dca0269f0310058013e66c984ed3.mugqic.done
)
macs2_callpeak_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: homer_annotate_peaks
#-------------------------------------------------------------------------------
STEP=homer_annotate_peaks
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_1_JOB_ID: homer_annotate_peaks.A549-DEX_NIPBL_MF1
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.A549-DEX_NIPBL_MF1
JOB_DEPENDENCIES=$macs2_callpeak_1_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.A549-DEX_NIPBL_MF1.de663d0290e352e6568ab4a577c949c4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.A549-DEX_NIPBL_MF1.de663d0290e352e6568ab4a577c949c4.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/A549-DEX_NIPBL_MF1/A549-DEX_NIPBL_MF1 && \
annotatePeaks.pl \
  peak_call/A549-DEX_NIPBL_MF1/A549-DEX_NIPBL_MF1_peaks.narrowPeak \
  mm10 \
  -gsize mm10 \
  -cons -CpG \
  -go annotation/A549-DEX_NIPBL_MF1/A549-DEX_NIPBL_MF1 \
  -genomeOntology annotation/A549-DEX_NIPBL_MF1/A549-DEX_NIPBL_MF1 \
  > annotation/A549-DEX_NIPBL_MF1/A549-DEX_NIPBL_MF1.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/A549-DEX_NIPBL_MF1/A549-DEX_NIPBL_MF1.annotated.csv",
  "annotation/A549-DEX_NIPBL_MF1/A549-DEX_NIPBL_MF1",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.A549-DEX_NIPBL_MF1.de663d0290e352e6568ab4a577c949c4.mugqic.done
)
homer_annotate_peaks_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_2_JOB_ID: homer_annotate_peaks.A549-DEX_MED1_MF1
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.A549-DEX_MED1_MF1
JOB_DEPENDENCIES=$macs2_callpeak_2_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.A549-DEX_MED1_MF1.383c49ae77e6c96bb214f331f51117dc.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.A549-DEX_MED1_MF1.383c49ae77e6c96bb214f331f51117dc.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/A549-DEX_MED1_MF1/A549-DEX_MED1_MF1 && \
annotatePeaks.pl \
  peak_call/A549-DEX_MED1_MF1/A549-DEX_MED1_MF1_peaks.narrowPeak \
  mm10 \
  -gsize mm10 \
  -cons -CpG \
  -go annotation/A549-DEX_MED1_MF1/A549-DEX_MED1_MF1 \
  -genomeOntology annotation/A549-DEX_MED1_MF1/A549-DEX_MED1_MF1 \
  > annotation/A549-DEX_MED1_MF1/A549-DEX_MED1_MF1.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/A549-DEX_MED1_MF1/A549-DEX_MED1_MF1.annotated.csv",
  "annotation/A549-DEX_MED1_MF1/A549-DEX_MED1_MF1",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.A549-DEX_MED1_MF1.383c49ae77e6c96bb214f331f51117dc.mugqic.done
)
homer_annotate_peaks_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_3_JOB_ID: homer_annotate_peaks.A549-DEX_SMC1A_MF1
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.A549-DEX_SMC1A_MF1
JOB_DEPENDENCIES=$macs2_callpeak_3_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.A549-DEX_SMC1A_MF1.77951c7de21f0a345992bdaa890b66a0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.A549-DEX_SMC1A_MF1.77951c7de21f0a345992bdaa890b66a0.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/A549-DEX_SMC1A_MF1/A549-DEX_SMC1A_MF1 && \
annotatePeaks.pl \
  peak_call/A549-DEX_SMC1A_MF1/A549-DEX_SMC1A_MF1_peaks.narrowPeak \
  mm10 \
  -gsize mm10 \
  -cons -CpG \
  -go annotation/A549-DEX_SMC1A_MF1/A549-DEX_SMC1A_MF1 \
  -genomeOntology annotation/A549-DEX_SMC1A_MF1/A549-DEX_SMC1A_MF1 \
  > annotation/A549-DEX_SMC1A_MF1/A549-DEX_SMC1A_MF1.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/A549-DEX_SMC1A_MF1/A549-DEX_SMC1A_MF1.annotated.csv",
  "annotation/A549-DEX_SMC1A_MF1/A549-DEX_SMC1A_MF1",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.A549-DEX_SMC1A_MF1.77951c7de21f0a345992bdaa890b66a0.mugqic.done
)
homer_annotate_peaks_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_4_JOB_ID: homer_annotate_peaks.A549-DEX_CDK9_MF1
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.A549-DEX_CDK9_MF1
JOB_DEPENDENCIES=$macs2_callpeak_4_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.A549-DEX_CDK9_MF1.6e768cfe3ed709c98b7e5fb3125c0621.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.A549-DEX_CDK9_MF1.6e768cfe3ed709c98b7e5fb3125c0621.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/A549-DEX_CDK9_MF1/A549-DEX_CDK9_MF1 && \
annotatePeaks.pl \
  peak_call/A549-DEX_CDK9_MF1/A549-DEX_CDK9_MF1_peaks.narrowPeak \
  mm10 \
  -gsize mm10 \
  -cons -CpG \
  -go annotation/A549-DEX_CDK9_MF1/A549-DEX_CDK9_MF1 \
  -genomeOntology annotation/A549-DEX_CDK9_MF1/A549-DEX_CDK9_MF1 \
  > annotation/A549-DEX_CDK9_MF1/A549-DEX_CDK9_MF1.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/A549-DEX_CDK9_MF1/A549-DEX_CDK9_MF1.annotated.csv",
  "annotation/A549-DEX_CDK9_MF1/A549-DEX_CDK9_MF1",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.A549-DEX_CDK9_MF1.6e768cfe3ed709c98b7e5fb3125c0621.mugqic.done
)
homer_annotate_peaks_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_5_JOB_ID: homer_annotate_peaks.A549-DEX_BRD4_MF1
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.A549-DEX_BRD4_MF1
JOB_DEPENDENCIES=$macs2_callpeak_5_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.A549-DEX_BRD4_MF1.204e2fe72906c6fa3a70a511d9039c8c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.A549-DEX_BRD4_MF1.204e2fe72906c6fa3a70a511d9039c8c.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/A549-DEX_BRD4_MF1/A549-DEX_BRD4_MF1 && \
annotatePeaks.pl \
  peak_call/A549-DEX_BRD4_MF1/A549-DEX_BRD4_MF1_peaks.narrowPeak \
  mm10 \
  -gsize mm10 \
  -cons -CpG \
  -go annotation/A549-DEX_BRD4_MF1/A549-DEX_BRD4_MF1 \
  -genomeOntology annotation/A549-DEX_BRD4_MF1/A549-DEX_BRD4_MF1 \
  > annotation/A549-DEX_BRD4_MF1/A549-DEX_BRD4_MF1.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/A549-DEX_BRD4_MF1/A549-DEX_BRD4_MF1.annotated.csv",
  "annotation/A549-DEX_BRD4_MF1/A549-DEX_BRD4_MF1",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.A549-DEX_BRD4_MF1.204e2fe72906c6fa3a70a511d9039c8c.mugqic.done
)
homer_annotate_peaks_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_6_JOB_ID: homer_annotate_peaks_report
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks_report
JOB_DEPENDENCIES=$homer_annotate_peaks_1_JOB_ID:$homer_annotate_peaks_2_JOB_ID:$homer_annotate_peaks_3_JOB_ID:$homer_annotate_peaks_4_JOB_ID:$homer_annotate_peaks_5_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks_report.e63d7dd43c5b55894bf2c7883059bdfc.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks_report.e63d7dd43c5b55894bf2c7883059bdfc.mugqic.done'
mkdir -p report/annotation/ && \
cp /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/ChipSeq.homer_annotate_peaks.md report/ && \
for contrast in A549-DEX_NIPBL_MF1 A549-DEX_MED1_MF1 A549-DEX_SMC1A_MF1 A549-DEX_CDK9_MF1 A549-DEX_BRD4_MF1
do
  rsync -avP annotation/$contrast report/annotation/ && \
  echo -e "* [Gene Annotations for Design $contrast](annotation/$contrast/${contrast}.annotated.csv)
* [HOMER Gene Ontology Annotations for Design $contrast](annotation/$contrast/$contrast/geneOntology.html)
* [HOMER Genome Ontology Annotations for Design $contrast](annotation/$contrast/$contrast/GenomeOntology.html)" \
  >> report/ChipSeq.homer_annotate_peaks.md
done
homer_annotate_peaks_report.e63d7dd43c5b55894bf2c7883059bdfc.mugqic.done
)
homer_annotate_peaks_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: homer_find_motifs_genome
#-------------------------------------------------------------------------------
STEP=homer_find_motifs_genome
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_1_JOB_ID: homer_find_motifs_genome.A549-DEX_NIPBL_MF1
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.A549-DEX_NIPBL_MF1
JOB_DEPENDENCIES=$macs2_callpeak_1_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.A549-DEX_NIPBL_MF1.cf7a8a5b6a84198d23bc947dff0fb275.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.A549-DEX_NIPBL_MF1.cf7a8a5b6a84198d23bc947dff0fb275.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/A549-DEX_NIPBL_MF1/A549-DEX_NIPBL_MF1 && \
findMotifsGenome.pl \
  peak_call/A549-DEX_NIPBL_MF1/A549-DEX_NIPBL_MF1_peaks.narrowPeak \
  mm10 \
  annotation/A549-DEX_NIPBL_MF1/A549-DEX_NIPBL_MF1 \
  -preparsedDir annotation/A549-DEX_NIPBL_MF1/A549-DEX_NIPBL_MF1/preparsed \
  -p 4
homer_find_motifs_genome.A549-DEX_NIPBL_MF1.cf7a8a5b6a84198d23bc947dff0fb275.mugqic.done
)
homer_find_motifs_genome_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_2_JOB_ID: homer_find_motifs_genome.A549-DEX_MED1_MF1
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.A549-DEX_MED1_MF1
JOB_DEPENDENCIES=$macs2_callpeak_2_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.A549-DEX_MED1_MF1.46597b9f96f2d68b0af2b7e0bdd31db7.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.A549-DEX_MED1_MF1.46597b9f96f2d68b0af2b7e0bdd31db7.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/A549-DEX_MED1_MF1/A549-DEX_MED1_MF1 && \
findMotifsGenome.pl \
  peak_call/A549-DEX_MED1_MF1/A549-DEX_MED1_MF1_peaks.narrowPeak \
  mm10 \
  annotation/A549-DEX_MED1_MF1/A549-DEX_MED1_MF1 \
  -preparsedDir annotation/A549-DEX_MED1_MF1/A549-DEX_MED1_MF1/preparsed \
  -p 4
homer_find_motifs_genome.A549-DEX_MED1_MF1.46597b9f96f2d68b0af2b7e0bdd31db7.mugqic.done
)
homer_find_motifs_genome_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_3_JOB_ID: homer_find_motifs_genome.A549-DEX_SMC1A_MF1
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.A549-DEX_SMC1A_MF1
JOB_DEPENDENCIES=$macs2_callpeak_3_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.A549-DEX_SMC1A_MF1.11e86d9efc9852af2826f9a0e8756685.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.A549-DEX_SMC1A_MF1.11e86d9efc9852af2826f9a0e8756685.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/A549-DEX_SMC1A_MF1/A549-DEX_SMC1A_MF1 && \
findMotifsGenome.pl \
  peak_call/A549-DEX_SMC1A_MF1/A549-DEX_SMC1A_MF1_peaks.narrowPeak \
  mm10 \
  annotation/A549-DEX_SMC1A_MF1/A549-DEX_SMC1A_MF1 \
  -preparsedDir annotation/A549-DEX_SMC1A_MF1/A549-DEX_SMC1A_MF1/preparsed \
  -p 4
homer_find_motifs_genome.A549-DEX_SMC1A_MF1.11e86d9efc9852af2826f9a0e8756685.mugqic.done
)
homer_find_motifs_genome_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_4_JOB_ID: homer_find_motifs_genome.A549-DEX_CDK9_MF1
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.A549-DEX_CDK9_MF1
JOB_DEPENDENCIES=$macs2_callpeak_4_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.A549-DEX_CDK9_MF1.6d43c3bc40c68d49b874af94904ee8b4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.A549-DEX_CDK9_MF1.6d43c3bc40c68d49b874af94904ee8b4.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/A549-DEX_CDK9_MF1/A549-DEX_CDK9_MF1 && \
findMotifsGenome.pl \
  peak_call/A549-DEX_CDK9_MF1/A549-DEX_CDK9_MF1_peaks.narrowPeak \
  mm10 \
  annotation/A549-DEX_CDK9_MF1/A549-DEX_CDK9_MF1 \
  -preparsedDir annotation/A549-DEX_CDK9_MF1/A549-DEX_CDK9_MF1/preparsed \
  -p 4
homer_find_motifs_genome.A549-DEX_CDK9_MF1.6d43c3bc40c68d49b874af94904ee8b4.mugqic.done
)
homer_find_motifs_genome_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_5_JOB_ID: homer_find_motifs_genome.A549-DEX_BRD4_MF1
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.A549-DEX_BRD4_MF1
JOB_DEPENDENCIES=$macs2_callpeak_5_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.A549-DEX_BRD4_MF1.b53175e152eab2718834ff143437625d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.A549-DEX_BRD4_MF1.b53175e152eab2718834ff143437625d.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/A549-DEX_BRD4_MF1/A549-DEX_BRD4_MF1 && \
findMotifsGenome.pl \
  peak_call/A549-DEX_BRD4_MF1/A549-DEX_BRD4_MF1_peaks.narrowPeak \
  mm10 \
  annotation/A549-DEX_BRD4_MF1/A549-DEX_BRD4_MF1 \
  -preparsedDir annotation/A549-DEX_BRD4_MF1/A549-DEX_BRD4_MF1/preparsed \
  -p 4
homer_find_motifs_genome.A549-DEX_BRD4_MF1.b53175e152eab2718834ff143437625d.mugqic.done
)
homer_find_motifs_genome_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_6_JOB_ID: homer_find_motifs_genome_report
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome_report
JOB_DEPENDENCIES=$homer_find_motifs_genome_1_JOB_ID:$homer_find_motifs_genome_2_JOB_ID:$homer_find_motifs_genome_3_JOB_ID:$homer_find_motifs_genome_4_JOB_ID:$homer_find_motifs_genome_5_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome_report.d146ab67505164435daa8f88068769c6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome_report.d146ab67505164435daa8f88068769c6.mugqic.done'
mkdir -p report/annotation/ && \
cp /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/ChipSeq.homer_find_motifs_genome.md report/ && \
for contrast in A549-DEX_NIPBL_MF1 A549-DEX_MED1_MF1 A549-DEX_SMC1A_MF1 A549-DEX_CDK9_MF1 A549-DEX_BRD4_MF1
do
  rsync -avP annotation/$contrast report/annotation/ && \
  echo -e "* [HOMER _De Novo_ Motif Results for Design $contrast](annotation/$contrast/$contrast/homerResults.html)
* [HOMER Known Motif Results for Design $contrast](annotation/$contrast/$contrast/knownResults.html)" \
  >> report/ChipSeq.homer_find_motifs_genome.md
done
homer_find_motifs_genome_report.d146ab67505164435daa8f88068769c6.mugqic.done
)
homer_find_motifs_genome_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: annotation_graphs
#-------------------------------------------------------------------------------
STEP=annotation_graphs
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: annotation_graphs_1_JOB_ID: annotation_graphs
#-------------------------------------------------------------------------------
JOB_NAME=annotation_graphs
JOB_DEPENDENCIES=$homer_annotate_peaks_1_JOB_ID:$homer_annotate_peaks_2_JOB_ID:$homer_annotate_peaks_3_JOB_ID:$homer_annotate_peaks_4_JOB_ID:$homer_annotate_peaks_5_JOB_ID
JOB_DONE=job_output/annotation_graphs/annotation_graphs.e606333b9aa2460216ff89706860a38b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'annotation_graphs.e606333b9aa2460216ff89706860a38b.mugqic.done'
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
for contrast in A549-DEX_NIPBL_MF1 A549-DEX_MED1_MF1 A549-DEX_SMC1A_MF1 A549-DEX_CDK9_MF1 A549-DEX_BRD4_MF1
do
  cp --parents graphs/${contrast}_Misc_Graphs.ps report/
  convert -rotate 90 graphs/${contrast}_Misc_Graphs.ps report/graphs/${contrast}_Misc_Graphs.png
  echo -e "----

![Annotation Statistics for Design $contrast ([download high-res image](graphs/${contrast}_Misc_Graphs.ps))](graphs/${contrast}_Misc_Graphs.png)
" \
  >> report/ChipSeq.annotation_graphs.md
done
annotation_graphs.e606333b9aa2460216ff89706860a38b.mugqic.done
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
wget "http://mugqic.hpc.mcgill.ca/cgi-bin/pipeline.cgi?hostname=lg-1r17-n03&ip=10.241.129.13&pipeline=ChipSeq&steps=picard_sam_to_fastq,trimmomatic,merge_trimmomatic_stats,bwa_mem_picard_sort_sam,samtools_view_filter,picard_merge_sam_files,picard_mark_duplicates,metrics,homer_make_tag_directory,qc_metrics,homer_make_ucsc_file,macs2_callpeak,homer_annotate_peaks,homer_find_motifs_genome,annotation_graphs&samples=6" --quiet --output-document=/dev/null

