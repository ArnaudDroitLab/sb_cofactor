#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

#-------------------------------------------------------------------------------
# ChipSeq PBSScheduler Job Submission Bash script
# Version: 2.2.0
# Created on: 2017-05-15T18:31:37
# Steps:
#   picard_sam_to_fastq: 0 job... skipping
#   trimmomatic: 5 jobs
#   merge_trimmomatic_stats: 1 job
#   bwa_mem_picard_sort_sam: 6 jobs
#   samtools_view_filter: 6 jobs
#   picard_merge_sam_files: 10 jobs
#   picard_mark_duplicates: 11 jobs
#   metrics: 2 jobs
#   homer_make_tag_directory: 10 jobs
#   qc_metrics: 1 job
#   homer_make_ucsc_file: 11 jobs
#   macs2_callpeak: 17 jobs
#   homer_annotate_peaks: 17 jobs
#   homer_find_motifs_genome: 17 jobs
#   annotation_graphs: 1 job
#   TOTAL: 115 jobs
#-------------------------------------------------------------------------------

OUTPUT_DIR=/gs/scratch/efournier/CofactorHR/A549/output/chip-pipeline
JOB_OUTPUT_DIR=$OUTPUT_DIR/job_output
TIMESTAMP=`date +%FT%H.%M.%S`
JOB_LIST=$JOB_OUTPUT_DIR/ChipSeq_job_list_$TIMESTAMP
mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR


#-------------------------------------------------------------------------------
# STEP: trimmomatic
#-------------------------------------------------------------------------------
STEP=trimmomatic
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: trimmomatic_1_JOB_ID: trimmomatic.HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam.ee945075f2fc4240080f95b7c450e687.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam.ee945075f2fc4240080f95b7c450e687.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/A549_DEX_SMC1A_rep1 && \
`cat > trim/A549_DEX_SMC1A_rep1/HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.single.fastq.gz \
  trim/A549_DEX_SMC1A_rep1/HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/A549_DEX_SMC1A_rep1/HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/A549_DEX_SMC1A_rep1/HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam.trim.log
trimmomatic.HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam.ee945075f2fc4240080f95b7c450e687.mugqic.done
)
trimmomatic_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 | grep "[0-9]")
echo "$trimmomatic_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_2_JOB_ID: trimmomatic.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.8d190831766d68969af4f2071bad05c8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.8d190831766d68969af4f2071bad05c8.mugqic.done'
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
  /gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.single.fastq.gz \
  trim/A549_CTRL_SMC1A_rep2/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/A549_CTRL_SMC1A_rep2/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/A549_CTRL_SMC1A_rep2/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.trim.log
trimmomatic.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.8d190831766d68969af4f2071bad05c8.mugqic.done
)
trimmomatic_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 | grep "[0-9]")
echo "$trimmomatic_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_3_JOB_ID: trimmomatic.HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam.de852473ebeff98a7e86542d31b4f758.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam.de852473ebeff98a7e86542d31b4f758.mugqic.done'
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
  /gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.single.fastq.gz \
  trim/A549_CTRL_WCE_rep2/HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/A549_CTRL_WCE_rep2/HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/A549_CTRL_WCE_rep2/HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam.trim.log
trimmomatic.HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam.de852473ebeff98a7e86542d31b4f758.mugqic.done
)
trimmomatic_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 | grep "[0-9]")
echo "$trimmomatic_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_4_JOB_ID: trimmomatic.HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.f8bb789a922339a4c8ea2e217327cefe.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.f8bb789a922339a4c8ea2e217327cefe.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/trimmomatic/0.35 && \
mkdir -p trim/A549_CTRL_NIPBL_rep2 && \
`cat > trim/A549_CTRL_NIPBL_rep2/HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 6 \
  -phred33 \
  /gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.single.fastq.gz \
  trim/A549_CTRL_NIPBL_rep2/HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/A549_CTRL_NIPBL_rep2/HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/A549_CTRL_NIPBL_rep2/HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.trim.log
trimmomatic.HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.f8bb789a922339a4c8ea2e217327cefe.mugqic.done
)
trimmomatic_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 | grep "[0-9]")
echo "$trimmomatic_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: trimmomatic_5_JOB_ID: trimmomatic.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.b40400c6033420347750cdea293bdb6a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.b40400c6033420347750cdea293bdb6a.mugqic.done'
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
  /gs/scratch/efournier/CofactorHR/A549/raw/chip-seq/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.single.fastq.gz \
  trim/A549_CTRL_MED1_rep2/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.trim.single.fastq.gz \
  ILLUMINACLIP:trim/A549_CTRL_MED1_rep2/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/A549_CTRL_MED1_rep2/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.trim.log
trimmomatic.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.b40400c6033420347750cdea293bdb6a.mugqic.done
)
trimmomatic_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=6 | grep "[0-9]")
echo "$trimmomatic_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: merge_trimmomatic_stats
#-------------------------------------------------------------------------------
STEP=merge_trimmomatic_stats
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: merge_trimmomatic_stats_1_JOB_ID: merge_trimmomatic_stats
#-------------------------------------------------------------------------------
JOB_NAME=merge_trimmomatic_stats
JOB_DEPENDENCIES=$trimmomatic_1_JOB_ID:$trimmomatic_2_JOB_ID:$trimmomatic_3_JOB_ID:$trimmomatic_4_JOB_ID:$trimmomatic_5_JOB_ID
JOB_DONE=job_output/merge_trimmomatic_stats/merge_trimmomatic_stats.3668826523ddcec8fcc526bd82ae80f3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'merge_trimmomatic_stats.3668826523ddcec8fcc526bd82ae80f3.mugqic.done'
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
grep ^Input trim/A549_DEX_SMC1A_rep1/HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/A549_DEX_SMC1A_rep1	HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/A549_CTRL_BRD4_rep1/HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/A549_CTRL_BRD4_rep1	HI.1997.006.Index_2.A549_CTRL_BRD4_MF_CHIP5_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/A549_CTRL_SMC1A_rep2/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/A549_CTRL_SMC1A_rep2	HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/A549_CTRL_WCE_rep2/HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/A549_CTRL_WCE_rep2	HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/A549_CTRL_NIPBL_rep2/HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/A549_CTRL_NIPBL_rep2	HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/A549_CTRL_MED1_rep2/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/A549_CTRL_MED1_rep2	HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam	\1	\2/' | \
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
merge_trimmomatic_stats.3668826523ddcec8fcc526bd82ae80f3.mugqic.done
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
# JOB: bwa_mem_picard_sort_sam_1_JOB_ID: bwa_mem_picard_sort_sam.HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_1_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam.9c4d21faeb4aeac669c7ac4b11638748.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam.9c4d21faeb4aeac669c7ac4b11638748.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/A549_DEX_SMC1A_rep1/HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam	SM:A549_DEX_SMC1A_rep1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/bwa_index/Homo_sapiens.hg19.fa \
  trim/A549_DEX_SMC1A_rep1/HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/A549_DEX_SMC1A_rep1/HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam/HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam.9c4d21faeb4aeac669c7ac4b11638748.mugqic.done
)
bwa_mem_picard_sort_sam_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_2_JOB_ID: bwa_mem_picard_sort_sam.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_2_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.236863e7fae9fcb4335f7a92f239409b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.236863e7fae9fcb4335f7a92f239409b.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/A549_CTRL_SMC1A_rep2/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam	SM:A549_CTRL_SMC1A_rep2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/bwa_index/Homo_sapiens.hg19.fa \
  trim/A549_CTRL_SMC1A_rep2/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/A549_CTRL_SMC1A_rep2/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.236863e7fae9fcb4335f7a92f239409b.mugqic.done
)
bwa_mem_picard_sort_sam_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_3_JOB_ID: bwa_mem_picard_sort_sam.HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_3_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam.8d887f912ec79c183ea51e5c0f7d74fb.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam.8d887f912ec79c183ea51e5c0f7d74fb.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/A549_CTRL_WCE_rep2/HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam	SM:A549_CTRL_WCE_rep2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/bwa_index/Homo_sapiens.hg19.fa \
  trim/A549_CTRL_WCE_rep2/HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/A549_CTRL_WCE_rep2/HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam/HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam.8d887f912ec79c183ea51e5c0f7d74fb.mugqic.done
)
bwa_mem_picard_sort_sam_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_4_JOB_ID: bwa_mem_picard_sort_sam.HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_4_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.1624faebae0d04775e0f2fe9386b34bf.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.1624faebae0d04775e0f2fe9386b34bf.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/A549_CTRL_NIPBL_rep2/HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam	SM:A549_CTRL_NIPBL_rep2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/bwa_index/Homo_sapiens.hg19.fa \
  trim/A549_CTRL_NIPBL_rep2/HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/A549_CTRL_NIPBL_rep2/HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam/HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.1624faebae0d04775e0f2fe9386b34bf.mugqic.done
)
bwa_mem_picard_sort_sam_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_5_JOB_ID: bwa_mem_picard_sort_sam.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam
JOB_DEPENDENCIES=$trimmomatic_5_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.a3dbc83af100541f628528ad4bec6dda.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.a3dbc83af100541f628528ad4bec6dda.mugqic.done'
module load mugqic/bwa/0.7.12 mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
mkdir -p alignment/A549_CTRL_MED1_rep2/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam && \
bwa mem  \
  -M -t 7 \
  -R '@RG	ID:HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam	SM:A549_CTRL_MED1_rep2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.hg19/genome/bwa_index/Homo_sapiens.hg19.fa \
  trim/A549_CTRL_MED1_rep2/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.trim.single.fastq.gz | \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx15G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=/dev/stdin \
  OUTPUT=alignment/A549_CTRL_MED1_rep2/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.sorted.bam \
  SORT_ORDER=coordinate \
  MAX_RECORDS_IN_RAM=3750000
bwa_mem_picard_sort_sam.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.a3dbc83af100541f628528ad4bec6dda.mugqic.done
)
bwa_mem_picard_sort_sam_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=12 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_6_JOB_ID: bwa_mem_picard_sort_sam_report
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam_report
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_1_JOB_ID:$bwa_mem_picard_sort_sam_2_JOB_ID:$bwa_mem_picard_sort_sam_3_JOB_ID:$bwa_mem_picard_sort_sam_4_JOB_ID:$bwa_mem_picard_sort_sam_5_JOB_ID
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
bwa_mem_picard_sort_sam_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$bwa_mem_picard_sort_sam_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: samtools_view_filter
#-------------------------------------------------------------------------------
STEP=samtools_view_filter
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_1_JOB_ID: samtools_view_filter.HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_1_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam.c8edb9e0b4f589a99b8c47fb4ebc63ee.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam.c8edb9e0b4f589a99b8c47fb4ebc63ee.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/A549_DEX_SMC1A_rep1/HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam/HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam.sorted.bam \
  > alignment/A549_DEX_SMC1A_rep1/HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam/HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam.c8edb9e0b4f589a99b8c47fb4ebc63ee.mugqic.done
)
samtools_view_filter_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_2_JOB_ID: samtools_view_filter.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_2_JOB_ID
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
samtools_view_filter_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_3_JOB_ID: samtools_view_filter.HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.2449.006.Index_19.A549_ctrl_WCE_MF6_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_3_JOB_ID
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
samtools_view_filter_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_4_JOB_ID: samtools_view_filter.HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_4_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.c2326c9f209ad2732f12b57240f3ffca.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.c2326c9f209ad2732f12b57240f3ffca.mugqic.done'
module load mugqic/samtools/1.3 && \
samtools view -b -F4 -q 20 \
  alignment/A549_CTRL_NIPBL_rep2/HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam/HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.sorted.bam \
  > alignment/A549_CTRL_NIPBL_rep2/HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam/HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.sorted.filtered.bam
samtools_view_filter.HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.c2326c9f209ad2732f12b57240f3ffca.mugqic.done
)
samtools_view_filter_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_5_JOB_ID: samtools_view_filter.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_5_JOB_ID
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
samtools_view_filter_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_6_JOB_ID: samtools_view_filter_report
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter_report
JOB_DEPENDENCIES=$samtools_view_filter_1_JOB_ID:$samtools_view_filter_2_JOB_ID:$samtools_view_filter_3_JOB_ID:$samtools_view_filter_4_JOB_ID:$samtools_view_filter_5_JOB_ID
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
samtools_view_filter_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$samtools_view_filter_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: picard_merge_sam_files
#-------------------------------------------------------------------------------
STEP=picard_merge_sam_files
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_1_JOB_ID: symlink_readset_sample_bam.A549_DEX_BRD4_rep1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.A549_DEX_BRD4_rep1
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.A549_DEX_BRD4_rep1.62fe20cb023c57aacf7b420810c59c8b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.A549_DEX_BRD4_rep1.62fe20cb023c57aacf7b420810c59c8b.mugqic.done'
mkdir -p alignment/A549_DEX_BRD4_rep1 && \
ln -s -f HI.1613.003.Index_19.A549_Dex_SMC1_MF_ChIP1_3_rep1.bam/HI.1613.003.Index_19.A549_Dex_SMC1_MF_ChIP1_3_rep1.bam.sorted.filtered.bam alignment/A549_DEX_BRD4_rep1/A549_DEX_BRD4_rep1.merged.bam
symlink_readset_sample_bam.A549_DEX_BRD4_rep1.62fe20cb023c57aacf7b420810c59c8b.mugqic.done
)
picard_merge_sam_files_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$picard_merge_sam_files_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_2_JOB_ID: symlink_readset_sample_bam.A549_CTRL_SMC1A_rep1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.A549_CTRL_SMC1A_rep1
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.A549_CTRL_SMC1A_rep1.301863dc2fbaabe3fa9463fd1a0ab344.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.A549_CTRL_SMC1A_rep1.301863dc2fbaabe3fa9463fd1a0ab344.mugqic.done'
mkdir -p alignment/A549_CTRL_SMC1A_rep1 && \
ln -s -f HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam/HI.1613.004.Index_14.A549_ctrl_SMC1_MF_ChIP1_3_rep1.bam.sorted.filtered.bam alignment/A549_CTRL_SMC1A_rep1/A549_CTRL_SMC1A_rep1.merged.bam
symlink_readset_sample_bam.A549_CTRL_SMC1A_rep1.301863dc2fbaabe3fa9463fd1a0ab344.mugqic.done
)
picard_merge_sam_files_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$picard_merge_sam_files_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_3_JOB_ID: symlink_readset_sample_bam.A549_CTRL_WCE_rep1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.A549_CTRL_WCE_rep1
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.A549_CTRL_WCE_rep1.55133cc6a5e4bf8d3079533ed1955e28.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.A549_CTRL_WCE_rep1.55133cc6a5e4bf8d3079533ed1955e28.mugqic.done'
mkdir -p alignment/A549_CTRL_WCE_rep1 && \
ln -s -f HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam/HI.1613.004.Index_15.A549_ctrl_WCE_MF_ChIP1_8_rep1.bam.sorted.filtered.bam alignment/A549_CTRL_WCE_rep1/A549_CTRL_WCE_rep1.merged.bam
symlink_readset_sample_bam.A549_CTRL_WCE_rep1.55133cc6a5e4bf8d3079533ed1955e28.mugqic.done
)
picard_merge_sam_files_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$picard_merge_sam_files_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_4_JOB_ID: symlink_readset_sample_bam.A549_CTRL_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.A549_CTRL_MED1_rep1
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.A549_CTRL_MED1_rep1.93b78a91a30a2b05a8ec84948942bd09.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.A549_CTRL_MED1_rep1.93b78a91a30a2b05a8ec84948942bd09.mugqic.done'
mkdir -p alignment/A549_CTRL_MED1_rep1 && \
ln -s -f HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam/HI.1613.004.Index_18.A549_ctrl_MED1_MF_ChIP1_2_rep1.bam.sorted.filtered.bam alignment/A549_CTRL_MED1_rep1/A549_CTRL_MED1_rep1.merged.bam
symlink_readset_sample_bam.A549_CTRL_MED1_rep1.93b78a91a30a2b05a8ec84948942bd09.mugqic.done
)
picard_merge_sam_files_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$picard_merge_sam_files_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_5_JOB_ID: symlink_readset_sample_bam.A549_CTRL_NIPBL_rep1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.A549_CTRL_NIPBL_rep1
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.A549_CTRL_NIPBL_rep1.3dbc49899185a184af13b1997d2223db.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.A549_CTRL_NIPBL_rep1.3dbc49899185a184af13b1997d2223db.mugqic.done'
mkdir -p alignment/A549_CTRL_NIPBL_rep1 && \
ln -s -f HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam/HI.1613.004.Index_7.A549_ctrl_NIPBL_MF_ChIP1_1_rep1.bam.sorted.filtered.bam alignment/A549_CTRL_NIPBL_rep1/A549_CTRL_NIPBL_rep1.merged.bam
symlink_readset_sample_bam.A549_CTRL_NIPBL_rep1.3dbc49899185a184af13b1997d2223db.mugqic.done
)
picard_merge_sam_files_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$picard_merge_sam_files_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_6_JOB_ID: symlink_readset_sample_bam.A549_DEX_SMC1A_rep1
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.A549_DEX_SMC1A_rep1
JOB_DEPENDENCIES=$samtools_view_filter_1_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.A549_DEX_SMC1A_rep1.681aeed27f676634b4d462f05b1c4a60.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.A549_DEX_SMC1A_rep1.681aeed27f676634b4d462f05b1c4a60.mugqic.done'
mkdir -p alignment/A549_DEX_SMC1A_rep1 && \
ln -s -f HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam/HI.1997.006.Index_13.A549_DEX_BRD4_MF_CHIP5_rep1.bam.sorted.filtered.bam alignment/A549_DEX_SMC1A_rep1/A549_DEX_SMC1A_rep1.merged.bam
symlink_readset_sample_bam.A549_DEX_SMC1A_rep1.681aeed27f676634b4d462f05b1c4a60.mugqic.done
)
picard_merge_sam_files_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_7_JOB_ID: symlink_readset_sample_bam.A549_CTRL_SMC1A_rep2
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.A549_CTRL_SMC1A_rep2
JOB_DEPENDENCIES=$samtools_view_filter_2_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.A549_CTRL_SMC1A_rep2.e14f5c2700d1c3b7065db2d93c834ed5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.A549_CTRL_SMC1A_rep2.e14f5c2700d1c3b7065db2d93c834ed5.mugqic.done'
mkdir -p alignment/A549_CTRL_SMC1A_rep2 && \
ln -s -f HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam/HI.2449.006.Index_12.A549_ctrl_SMC1A_MF6_rep1.bam.sorted.filtered.bam alignment/A549_CTRL_SMC1A_rep2/A549_CTRL_SMC1A_rep2.merged.bam
symlink_readset_sample_bam.A549_CTRL_SMC1A_rep2.e14f5c2700d1c3b7065db2d93c834ed5.mugqic.done
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
JOB_DEPENDENCIES=$samtools_view_filter_3_JOB_ID
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
# JOB: picard_merge_sam_files_9_JOB_ID: symlink_readset_sample_bam.A549_CTRL_NIPBL_rep2
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.A549_CTRL_NIPBL_rep2
JOB_DEPENDENCIES=$samtools_view_filter_4_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.A549_CTRL_NIPBL_rep2.0629d48f508e69caf2bec0ebfe1bdc1a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.A549_CTRL_NIPBL_rep2.0629d48f508e69caf2bec0ebfe1bdc1a.mugqic.done'
mkdir -p alignment/A549_CTRL_NIPBL_rep2 && \
ln -s -f HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam/HI.2449.006.Index_4.A549_ctrl_NIPBL_MF6_rep1.bam.sorted.filtered.bam alignment/A549_CTRL_NIPBL_rep2/A549_CTRL_NIPBL_rep2.merged.bam
symlink_readset_sample_bam.A549_CTRL_NIPBL_rep2.0629d48f508e69caf2bec0ebfe1bdc1a.mugqic.done
)
picard_merge_sam_files_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_10_JOB_ID: symlink_readset_sample_bam.A549_CTRL_MED1_rep2
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.A549_CTRL_MED1_rep2
JOB_DEPENDENCIES=$samtools_view_filter_5_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.A549_CTRL_MED1_rep2.0dd9318287ce9f4696134ffa8bca3173.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.A549_CTRL_MED1_rep2.0dd9318287ce9f4696134ffa8bca3173.mugqic.done'
mkdir -p alignment/A549_CTRL_MED1_rep2 && \
ln -s -f HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam/HI.2449.006.Index_5.A549_ctrl_MED1_MF6_rep1.bam.sorted.filtered.bam alignment/A549_CTRL_MED1_rep2/A549_CTRL_MED1_rep2.merged.bam
symlink_readset_sample_bam.A549_CTRL_MED1_rep2.0dd9318287ce9f4696134ffa8bca3173.mugqic.done
)
picard_merge_sam_files_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_merge_sam_files_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: picard_mark_duplicates
#-------------------------------------------------------------------------------
STEP=picard_mark_duplicates
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_1_JOB_ID: picard_mark_duplicates.A549_DEX_BRD4_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.A549_DEX_BRD4_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_1_JOB_ID
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
picard_mark_duplicates_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_2_JOB_ID: picard_mark_duplicates.A549_CTRL_SMC1A_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.A549_CTRL_SMC1A_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_2_JOB_ID
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
picard_mark_duplicates_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_3_JOB_ID: picard_mark_duplicates.A549_CTRL_WCE_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.A549_CTRL_WCE_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_3_JOB_ID
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
picard_mark_duplicates_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_4_JOB_ID: picard_mark_duplicates.A549_CTRL_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.A549_CTRL_MED1_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_4_JOB_ID
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
picard_mark_duplicates_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_5_JOB_ID: picard_mark_duplicates.A549_CTRL_NIPBL_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.A549_CTRL_NIPBL_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_5_JOB_ID
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
picard_mark_duplicates_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_6_JOB_ID: picard_mark_duplicates.A549_DEX_SMC1A_rep1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.A549_DEX_SMC1A_rep1
JOB_DEPENDENCIES=$picard_merge_sam_files_6_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.A549_DEX_SMC1A_rep1.5407beeece1d57f26c8cca871dd08f80.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.A549_DEX_SMC1A_rep1.5407beeece1d57f26c8cca871dd08f80.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/A549_DEX_SMC1A_rep1/A549_DEX_SMC1A_rep1.merged.bam \
  OUTPUT=alignment/A549_DEX_SMC1A_rep1/A549_DEX_SMC1A_rep1.sorted.dup.bam \
  METRICS_FILE=alignment/A549_DEX_SMC1A_rep1/A549_DEX_SMC1A_rep1.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.A549_DEX_SMC1A_rep1.5407beeece1d57f26c8cca871dd08f80.mugqic.done
)
picard_mark_duplicates_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_7_JOB_ID: picard_mark_duplicates.A549_CTRL_SMC1A_rep2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.A549_CTRL_SMC1A_rep2
JOB_DEPENDENCIES=$picard_merge_sam_files_7_JOB_ID
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
# JOB: picard_mark_duplicates_9_JOB_ID: picard_mark_duplicates.A549_CTRL_NIPBL_rep2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.A549_CTRL_NIPBL_rep2
JOB_DEPENDENCIES=$picard_merge_sam_files_9_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.A549_CTRL_NIPBL_rep2.b9c62bfc7b0940d5e8576775d7f44b25.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.A549_CTRL_NIPBL_rep2.b9c62bfc7b0940d5e8576775d7f44b25.mugqic.done'
module load mugqic/java/openjdk-jdk1.8.0_72 mugqic/picard/1.123 && \
java -Djava.io.tmpdir=/localscratch/ -XX:ParallelGCThreads=1 -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  TMP_DIR=/localscratch/ \
  INPUT=alignment/A549_CTRL_NIPBL_rep2/A549_CTRL_NIPBL_rep2.merged.bam \
  OUTPUT=alignment/A549_CTRL_NIPBL_rep2/A549_CTRL_NIPBL_rep2.sorted.dup.bam \
  METRICS_FILE=alignment/A549_CTRL_NIPBL_rep2/A549_CTRL_NIPBL_rep2.sorted.dup.metrics \
  MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.A549_CTRL_NIPBL_rep2.b9c62bfc7b0940d5e8576775d7f44b25.mugqic.done
)
picard_mark_duplicates_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_10_JOB_ID: picard_mark_duplicates.A549_CTRL_MED1_rep2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.A549_CTRL_MED1_rep2
JOB_DEPENDENCIES=$picard_merge_sam_files_10_JOB_ID
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
picard_mark_duplicates_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=48:00:0 -q metaq -l nodes=1:ppn=2 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_11_JOB_ID: picard_mark_duplicates_report
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates_report
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_2_JOB_ID:$picard_mark_duplicates_3_JOB_ID:$picard_mark_duplicates_4_JOB_ID:$picard_mark_duplicates_5_JOB_ID:$picard_mark_duplicates_6_JOB_ID:$picard_mark_duplicates_7_JOB_ID:$picard_mark_duplicates_8_JOB_ID:$picard_mark_duplicates_9_JOB_ID:$picard_mark_duplicates_10_JOB_ID
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
picard_mark_duplicates_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$picard_mark_duplicates_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: metrics
#-------------------------------------------------------------------------------
STEP=metrics
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: metrics_1_JOB_ID: metrics.flagstat
#-------------------------------------------------------------------------------
JOB_NAME=metrics.flagstat
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_2_JOB_ID:$picard_mark_duplicates_3_JOB_ID:$picard_mark_duplicates_4_JOB_ID:$picard_mark_duplicates_5_JOB_ID:$picard_mark_duplicates_6_JOB_ID:$picard_mark_duplicates_7_JOB_ID:$picard_mark_duplicates_8_JOB_ID:$picard_mark_duplicates_9_JOB_ID:$picard_mark_duplicates_10_JOB_ID
JOB_DONE=job_output/metrics/metrics.flagstat.fa22288a9e6026177aa55876d23cd909.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'metrics.flagstat.fa22288a9e6026177aa55876d23cd909.mugqic.done'
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
  alignment/A549_DEX_SMC1A_rep1/A549_DEX_SMC1A_rep1.sorted.dup.bam \
  > alignment/A549_DEX_SMC1A_rep1/A549_DEX_SMC1A_rep1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/A549_CTRL_BRD4_rep1/A549_CTRL_BRD4_rep1.sorted.dup.bam \
  > alignment/A549_CTRL_BRD4_rep1/A549_CTRL_BRD4_rep1.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/A549_CTRL_SMC1A_rep2/A549_CTRL_SMC1A_rep2.sorted.dup.bam \
  > alignment/A549_CTRL_SMC1A_rep2/A549_CTRL_SMC1A_rep2.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/A549_CTRL_WCE_rep2/A549_CTRL_WCE_rep2.sorted.dup.bam \
  > alignment/A549_CTRL_WCE_rep2/A549_CTRL_WCE_rep2.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/A549_CTRL_NIPBL_rep2/A549_CTRL_NIPBL_rep2.sorted.dup.bam \
  > alignment/A549_CTRL_NIPBL_rep2/A549_CTRL_NIPBL_rep2.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/A549_CTRL_MED1_rep2/A549_CTRL_MED1_rep2.sorted.dup.bam \
  > alignment/A549_CTRL_MED1_rep2/A549_CTRL_MED1_rep2.sorted.dup.bam.flagstat
metrics.flagstat.fa22288a9e6026177aa55876d23cd909.mugqic.done
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
JOB_DONE=job_output/metrics/metrics_report.67cc07ddc444824bc219293d8e97a0c6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'metrics_report.67cc07ddc444824bc219293d8e97a0c6.mugqic.done'
module load mugqic/pandoc/1.15.2 && \
for sample in A549_DEX_CDK9_rep1 A549_DEX_MED1_rep1 A549_DEX_BRD4_rep1 A549_DEX_WCE_rep1 A549_DEX_NIPBL_rep1 A549_CTRL_SMC1A_rep1 A549_CTRL_WCE_rep1 A549_CTRL_CDK9_rep1 A549_CTRL_MED1_rep1 A549_CTRL_NIPBL_rep1 A549_DEX_SMC1A_rep1 A549_CTRL_BRD4_rep1 A549_CTRL_SMC1A_rep2 A549_CTRL_WCE_rep2 A549_CTRL_NIPBL_rep2 A549_CTRL_MED1_rep2
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

metrics_report.67cc07ddc444824bc219293d8e97a0c6.mugqic.done
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
# JOB: homer_make_tag_directory_1_JOB_ID: homer_make_tag_directory.A549_DEX_BRD4_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.A549_DEX_BRD4_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID
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
homer_make_tag_directory_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_2_JOB_ID: homer_make_tag_directory.A549_CTRL_SMC1A_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.A549_CTRL_SMC1A_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_2_JOB_ID
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
homer_make_tag_directory_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_3_JOB_ID: homer_make_tag_directory.A549_CTRL_WCE_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.A549_CTRL_WCE_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_3_JOB_ID
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
homer_make_tag_directory_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_4_JOB_ID: homer_make_tag_directory.A549_CTRL_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.A549_CTRL_MED1_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_4_JOB_ID
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
homer_make_tag_directory_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_5_JOB_ID: homer_make_tag_directory.A549_CTRL_NIPBL_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.A549_CTRL_NIPBL_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_5_JOB_ID
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
homer_make_tag_directory_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_6_JOB_ID: homer_make_tag_directory.A549_DEX_SMC1A_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.A549_DEX_SMC1A_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_6_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.A549_DEX_SMC1A_rep1.b32e4a3e11b54e598e1684cd2d93c2ea.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.A549_DEX_SMC1A_rep1.b32e4a3e11b54e598e1684cd2d93c2ea.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/A549_DEX_SMC1A_rep1 \
  alignment/A549_DEX_SMC1A_rep1/A549_DEX_SMC1A_rep1.sorted.dup.bam \
  -checkGC -genome hg19
homer_make_tag_directory.A549_DEX_SMC1A_rep1.b32e4a3e11b54e598e1684cd2d93c2ea.mugqic.done
)
homer_make_tag_directory_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_7_JOB_ID: homer_make_tag_directory.A549_CTRL_SMC1A_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.A549_CTRL_SMC1A_rep2
JOB_DEPENDENCIES=$picard_mark_duplicates_7_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.A549_CTRL_SMC1A_rep2.b156f9636c98fe7e868ab2e811944176.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.A549_CTRL_SMC1A_rep2.b156f9636c98fe7e868ab2e811944176.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/A549_CTRL_SMC1A_rep2 \
  alignment/A549_CTRL_SMC1A_rep2/A549_CTRL_SMC1A_rep2.sorted.dup.bam \
  -checkGC -genome hg19
homer_make_tag_directory.A549_CTRL_SMC1A_rep2.b156f9636c98fe7e868ab2e811944176.mugqic.done
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
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.A549_CTRL_WCE_rep2.9ae0c01111c6204abbb268573649e57b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.A549_CTRL_WCE_rep2.9ae0c01111c6204abbb268573649e57b.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/A549_CTRL_WCE_rep2 \
  alignment/A549_CTRL_WCE_rep2/A549_CTRL_WCE_rep2.sorted.dup.bam \
  -checkGC -genome hg19
homer_make_tag_directory.A549_CTRL_WCE_rep2.9ae0c01111c6204abbb268573649e57b.mugqic.done
)
homer_make_tag_directory_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_9_JOB_ID: homer_make_tag_directory.A549_CTRL_NIPBL_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.A549_CTRL_NIPBL_rep2
JOB_DEPENDENCIES=$picard_mark_duplicates_9_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.A549_CTRL_NIPBL_rep2.73a72f509e0b3fd66058c8f9101a72a4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.A549_CTRL_NIPBL_rep2.73a72f509e0b3fd66058c8f9101a72a4.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/A549_CTRL_NIPBL_rep2 \
  alignment/A549_CTRL_NIPBL_rep2/A549_CTRL_NIPBL_rep2.sorted.dup.bam \
  -checkGC -genome hg19
homer_make_tag_directory.A549_CTRL_NIPBL_rep2.73a72f509e0b3fd66058c8f9101a72a4.mugqic.done
)
homer_make_tag_directory_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_10_JOB_ID: homer_make_tag_directory.A549_CTRL_MED1_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.A549_CTRL_MED1_rep2
JOB_DEPENDENCIES=$picard_mark_duplicates_10_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.A549_CTRL_MED1_rep2.c9481b4a02810375ec2375452d8f1af2.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.A549_CTRL_MED1_rep2.c9481b4a02810375ec2375452d8f1af2.mugqic.done'
module load mugqic/samtools/1.3 mugqic/homer/4.7 && \
makeTagDirectory \
  tags/A549_CTRL_MED1_rep2 \
  alignment/A549_CTRL_MED1_rep2/A549_CTRL_MED1_rep2.sorted.dup.bam \
  -checkGC -genome hg19
homer_make_tag_directory.A549_CTRL_MED1_rep2.c9481b4a02810375ec2375452d8f1af2.mugqic.done
)
homer_make_tag_directory_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_tag_directory_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: qc_metrics
#-------------------------------------------------------------------------------
STEP=qc_metrics
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: qc_metrics_1_JOB_ID: qc_plots_R
#-------------------------------------------------------------------------------
JOB_NAME=qc_plots_R
JOB_DEPENDENCIES=$homer_make_tag_directory_1_JOB_ID:$homer_make_tag_directory_2_JOB_ID:$homer_make_tag_directory_3_JOB_ID:$homer_make_tag_directory_4_JOB_ID:$homer_make_tag_directory_5_JOB_ID:$homer_make_tag_directory_6_JOB_ID:$homer_make_tag_directory_7_JOB_ID:$homer_make_tag_directory_8_JOB_ID:$homer_make_tag_directory_9_JOB_ID:$homer_make_tag_directory_10_JOB_ID
JOB_DONE=job_output/qc_metrics/qc_plots_R.d797c6e24129c9730fb6776cd7690621.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'qc_plots_R.d797c6e24129c9730fb6776cd7690621.mugqic.done'
module load mugqic/mugqic_tools/2.1.5 mugqic/R_Bioconductor/3.2.3_3.2 && \
mkdir -p graphs && \
Rscript $R_TOOLS/chipSeqGenerateQCMetrics.R \
  ../../raw/chip-seq/design.txt \
  /gs/scratch/efournier/CofactorHR/A549/output/chip-pipeline && \
cp /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/ChipSeq.qc_metrics.md report/ChipSeq.qc_metrics.md && \
for sample in A549_DEX_CDK9_rep1 A549_DEX_MED1_rep1 A549_DEX_BRD4_rep1 A549_DEX_WCE_rep1 A549_DEX_NIPBL_rep1 A549_CTRL_SMC1A_rep1 A549_CTRL_WCE_rep1 A549_CTRL_CDK9_rep1 A549_CTRL_MED1_rep1 A549_CTRL_NIPBL_rep1 A549_DEX_SMC1A_rep1 A549_CTRL_BRD4_rep1 A549_CTRL_SMC1A_rep2 A549_CTRL_WCE_rep2 A549_CTRL_NIPBL_rep2 A549_CTRL_MED1_rep2
do
  cp --parents graphs/${sample}_QC_Metrics.ps report/
  convert -rotate 90 graphs/${sample}_QC_Metrics.ps report/graphs/${sample}_QC_Metrics.png
  echo -e "----

![QC Metrics for Sample $sample ([download high-res image](graphs/${sample}_QC_Metrics.ps))](graphs/${sample}_QC_Metrics.png)
" \
  >> report/ChipSeq.qc_metrics.md
done
qc_plots_R.d797c6e24129c9730fb6776cd7690621.mugqic.done
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
# JOB: homer_make_ucsc_file_1_JOB_ID: homer_make_ucsc_file.A549_DEX_BRD4_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.A549_DEX_BRD4_rep1
JOB_DEPENDENCIES=$homer_make_tag_directory_1_JOB_ID
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
homer_make_ucsc_file_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_2_JOB_ID: homer_make_ucsc_file.A549_CTRL_SMC1A_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.A549_CTRL_SMC1A_rep1
JOB_DEPENDENCIES=$homer_make_tag_directory_2_JOB_ID
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
homer_make_ucsc_file_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_3_JOB_ID: homer_make_ucsc_file.A549_CTRL_WCE_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.A549_CTRL_WCE_rep1
JOB_DEPENDENCIES=$homer_make_tag_directory_3_JOB_ID
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
homer_make_ucsc_file_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_4_JOB_ID: homer_make_ucsc_file.A549_CTRL_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.A549_CTRL_MED1_rep1
JOB_DEPENDENCIES=$homer_make_tag_directory_4_JOB_ID
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
homer_make_ucsc_file_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_5_JOB_ID: homer_make_ucsc_file.A549_CTRL_NIPBL_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.A549_CTRL_NIPBL_rep1
JOB_DEPENDENCIES=$homer_make_tag_directory_5_JOB_ID
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
homer_make_ucsc_file_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_6_JOB_ID: homer_make_ucsc_file.A549_DEX_SMC1A_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.A549_DEX_SMC1A_rep1
JOB_DEPENDENCIES=$homer_make_tag_directory_6_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.A549_DEX_SMC1A_rep1.6c66c25a22b612f13fc5a4a7c514efab.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.A549_DEX_SMC1A_rep1.6c66c25a22b612f13fc5a4a7c514efab.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/A549_DEX_SMC1A_rep1 && \
makeUCSCfile \
  tags/A549_DEX_SMC1A_rep1 | \
gzip -1 -c > tracks/A549_DEX_SMC1A_rep1/A549_DEX_SMC1A_rep1.ucsc.bedGraph.gz
homer_make_ucsc_file.A549_DEX_SMC1A_rep1.6c66c25a22b612f13fc5a4a7c514efab.mugqic.done
)
homer_make_ucsc_file_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_7_JOB_ID: homer_make_ucsc_file.A549_CTRL_SMC1A_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.A549_CTRL_SMC1A_rep2
JOB_DEPENDENCIES=$homer_make_tag_directory_7_JOB_ID
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
# JOB: homer_make_ucsc_file_9_JOB_ID: homer_make_ucsc_file.A549_CTRL_NIPBL_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.A549_CTRL_NIPBL_rep2
JOB_DEPENDENCIES=$homer_make_tag_directory_9_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.A549_CTRL_NIPBL_rep2.ad4198c71841936cbb30a57b9cdcb26c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.A549_CTRL_NIPBL_rep2.ad4198c71841936cbb30a57b9cdcb26c.mugqic.done'
module load mugqic/homer/4.7 && \
mkdir -p tracks/A549_CTRL_NIPBL_rep2 && \
makeUCSCfile \
  tags/A549_CTRL_NIPBL_rep2 | \
gzip -1 -c > tracks/A549_CTRL_NIPBL_rep2/A549_CTRL_NIPBL_rep2.ucsc.bedGraph.gz
homer_make_ucsc_file.A549_CTRL_NIPBL_rep2.ad4198c71841936cbb30a57b9cdcb26c.mugqic.done
)
homer_make_ucsc_file_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_10_JOB_ID: homer_make_ucsc_file.A549_CTRL_MED1_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.A549_CTRL_MED1_rep2
JOB_DEPENDENCIES=$homer_make_tag_directory_10_JOB_ID
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
homer_make_ucsc_file_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_11_JOB_ID: homer_make_ucsc_file_report
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file_report
JOB_DEPENDENCIES=$homer_make_ucsc_file_10_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file_report.c9cac04b3bc71f34cbcde042cf5b419d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file_report.c9cac04b3bc71f34cbcde042cf5b419d.mugqic.done'
mkdir -p report && \
zip -r report/tracks.zip tracks/*/*.ucsc.bedGraph.gz && \
cp /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/ChipSeq.homer_make_ucsc_file.md report/
homer_make_ucsc_file_report.c9cac04b3bc71f34cbcde042cf5b419d.mugqic.done
)
homer_make_ucsc_file_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_make_ucsc_file_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: macs2_callpeak
#-------------------------------------------------------------------------------
STEP=macs2_callpeak
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_1_JOB_ID: macs2_callpeak.A549_CTRL_BRD4
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_CTRL_BRD4
JOB_DEPENDENCIES=$picard_mark_duplicates_3_JOB_ID:$picard_mark_duplicates_8_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_CTRL_BRD4.00063051b17e41fb5400ae9f091c430a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_CTRL_BRD4.00063051b17e41fb5400ae9f091c430a.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/A549_CTRL_BRD4 && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2509729011.2 \
  --treatment \
  alignment/A549_CTRL_BRD4_rep1/A549_CTRL_BRD4_rep1.sorted.dup.bam \
  --control \
  alignment/A549_CTRL_WCE_rep1/A549_CTRL_WCE_rep1.sorted.dup.bam \
  alignment/A549_CTRL_WCE_rep2/A549_CTRL_WCE_rep2.sorted.dup.bam \
  --name peak_call/A549_CTRL_BRD4/A549_CTRL_BRD4 \
  >& peak_call/A549_CTRL_BRD4/A549_CTRL_BRD4.diag.macs.out
macs2_callpeak.A549_CTRL_BRD4.00063051b17e41fb5400ae9f091c430a.mugqic.done
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
JOB_DEPENDENCIES=$picard_mark_duplicates_3_JOB_ID:$picard_mark_duplicates_8_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_CTRL_CDK9.83beb553134d44e3031b25670cfeeb52.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_CTRL_CDK9.83beb553134d44e3031b25670cfeeb52.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/A549_CTRL_CDK9 && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2509729011.2 \
  --treatment \
  alignment/A549_CTRL_CDK9_rep1/A549_CTRL_CDK9_rep1.sorted.dup.bam \
  --control \
  alignment/A549_CTRL_WCE_rep1/A549_CTRL_WCE_rep1.sorted.dup.bam \
  alignment/A549_CTRL_WCE_rep2/A549_CTRL_WCE_rep2.sorted.dup.bam \
  --name peak_call/A549_CTRL_CDK9/A549_CTRL_CDK9 \
  >& peak_call/A549_CTRL_CDK9/A549_CTRL_CDK9.diag.macs.out
macs2_callpeak.A549_CTRL_CDK9.83beb553134d44e3031b25670cfeeb52.mugqic.done
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
JOB_DEPENDENCIES=$picard_mark_duplicates_3_JOB_ID:$picard_mark_duplicates_4_JOB_ID:$picard_mark_duplicates_8_JOB_ID:$picard_mark_duplicates_10_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_CTRL_MED1.0c0b07656529074d3531864287d487fd.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_CTRL_MED1.0c0b07656529074d3531864287d487fd.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/A549_CTRL_MED1 && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2509729011.2 \
  --treatment \
  alignment/A549_CTRL_MED1_rep1/A549_CTRL_MED1_rep1.sorted.dup.bam \
  alignment/A549_CTRL_MED1_rep2/A549_CTRL_MED1_rep2.sorted.dup.bam \
  --control \
  alignment/A549_CTRL_WCE_rep1/A549_CTRL_WCE_rep1.sorted.dup.bam \
  alignment/A549_CTRL_WCE_rep2/A549_CTRL_WCE_rep2.sorted.dup.bam \
  --name peak_call/A549_CTRL_MED1/A549_CTRL_MED1 \
  >& peak_call/A549_CTRL_MED1/A549_CTRL_MED1.diag.macs.out
macs2_callpeak.A549_CTRL_MED1.0c0b07656529074d3531864287d487fd.mugqic.done
)
macs2_callpeak_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_4_JOB_ID: macs2_callpeak.A549_CTRL_NIPBL
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_CTRL_NIPBL
JOB_DEPENDENCIES=$picard_mark_duplicates_3_JOB_ID:$picard_mark_duplicates_5_JOB_ID:$picard_mark_duplicates_8_JOB_ID:$picard_mark_duplicates_9_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_CTRL_NIPBL.32bdbf16e5347684abc0921837f89437.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_CTRL_NIPBL.32bdbf16e5347684abc0921837f89437.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/A549_CTRL_NIPBL && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2509729011.2 \
  --treatment \
  alignment/A549_CTRL_NIPBL_rep1/A549_CTRL_NIPBL_rep1.sorted.dup.bam \
  alignment/A549_CTRL_NIPBL_rep2/A549_CTRL_NIPBL_rep2.sorted.dup.bam \
  --control \
  alignment/A549_CTRL_WCE_rep1/A549_CTRL_WCE_rep1.sorted.dup.bam \
  alignment/A549_CTRL_WCE_rep2/A549_CTRL_WCE_rep2.sorted.dup.bam \
  --name peak_call/A549_CTRL_NIPBL/A549_CTRL_NIPBL \
  >& peak_call/A549_CTRL_NIPBL/A549_CTRL_NIPBL.diag.macs.out
macs2_callpeak.A549_CTRL_NIPBL.32bdbf16e5347684abc0921837f89437.mugqic.done
)
macs2_callpeak_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_5_JOB_ID: macs2_callpeak.A549_CTRL_SMC1A
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_CTRL_SMC1A
JOB_DEPENDENCIES=$picard_mark_duplicates_2_JOB_ID:$picard_mark_duplicates_3_JOB_ID:$picard_mark_duplicates_7_JOB_ID:$picard_mark_duplicates_8_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_CTRL_SMC1A.b815f80a367d46e2d9366a3dede9f70f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_CTRL_SMC1A.b815f80a367d46e2d9366a3dede9f70f.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/A549_CTRL_SMC1A && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2509729011.2 \
  --treatment \
  alignment/A549_CTRL_SMC1A_rep1/A549_CTRL_SMC1A_rep1.sorted.dup.bam \
  alignment/A549_CTRL_SMC1A_rep2/A549_CTRL_SMC1A_rep2.sorted.dup.bam \
  --control \
  alignment/A549_CTRL_WCE_rep1/A549_CTRL_WCE_rep1.sorted.dup.bam \
  alignment/A549_CTRL_WCE_rep2/A549_CTRL_WCE_rep2.sorted.dup.bam \
  --name peak_call/A549_CTRL_SMC1A/A549_CTRL_SMC1A \
  >& peak_call/A549_CTRL_SMC1A/A549_CTRL_SMC1A.diag.macs.out
macs2_callpeak.A549_CTRL_SMC1A.b815f80a367d46e2d9366a3dede9f70f.mugqic.done
)
macs2_callpeak_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_6_JOB_ID: macs2_callpeak.A549_DEX_BRD4
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_DEX_BRD4
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_DEX_BRD4.df6e3990741b83246e4e30a7be9ec96d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_DEX_BRD4.df6e3990741b83246e4e30a7be9ec96d.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/A549_DEX_BRD4 && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2509729011.2 \
  --treatment \
  alignment/A549_DEX_BRD4_rep1/A549_DEX_BRD4_rep1.sorted.dup.bam \
  --control \
  alignment/A549_DEX_WCE_rep1/A549_DEX_WCE_rep1.sorted.dup.bam \
  --name peak_call/A549_DEX_BRD4/A549_DEX_BRD4 \
  >& peak_call/A549_DEX_BRD4/A549_DEX_BRD4.diag.macs.out
macs2_callpeak.A549_DEX_BRD4.df6e3990741b83246e4e30a7be9ec96d.mugqic.done
)
macs2_callpeak_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_7_JOB_ID: macs2_callpeak.A549_DEX_CDK9
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_DEX_CDK9
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_DEX_CDK9.90cfc14de97d010780f9e057c76333b0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_DEX_CDK9.90cfc14de97d010780f9e057c76333b0.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/A549_DEX_CDK9 && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2509729011.2 \
  --treatment \
  alignment/A549_DEX_CDK9_rep1/A549_DEX_CDK9_rep1.sorted.dup.bam \
  --control \
  alignment/A549_DEX_WCE_rep1/A549_DEX_WCE_rep1.sorted.dup.bam \
  --name peak_call/A549_DEX_CDK9/A549_DEX_CDK9 \
  >& peak_call/A549_DEX_CDK9/A549_DEX_CDK9.diag.macs.out
macs2_callpeak.A549_DEX_CDK9.90cfc14de97d010780f9e057c76333b0.mugqic.done
)
macs2_callpeak_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$macs2_callpeak_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_8_JOB_ID: macs2_callpeak.A549_DEX_MED1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_DEX_MED1
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_DEX_MED1.a02f22993dcaed7dee50bba10abb2031.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_DEX_MED1.a02f22993dcaed7dee50bba10abb2031.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/A549_DEX_MED1 && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2509729011.2 \
  --treatment \
  alignment/A549_DEX_MED1_rep1/A549_DEX_MED1_rep1.sorted.dup.bam \
  --control \
  alignment/A549_DEX_WCE_rep1/A549_DEX_WCE_rep1.sorted.dup.bam \
  --name peak_call/A549_DEX_MED1/A549_DEX_MED1 \
  >& peak_call/A549_DEX_MED1/A549_DEX_MED1.diag.macs.out
macs2_callpeak.A549_DEX_MED1.a02f22993dcaed7dee50bba10abb2031.mugqic.done
)
macs2_callpeak_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$macs2_callpeak_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_9_JOB_ID: macs2_callpeak.A549_DEX_NIPBL
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_DEX_NIPBL
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_DEX_NIPBL.3adce8c7906b7545d5e9dab369af711c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_DEX_NIPBL.3adce8c7906b7545d5e9dab369af711c.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/A549_DEX_NIPBL && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2509729011.2 \
  --treatment \
  alignment/A549_DEX_NIPBL_rep1/A549_DEX_NIPBL_rep1.sorted.dup.bam \
  --control \
  alignment/A549_DEX_WCE_rep1/A549_DEX_WCE_rep1.sorted.dup.bam \
  --name peak_call/A549_DEX_NIPBL/A549_DEX_NIPBL \
  >& peak_call/A549_DEX_NIPBL/A549_DEX_NIPBL.diag.macs.out
macs2_callpeak.A549_DEX_NIPBL.3adce8c7906b7545d5e9dab369af711c.mugqic.done
)
macs2_callpeak_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$macs2_callpeak_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_10_JOB_ID: macs2_callpeak.A549_DEX_SMC1A_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_DEX_SMC1A_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_6_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_DEX_SMC1A_rep1.91d061fe0715ed9ab5a37ff069c8a69f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_DEX_SMC1A_rep1.91d061fe0715ed9ab5a37ff069c8a69f.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/A549_DEX_SMC1A_rep1 && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2509729011.2 \
  --treatment \
  alignment/A549_DEX_SMC1A_rep1/A549_DEX_SMC1A_rep1.sorted.dup.bam \
  --control \
  alignment/A549_DEX_WCE_rep1/A549_DEX_WCE_rep1.sorted.dup.bam \
  --name peak_call/A549_DEX_SMC1A_rep1/A549_DEX_SMC1A_rep1 \
  >& peak_call/A549_DEX_SMC1A_rep1/A549_DEX_SMC1A_rep1.diag.macs.out
macs2_callpeak.A549_DEX_SMC1A_rep1.91d061fe0715ed9ab5a37ff069c8a69f.mugqic.done
)
macs2_callpeak_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_11_JOB_ID: macs2_callpeak.A549_CTRL_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_CTRL_MED1_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_3_JOB_ID:$picard_mark_duplicates_4_JOB_ID:$picard_mark_duplicates_8_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_CTRL_MED1_rep1.5fa8c840ad9d051aa3dd52f8cf18c15e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_CTRL_MED1_rep1.5fa8c840ad9d051aa3dd52f8cf18c15e.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/A549_CTRL_MED1_rep1 && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2509729011.2 \
  --treatment \
  alignment/A549_CTRL_MED1_rep1/A549_CTRL_MED1_rep1.sorted.dup.bam \
  --control \
  alignment/A549_CTRL_WCE_rep1/A549_CTRL_WCE_rep1.sorted.dup.bam \
  alignment/A549_CTRL_WCE_rep2/A549_CTRL_WCE_rep2.sorted.dup.bam \
  --name peak_call/A549_CTRL_MED1_rep1/A549_CTRL_MED1_rep1 \
  >& peak_call/A549_CTRL_MED1_rep1/A549_CTRL_MED1_rep1.diag.macs.out
macs2_callpeak.A549_CTRL_MED1_rep1.5fa8c840ad9d051aa3dd52f8cf18c15e.mugqic.done
)
macs2_callpeak_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_12_JOB_ID: macs2_callpeak.A549_CTRL_NIPBL_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_CTRL_NIPBL_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_3_JOB_ID:$picard_mark_duplicates_5_JOB_ID:$picard_mark_duplicates_8_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_CTRL_NIPBL_rep1.9f777f83442079bb8eb0605be0c7a33f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_CTRL_NIPBL_rep1.9f777f83442079bb8eb0605be0c7a33f.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/A549_CTRL_NIPBL_rep1 && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2509729011.2 \
  --treatment \
  alignment/A549_CTRL_NIPBL_rep1/A549_CTRL_NIPBL_rep1.sorted.dup.bam \
  --control \
  alignment/A549_CTRL_WCE_rep1/A549_CTRL_WCE_rep1.sorted.dup.bam \
  alignment/A549_CTRL_WCE_rep2/A549_CTRL_WCE_rep2.sorted.dup.bam \
  --name peak_call/A549_CTRL_NIPBL_rep1/A549_CTRL_NIPBL_rep1 \
  >& peak_call/A549_CTRL_NIPBL_rep1/A549_CTRL_NIPBL_rep1.diag.macs.out
macs2_callpeak.A549_CTRL_NIPBL_rep1.9f777f83442079bb8eb0605be0c7a33f.mugqic.done
)
macs2_callpeak_12_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_13_JOB_ID: macs2_callpeak.A549_CTRL_SMC1A_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_CTRL_SMC1A_rep1
JOB_DEPENDENCIES=$picard_mark_duplicates_2_JOB_ID:$picard_mark_duplicates_3_JOB_ID:$picard_mark_duplicates_8_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_CTRL_SMC1A_rep1.a3fa6a5c579ba8d18b5b5b3b93ac3a3b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_CTRL_SMC1A_rep1.a3fa6a5c579ba8d18b5b5b3b93ac3a3b.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/A549_CTRL_SMC1A_rep1 && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2509729011.2 \
  --treatment \
  alignment/A549_CTRL_SMC1A_rep1/A549_CTRL_SMC1A_rep1.sorted.dup.bam \
  --control \
  alignment/A549_CTRL_WCE_rep1/A549_CTRL_WCE_rep1.sorted.dup.bam \
  alignment/A549_CTRL_WCE_rep2/A549_CTRL_WCE_rep2.sorted.dup.bam \
  --name peak_call/A549_CTRL_SMC1A_rep1/A549_CTRL_SMC1A_rep1 \
  >& peak_call/A549_CTRL_SMC1A_rep1/A549_CTRL_SMC1A_rep1.diag.macs.out
macs2_callpeak.A549_CTRL_SMC1A_rep1.a3fa6a5c579ba8d18b5b5b3b93ac3a3b.mugqic.done
)
macs2_callpeak_13_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_13_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_14_JOB_ID: macs2_callpeak.A549_CTRL_MED1_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_CTRL_MED1_rep2
JOB_DEPENDENCIES=$picard_mark_duplicates_3_JOB_ID:$picard_mark_duplicates_8_JOB_ID:$picard_mark_duplicates_10_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_CTRL_MED1_rep2.6f9a8812035287afce48363582ec8240.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_CTRL_MED1_rep2.6f9a8812035287afce48363582ec8240.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/A549_CTRL_MED1_rep2 && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2509729011.2 \
  --treatment \
  alignment/A549_CTRL_MED1_rep2/A549_CTRL_MED1_rep2.sorted.dup.bam \
  --control \
  alignment/A549_CTRL_WCE_rep1/A549_CTRL_WCE_rep1.sorted.dup.bam \
  alignment/A549_CTRL_WCE_rep2/A549_CTRL_WCE_rep2.sorted.dup.bam \
  --name peak_call/A549_CTRL_MED1_rep2/A549_CTRL_MED1_rep2 \
  >& peak_call/A549_CTRL_MED1_rep2/A549_CTRL_MED1_rep2.diag.macs.out
macs2_callpeak.A549_CTRL_MED1_rep2.6f9a8812035287afce48363582ec8240.mugqic.done
)
macs2_callpeak_14_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_14_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_15_JOB_ID: macs2_callpeak.A549_CTRL_NIPBL_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_CTRL_NIPBL_rep2
JOB_DEPENDENCIES=$picard_mark_duplicates_3_JOB_ID:$picard_mark_duplicates_8_JOB_ID:$picard_mark_duplicates_9_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_CTRL_NIPBL_rep2.9204000ce1d7e652e83cb0410938701e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_CTRL_NIPBL_rep2.9204000ce1d7e652e83cb0410938701e.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/A549_CTRL_NIPBL_rep2 && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2509729011.2 \
  --treatment \
  alignment/A549_CTRL_NIPBL_rep2/A549_CTRL_NIPBL_rep2.sorted.dup.bam \
  --control \
  alignment/A549_CTRL_WCE_rep1/A549_CTRL_WCE_rep1.sorted.dup.bam \
  alignment/A549_CTRL_WCE_rep2/A549_CTRL_WCE_rep2.sorted.dup.bam \
  --name peak_call/A549_CTRL_NIPBL_rep2/A549_CTRL_NIPBL_rep2 \
  >& peak_call/A549_CTRL_NIPBL_rep2/A549_CTRL_NIPBL_rep2.diag.macs.out
macs2_callpeak.A549_CTRL_NIPBL_rep2.9204000ce1d7e652e83cb0410938701e.mugqic.done
)
macs2_callpeak_15_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_15_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_16_JOB_ID: macs2_callpeak.A549_CTRL_SMC1A_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_CTRL_SMC1A_rep2
JOB_DEPENDENCIES=$picard_mark_duplicates_3_JOB_ID:$picard_mark_duplicates_7_JOB_ID:$picard_mark_duplicates_8_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_CTRL_SMC1A_rep2.990bfe6a0f6d544516580fbb8d0b2a0e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_CTRL_SMC1A_rep2.990bfe6a0f6d544516580fbb8d0b2a0e.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/A549_CTRL_SMC1A_rep2 && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2509729011.2 \
  --treatment \
  alignment/A549_CTRL_SMC1A_rep2/A549_CTRL_SMC1A_rep2.sorted.dup.bam \
  --control \
  alignment/A549_CTRL_WCE_rep1/A549_CTRL_WCE_rep1.sorted.dup.bam \
  alignment/A549_CTRL_WCE_rep2/A549_CTRL_WCE_rep2.sorted.dup.bam \
  --name peak_call/A549_CTRL_SMC1A_rep2/A549_CTRL_SMC1A_rep2 \
  >& peak_call/A549_CTRL_SMC1A_rep2/A549_CTRL_SMC1A_rep2.diag.macs.out
macs2_callpeak.A549_CTRL_SMC1A_rep2.990bfe6a0f6d544516580fbb8d0b2a0e.mugqic.done
)
macs2_callpeak_16_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_16_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_17_JOB_ID: macs2_callpeak_report
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_report
JOB_DEPENDENCIES=$macs2_callpeak_1_JOB_ID:$macs2_callpeak_2_JOB_ID:$macs2_callpeak_3_JOB_ID:$macs2_callpeak_4_JOB_ID:$macs2_callpeak_5_JOB_ID:$macs2_callpeak_6_JOB_ID:$macs2_callpeak_7_JOB_ID:$macs2_callpeak_8_JOB_ID:$macs2_callpeak_9_JOB_ID:$macs2_callpeak_10_JOB_ID:$macs2_callpeak_11_JOB_ID:$macs2_callpeak_12_JOB_ID:$macs2_callpeak_13_JOB_ID:$macs2_callpeak_14_JOB_ID:$macs2_callpeak_15_JOB_ID:$macs2_callpeak_16_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_report.99b505db119be7008d733d72818c4b4a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_report.99b505db119be7008d733d72818c4b4a.mugqic.done'
mkdir -p report && \
cp /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/ChipSeq.macs2_callpeak.md report/ && \
for contrast in A549_CTRL_BRD4 A549_CTRL_CDK9 A549_CTRL_MED1 A549_CTRL_NIPBL A549_CTRL_SMC1A A549_DEX_BRD4 A549_DEX_CDK9 A549_DEX_MED1 A549_DEX_NIPBL A549_DEX_SMC1A_rep1 A549_CTRL_MED1_rep1 A549_CTRL_NIPBL_rep1 A549_CTRL_SMC1A_rep1 A549_CTRL_MED1_rep2 A549_CTRL_NIPBL_rep2 A549_CTRL_SMC1A_rep2
do
  cp -a --parents peak_call/$contrast/ report/ && \
  echo -e "* [Peak Calls File for Design $contrast](peak_call/$contrast/${contrast}_peaks.xls)" \
  >> report/ChipSeq.macs2_callpeak.md
done
macs2_callpeak_report.99b505db119be7008d733d72818c4b4a.mugqic.done
)
macs2_callpeak_17_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$macs2_callpeak_17_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


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
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.A549_CTRL_BRD4.b6fcf04a12d9d7d8bb983750420b036b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.A549_CTRL_BRD4.b6fcf04a12d9d7d8bb983750420b036b.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/A549_CTRL_BRD4/A549_CTRL_BRD4 && \
annotatePeaks.pl \
  peak_call/A549_CTRL_BRD4/A549_CTRL_BRD4_peaks.narrowPeak \
  hg19 \
  -gsize hg19 \
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
homer_annotate_peaks.A549_CTRL_BRD4.b6fcf04a12d9d7d8bb983750420b036b.mugqic.done
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
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.A549_CTRL_CDK9.25b46eaf0a898fdd86fb2435cf69ccff.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.A549_CTRL_CDK9.25b46eaf0a898fdd86fb2435cf69ccff.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/A549_CTRL_CDK9/A549_CTRL_CDK9 && \
annotatePeaks.pl \
  peak_call/A549_CTRL_CDK9/A549_CTRL_CDK9_peaks.narrowPeak \
  hg19 \
  -gsize hg19 \
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
homer_annotate_peaks.A549_CTRL_CDK9.25b46eaf0a898fdd86fb2435cf69ccff.mugqic.done
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
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.A549_CTRL_MED1.48541b8bda17425ee2e7faccc90315a5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.A549_CTRL_MED1.48541b8bda17425ee2e7faccc90315a5.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/A549_CTRL_MED1/A549_CTRL_MED1 && \
annotatePeaks.pl \
  peak_call/A549_CTRL_MED1/A549_CTRL_MED1_peaks.narrowPeak \
  hg19 \
  -gsize hg19 \
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
homer_annotate_peaks.A549_CTRL_MED1.48541b8bda17425ee2e7faccc90315a5.mugqic.done
)
homer_annotate_peaks_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_4_JOB_ID: homer_annotate_peaks.A549_CTRL_NIPBL
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.A549_CTRL_NIPBL
JOB_DEPENDENCIES=$macs2_callpeak_4_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.A549_CTRL_NIPBL.8e353268a56697c9e0817c5e12297df3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.A549_CTRL_NIPBL.8e353268a56697c9e0817c5e12297df3.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/A549_CTRL_NIPBL/A549_CTRL_NIPBL && \
annotatePeaks.pl \
  peak_call/A549_CTRL_NIPBL/A549_CTRL_NIPBL_peaks.narrowPeak \
  hg19 \
  -gsize hg19 \
  -cons -CpG \
  -go annotation/A549_CTRL_NIPBL/A549_CTRL_NIPBL \
  -genomeOntology annotation/A549_CTRL_NIPBL/A549_CTRL_NIPBL \
  > annotation/A549_CTRL_NIPBL/A549_CTRL_NIPBL.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/A549_CTRL_NIPBL/A549_CTRL_NIPBL.annotated.csv",
  "annotation/A549_CTRL_NIPBL/A549_CTRL_NIPBL",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.A549_CTRL_NIPBL.8e353268a56697c9e0817c5e12297df3.mugqic.done
)
homer_annotate_peaks_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_5_JOB_ID: homer_annotate_peaks.A549_CTRL_SMC1A
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.A549_CTRL_SMC1A
JOB_DEPENDENCIES=$macs2_callpeak_5_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.A549_CTRL_SMC1A.d1da66e913132bd9703b2407b6ddff9a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.A549_CTRL_SMC1A.d1da66e913132bd9703b2407b6ddff9a.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/A549_CTRL_SMC1A/A549_CTRL_SMC1A && \
annotatePeaks.pl \
  peak_call/A549_CTRL_SMC1A/A549_CTRL_SMC1A_peaks.narrowPeak \
  hg19 \
  -gsize hg19 \
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
homer_annotate_peaks.A549_CTRL_SMC1A.d1da66e913132bd9703b2407b6ddff9a.mugqic.done
)
homer_annotate_peaks_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_6_JOB_ID: homer_annotate_peaks.A549_DEX_BRD4
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.A549_DEX_BRD4
JOB_DEPENDENCIES=$macs2_callpeak_6_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.A549_DEX_BRD4.345543c9fdc20c19a9113dc9178d31d8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.A549_DEX_BRD4.345543c9fdc20c19a9113dc9178d31d8.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/A549_DEX_BRD4/A549_DEX_BRD4 && \
annotatePeaks.pl \
  peak_call/A549_DEX_BRD4/A549_DEX_BRD4_peaks.narrowPeak \
  hg19 \
  -gsize hg19 \
  -cons -CpG \
  -go annotation/A549_DEX_BRD4/A549_DEX_BRD4 \
  -genomeOntology annotation/A549_DEX_BRD4/A549_DEX_BRD4 \
  > annotation/A549_DEX_BRD4/A549_DEX_BRD4.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/A549_DEX_BRD4/A549_DEX_BRD4.annotated.csv",
  "annotation/A549_DEX_BRD4/A549_DEX_BRD4",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.A549_DEX_BRD4.345543c9fdc20c19a9113dc9178d31d8.mugqic.done
)
homer_annotate_peaks_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_7_JOB_ID: homer_annotate_peaks.A549_DEX_CDK9
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.A549_DEX_CDK9
JOB_DEPENDENCIES=$macs2_callpeak_7_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.A549_DEX_CDK9.d5ed003446f13fc66066d838db734154.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.A549_DEX_CDK9.d5ed003446f13fc66066d838db734154.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/A549_DEX_CDK9/A549_DEX_CDK9 && \
annotatePeaks.pl \
  peak_call/A549_DEX_CDK9/A549_DEX_CDK9_peaks.narrowPeak \
  hg19 \
  -gsize hg19 \
  -cons -CpG \
  -go annotation/A549_DEX_CDK9/A549_DEX_CDK9 \
  -genomeOntology annotation/A549_DEX_CDK9/A549_DEX_CDK9 \
  > annotation/A549_DEX_CDK9/A549_DEX_CDK9.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/A549_DEX_CDK9/A549_DEX_CDK9.annotated.csv",
  "annotation/A549_DEX_CDK9/A549_DEX_CDK9",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.A549_DEX_CDK9.d5ed003446f13fc66066d838db734154.mugqic.done
)
homer_annotate_peaks_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_8_JOB_ID: homer_annotate_peaks.A549_DEX_MED1
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.A549_DEX_MED1
JOB_DEPENDENCIES=$macs2_callpeak_8_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.A549_DEX_MED1.9e2bf672ce34d020c48d84e9216ce784.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.A549_DEX_MED1.9e2bf672ce34d020c48d84e9216ce784.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/A549_DEX_MED1/A549_DEX_MED1 && \
annotatePeaks.pl \
  peak_call/A549_DEX_MED1/A549_DEX_MED1_peaks.narrowPeak \
  hg19 \
  -gsize hg19 \
  -cons -CpG \
  -go annotation/A549_DEX_MED1/A549_DEX_MED1 \
  -genomeOntology annotation/A549_DEX_MED1/A549_DEX_MED1 \
  > annotation/A549_DEX_MED1/A549_DEX_MED1.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/A549_DEX_MED1/A549_DEX_MED1.annotated.csv",
  "annotation/A549_DEX_MED1/A549_DEX_MED1",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.A549_DEX_MED1.9e2bf672ce34d020c48d84e9216ce784.mugqic.done
)
homer_annotate_peaks_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_9_JOB_ID: homer_annotate_peaks.A549_DEX_NIPBL
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.A549_DEX_NIPBL
JOB_DEPENDENCIES=$macs2_callpeak_9_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.A549_DEX_NIPBL.f140f66c08d34af31b1e88570f8ddd58.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.A549_DEX_NIPBL.f140f66c08d34af31b1e88570f8ddd58.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/A549_DEX_NIPBL/A549_DEX_NIPBL && \
annotatePeaks.pl \
  peak_call/A549_DEX_NIPBL/A549_DEX_NIPBL_peaks.narrowPeak \
  hg19 \
  -gsize hg19 \
  -cons -CpG \
  -go annotation/A549_DEX_NIPBL/A549_DEX_NIPBL \
  -genomeOntology annotation/A549_DEX_NIPBL/A549_DEX_NIPBL \
  > annotation/A549_DEX_NIPBL/A549_DEX_NIPBL.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/A549_DEX_NIPBL/A549_DEX_NIPBL.annotated.csv",
  "annotation/A549_DEX_NIPBL/A549_DEX_NIPBL",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.A549_DEX_NIPBL.f140f66c08d34af31b1e88570f8ddd58.mugqic.done
)
homer_annotate_peaks_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_10_JOB_ID: homer_annotate_peaks.A549_DEX_SMC1A_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.A549_DEX_SMC1A_rep1
JOB_DEPENDENCIES=$macs2_callpeak_10_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.A549_DEX_SMC1A_rep1.646178cd746142e9024364b7f0efd5ae.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.A549_DEX_SMC1A_rep1.646178cd746142e9024364b7f0efd5ae.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/A549_DEX_SMC1A_rep1/A549_DEX_SMC1A_rep1 && \
annotatePeaks.pl \
  peak_call/A549_DEX_SMC1A_rep1/A549_DEX_SMC1A_rep1_peaks.narrowPeak \
  hg19 \
  -gsize hg19 \
  -cons -CpG \
  -go annotation/A549_DEX_SMC1A_rep1/A549_DEX_SMC1A_rep1 \
  -genomeOntology annotation/A549_DEX_SMC1A_rep1/A549_DEX_SMC1A_rep1 \
  > annotation/A549_DEX_SMC1A_rep1/A549_DEX_SMC1A_rep1.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/A549_DEX_SMC1A_rep1/A549_DEX_SMC1A_rep1.annotated.csv",
  "annotation/A549_DEX_SMC1A_rep1/A549_DEX_SMC1A_rep1",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.A549_DEX_SMC1A_rep1.646178cd746142e9024364b7f0efd5ae.mugqic.done
)
homer_annotate_peaks_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_11_JOB_ID: homer_annotate_peaks.A549_CTRL_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.A549_CTRL_MED1_rep1
JOB_DEPENDENCIES=$macs2_callpeak_11_JOB_ID
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
homer_annotate_peaks_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_12_JOB_ID: homer_annotate_peaks.A549_CTRL_NIPBL_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.A549_CTRL_NIPBL_rep1
JOB_DEPENDENCIES=$macs2_callpeak_12_JOB_ID
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
homer_annotate_peaks_12_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_13_JOB_ID: homer_annotate_peaks.A549_CTRL_SMC1A_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.A549_CTRL_SMC1A_rep1
JOB_DEPENDENCIES=$macs2_callpeak_13_JOB_ID
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
homer_annotate_peaks_13_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_13_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_14_JOB_ID: homer_annotate_peaks.A549_CTRL_MED1_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.A549_CTRL_MED1_rep2
JOB_DEPENDENCIES=$macs2_callpeak_14_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.A549_CTRL_MED1_rep2.e9f949810e69a9a94a3d7cce638c4766.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.A549_CTRL_MED1_rep2.e9f949810e69a9a94a3d7cce638c4766.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/A549_CTRL_MED1_rep2/A549_CTRL_MED1_rep2 && \
annotatePeaks.pl \
  peak_call/A549_CTRL_MED1_rep2/A549_CTRL_MED1_rep2_peaks.narrowPeak \
  hg19 \
  -gsize hg19 \
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
homer_annotate_peaks.A549_CTRL_MED1_rep2.e9f949810e69a9a94a3d7cce638c4766.mugqic.done
)
homer_annotate_peaks_14_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_14_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_15_JOB_ID: homer_annotate_peaks.A549_CTRL_NIPBL_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.A549_CTRL_NIPBL_rep2
JOB_DEPENDENCIES=$macs2_callpeak_15_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.A549_CTRL_NIPBL_rep2.80c9e8f847b4f5dfad81f6320c5377d8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.A549_CTRL_NIPBL_rep2.80c9e8f847b4f5dfad81f6320c5377d8.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/A549_CTRL_NIPBL_rep2/A549_CTRL_NIPBL_rep2 && \
annotatePeaks.pl \
  peak_call/A549_CTRL_NIPBL_rep2/A549_CTRL_NIPBL_rep2_peaks.narrowPeak \
  hg19 \
  -gsize hg19 \
  -cons -CpG \
  -go annotation/A549_CTRL_NIPBL_rep2/A549_CTRL_NIPBL_rep2 \
  -genomeOntology annotation/A549_CTRL_NIPBL_rep2/A549_CTRL_NIPBL_rep2 \
  > annotation/A549_CTRL_NIPBL_rep2/A549_CTRL_NIPBL_rep2.annotated.csv && \
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "annotation/A549_CTRL_NIPBL_rep2/A549_CTRL_NIPBL_rep2.annotated.csv",
  "annotation/A549_CTRL_NIPBL_rep2/A549_CTRL_NIPBL_rep2",
  -2000,
  -10000,
  -10000,
  -100000,
  100000
)'
homer_annotate_peaks.A549_CTRL_NIPBL_rep2.80c9e8f847b4f5dfad81f6320c5377d8.mugqic.done
)
homer_annotate_peaks_15_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_15_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_16_JOB_ID: homer_annotate_peaks.A549_CTRL_SMC1A_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.A549_CTRL_SMC1A_rep2
JOB_DEPENDENCIES=$macs2_callpeak_16_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks.A549_CTRL_SMC1A_rep2.84906e93638e466c7be34d5b287638c4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks.A549_CTRL_SMC1A_rep2.84906e93638e466c7be34d5b287638c4.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.7 mugqic/mugqic_tools/2.1.5 && \
mkdir -p annotation/A549_CTRL_SMC1A_rep2/A549_CTRL_SMC1A_rep2 && \
annotatePeaks.pl \
  peak_call/A549_CTRL_SMC1A_rep2/A549_CTRL_SMC1A_rep2_peaks.narrowPeak \
  hg19 \
  -gsize hg19 \
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
homer_annotate_peaks.A549_CTRL_SMC1A_rep2.84906e93638e466c7be34d5b287638c4.mugqic.done
)
homer_annotate_peaks_16_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_16_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_17_JOB_ID: homer_annotate_peaks_report
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks_report
JOB_DEPENDENCIES=$homer_annotate_peaks_1_JOB_ID:$homer_annotate_peaks_2_JOB_ID:$homer_annotate_peaks_3_JOB_ID:$homer_annotate_peaks_4_JOB_ID:$homer_annotate_peaks_5_JOB_ID:$homer_annotate_peaks_6_JOB_ID:$homer_annotate_peaks_7_JOB_ID:$homer_annotate_peaks_8_JOB_ID:$homer_annotate_peaks_9_JOB_ID:$homer_annotate_peaks_10_JOB_ID:$homer_annotate_peaks_11_JOB_ID:$homer_annotate_peaks_12_JOB_ID:$homer_annotate_peaks_13_JOB_ID:$homer_annotate_peaks_14_JOB_ID:$homer_annotate_peaks_15_JOB_ID:$homer_annotate_peaks_16_JOB_ID
JOB_DONE=job_output/homer_annotate_peaks/homer_annotate_peaks_report.9e9da26c00b01001f54b1ba97eac169e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_annotate_peaks_report.9e9da26c00b01001f54b1ba97eac169e.mugqic.done'
mkdir -p report/annotation/ && \
cp /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/ChipSeq.homer_annotate_peaks.md report/ && \
for contrast in A549_CTRL_BRD4 A549_CTRL_CDK9 A549_CTRL_MED1 A549_CTRL_NIPBL A549_CTRL_SMC1A A549_DEX_BRD4 A549_DEX_CDK9 A549_DEX_MED1 A549_DEX_NIPBL A549_DEX_SMC1A_rep1 A549_CTRL_MED1_rep1 A549_CTRL_NIPBL_rep1 A549_CTRL_SMC1A_rep1 A549_CTRL_MED1_rep2 A549_CTRL_NIPBL_rep2 A549_CTRL_SMC1A_rep2
do
  rsync -avP annotation/$contrast report/annotation/ && \
  echo -e "* [Gene Annotations for Design $contrast](annotation/$contrast/${contrast}.annotated.csv)
* [HOMER Gene Ontology Annotations for Design $contrast](annotation/$contrast/$contrast/geneOntology.html)
* [HOMER Genome Ontology Annotations for Design $contrast](annotation/$contrast/$contrast/GenomeOntology.html)" \
  >> report/ChipSeq.homer_annotate_peaks.md
done
homer_annotate_peaks_report.9e9da26c00b01001f54b1ba97eac169e.mugqic.done
)
homer_annotate_peaks_17_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_annotate_peaks_17_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


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
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.A549_CTRL_BRD4.8df8c3a7aeea2da881d8d37c7c9ae6de.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.A549_CTRL_BRD4.8df8c3a7aeea2da881d8d37c7c9ae6de.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/A549_CTRL_BRD4/A549_CTRL_BRD4 && \
findMotifsGenome.pl \
  peak_call/A549_CTRL_BRD4/A549_CTRL_BRD4_peaks.narrowPeak \
  hg19 \
  annotation/A549_CTRL_BRD4/A549_CTRL_BRD4 \
  -preparsedDir annotation/A549_CTRL_BRD4/A549_CTRL_BRD4/preparsed \
  -p 4
homer_find_motifs_genome.A549_CTRL_BRD4.8df8c3a7aeea2da881d8d37c7c9ae6de.mugqic.done
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
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.A549_CTRL_CDK9.4753a0cabb3f83f63732264ce86d97b9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.A549_CTRL_CDK9.4753a0cabb3f83f63732264ce86d97b9.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/A549_CTRL_CDK9/A549_CTRL_CDK9 && \
findMotifsGenome.pl \
  peak_call/A549_CTRL_CDK9/A549_CTRL_CDK9_peaks.narrowPeak \
  hg19 \
  annotation/A549_CTRL_CDK9/A549_CTRL_CDK9 \
  -preparsedDir annotation/A549_CTRL_CDK9/A549_CTRL_CDK9/preparsed \
  -p 4
homer_find_motifs_genome.A549_CTRL_CDK9.4753a0cabb3f83f63732264ce86d97b9.mugqic.done
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
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.A549_CTRL_MED1.37833b73eead4192681bc607be4c50b2.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.A549_CTRL_MED1.37833b73eead4192681bc607be4c50b2.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/A549_CTRL_MED1/A549_CTRL_MED1 && \
findMotifsGenome.pl \
  peak_call/A549_CTRL_MED1/A549_CTRL_MED1_peaks.narrowPeak \
  hg19 \
  annotation/A549_CTRL_MED1/A549_CTRL_MED1 \
  -preparsedDir annotation/A549_CTRL_MED1/A549_CTRL_MED1/preparsed \
  -p 4
homer_find_motifs_genome.A549_CTRL_MED1.37833b73eead4192681bc607be4c50b2.mugqic.done
)
homer_find_motifs_genome_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_4_JOB_ID: homer_find_motifs_genome.A549_CTRL_NIPBL
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.A549_CTRL_NIPBL
JOB_DEPENDENCIES=$macs2_callpeak_4_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.A549_CTRL_NIPBL.5c6f5c80de76f28bcb2649f2f1637e47.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.A549_CTRL_NIPBL.5c6f5c80de76f28bcb2649f2f1637e47.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/A549_CTRL_NIPBL/A549_CTRL_NIPBL && \
findMotifsGenome.pl \
  peak_call/A549_CTRL_NIPBL/A549_CTRL_NIPBL_peaks.narrowPeak \
  hg19 \
  annotation/A549_CTRL_NIPBL/A549_CTRL_NIPBL \
  -preparsedDir annotation/A549_CTRL_NIPBL/A549_CTRL_NIPBL/preparsed \
  -p 4
homer_find_motifs_genome.A549_CTRL_NIPBL.5c6f5c80de76f28bcb2649f2f1637e47.mugqic.done
)
homer_find_motifs_genome_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_5_JOB_ID: homer_find_motifs_genome.A549_CTRL_SMC1A
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.A549_CTRL_SMC1A
JOB_DEPENDENCIES=$macs2_callpeak_5_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.A549_CTRL_SMC1A.d15d78104735482e14c3b97988cd72d9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.A549_CTRL_SMC1A.d15d78104735482e14c3b97988cd72d9.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/A549_CTRL_SMC1A/A549_CTRL_SMC1A && \
findMotifsGenome.pl \
  peak_call/A549_CTRL_SMC1A/A549_CTRL_SMC1A_peaks.narrowPeak \
  hg19 \
  annotation/A549_CTRL_SMC1A/A549_CTRL_SMC1A \
  -preparsedDir annotation/A549_CTRL_SMC1A/A549_CTRL_SMC1A/preparsed \
  -p 4
homer_find_motifs_genome.A549_CTRL_SMC1A.d15d78104735482e14c3b97988cd72d9.mugqic.done
)
homer_find_motifs_genome_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_6_JOB_ID: homer_find_motifs_genome.A549_DEX_BRD4
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.A549_DEX_BRD4
JOB_DEPENDENCIES=$macs2_callpeak_6_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.A549_DEX_BRD4.410e1af67b1642071ed9062c7befff58.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.A549_DEX_BRD4.410e1af67b1642071ed9062c7befff58.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/A549_DEX_BRD4/A549_DEX_BRD4 && \
findMotifsGenome.pl \
  peak_call/A549_DEX_BRD4/A549_DEX_BRD4_peaks.narrowPeak \
  hg19 \
  annotation/A549_DEX_BRD4/A549_DEX_BRD4 \
  -preparsedDir annotation/A549_DEX_BRD4/A549_DEX_BRD4/preparsed \
  -p 4
homer_find_motifs_genome.A549_DEX_BRD4.410e1af67b1642071ed9062c7befff58.mugqic.done
)
homer_find_motifs_genome_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_7_JOB_ID: homer_find_motifs_genome.A549_DEX_CDK9
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.A549_DEX_CDK9
JOB_DEPENDENCIES=$macs2_callpeak_7_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.A549_DEX_CDK9.352a96ad9419a0f240e38a77f02a02c6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.A549_DEX_CDK9.352a96ad9419a0f240e38a77f02a02c6.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/A549_DEX_CDK9/A549_DEX_CDK9 && \
findMotifsGenome.pl \
  peak_call/A549_DEX_CDK9/A549_DEX_CDK9_peaks.narrowPeak \
  hg19 \
  annotation/A549_DEX_CDK9/A549_DEX_CDK9 \
  -preparsedDir annotation/A549_DEX_CDK9/A549_DEX_CDK9/preparsed \
  -p 4
homer_find_motifs_genome.A549_DEX_CDK9.352a96ad9419a0f240e38a77f02a02c6.mugqic.done
)
homer_find_motifs_genome_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_8_JOB_ID: homer_find_motifs_genome.A549_DEX_MED1
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.A549_DEX_MED1
JOB_DEPENDENCIES=$macs2_callpeak_8_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.A549_DEX_MED1.33ef9be918408c32f5a4a36ef1d0ff66.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.A549_DEX_MED1.33ef9be918408c32f5a4a36ef1d0ff66.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/A549_DEX_MED1/A549_DEX_MED1 && \
findMotifsGenome.pl \
  peak_call/A549_DEX_MED1/A549_DEX_MED1_peaks.narrowPeak \
  hg19 \
  annotation/A549_DEX_MED1/A549_DEX_MED1 \
  -preparsedDir annotation/A549_DEX_MED1/A549_DEX_MED1/preparsed \
  -p 4
homer_find_motifs_genome.A549_DEX_MED1.33ef9be918408c32f5a4a36ef1d0ff66.mugqic.done
)
homer_find_motifs_genome_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_9_JOB_ID: homer_find_motifs_genome.A549_DEX_NIPBL
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.A549_DEX_NIPBL
JOB_DEPENDENCIES=$macs2_callpeak_9_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.A549_DEX_NIPBL.bf77c956127575c962e1480c413607fb.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.A549_DEX_NIPBL.bf77c956127575c962e1480c413607fb.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/A549_DEX_NIPBL/A549_DEX_NIPBL && \
findMotifsGenome.pl \
  peak_call/A549_DEX_NIPBL/A549_DEX_NIPBL_peaks.narrowPeak \
  hg19 \
  annotation/A549_DEX_NIPBL/A549_DEX_NIPBL \
  -preparsedDir annotation/A549_DEX_NIPBL/A549_DEX_NIPBL/preparsed \
  -p 4
homer_find_motifs_genome.A549_DEX_NIPBL.bf77c956127575c962e1480c413607fb.mugqic.done
)
homer_find_motifs_genome_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_10_JOB_ID: homer_find_motifs_genome.A549_DEX_SMC1A_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.A549_DEX_SMC1A_rep1
JOB_DEPENDENCIES=$macs2_callpeak_10_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.A549_DEX_SMC1A_rep1.4cda0a44dbcfa24277c296bf87ed40a1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.A549_DEX_SMC1A_rep1.4cda0a44dbcfa24277c296bf87ed40a1.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/A549_DEX_SMC1A_rep1/A549_DEX_SMC1A_rep1 && \
findMotifsGenome.pl \
  peak_call/A549_DEX_SMC1A_rep1/A549_DEX_SMC1A_rep1_peaks.narrowPeak \
  hg19 \
  annotation/A549_DEX_SMC1A_rep1/A549_DEX_SMC1A_rep1 \
  -preparsedDir annotation/A549_DEX_SMC1A_rep1/A549_DEX_SMC1A_rep1/preparsed \
  -p 4
homer_find_motifs_genome.A549_DEX_SMC1A_rep1.4cda0a44dbcfa24277c296bf87ed40a1.mugqic.done
)
homer_find_motifs_genome_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_11_JOB_ID: homer_find_motifs_genome.A549_CTRL_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.A549_CTRL_MED1_rep1
JOB_DEPENDENCIES=$macs2_callpeak_11_JOB_ID
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
homer_find_motifs_genome_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_12_JOB_ID: homer_find_motifs_genome.A549_CTRL_NIPBL_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.A549_CTRL_NIPBL_rep1
JOB_DEPENDENCIES=$macs2_callpeak_12_JOB_ID
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
homer_find_motifs_genome_12_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_13_JOB_ID: homer_find_motifs_genome.A549_CTRL_SMC1A_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.A549_CTRL_SMC1A_rep1
JOB_DEPENDENCIES=$macs2_callpeak_13_JOB_ID
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
homer_find_motifs_genome_13_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_13_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_14_JOB_ID: homer_find_motifs_genome.A549_CTRL_MED1_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.A549_CTRL_MED1_rep2
JOB_DEPENDENCIES=$macs2_callpeak_14_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.A549_CTRL_MED1_rep2.8ba725ade5fe942cbea7e75cd6c2dae0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.A549_CTRL_MED1_rep2.8ba725ade5fe942cbea7e75cd6c2dae0.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/A549_CTRL_MED1_rep2/A549_CTRL_MED1_rep2 && \
findMotifsGenome.pl \
  peak_call/A549_CTRL_MED1_rep2/A549_CTRL_MED1_rep2_peaks.narrowPeak \
  hg19 \
  annotation/A549_CTRL_MED1_rep2/A549_CTRL_MED1_rep2 \
  -preparsedDir annotation/A549_CTRL_MED1_rep2/A549_CTRL_MED1_rep2/preparsed \
  -p 4
homer_find_motifs_genome.A549_CTRL_MED1_rep2.8ba725ade5fe942cbea7e75cd6c2dae0.mugqic.done
)
homer_find_motifs_genome_14_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_14_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_15_JOB_ID: homer_find_motifs_genome.A549_CTRL_NIPBL_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.A549_CTRL_NIPBL_rep2
JOB_DEPENDENCIES=$macs2_callpeak_15_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.A549_CTRL_NIPBL_rep2.8b7c3efe1749e6afc440a231e841cf93.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.A549_CTRL_NIPBL_rep2.8b7c3efe1749e6afc440a231e841cf93.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/A549_CTRL_NIPBL_rep2/A549_CTRL_NIPBL_rep2 && \
findMotifsGenome.pl \
  peak_call/A549_CTRL_NIPBL_rep2/A549_CTRL_NIPBL_rep2_peaks.narrowPeak \
  hg19 \
  annotation/A549_CTRL_NIPBL_rep2/A549_CTRL_NIPBL_rep2 \
  -preparsedDir annotation/A549_CTRL_NIPBL_rep2/A549_CTRL_NIPBL_rep2/preparsed \
  -p 4
homer_find_motifs_genome.A549_CTRL_NIPBL_rep2.8b7c3efe1749e6afc440a231e841cf93.mugqic.done
)
homer_find_motifs_genome_15_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_15_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_16_JOB_ID: homer_find_motifs_genome.A549_CTRL_SMC1A_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.A549_CTRL_SMC1A_rep2
JOB_DEPENDENCIES=$macs2_callpeak_16_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome.A549_CTRL_SMC1A_rep2.04b8bde5e3fdbcfa15c263de09071fce.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome.A549_CTRL_SMC1A_rep2.04b8bde5e3fdbcfa15c263de09071fce.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/weblogo/3.3 mugqic/homer/4.7 && \
mkdir -p annotation/A549_CTRL_SMC1A_rep2/A549_CTRL_SMC1A_rep2 && \
findMotifsGenome.pl \
  peak_call/A549_CTRL_SMC1A_rep2/A549_CTRL_SMC1A_rep2_peaks.narrowPeak \
  hg19 \
  annotation/A549_CTRL_SMC1A_rep2/A549_CTRL_SMC1A_rep2 \
  -preparsedDir annotation/A549_CTRL_SMC1A_rep2/A549_CTRL_SMC1A_rep2/preparsed \
  -p 4
homer_find_motifs_genome.A549_CTRL_SMC1A_rep2.04b8bde5e3fdbcfa15c263de09071fce.mugqic.done
)
homer_find_motifs_genome_16_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_16_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_17_JOB_ID: homer_find_motifs_genome_report
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome_report
JOB_DEPENDENCIES=$homer_find_motifs_genome_1_JOB_ID:$homer_find_motifs_genome_2_JOB_ID:$homer_find_motifs_genome_3_JOB_ID:$homer_find_motifs_genome_4_JOB_ID:$homer_find_motifs_genome_5_JOB_ID:$homer_find_motifs_genome_6_JOB_ID:$homer_find_motifs_genome_7_JOB_ID:$homer_find_motifs_genome_8_JOB_ID:$homer_find_motifs_genome_9_JOB_ID:$homer_find_motifs_genome_10_JOB_ID:$homer_find_motifs_genome_11_JOB_ID:$homer_find_motifs_genome_12_JOB_ID:$homer_find_motifs_genome_13_JOB_ID:$homer_find_motifs_genome_14_JOB_ID:$homer_find_motifs_genome_15_JOB_ID:$homer_find_motifs_genome_16_JOB_ID
JOB_DONE=job_output/homer_find_motifs_genome/homer_find_motifs_genome_report.9905966bef4fd299070c52342290bbfa.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_find_motifs_genome_report.9905966bef4fd299070c52342290bbfa.mugqic.done'
mkdir -p report/annotation/ && \
cp /cvmfs/soft.mugqic/CentOS6/software/mugqic_pipelines/mugqic_pipelines-2.2.0/bfx/report/ChipSeq.homer_find_motifs_genome.md report/ && \
for contrast in A549_CTRL_BRD4 A549_CTRL_CDK9 A549_CTRL_MED1 A549_CTRL_NIPBL A549_CTRL_SMC1A A549_DEX_BRD4 A549_DEX_CDK9 A549_DEX_MED1 A549_DEX_NIPBL A549_DEX_SMC1A_rep1 A549_CTRL_MED1_rep1 A549_CTRL_NIPBL_rep1 A549_CTRL_SMC1A_rep1 A549_CTRL_MED1_rep2 A549_CTRL_NIPBL_rep2 A549_CTRL_SMC1A_rep2
do
  rsync -avP annotation/$contrast report/annotation/ && \
  echo -e "* [HOMER _De Novo_ Motif Results for Design $contrast](annotation/$contrast/$contrast/homerResults.html)
* [HOMER Known Motif Results for Design $contrast](annotation/$contrast/$contrast/knownResults.html)" \
  >> report/ChipSeq.homer_find_motifs_genome.md
done
homer_find_motifs_genome_report.9905966bef4fd299070c52342290bbfa.mugqic.done
)
homer_find_motifs_genome_17_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$homer_find_motifs_genome_17_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: annotation_graphs
#-------------------------------------------------------------------------------
STEP=annotation_graphs
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: annotation_graphs_1_JOB_ID: annotation_graphs
#-------------------------------------------------------------------------------
JOB_NAME=annotation_graphs
JOB_DEPENDENCIES=$homer_annotate_peaks_1_JOB_ID:$homer_annotate_peaks_2_JOB_ID:$homer_annotate_peaks_3_JOB_ID:$homer_annotate_peaks_4_JOB_ID:$homer_annotate_peaks_5_JOB_ID:$homer_annotate_peaks_6_JOB_ID:$homer_annotate_peaks_7_JOB_ID:$homer_annotate_peaks_8_JOB_ID:$homer_annotate_peaks_9_JOB_ID:$homer_annotate_peaks_10_JOB_ID:$homer_annotate_peaks_11_JOB_ID:$homer_annotate_peaks_12_JOB_ID:$homer_annotate_peaks_13_JOB_ID:$homer_annotate_peaks_14_JOB_ID:$homer_annotate_peaks_15_JOB_ID:$homer_annotate_peaks_16_JOB_ID
JOB_DONE=job_output/annotation_graphs/annotation_graphs.c5fda80d0136f850c16d34e658b91977.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'annotation_graphs.c5fda80d0136f850c16d34e658b91977.mugqic.done'
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
for contrast in A549_CTRL_BRD4 A549_CTRL_CDK9 A549_CTRL_MED1 A549_CTRL_NIPBL A549_CTRL_SMC1A A549_DEX_BRD4 A549_DEX_CDK9 A549_DEX_MED1 A549_DEX_NIPBL A549_DEX_SMC1A_rep1 A549_CTRL_MED1_rep1 A549_CTRL_NIPBL_rep1 A549_CTRL_SMC1A_rep1 A549_CTRL_MED1_rep2 A549_CTRL_NIPBL_rep2 A549_CTRL_SMC1A_rep2
do
  cp --parents graphs/${contrast}_Misc_Graphs.ps report/
  convert -rotate 90 graphs/${contrast}_Misc_Graphs.ps report/graphs/${contrast}_Misc_Graphs.png
  echo -e "----

![Annotation Statistics for Design $contrast ([download high-res image](graphs/${contrast}_Misc_Graphs.ps))](graphs/${contrast}_Misc_Graphs.png)
" \
  >> report/ChipSeq.annotation_graphs.md
done
annotation_graphs.c5fda80d0136f850c16d34e658b91977.mugqic.done
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
wget "http://mugqic.hpc.mcgill.ca/cgi-bin/pipeline.cgi?hostname=lg-1r17-n04&ip=10.241.129.14&pipeline=ChipSeq&steps=picard_sam_to_fastq,trimmomatic,merge_trimmomatic_stats,bwa_mem_picard_sort_sam,samtools_view_filter,picard_merge_sam_files,picard_mark_duplicates,metrics,homer_make_tag_directory,qc_metrics,homer_make_ucsc_file,macs2_callpeak,homer_annotate_peaks,homer_find_motifs_genome,annotation_graphs&samples=16" --quiet --output-document=/dev/null

