#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

#-------------------------------------------------------------------------------
# ChipSeq PBSScheduler Job Submission Bash script
# Version: 2.2.0
# Created on: 2017-08-30T18:30:31
# Steps:
#   picard_sam_to_fastq: 0 job... skipping
#   trimmomatic: 0 job... skipping
#   merge_trimmomatic_stats: 0 job... skipping
#   bwa_mem_picard_sort_sam: 0 job... skipping
#   samtools_view_filter: 0 job... skipping
#   picard_merge_sam_files: 0 job... skipping
#   picard_mark_duplicates: 0 job... skipping
#   metrics: 0 job... skipping
#   homer_make_tag_directory: 8 jobs
#   qc_metrics: 1 job
#   homer_make_ucsc_file: 9 jobs
#   macs2_callpeak: 0 job... skipping
#   homer_annotate_peaks: 9 jobs
#   homer_find_motifs_genome: 9 jobs
#   annotation_graphs: 1 job
#   TOTAL: 37 jobs
#-------------------------------------------------------------------------------

OUTPUT_DIR=/gs/scratch/efournier/sb_cofactor_hr/A549/output/chip-pipeline-GRCh38
JOB_OUTPUT_DIR=$OUTPUT_DIR/job_output
TIMESTAMP=`date +%FT%H.%M.%S`
JOB_LIST=$JOB_OUTPUT_DIR/ChipSeq_job_list_$TIMESTAMP
mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR


#-------------------------------------------------------------------------------
# STEP: homer_make_tag_directory
#-------------------------------------------------------------------------------
STEP=homer_make_tag_directory
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_1_JOB_ID: homer_make_tag_directory.A549_CTRL_BRD4_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.A549_CTRL_BRD4_rep1
JOB_DEPENDENCIES=
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$homer_make_tag_directory_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_2_JOB_ID: homer_make_tag_directory.A549_CTRL_CDK9_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.A549_CTRL_CDK9_rep1
JOB_DEPENDENCIES=
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$homer_make_tag_directory_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_3_JOB_ID: homer_make_tag_directory.A549_CTRL_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.A549_CTRL_MED1_rep1
JOB_DEPENDENCIES=
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$homer_make_tag_directory_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_4_JOB_ID: homer_make_tag_directory.A549_CTRL_MED1_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.A549_CTRL_MED1_rep2
JOB_DEPENDENCIES=
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$homer_make_tag_directory_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_5_JOB_ID: homer_make_tag_directory.A549_CTRL_SMC1A_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.A549_CTRL_SMC1A_rep1
JOB_DEPENDENCIES=
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$homer_make_tag_directory_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_6_JOB_ID: homer_make_tag_directory.A549_CTRL_SMC1A_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.A549_CTRL_SMC1A_rep2
JOB_DEPENDENCIES=
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$homer_make_tag_directory_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_7_JOB_ID: homer_make_tag_directory.A549_CTRL_WCE_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.A549_CTRL_WCE_rep1
JOB_DEPENDENCIES=
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$homer_make_tag_directory_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_8_JOB_ID: homer_make_tag_directory.A549_CTRL_WCE_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.A549_CTRL_WCE_rep2
JOB_DEPENDENCIES=
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
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
# STEP: homer_annotate_peaks
#-------------------------------------------------------------------------------
STEP=homer_annotate_peaks
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_1_JOB_ID: homer_annotate_peaks.A549_CTRL_BRD4
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.A549_CTRL_BRD4
JOB_DEPENDENCIES=
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m | grep "[0-9]")
echo "$homer_annotate_peaks_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_2_JOB_ID: homer_annotate_peaks.A549_CTRL_CDK9
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.A549_CTRL_CDK9
JOB_DEPENDENCIES=
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m | grep "[0-9]")
echo "$homer_annotate_peaks_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_3_JOB_ID: homer_annotate_peaks.A549_CTRL_MED1
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.A549_CTRL_MED1
JOB_DEPENDENCIES=
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m | grep "[0-9]")
echo "$homer_annotate_peaks_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_4_JOB_ID: homer_annotate_peaks.A549_CTRL_SMC1A
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.A549_CTRL_SMC1A
JOB_DEPENDENCIES=
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m | grep "[0-9]")
echo "$homer_annotate_peaks_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_5_JOB_ID: homer_annotate_peaks.A549_CTRL_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.A549_CTRL_MED1_rep1
JOB_DEPENDENCIES=
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m | grep "[0-9]")
echo "$homer_annotate_peaks_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_6_JOB_ID: homer_annotate_peaks.A549_CTRL_SMC1A_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.A549_CTRL_SMC1A_rep1
JOB_DEPENDENCIES=
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m | grep "[0-9]")
echo "$homer_annotate_peaks_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_7_JOB_ID: homer_annotate_peaks.A549_CTRL_MED1_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.A549_CTRL_MED1_rep2
JOB_DEPENDENCIES=
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m | grep "[0-9]")
echo "$homer_annotate_peaks_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_annotate_peaks_8_JOB_ID: homer_annotate_peaks.A549_CTRL_SMC1A_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_annotate_peaks.A549_CTRL_SMC1A_rep2
JOB_DEPENDENCIES=
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=2 -l pmem=2700m | grep "[0-9]")
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
JOB_DEPENDENCIES=
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 | grep "[0-9]")
echo "$homer_find_motifs_genome_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_2_JOB_ID: homer_find_motifs_genome.A549_CTRL_CDK9
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.A549_CTRL_CDK9
JOB_DEPENDENCIES=
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 | grep "[0-9]")
echo "$homer_find_motifs_genome_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_3_JOB_ID: homer_find_motifs_genome.A549_CTRL_MED1
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.A549_CTRL_MED1
JOB_DEPENDENCIES=
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 | grep "[0-9]")
echo "$homer_find_motifs_genome_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_4_JOB_ID: homer_find_motifs_genome.A549_CTRL_SMC1A
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.A549_CTRL_SMC1A
JOB_DEPENDENCIES=
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 | grep "[0-9]")
echo "$homer_find_motifs_genome_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_5_JOB_ID: homer_find_motifs_genome.A549_CTRL_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.A549_CTRL_MED1_rep1
JOB_DEPENDENCIES=
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 | grep "[0-9]")
echo "$homer_find_motifs_genome_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_6_JOB_ID: homer_find_motifs_genome.A549_CTRL_SMC1A_rep1
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.A549_CTRL_SMC1A_rep1
JOB_DEPENDENCIES=
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 | grep "[0-9]")
echo "$homer_find_motifs_genome_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_7_JOB_ID: homer_find_motifs_genome.A549_CTRL_MED1_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.A549_CTRL_MED1_rep2
JOB_DEPENDENCIES=
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 | grep "[0-9]")
echo "$homer_find_motifs_genome_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_find_motifs_genome_8_JOB_ID: homer_find_motifs_genome.A549_CTRL_SMC1A_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_find_motifs_genome.A549_CTRL_SMC1A_rep2
JOB_DEPENDENCIES=
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=4 | grep "[0-9]")
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
wget "http://mugqic.hpc.mcgill.ca/cgi-bin/pipeline.cgi?hostname=lg-1r14-n04&ip=10.241.129.4&pipeline=ChipSeq&steps=picard_sam_to_fastq,trimmomatic,merge_trimmomatic_stats,bwa_mem_picard_sort_sam,samtools_view_filter,picard_merge_sam_files,picard_mark_duplicates,metrics,homer_make_tag_directory,qc_metrics,homer_make_ucsc_file,macs2_callpeak,homer_annotate_peaks,homer_find_motifs_genome,annotation_graphs&samples=8" --quiet --output-document=/dev/null

