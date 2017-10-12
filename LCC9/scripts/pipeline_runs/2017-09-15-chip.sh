#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

#-------------------------------------------------------------------------------
# ChipSeq PBSScheduler Job Submission Bash script
# Version: 2.2.0
# Created on: 2017-09-15T18:08:19
# Steps:
#   picard_sam_to_fastq: 0 job... skipping
#   trimmomatic: 0 job... skipping
#   merge_trimmomatic_stats: 0 job... skipping
#   bwa_mem_picard_sort_sam: 0 job... skipping
#   samtools_view_filter: 0 job... skipping
#   picard_merge_sam_files: 0 job... skipping
#   picard_mark_duplicates: 0 job... skipping
#   metrics: 0 job... skipping
#   homer_make_tag_directory: 12 jobs
#   qc_metrics: 1 job
#   homer_make_ucsc_file: 13 jobs
#   macs2_callpeak: 11 jobs
#   homer_annotate_peaks: 11 jobs
#   homer_find_motifs_genome: 9 jobs
#   annotation_graphs: 1 job
#   TOTAL: 58 jobs
#-------------------------------------------------------------------------------

OUTPUT_DIR=/gs/scratch/efournier/sb_cofactor_hr/LCC9/output/chip-pipeline-GRCh38
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
# JOB: homer_make_tag_directory_1_JOB_ID: homer_make_tag_directory.LCC9_CTRL_MED1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.LCC9_CTRL_MED1
JOB_DEPENDENCIES=
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$homer_make_tag_directory_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_2_JOB_ID: homer_make_tag_directory.LCC9_E2_WCE
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.LCC9_E2_WCE
JOB_DEPENDENCIES=
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$homer_make_tag_directory_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_3_JOB_ID: homer_make_tag_directory.LCC9_E2_MED1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.LCC9_E2_MED1
JOB_DEPENDENCIES=
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$homer_make_tag_directory_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_4_JOB_ID: homer_make_tag_directory.LCC9_CTRL_WCE
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.LCC9_CTRL_WCE
JOB_DEPENDENCIES=
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$homer_make_tag_directory_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_5_JOB_ID: homer_make_tag_directory.LCC9_CTRL_SMC1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.LCC9_CTRL_SMC1
JOB_DEPENDENCIES=
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$homer_make_tag_directory_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_6_JOB_ID: homer_make_tag_directory.LCC9_E2_SMC1
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.LCC9_E2_SMC1
JOB_DEPENDENCIES=
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$homer_make_tag_directory_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_7_JOB_ID: homer_make_tag_directory.LCC9_E2_NIPBL
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.LCC9_E2_NIPBL
JOB_DEPENDENCIES=
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$homer_make_tag_directory_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_8_JOB_ID: homer_make_tag_directory.LCC9_CTRL_NIPBL
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.LCC9_CTRL_NIPBL
JOB_DEPENDENCIES=
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$homer_make_tag_directory_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_9_JOB_ID: homer_make_tag_directory.LCC9_CTRL_POL2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.LCC9_CTRL_POL2
JOB_DEPENDENCIES=
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$homer_make_tag_directory_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_10_JOB_ID: homer_make_tag_directory.LCC9_E2_POL2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.LCC9_E2_POL2
JOB_DEPENDENCIES=
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$homer_make_tag_directory_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_11_JOB_ID: homer_make_tag_directory.LCC9_CTRL_ERA
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.LCC9_CTRL_ERA
JOB_DEPENDENCIES=
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$homer_make_tag_directory_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_12_JOB_ID: homer_make_tag_directory.LCC9_E2_ERA
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.LCC9_E2_ERA
JOB_DEPENDENCIES=
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
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
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
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.LCC9_CTRL_MED1.1bb2ceb09a4a548959c6129ddb00d0fe.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.LCC9_CTRL_MED1.1bb2ceb09a4a548959c6129ddb00d0fe.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/LCC9_CTRL_MED1 && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2479938032.8 \
  --treatment \
  alignment/LCC9_CTRL_MED1/LCC9_CTRL_MED1.sorted.dup.bam \
  --control \
  alignment/LCC9_CTRL_WCE/LCC9_CTRL_WCE.sorted.dup.bam \
  --name peak_call/LCC9_CTRL_MED1/LCC9_CTRL_MED1 \
  >& peak_call/LCC9_CTRL_MED1/LCC9_CTRL_MED1.diag.macs.out
macs2_callpeak.LCC9_CTRL_MED1.1bb2ceb09a4a548959c6129ddb00d0fe.mugqic.done
)
macs2_callpeak_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$macs2_callpeak_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_2_JOB_ID: macs2_callpeak.LCC9_CTRL_NIPBL
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.LCC9_CTRL_NIPBL
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.LCC9_CTRL_NIPBL.33ad8664f4069c56242ccadabfa2ceaf.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.LCC9_CTRL_NIPBL.33ad8664f4069c56242ccadabfa2ceaf.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/LCC9_CTRL_NIPBL && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2479938032.8 \
  --treatment \
  alignment/LCC9_CTRL_NIPBL/LCC9_CTRL_NIPBL.sorted.dup.bam \
  --control \
  alignment/LCC9_CTRL_WCE/LCC9_CTRL_WCE.sorted.dup.bam \
  --name peak_call/LCC9_CTRL_NIPBL/LCC9_CTRL_NIPBL \
  >& peak_call/LCC9_CTRL_NIPBL/LCC9_CTRL_NIPBL.diag.macs.out
macs2_callpeak.LCC9_CTRL_NIPBL.33ad8664f4069c56242ccadabfa2ceaf.mugqic.done
)
macs2_callpeak_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$macs2_callpeak_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_3_JOB_ID: macs2_callpeak.LCC9_CTRL_SMC1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.LCC9_CTRL_SMC1
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.LCC9_CTRL_SMC1.8b89b1785f340c2924ba7c8cce1a904b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.LCC9_CTRL_SMC1.8b89b1785f340c2924ba7c8cce1a904b.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/LCC9_CTRL_SMC1 && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2479938032.8 \
  --treatment \
  alignment/LCC9_CTRL_SMC1/LCC9_CTRL_SMC1.sorted.dup.bam \
  --control \
  alignment/LCC9_CTRL_WCE/LCC9_CTRL_WCE.sorted.dup.bam \
  --name peak_call/LCC9_CTRL_SMC1/LCC9_CTRL_SMC1 \
  >& peak_call/LCC9_CTRL_SMC1/LCC9_CTRL_SMC1.diag.macs.out
macs2_callpeak.LCC9_CTRL_SMC1.8b89b1785f340c2924ba7c8cce1a904b.mugqic.done
)
macs2_callpeak_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$macs2_callpeak_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_4_JOB_ID: macs2_callpeak.LCC9_CTRL_POL2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.LCC9_CTRL_POL2
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.LCC9_CTRL_POL2.6ad259bb2e50b6fb74dcf28b43ec77e4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.LCC9_CTRL_POL2.6ad259bb2e50b6fb74dcf28b43ec77e4.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/LCC9_CTRL_POL2 && \
macs2 callpeak --format BAM --broad --nomodel \
  --gsize 2479938032.8 \
  --treatment \
  alignment/LCC9_CTRL_POL2/LCC9_CTRL_POL2.sorted.dup.bam \
  --control \
  alignment/LCC9_CTRL_WCE/LCC9_CTRL_WCE.sorted.dup.bam \
  --name peak_call/LCC9_CTRL_POL2/LCC9_CTRL_POL2 \
  >& peak_call/LCC9_CTRL_POL2/LCC9_CTRL_POL2.diag.macs.out
macs2_callpeak.LCC9_CTRL_POL2.6ad259bb2e50b6fb74dcf28b43ec77e4.mugqic.done
)
macs2_callpeak_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$macs2_callpeak_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_5_JOB_ID: macs2_callpeak.LCC9_CTRL_ERA
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.LCC9_CTRL_ERA
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.LCC9_CTRL_ERA.6afe33ef2434444c1c1eb7d6cfc116fa.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.LCC9_CTRL_ERA.6afe33ef2434444c1c1eb7d6cfc116fa.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/LCC9_CTRL_ERA && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2479938032.8 \
  --treatment \
  alignment/LCC9_CTRL_ERA/LCC9_CTRL_ERA.sorted.dup.bam \
  alignment/LCC9_E2_ERA/LCC9_E2_ERA.sorted.dup.bam \
  --control \
  alignment/LCC9_CTRL_WCE/LCC9_CTRL_WCE.sorted.dup.bam \
  --name peak_call/LCC9_CTRL_ERA/LCC9_CTRL_ERA \
  >& peak_call/LCC9_CTRL_ERA/LCC9_CTRL_ERA.diag.macs.out
macs2_callpeak.LCC9_CTRL_ERA.6afe33ef2434444c1c1eb7d6cfc116fa.mugqic.done
)
macs2_callpeak_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$macs2_callpeak_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_6_JOB_ID: macs2_callpeak.LCC9_E2_MED1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.LCC9_E2_MED1
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.LCC9_E2_MED1.5d0b7ee5db6c915a720134cfc1234fa9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.LCC9_E2_MED1.5d0b7ee5db6c915a720134cfc1234fa9.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/LCC9_E2_MED1 && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2479938032.8 \
  --treatment \
  alignment/LCC9_E2_MED1/LCC9_E2_MED1.sorted.dup.bam \
  --control \
  alignment/LCC9_E2_WCE/LCC9_E2_WCE.sorted.dup.bam \
  --name peak_call/LCC9_E2_MED1/LCC9_E2_MED1 \
  >& peak_call/LCC9_E2_MED1/LCC9_E2_MED1.diag.macs.out
macs2_callpeak.LCC9_E2_MED1.5d0b7ee5db6c915a720134cfc1234fa9.mugqic.done
)
macs2_callpeak_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$macs2_callpeak_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_7_JOB_ID: macs2_callpeak.LCC9_E2_NIPBL
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.LCC9_E2_NIPBL
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.LCC9_E2_NIPBL.1ad07ac3c1e0cff2e9b8764bf0ab83fc.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.LCC9_E2_NIPBL.1ad07ac3c1e0cff2e9b8764bf0ab83fc.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/LCC9_E2_NIPBL && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2479938032.8 \
  --treatment \
  alignment/LCC9_E2_NIPBL/LCC9_E2_NIPBL.sorted.dup.bam \
  --control \
  alignment/LCC9_E2_WCE/LCC9_E2_WCE.sorted.dup.bam \
  --name peak_call/LCC9_E2_NIPBL/LCC9_E2_NIPBL \
  >& peak_call/LCC9_E2_NIPBL/LCC9_E2_NIPBL.diag.macs.out
macs2_callpeak.LCC9_E2_NIPBL.1ad07ac3c1e0cff2e9b8764bf0ab83fc.mugqic.done
)
macs2_callpeak_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$macs2_callpeak_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_8_JOB_ID: macs2_callpeak.LCC9_E2_SMC1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.LCC9_E2_SMC1
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.LCC9_E2_SMC1.56493423bbd8ce18df6273b79db9e744.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.LCC9_E2_SMC1.56493423bbd8ce18df6273b79db9e744.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/LCC9_E2_SMC1 && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2479938032.8 \
  --treatment \
  alignment/LCC9_E2_SMC1/LCC9_E2_SMC1.sorted.dup.bam \
  --control \
  alignment/LCC9_E2_WCE/LCC9_E2_WCE.sorted.dup.bam \
  --name peak_call/LCC9_E2_SMC1/LCC9_E2_SMC1 \
  >& peak_call/LCC9_E2_SMC1/LCC9_E2_SMC1.diag.macs.out
macs2_callpeak.LCC9_E2_SMC1.56493423bbd8ce18df6273b79db9e744.mugqic.done
)
macs2_callpeak_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$macs2_callpeak_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_9_JOB_ID: macs2_callpeak.LCC9_E2_POL2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.LCC9_E2_POL2
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.LCC9_E2_POL2.3e716fb30204189a14db55f0f3707e0f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.LCC9_E2_POL2.3e716fb30204189a14db55f0f3707e0f.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/LCC9_E2_POL2 && \
macs2 callpeak --format BAM --broad --nomodel \
  --gsize 2479938032.8 \
  --treatment \
  alignment/LCC9_E2_POL2/LCC9_E2_POL2.sorted.dup.bam \
  --control \
  alignment/LCC9_E2_WCE/LCC9_E2_WCE.sorted.dup.bam \
  --name peak_call/LCC9_E2_POL2/LCC9_E2_POL2 \
  >& peak_call/LCC9_E2_POL2/LCC9_E2_POL2.diag.macs.out
macs2_callpeak.LCC9_E2_POL2.3e716fb30204189a14db55f0f3707e0f.mugqic.done
)
macs2_callpeak_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
echo "$macs2_callpeak_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_10_JOB_ID: macs2_callpeak.LCC9_E2_ERA
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.LCC9_E2_ERA
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.LCC9_E2_ERA.4b0242b84a5ef961cff9ab30c74cde4c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.LCC9_E2_ERA.4b0242b84a5ef961cff9ab30c74cde4c.mugqic.done'
module load mugqic/python/2.7.8 mugqic/MACS2/2.1.0.20140616 && \
mkdir -p peak_call/LCC9_E2_ERA && \
macs2 callpeak --format BAM --nomodel \
  --gsize 2479938032.8 \
  --treatment \
  alignment/LCC9_E2_ERA/LCC9_E2_ERA.sorted.dup.bam \
  --control \
  alignment/LCC9_E2_WCE/LCC9_E2_WCE.sorted.dup.bam \
  --name peak_call/LCC9_E2_ERA/LCC9_E2_ERA \
  >& peak_call/LCC9_E2_ERA/LCC9_E2_ERA.diag.macs.out
macs2_callpeak.LCC9_E2_ERA.4b0242b84a5ef961cff9ab30c74cde4c.mugqic.done
)
macs2_callpeak_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 | grep "[0-9]")
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
wget "http://mugqic.hpc.mcgill.ca/cgi-bin/pipeline.cgi?hostname=lg-1r14-n04&ip=10.241.129.4&pipeline=ChipSeq&steps=picard_sam_to_fastq,trimmomatic,merge_trimmomatic_stats,bwa_mem_picard_sort_sam,samtools_view_filter,picard_merge_sam_files,picard_mark_duplicates,metrics,homer_make_tag_directory,qc_metrics,homer_make_ucsc_file,macs2_callpeak,homer_annotate_peaks,homer_find_motifs_genome,annotation_graphs&samples=12" --quiet --output-document=/dev/null

