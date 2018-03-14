#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

#-------------------------------------------------------------------------------
# ChipSeq SlurmScheduler Job Submission Bash script
# Version: 3.0.1-beta
# Created on: 2018-03-13T19:42:41
# Steps:
#   macs2_callpeak: 49 jobs
#   TOTAL: 49 jobs
#-------------------------------------------------------------------------------

OUTPUT_DIR=/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38
JOB_OUTPUT_DIR=$OUTPUT_DIR/job_output
TIMESTAMP=`date +%FT%H.%M.%S`
JOB_LIST=$JOB_OUTPUT_DIR/ChipSeq_job_list_$TIMESTAMP
export CONFIG_FILES="/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini,input/chipseq.numpy.bug.ini"
mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR


#-------------------------------------------------------------------------------
# STEP: macs2_callpeak
#-------------------------------------------------------------------------------
STEP=macs2_callpeak
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_1_JOB_ID: macs2_callpeak.A549_CTRL_BRD4
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_CTRL_BRD4
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_CTRL_BRD4.4a5078ef8dcc92af38824e86f5c42861.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_CTRL_BRD4.4a5078ef8dcc92af38824e86f5c42861.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/A549_CTRL_BRD4 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/A549_CTRL_BRD4_rep1/A549_CTRL_BRD4_rep1.sorted.dup.bam \
  alignment/A549_CTRL_BRD4_rep2/A549_CTRL_BRD4_rep2.sorted.dup.bam \
  --control \
  alignment/A549_CTRL_WCE_rep1/A549_CTRL_WCE_rep1.sorted.dup.bam \
  alignment/A549_CTRL_WCE_rep2/A549_CTRL_WCE_rep2.sorted.dup.bam \
  alignment/A549_CTRL_WCE_rep3/A549_CTRL_WCE_rep3.sorted.dup.bam \
  --name peak_call/A549_CTRL_BRD4/A549_CTRL_BRD4 \
  >& peak_call/A549_CTRL_BRD4/A549_CTRL_BRD4.diag.macs.out
macs2_callpeak.A549_CTRL_BRD4.4a5078ef8dcc92af38824e86f5c42861.mugqic.done
)
macs2_callpeak_1_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE &&  $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_2_JOB_ID: macs2_callpeak_bigBed.A549_CTRL_BRD4
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.A549_CTRL_BRD4
JOB_DEPENDENCIES=$macs2_callpeak_1_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.A549_CTRL_BRD4.e9d70fb51085c807e53b26bffe02bc32.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.A549_CTRL_BRD4.e9d70fb51085c807e53b26bffe02bc32.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/A549_CTRL_BRD4/A549_CTRL_BRD4_peaks.narrowPeak > peak_call/A549_CTRL_BRD4/A549_CTRL_BRD4_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/A549_CTRL_BRD4/A549_CTRL_BRD4_peaks.narrowPeak.bed \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/A549_CTRL_BRD4/A549_CTRL_BRD4_peaks.narrowPeak.bb
macs2_callpeak_bigBed.A549_CTRL_BRD4.e9d70fb51085c807e53b26bffe02bc32.mugqic.done
)
macs2_callpeak_2_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE &&  $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_3_JOB_ID: macs2_callpeak.A549_CTRL_CDK9
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_CTRL_CDK9
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_CTRL_CDK9.a1c79bf2a648ace6a7594403b8615a68.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_CTRL_CDK9.a1c79bf2a648ace6a7594403b8615a68.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/A549_CTRL_CDK9 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/A549_CTRL_CDK9_rep1/A549_CTRL_CDK9_rep1.sorted.dup.bam \
  alignment/A549_CTRL_CDK9_rep2/A549_CTRL_CDK9_rep2.sorted.dup.bam \
  --control \
  alignment/A549_CTRL_WCE_rep1/A549_CTRL_WCE_rep1.sorted.dup.bam \
  alignment/A549_CTRL_WCE_rep2/A549_CTRL_WCE_rep2.sorted.dup.bam \
  alignment/A549_CTRL_WCE_rep3/A549_CTRL_WCE_rep3.sorted.dup.bam \
  --name peak_call/A549_CTRL_CDK9/A549_CTRL_CDK9 \
  >& peak_call/A549_CTRL_CDK9/A549_CTRL_CDK9.diag.macs.out
macs2_callpeak.A549_CTRL_CDK9.a1c79bf2a648ace6a7594403b8615a68.mugqic.done
)
macs2_callpeak_3_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE &&  $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_4_JOB_ID: macs2_callpeak_bigBed.A549_CTRL_CDK9
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.A549_CTRL_CDK9
JOB_DEPENDENCIES=$macs2_callpeak_3_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.A549_CTRL_CDK9.58fc413bc7ab3c17392537d847754d73.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.A549_CTRL_CDK9.58fc413bc7ab3c17392537d847754d73.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/A549_CTRL_CDK9/A549_CTRL_CDK9_peaks.narrowPeak > peak_call/A549_CTRL_CDK9/A549_CTRL_CDK9_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/A549_CTRL_CDK9/A549_CTRL_CDK9_peaks.narrowPeak.bed \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/A549_CTRL_CDK9/A549_CTRL_CDK9_peaks.narrowPeak.bb
macs2_callpeak_bigBed.A549_CTRL_CDK9.58fc413bc7ab3c17392537d847754d73.mugqic.done
)
macs2_callpeak_4_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE &&  $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_5_JOB_ID: macs2_callpeak.A549_CTRL_MED1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_CTRL_MED1
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_CTRL_MED1.a907f53bff399e03691d20d8ba8233b5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_CTRL_MED1.a907f53bff399e03691d20d8ba8233b5.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/A549_CTRL_MED1 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/A549_CTRL_MED1_rep1/A549_CTRL_MED1_rep1.sorted.dup.bam \
  alignment/A549_CTRL_MED1_rep2/A549_CTRL_MED1_rep2.sorted.dup.bam \
  --control \
  alignment/A549_CTRL_WCE_rep1/A549_CTRL_WCE_rep1.sorted.dup.bam \
  alignment/A549_CTRL_WCE_rep2/A549_CTRL_WCE_rep2.sorted.dup.bam \
  alignment/A549_CTRL_WCE_rep3/A549_CTRL_WCE_rep3.sorted.dup.bam \
  --name peak_call/A549_CTRL_MED1/A549_CTRL_MED1 \
  >& peak_call/A549_CTRL_MED1/A549_CTRL_MED1.diag.macs.out
macs2_callpeak.A549_CTRL_MED1.a907f53bff399e03691d20d8ba8233b5.mugqic.done
)
macs2_callpeak_5_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE &&  $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_6_JOB_ID: macs2_callpeak_bigBed.A549_CTRL_MED1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.A549_CTRL_MED1
JOB_DEPENDENCIES=$macs2_callpeak_5_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.A549_CTRL_MED1.3d902f3b4b055db59fd6c26b6c062553.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.A549_CTRL_MED1.3d902f3b4b055db59fd6c26b6c062553.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/A549_CTRL_MED1/A549_CTRL_MED1_peaks.narrowPeak > peak_call/A549_CTRL_MED1/A549_CTRL_MED1_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/A549_CTRL_MED1/A549_CTRL_MED1_peaks.narrowPeak.bed \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/A549_CTRL_MED1/A549_CTRL_MED1_peaks.narrowPeak.bb
macs2_callpeak_bigBed.A549_CTRL_MED1.3d902f3b4b055db59fd6c26b6c062553.mugqic.done
)
macs2_callpeak_6_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE &&  $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_7_JOB_ID: macs2_callpeak.A549_CTRL_NIPBL
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_CTRL_NIPBL
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_CTRL_NIPBL.51cb8030b6878e4eec9a5dcfdd758694.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_CTRL_NIPBL.51cb8030b6878e4eec9a5dcfdd758694.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/A549_CTRL_NIPBL && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/A549_CTRL_NIPBL_rep1/A549_CTRL_NIPBL_rep1.sorted.dup.bam \
  alignment/A549_CTRL_NIPBL_rep2/A549_CTRL_NIPBL_rep2.sorted.dup.bam \
  --control \
  alignment/A549_CTRL_WCE_rep1/A549_CTRL_WCE_rep1.sorted.dup.bam \
  alignment/A549_CTRL_WCE_rep2/A549_CTRL_WCE_rep2.sorted.dup.bam \
  alignment/A549_CTRL_WCE_rep3/A549_CTRL_WCE_rep3.sorted.dup.bam \
  --name peak_call/A549_CTRL_NIPBL/A549_CTRL_NIPBL \
  >& peak_call/A549_CTRL_NIPBL/A549_CTRL_NIPBL.diag.macs.out
macs2_callpeak.A549_CTRL_NIPBL.51cb8030b6878e4eec9a5dcfdd758694.mugqic.done
)
macs2_callpeak_7_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE &&  $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_8_JOB_ID: macs2_callpeak_bigBed.A549_CTRL_NIPBL
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.A549_CTRL_NIPBL
JOB_DEPENDENCIES=$macs2_callpeak_7_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.A549_CTRL_NIPBL.24e32f61027bd077ecf2d10aa62e6fba.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.A549_CTRL_NIPBL.24e32f61027bd077ecf2d10aa62e6fba.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/A549_CTRL_NIPBL/A549_CTRL_NIPBL_peaks.narrowPeak > peak_call/A549_CTRL_NIPBL/A549_CTRL_NIPBL_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/A549_CTRL_NIPBL/A549_CTRL_NIPBL_peaks.narrowPeak.bed \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/A549_CTRL_NIPBL/A549_CTRL_NIPBL_peaks.narrowPeak.bb
macs2_callpeak_bigBed.A549_CTRL_NIPBL.24e32f61027bd077ecf2d10aa62e6fba.mugqic.done
)
macs2_callpeak_8_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE &&  $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_9_JOB_ID: macs2_callpeak.A549_CTRL_SMC1A
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_CTRL_SMC1A
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_CTRL_SMC1A.4dc505f26df06a74e7a97d918efabbcf.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_CTRL_SMC1A.4dc505f26df06a74e7a97d918efabbcf.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/A549_CTRL_SMC1A && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/A549_CTRL_SMC1A_rep1/A549_CTRL_SMC1A_rep1.sorted.dup.bam \
  alignment/A549_CTRL_SMC1A_rep2/A549_CTRL_SMC1A_rep2.sorted.dup.bam \
  --control \
  alignment/A549_CTRL_WCE_rep1/A549_CTRL_WCE_rep1.sorted.dup.bam \
  alignment/A549_CTRL_WCE_rep2/A549_CTRL_WCE_rep2.sorted.dup.bam \
  alignment/A549_CTRL_WCE_rep3/A549_CTRL_WCE_rep3.sorted.dup.bam \
  --name peak_call/A549_CTRL_SMC1A/A549_CTRL_SMC1A \
  >& peak_call/A549_CTRL_SMC1A/A549_CTRL_SMC1A.diag.macs.out
macs2_callpeak.A549_CTRL_SMC1A.4dc505f26df06a74e7a97d918efabbcf.mugqic.done
)
macs2_callpeak_9_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE &&  $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_10_JOB_ID: macs2_callpeak_bigBed.A549_CTRL_SMC1A
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.A549_CTRL_SMC1A
JOB_DEPENDENCIES=$macs2_callpeak_9_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.A549_CTRL_SMC1A.2baebbf8dca0c60c1536501694f19c5d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.A549_CTRL_SMC1A.2baebbf8dca0c60c1536501694f19c5d.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/A549_CTRL_SMC1A/A549_CTRL_SMC1A_peaks.narrowPeak > peak_call/A549_CTRL_SMC1A/A549_CTRL_SMC1A_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/A549_CTRL_SMC1A/A549_CTRL_SMC1A_peaks.narrowPeak.bed \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/A549_CTRL_SMC1A/A549_CTRL_SMC1A_peaks.narrowPeak.bb
macs2_callpeak_bigBed.A549_CTRL_SMC1A.2baebbf8dca0c60c1536501694f19c5d.mugqic.done
)
macs2_callpeak_10_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE &&  $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_11_JOB_ID: macs2_callpeak.A549_DEX_BRD4
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_DEX_BRD4
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_DEX_BRD4.d649030be42e7aa63db0e095e4e102c8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_DEX_BRD4.d649030be42e7aa63db0e095e4e102c8.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/A549_DEX_BRD4 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/A549_DEX_BRD4_rep1/A549_DEX_BRD4_rep1.sorted.dup.bam \
  alignment/A549_DEX_BRD4_rep2/A549_DEX_BRD4_rep2.sorted.dup.bam \
  --control \
  alignment/A549_DEX_WCE_rep1/A549_DEX_WCE_rep1.sorted.dup.bam \
  alignment/A549_DEX_WCE_rep2/A549_DEX_WCE_rep2.sorted.dup.bam \
  --name peak_call/A549_DEX_BRD4/A549_DEX_BRD4 \
  >& peak_call/A549_DEX_BRD4/A549_DEX_BRD4.diag.macs.out
macs2_callpeak.A549_DEX_BRD4.d649030be42e7aa63db0e095e4e102c8.mugqic.done
)
macs2_callpeak_11_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE &&  $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_12_JOB_ID: macs2_callpeak_bigBed.A549_DEX_BRD4
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.A549_DEX_BRD4
JOB_DEPENDENCIES=$macs2_callpeak_11_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.A549_DEX_BRD4.86a647ecd2bbea4738f39d87d93ffbf3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.A549_DEX_BRD4.86a647ecd2bbea4738f39d87d93ffbf3.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/A549_DEX_BRD4/A549_DEX_BRD4_peaks.narrowPeak > peak_call/A549_DEX_BRD4/A549_DEX_BRD4_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/A549_DEX_BRD4/A549_DEX_BRD4_peaks.narrowPeak.bed \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/A549_DEX_BRD4/A549_DEX_BRD4_peaks.narrowPeak.bb
macs2_callpeak_bigBed.A549_DEX_BRD4.86a647ecd2bbea4738f39d87d93ffbf3.mugqic.done
)
macs2_callpeak_12_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE &&  $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_13_JOB_ID: macs2_callpeak.A549_DEX_CDK9
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_DEX_CDK9
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_DEX_CDK9.e436040559f91de0b1e95d46c85f51d8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_DEX_CDK9.e436040559f91de0b1e95d46c85f51d8.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/A549_DEX_CDK9 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/A549_DEX_CDK9_rep1/A549_DEX_CDK9_rep1.sorted.dup.bam \
  alignment/A549_DEX_CDK9_rep2/A549_DEX_CDK9_rep2.sorted.dup.bam \
  --control \
  alignment/A549_DEX_WCE_rep1/A549_DEX_WCE_rep1.sorted.dup.bam \
  alignment/A549_DEX_WCE_rep2/A549_DEX_WCE_rep2.sorted.dup.bam \
  --name peak_call/A549_DEX_CDK9/A549_DEX_CDK9 \
  >& peak_call/A549_DEX_CDK9/A549_DEX_CDK9.diag.macs.out
macs2_callpeak.A549_DEX_CDK9.e436040559f91de0b1e95d46c85f51d8.mugqic.done
)
macs2_callpeak_13_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE &&  $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_13_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_14_JOB_ID: macs2_callpeak_bigBed.A549_DEX_CDK9
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.A549_DEX_CDK9
JOB_DEPENDENCIES=$macs2_callpeak_13_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.A549_DEX_CDK9.4477947b8f7ee67a583182dce604bc7d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.A549_DEX_CDK9.4477947b8f7ee67a583182dce604bc7d.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/A549_DEX_CDK9/A549_DEX_CDK9_peaks.narrowPeak > peak_call/A549_DEX_CDK9/A549_DEX_CDK9_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/A549_DEX_CDK9/A549_DEX_CDK9_peaks.narrowPeak.bed \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/A549_DEX_CDK9/A549_DEX_CDK9_peaks.narrowPeak.bb
macs2_callpeak_bigBed.A549_DEX_CDK9.4477947b8f7ee67a583182dce604bc7d.mugqic.done
)
macs2_callpeak_14_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE &&  $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_14_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_15_JOB_ID: macs2_callpeak.A549_DEX_MED1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_DEX_MED1
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_DEX_MED1.4c2cafe33615580105a8b0b725292e60.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_DEX_MED1.4c2cafe33615580105a8b0b725292e60.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/A549_DEX_MED1 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/A549_DEX_MED1_rep1/A549_DEX_MED1_rep1.sorted.dup.bam \
  --control \
  alignment/A549_DEX_WCE_rep1/A549_DEX_WCE_rep1.sorted.dup.bam \
  alignment/A549_DEX_WCE_rep2/A549_DEX_WCE_rep2.sorted.dup.bam \
  --name peak_call/A549_DEX_MED1/A549_DEX_MED1 \
  >& peak_call/A549_DEX_MED1/A549_DEX_MED1.diag.macs.out
macs2_callpeak.A549_DEX_MED1.4c2cafe33615580105a8b0b725292e60.mugqic.done
)
macs2_callpeak_15_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE &&  $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_15_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_16_JOB_ID: macs2_callpeak_bigBed.A549_DEX_MED1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.A549_DEX_MED1
JOB_DEPENDENCIES=$macs2_callpeak_15_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.A549_DEX_MED1.ae3041fa08e626bd0a912ed2f1163a8d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.A549_DEX_MED1.ae3041fa08e626bd0a912ed2f1163a8d.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/A549_DEX_MED1/A549_DEX_MED1_peaks.narrowPeak > peak_call/A549_DEX_MED1/A549_DEX_MED1_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/A549_DEX_MED1/A549_DEX_MED1_peaks.narrowPeak.bed \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/A549_DEX_MED1/A549_DEX_MED1_peaks.narrowPeak.bb
macs2_callpeak_bigBed.A549_DEX_MED1.ae3041fa08e626bd0a912ed2f1163a8d.mugqic.done
)
macs2_callpeak_16_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE &&  $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_16_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_17_JOB_ID: macs2_callpeak.A549_DEX_NIPBL
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_DEX_NIPBL
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_DEX_NIPBL.ae2f494b5b1b83f51afd3b4e7e556cf8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_DEX_NIPBL.ae2f494b5b1b83f51afd3b4e7e556cf8.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/A549_DEX_NIPBL && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/A549_DEX_NIPBL_rep1/A549_DEX_NIPBL_rep1.sorted.dup.bam \
  --control \
  alignment/A549_DEX_WCE_rep1/A549_DEX_WCE_rep1.sorted.dup.bam \
  alignment/A549_DEX_WCE_rep2/A549_DEX_WCE_rep2.sorted.dup.bam \
  --name peak_call/A549_DEX_NIPBL/A549_DEX_NIPBL \
  >& peak_call/A549_DEX_NIPBL/A549_DEX_NIPBL.diag.macs.out
macs2_callpeak.A549_DEX_NIPBL.ae2f494b5b1b83f51afd3b4e7e556cf8.mugqic.done
)
macs2_callpeak_17_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE &&  $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_17_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_18_JOB_ID: macs2_callpeak_bigBed.A549_DEX_NIPBL
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.A549_DEX_NIPBL
JOB_DEPENDENCIES=$macs2_callpeak_17_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.A549_DEX_NIPBL.ce3ba7f3d4a6d05d29e74fcb1ed8d27f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.A549_DEX_NIPBL.ce3ba7f3d4a6d05d29e74fcb1ed8d27f.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/A549_DEX_NIPBL/A549_DEX_NIPBL_peaks.narrowPeak > peak_call/A549_DEX_NIPBL/A549_DEX_NIPBL_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/A549_DEX_NIPBL/A549_DEX_NIPBL_peaks.narrowPeak.bed \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/A549_DEX_NIPBL/A549_DEX_NIPBL_peaks.narrowPeak.bb
macs2_callpeak_bigBed.A549_DEX_NIPBL.ce3ba7f3d4a6d05d29e74fcb1ed8d27f.mugqic.done
)
macs2_callpeak_18_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE &&  $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_18_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_19_JOB_ID: macs2_callpeak.A549_DEX_SMC1A
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_DEX_SMC1A
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_DEX_SMC1A.5cb813179828fabcd92c4c70aa5bf896.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_DEX_SMC1A.5cb813179828fabcd92c4c70aa5bf896.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/A549_DEX_SMC1A && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/A549_DEX_SMC1A_rep1/A549_DEX_SMC1A_rep1.sorted.dup.bam \
  --control \
  alignment/A549_DEX_WCE_rep1/A549_DEX_WCE_rep1.sorted.dup.bam \
  alignment/A549_DEX_WCE_rep2/A549_DEX_WCE_rep2.sorted.dup.bam \
  --name peak_call/A549_DEX_SMC1A/A549_DEX_SMC1A \
  >& peak_call/A549_DEX_SMC1A/A549_DEX_SMC1A.diag.macs.out
macs2_callpeak.A549_DEX_SMC1A.5cb813179828fabcd92c4c70aa5bf896.mugqic.done
)
macs2_callpeak_19_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE &&  $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_19_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_20_JOB_ID: macs2_callpeak_bigBed.A549_DEX_SMC1A
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.A549_DEX_SMC1A
JOB_DEPENDENCIES=$macs2_callpeak_19_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.A549_DEX_SMC1A.70162506358a513207210409c88aafe6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.A549_DEX_SMC1A.70162506358a513207210409c88aafe6.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/A549_DEX_SMC1A/A549_DEX_SMC1A_peaks.narrowPeak > peak_call/A549_DEX_SMC1A/A549_DEX_SMC1A_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/A549_DEX_SMC1A/A549_DEX_SMC1A_peaks.narrowPeak.bed \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/A549_DEX_SMC1A/A549_DEX_SMC1A_peaks.narrowPeak.bb
macs2_callpeak_bigBed.A549_DEX_SMC1A.70162506358a513207210409c88aafe6.mugqic.done
)
macs2_callpeak_20_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE &&  $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_20_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_21_JOB_ID: macs2_callpeak.A549_CTRL_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_CTRL_MED1_rep1
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_CTRL_MED1_rep1.0549f1498ae09c43d1236d4fd3e6421a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_CTRL_MED1_rep1.0549f1498ae09c43d1236d4fd3e6421a.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/A549_CTRL_MED1_rep1 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/A549_CTRL_MED1_rep1/A549_CTRL_MED1_rep1.sorted.dup.bam \
  --control \
  alignment/A549_CTRL_WCE_rep1/A549_CTRL_WCE_rep1.sorted.dup.bam \
  alignment/A549_CTRL_WCE_rep2/A549_CTRL_WCE_rep2.sorted.dup.bam \
  --name peak_call/A549_CTRL_MED1_rep1/A549_CTRL_MED1_rep1 \
  >& peak_call/A549_CTRL_MED1_rep1/A549_CTRL_MED1_rep1.diag.macs.out
macs2_callpeak.A549_CTRL_MED1_rep1.0549f1498ae09c43d1236d4fd3e6421a.mugqic.done
)
macs2_callpeak_21_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE &&  $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_21_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_22_JOB_ID: macs2_callpeak_bigBed.A549_CTRL_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.A549_CTRL_MED1_rep1
JOB_DEPENDENCIES=$macs2_callpeak_21_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.A549_CTRL_MED1_rep1.b8f6e6ae25afe3bff39c616e618c65f5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.A549_CTRL_MED1_rep1.b8f6e6ae25afe3bff39c616e618c65f5.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/A549_CTRL_MED1_rep1/A549_CTRL_MED1_rep1_peaks.narrowPeak > peak_call/A549_CTRL_MED1_rep1/A549_CTRL_MED1_rep1_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/A549_CTRL_MED1_rep1/A549_CTRL_MED1_rep1_peaks.narrowPeak.bed \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/A549_CTRL_MED1_rep1/A549_CTRL_MED1_rep1_peaks.narrowPeak.bb
macs2_callpeak_bigBed.A549_CTRL_MED1_rep1.b8f6e6ae25afe3bff39c616e618c65f5.mugqic.done
)
macs2_callpeak_22_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE &&  $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_22_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_23_JOB_ID: macs2_callpeak.A549_CTRL_NIPBL_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_CTRL_NIPBL_rep1
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_CTRL_NIPBL_rep1.36c22c3363a1c8751aac406e6ce48766.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_CTRL_NIPBL_rep1.36c22c3363a1c8751aac406e6ce48766.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/A549_CTRL_NIPBL_rep1 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/A549_CTRL_NIPBL_rep1/A549_CTRL_NIPBL_rep1.sorted.dup.bam \
  --control \
  alignment/A549_CTRL_WCE_rep1/A549_CTRL_WCE_rep1.sorted.dup.bam \
  alignment/A549_CTRL_WCE_rep2/A549_CTRL_WCE_rep2.sorted.dup.bam \
  --name peak_call/A549_CTRL_NIPBL_rep1/A549_CTRL_NIPBL_rep1 \
  >& peak_call/A549_CTRL_NIPBL_rep1/A549_CTRL_NIPBL_rep1.diag.macs.out
macs2_callpeak.A549_CTRL_NIPBL_rep1.36c22c3363a1c8751aac406e6ce48766.mugqic.done
)
macs2_callpeak_23_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE &&  $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_23_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_24_JOB_ID: macs2_callpeak_bigBed.A549_CTRL_NIPBL_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.A549_CTRL_NIPBL_rep1
JOB_DEPENDENCIES=$macs2_callpeak_23_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.A549_CTRL_NIPBL_rep1.2ad6c1a9050754bf94b494e19a80b940.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.A549_CTRL_NIPBL_rep1.2ad6c1a9050754bf94b494e19a80b940.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/A549_CTRL_NIPBL_rep1/A549_CTRL_NIPBL_rep1_peaks.narrowPeak > peak_call/A549_CTRL_NIPBL_rep1/A549_CTRL_NIPBL_rep1_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/A549_CTRL_NIPBL_rep1/A549_CTRL_NIPBL_rep1_peaks.narrowPeak.bed \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/A549_CTRL_NIPBL_rep1/A549_CTRL_NIPBL_rep1_peaks.narrowPeak.bb
macs2_callpeak_bigBed.A549_CTRL_NIPBL_rep1.2ad6c1a9050754bf94b494e19a80b940.mugqic.done
)
macs2_callpeak_24_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE &&  $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_24_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_25_JOB_ID: macs2_callpeak.A549_CTRL_SMC1A_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_CTRL_SMC1A_rep1
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_CTRL_SMC1A_rep1.61787260f5d087cd990c700bf9c1ed3b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_CTRL_SMC1A_rep1.61787260f5d087cd990c700bf9c1ed3b.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/A549_CTRL_SMC1A_rep1 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/A549_CTRL_SMC1A_rep1/A549_CTRL_SMC1A_rep1.sorted.dup.bam \
  --control \
  alignment/A549_CTRL_WCE_rep1/A549_CTRL_WCE_rep1.sorted.dup.bam \
  alignment/A549_CTRL_WCE_rep2/A549_CTRL_WCE_rep2.sorted.dup.bam \
  --name peak_call/A549_CTRL_SMC1A_rep1/A549_CTRL_SMC1A_rep1 \
  >& peak_call/A549_CTRL_SMC1A_rep1/A549_CTRL_SMC1A_rep1.diag.macs.out
macs2_callpeak.A549_CTRL_SMC1A_rep1.61787260f5d087cd990c700bf9c1ed3b.mugqic.done
)
macs2_callpeak_25_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE &&  $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_25_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_26_JOB_ID: macs2_callpeak_bigBed.A549_CTRL_SMC1A_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.A549_CTRL_SMC1A_rep1
JOB_DEPENDENCIES=$macs2_callpeak_25_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.A549_CTRL_SMC1A_rep1.67a2fedaebaabda32a884087783dd23a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.A549_CTRL_SMC1A_rep1.67a2fedaebaabda32a884087783dd23a.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/A549_CTRL_SMC1A_rep1/A549_CTRL_SMC1A_rep1_peaks.narrowPeak > peak_call/A549_CTRL_SMC1A_rep1/A549_CTRL_SMC1A_rep1_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/A549_CTRL_SMC1A_rep1/A549_CTRL_SMC1A_rep1_peaks.narrowPeak.bed \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/A549_CTRL_SMC1A_rep1/A549_CTRL_SMC1A_rep1_peaks.narrowPeak.bb
macs2_callpeak_bigBed.A549_CTRL_SMC1A_rep1.67a2fedaebaabda32a884087783dd23a.mugqic.done
)
macs2_callpeak_26_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE &&  $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_26_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_27_JOB_ID: macs2_callpeak.A549_CTRL_MED1_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_CTRL_MED1_rep2
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_CTRL_MED1_rep2.585e3b1cb68a806d07b58252020db52d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_CTRL_MED1_rep2.585e3b1cb68a806d07b58252020db52d.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/A549_CTRL_MED1_rep2 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/A549_CTRL_MED1_rep2/A549_CTRL_MED1_rep2.sorted.dup.bam \
  --control \
  alignment/A549_CTRL_WCE_rep1/A549_CTRL_WCE_rep1.sorted.dup.bam \
  alignment/A549_CTRL_WCE_rep2/A549_CTRL_WCE_rep2.sorted.dup.bam \
  --name peak_call/A549_CTRL_MED1_rep2/A549_CTRL_MED1_rep2 \
  >& peak_call/A549_CTRL_MED1_rep2/A549_CTRL_MED1_rep2.diag.macs.out
macs2_callpeak.A549_CTRL_MED1_rep2.585e3b1cb68a806d07b58252020db52d.mugqic.done
)
macs2_callpeak_27_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE &&  $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_27_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_28_JOB_ID: macs2_callpeak_bigBed.A549_CTRL_MED1_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.A549_CTRL_MED1_rep2
JOB_DEPENDENCIES=$macs2_callpeak_27_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.A549_CTRL_MED1_rep2.6483961365bbb5fa0bce7d3309386f4b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.A549_CTRL_MED1_rep2.6483961365bbb5fa0bce7d3309386f4b.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/A549_CTRL_MED1_rep2/A549_CTRL_MED1_rep2_peaks.narrowPeak > peak_call/A549_CTRL_MED1_rep2/A549_CTRL_MED1_rep2_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/A549_CTRL_MED1_rep2/A549_CTRL_MED1_rep2_peaks.narrowPeak.bed \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/A549_CTRL_MED1_rep2/A549_CTRL_MED1_rep2_peaks.narrowPeak.bb
macs2_callpeak_bigBed.A549_CTRL_MED1_rep2.6483961365bbb5fa0bce7d3309386f4b.mugqic.done
)
macs2_callpeak_28_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE &&  $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_28_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_29_JOB_ID: macs2_callpeak.A549_CTRL_NIPBL_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_CTRL_NIPBL_rep2
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_CTRL_NIPBL_rep2.6a4f07ea48bcf562f192d08db803886b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_CTRL_NIPBL_rep2.6a4f07ea48bcf562f192d08db803886b.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/A549_CTRL_NIPBL_rep2 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/A549_CTRL_NIPBL_rep2/A549_CTRL_NIPBL_rep2.sorted.dup.bam \
  --control \
  alignment/A549_CTRL_WCE_rep1/A549_CTRL_WCE_rep1.sorted.dup.bam \
  alignment/A549_CTRL_WCE_rep2/A549_CTRL_WCE_rep2.sorted.dup.bam \
  --name peak_call/A549_CTRL_NIPBL_rep2/A549_CTRL_NIPBL_rep2 \
  >& peak_call/A549_CTRL_NIPBL_rep2/A549_CTRL_NIPBL_rep2.diag.macs.out
macs2_callpeak.A549_CTRL_NIPBL_rep2.6a4f07ea48bcf562f192d08db803886b.mugqic.done
)
macs2_callpeak_29_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE &&  $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_29_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_30_JOB_ID: macs2_callpeak_bigBed.A549_CTRL_NIPBL_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.A549_CTRL_NIPBL_rep2
JOB_DEPENDENCIES=$macs2_callpeak_29_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.A549_CTRL_NIPBL_rep2.679e611c64d62d2f5cd482c36f299fda.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.A549_CTRL_NIPBL_rep2.679e611c64d62d2f5cd482c36f299fda.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/A549_CTRL_NIPBL_rep2/A549_CTRL_NIPBL_rep2_peaks.narrowPeak > peak_call/A549_CTRL_NIPBL_rep2/A549_CTRL_NIPBL_rep2_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/A549_CTRL_NIPBL_rep2/A549_CTRL_NIPBL_rep2_peaks.narrowPeak.bed \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/A549_CTRL_NIPBL_rep2/A549_CTRL_NIPBL_rep2_peaks.narrowPeak.bb
macs2_callpeak_bigBed.A549_CTRL_NIPBL_rep2.679e611c64d62d2f5cd482c36f299fda.mugqic.done
)
macs2_callpeak_30_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE &&  $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_30_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_31_JOB_ID: macs2_callpeak.A549_CTRL_SMC1A_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_CTRL_SMC1A_rep2
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_CTRL_SMC1A_rep2.b934405ad28fe18500623b7af48a2c50.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_CTRL_SMC1A_rep2.b934405ad28fe18500623b7af48a2c50.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/A549_CTRL_SMC1A_rep2 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/A549_CTRL_SMC1A_rep2/A549_CTRL_SMC1A_rep2.sorted.dup.bam \
  --control \
  alignment/A549_CTRL_WCE_rep1/A549_CTRL_WCE_rep1.sorted.dup.bam \
  alignment/A549_CTRL_WCE_rep2/A549_CTRL_WCE_rep2.sorted.dup.bam \
  --name peak_call/A549_CTRL_SMC1A_rep2/A549_CTRL_SMC1A_rep2 \
  >& peak_call/A549_CTRL_SMC1A_rep2/A549_CTRL_SMC1A_rep2.diag.macs.out
macs2_callpeak.A549_CTRL_SMC1A_rep2.b934405ad28fe18500623b7af48a2c50.mugqic.done
)
macs2_callpeak_31_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE &&  $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_31_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_32_JOB_ID: macs2_callpeak_bigBed.A549_CTRL_SMC1A_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.A549_CTRL_SMC1A_rep2
JOB_DEPENDENCIES=$macs2_callpeak_31_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.A549_CTRL_SMC1A_rep2.71f2b4437dc04324b5014056528cbe9d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.A549_CTRL_SMC1A_rep2.71f2b4437dc04324b5014056528cbe9d.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/A549_CTRL_SMC1A_rep2/A549_CTRL_SMC1A_rep2_peaks.narrowPeak > peak_call/A549_CTRL_SMC1A_rep2/A549_CTRL_SMC1A_rep2_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/A549_CTRL_SMC1A_rep2/A549_CTRL_SMC1A_rep2_peaks.narrowPeak.bed \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/A549_CTRL_SMC1A_rep2/A549_CTRL_SMC1A_rep2_peaks.narrowPeak.bb
macs2_callpeak_bigBed.A549_CTRL_SMC1A_rep2.71f2b4437dc04324b5014056528cbe9d.mugqic.done
)
macs2_callpeak_32_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE &&  $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_32_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_33_JOB_ID: macs2_callpeak.A549_CTRL_BRD4_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_CTRL_BRD4_rep1
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_CTRL_BRD4_rep1.7896ea17b4595ac6d13192f7b9e93089.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_CTRL_BRD4_rep1.7896ea17b4595ac6d13192f7b9e93089.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/A549_CTRL_BRD4_rep1 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/A549_CTRL_BRD4_rep1/A549_CTRL_BRD4_rep1.sorted.dup.bam \
  --control \
  alignment/A549_CTRL_WCE_rep3/A549_CTRL_WCE_rep3.sorted.dup.bam \
  --name peak_call/A549_CTRL_BRD4_rep1/A549_CTRL_BRD4_rep1 \
  >& peak_call/A549_CTRL_BRD4_rep1/A549_CTRL_BRD4_rep1.diag.macs.out
macs2_callpeak.A549_CTRL_BRD4_rep1.7896ea17b4595ac6d13192f7b9e93089.mugqic.done
)
macs2_callpeak_33_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE &&  $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_33_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_34_JOB_ID: macs2_callpeak_bigBed.A549_CTRL_BRD4_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.A549_CTRL_BRD4_rep1
JOB_DEPENDENCIES=$macs2_callpeak_33_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.A549_CTRL_BRD4_rep1.260d3f788b99c163810da94b9f9fbc9d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.A549_CTRL_BRD4_rep1.260d3f788b99c163810da94b9f9fbc9d.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/A549_CTRL_BRD4_rep1/A549_CTRL_BRD4_rep1_peaks.narrowPeak > peak_call/A549_CTRL_BRD4_rep1/A549_CTRL_BRD4_rep1_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/A549_CTRL_BRD4_rep1/A549_CTRL_BRD4_rep1_peaks.narrowPeak.bed \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/A549_CTRL_BRD4_rep1/A549_CTRL_BRD4_rep1_peaks.narrowPeak.bb
macs2_callpeak_bigBed.A549_CTRL_BRD4_rep1.260d3f788b99c163810da94b9f9fbc9d.mugqic.done
)
macs2_callpeak_34_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE &&  $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_34_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_35_JOB_ID: macs2_callpeak.A549_CTRL_BRD4_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_CTRL_BRD4_rep2
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_CTRL_BRD4_rep2.7a6253bf658306201bf0787289ad4d49.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_CTRL_BRD4_rep2.7a6253bf658306201bf0787289ad4d49.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/A549_CTRL_BRD4_rep2 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/A549_CTRL_BRD4_rep2/A549_CTRL_BRD4_rep2.sorted.dup.bam \
  --control \
  alignment/A549_CTRL_WCE_rep3/A549_CTRL_WCE_rep3.sorted.dup.bam \
  --name peak_call/A549_CTRL_BRD4_rep2/A549_CTRL_BRD4_rep2 \
  >& peak_call/A549_CTRL_BRD4_rep2/A549_CTRL_BRD4_rep2.diag.macs.out
macs2_callpeak.A549_CTRL_BRD4_rep2.7a6253bf658306201bf0787289ad4d49.mugqic.done
)
macs2_callpeak_35_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE &&  $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_35_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_36_JOB_ID: macs2_callpeak_bigBed.A549_CTRL_BRD4_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.A549_CTRL_BRD4_rep2
JOB_DEPENDENCIES=$macs2_callpeak_35_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.A549_CTRL_BRD4_rep2.7196485f9d1352794c11a9e3ffc882fc.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.A549_CTRL_BRD4_rep2.7196485f9d1352794c11a9e3ffc882fc.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/A549_CTRL_BRD4_rep2/A549_CTRL_BRD4_rep2_peaks.narrowPeak > peak_call/A549_CTRL_BRD4_rep2/A549_CTRL_BRD4_rep2_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/A549_CTRL_BRD4_rep2/A549_CTRL_BRD4_rep2_peaks.narrowPeak.bed \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/A549_CTRL_BRD4_rep2/A549_CTRL_BRD4_rep2_peaks.narrowPeak.bb
macs2_callpeak_bigBed.A549_CTRL_BRD4_rep2.7196485f9d1352794c11a9e3ffc882fc.mugqic.done
)
macs2_callpeak_36_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE &&  $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_36_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_37_JOB_ID: macs2_callpeak.A549_CTRL_CDK9_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_CTRL_CDK9_rep1
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_CTRL_CDK9_rep1.b7fb082777f463a8f78f98cd5c92fd82.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_CTRL_CDK9_rep1.b7fb082777f463a8f78f98cd5c92fd82.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/A549_CTRL_CDK9_rep1 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/A549_CTRL_CDK9_rep1/A549_CTRL_CDK9_rep1.sorted.dup.bam \
  --control \
  alignment/A549_CTRL_WCE_rep3/A549_CTRL_WCE_rep3.sorted.dup.bam \
  --name peak_call/A549_CTRL_CDK9_rep1/A549_CTRL_CDK9_rep1 \
  >& peak_call/A549_CTRL_CDK9_rep1/A549_CTRL_CDK9_rep1.diag.macs.out
macs2_callpeak.A549_CTRL_CDK9_rep1.b7fb082777f463a8f78f98cd5c92fd82.mugqic.done
)
macs2_callpeak_37_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE &&  $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_37_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_38_JOB_ID: macs2_callpeak_bigBed.A549_CTRL_CDK9_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.A549_CTRL_CDK9_rep1
JOB_DEPENDENCIES=$macs2_callpeak_37_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.A549_CTRL_CDK9_rep1.1a59d051bd370548fe5a9f9e0605071d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.A549_CTRL_CDK9_rep1.1a59d051bd370548fe5a9f9e0605071d.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/A549_CTRL_CDK9_rep1/A549_CTRL_CDK9_rep1_peaks.narrowPeak > peak_call/A549_CTRL_CDK9_rep1/A549_CTRL_CDK9_rep1_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/A549_CTRL_CDK9_rep1/A549_CTRL_CDK9_rep1_peaks.narrowPeak.bed \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/A549_CTRL_CDK9_rep1/A549_CTRL_CDK9_rep1_peaks.narrowPeak.bb
macs2_callpeak_bigBed.A549_CTRL_CDK9_rep1.1a59d051bd370548fe5a9f9e0605071d.mugqic.done
)
macs2_callpeak_38_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE &&  $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_38_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_39_JOB_ID: macs2_callpeak.A549_CTRL_CDK9_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_CTRL_CDK9_rep2
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_CTRL_CDK9_rep2.32e304aafbdf97bddf4d2df064f082b6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_CTRL_CDK9_rep2.32e304aafbdf97bddf4d2df064f082b6.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/A549_CTRL_CDK9_rep2 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/A549_CTRL_CDK9_rep2/A549_CTRL_CDK9_rep2.sorted.dup.bam \
  --control \
  alignment/A549_CTRL_WCE_rep3/A549_CTRL_WCE_rep3.sorted.dup.bam \
  --name peak_call/A549_CTRL_CDK9_rep2/A549_CTRL_CDK9_rep2 \
  >& peak_call/A549_CTRL_CDK9_rep2/A549_CTRL_CDK9_rep2.diag.macs.out
macs2_callpeak.A549_CTRL_CDK9_rep2.32e304aafbdf97bddf4d2df064f082b6.mugqic.done
)
macs2_callpeak_39_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE &&  $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_39_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_40_JOB_ID: macs2_callpeak_bigBed.A549_CTRL_CDK9_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.A549_CTRL_CDK9_rep2
JOB_DEPENDENCIES=$macs2_callpeak_39_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.A549_CTRL_CDK9_rep2.5cce1ff10bd3b62d09b2a7f560ba3580.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.A549_CTRL_CDK9_rep2.5cce1ff10bd3b62d09b2a7f560ba3580.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/A549_CTRL_CDK9_rep2/A549_CTRL_CDK9_rep2_peaks.narrowPeak > peak_call/A549_CTRL_CDK9_rep2/A549_CTRL_CDK9_rep2_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/A549_CTRL_CDK9_rep2/A549_CTRL_CDK9_rep2_peaks.narrowPeak.bed \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/A549_CTRL_CDK9_rep2/A549_CTRL_CDK9_rep2_peaks.narrowPeak.bb
macs2_callpeak_bigBed.A549_CTRL_CDK9_rep2.5cce1ff10bd3b62d09b2a7f560ba3580.mugqic.done
)
macs2_callpeak_40_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE &&  $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_40_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_41_JOB_ID: macs2_callpeak.A549_DEX_BRD4_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_DEX_BRD4_rep1
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_DEX_BRD4_rep1.559290731c467ac03fc5dae128bd0a73.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_DEX_BRD4_rep1.559290731c467ac03fc5dae128bd0a73.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/A549_DEX_BRD4_rep1 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/A549_DEX_BRD4_rep1/A549_DEX_BRD4_rep1.sorted.dup.bam \
  --control \
  alignment/A549_DEX_WCE_rep1/A549_DEX_WCE_rep1.sorted.dup.bam \
  --name peak_call/A549_DEX_BRD4_rep1/A549_DEX_BRD4_rep1 \
  >& peak_call/A549_DEX_BRD4_rep1/A549_DEX_BRD4_rep1.diag.macs.out
macs2_callpeak.A549_DEX_BRD4_rep1.559290731c467ac03fc5dae128bd0a73.mugqic.done
)
macs2_callpeak_41_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE &&  $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_41_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_42_JOB_ID: macs2_callpeak_bigBed.A549_DEX_BRD4_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.A549_DEX_BRD4_rep1
JOB_DEPENDENCIES=$macs2_callpeak_41_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.A549_DEX_BRD4_rep1.32a3b6897e67c4273f332429e57b12d1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.A549_DEX_BRD4_rep1.32a3b6897e67c4273f332429e57b12d1.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/A549_DEX_BRD4_rep1/A549_DEX_BRD4_rep1_peaks.narrowPeak > peak_call/A549_DEX_BRD4_rep1/A549_DEX_BRD4_rep1_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/A549_DEX_BRD4_rep1/A549_DEX_BRD4_rep1_peaks.narrowPeak.bed \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/A549_DEX_BRD4_rep1/A549_DEX_BRD4_rep1_peaks.narrowPeak.bb
macs2_callpeak_bigBed.A549_DEX_BRD4_rep1.32a3b6897e67c4273f332429e57b12d1.mugqic.done
)
macs2_callpeak_42_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE &&  $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_42_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_43_JOB_ID: macs2_callpeak.A549_DEX_BRD4_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_DEX_BRD4_rep2
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_DEX_BRD4_rep2.fce237994323868b5148ea97dbbd8ec7.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_DEX_BRD4_rep2.fce237994323868b5148ea97dbbd8ec7.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/A549_DEX_BRD4_rep2 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/A549_DEX_BRD4_rep2/A549_DEX_BRD4_rep2.sorted.dup.bam \
  --control \
  alignment/A549_DEX_WCE_rep2/A549_DEX_WCE_rep2.sorted.dup.bam \
  --name peak_call/A549_DEX_BRD4_rep2/A549_DEX_BRD4_rep2 \
  >& peak_call/A549_DEX_BRD4_rep2/A549_DEX_BRD4_rep2.diag.macs.out
macs2_callpeak.A549_DEX_BRD4_rep2.fce237994323868b5148ea97dbbd8ec7.mugqic.done
)
macs2_callpeak_43_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE &&  $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_43_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_44_JOB_ID: macs2_callpeak_bigBed.A549_DEX_BRD4_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.A549_DEX_BRD4_rep2
JOB_DEPENDENCIES=$macs2_callpeak_43_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.A549_DEX_BRD4_rep2.4d357e008257f521d03fa008697c7da5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.A549_DEX_BRD4_rep2.4d357e008257f521d03fa008697c7da5.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/A549_DEX_BRD4_rep2/A549_DEX_BRD4_rep2_peaks.narrowPeak > peak_call/A549_DEX_BRD4_rep2/A549_DEX_BRD4_rep2_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/A549_DEX_BRD4_rep2/A549_DEX_BRD4_rep2_peaks.narrowPeak.bed \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/A549_DEX_BRD4_rep2/A549_DEX_BRD4_rep2_peaks.narrowPeak.bb
macs2_callpeak_bigBed.A549_DEX_BRD4_rep2.4d357e008257f521d03fa008697c7da5.mugqic.done
)
macs2_callpeak_44_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE &&  $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_44_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_45_JOB_ID: macs2_callpeak.A549_DEX_CDK9_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_DEX_CDK9_rep1
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_DEX_CDK9_rep1.d97064da80fd1a3e55ff8d9dfb3102d5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_DEX_CDK9_rep1.d97064da80fd1a3e55ff8d9dfb3102d5.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/A549_DEX_CDK9_rep1 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/A549_DEX_CDK9_rep1/A549_DEX_CDK9_rep1.sorted.dup.bam \
  --control \
  alignment/A549_DEX_WCE_rep1/A549_DEX_WCE_rep1.sorted.dup.bam \
  --name peak_call/A549_DEX_CDK9_rep1/A549_DEX_CDK9_rep1 \
  >& peak_call/A549_DEX_CDK9_rep1/A549_DEX_CDK9_rep1.diag.macs.out
macs2_callpeak.A549_DEX_CDK9_rep1.d97064da80fd1a3e55ff8d9dfb3102d5.mugqic.done
)
macs2_callpeak_45_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE &&  $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_45_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_46_JOB_ID: macs2_callpeak_bigBed.A549_DEX_CDK9_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.A549_DEX_CDK9_rep1
JOB_DEPENDENCIES=$macs2_callpeak_45_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.A549_DEX_CDK9_rep1.181989cfafdaa3aed08167070293e989.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.A549_DEX_CDK9_rep1.181989cfafdaa3aed08167070293e989.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/A549_DEX_CDK9_rep1/A549_DEX_CDK9_rep1_peaks.narrowPeak > peak_call/A549_DEX_CDK9_rep1/A549_DEX_CDK9_rep1_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/A549_DEX_CDK9_rep1/A549_DEX_CDK9_rep1_peaks.narrowPeak.bed \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/A549_DEX_CDK9_rep1/A549_DEX_CDK9_rep1_peaks.narrowPeak.bb
macs2_callpeak_bigBed.A549_DEX_CDK9_rep1.181989cfafdaa3aed08167070293e989.mugqic.done
)
macs2_callpeak_46_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE &&  $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_46_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_47_JOB_ID: macs2_callpeak.A549_DEX_CDK9_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_DEX_CDK9_rep2
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_DEX_CDK9_rep2.112fa45931368f0da1f04d2244e4dd25.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_DEX_CDK9_rep2.112fa45931368f0da1f04d2244e4dd25.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/A549_DEX_CDK9_rep2 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/A549_DEX_CDK9_rep2/A549_DEX_CDK9_rep2.sorted.dup.bam \
  --control \
  alignment/A549_DEX_WCE_rep2/A549_DEX_WCE_rep2.sorted.dup.bam \
  --name peak_call/A549_DEX_CDK9_rep2/A549_DEX_CDK9_rep2 \
  >& peak_call/A549_DEX_CDK9_rep2/A549_DEX_CDK9_rep2.diag.macs.out
macs2_callpeak.A549_DEX_CDK9_rep2.112fa45931368f0da1f04d2244e4dd25.mugqic.done
)
macs2_callpeak_47_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE &&  $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_47_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_48_JOB_ID: macs2_callpeak_bigBed.A549_DEX_CDK9_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.A549_DEX_CDK9_rep2
JOB_DEPENDENCIES=$macs2_callpeak_47_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.A549_DEX_CDK9_rep2.0494953ef0ec511912857fa9ef81833c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.A549_DEX_CDK9_rep2.0494953ef0ec511912857fa9ef81833c.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/A549_DEX_CDK9_rep2/A549_DEX_CDK9_rep2_peaks.narrowPeak > peak_call/A549_DEX_CDK9_rep2/A549_DEX_CDK9_rep2_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/A549_DEX_CDK9_rep2/A549_DEX_CDK9_rep2_peaks.narrowPeak.bed \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/A549_DEX_CDK9_rep2/A549_DEX_CDK9_rep2_peaks.narrowPeak.bb
macs2_callpeak_bigBed.A549_DEX_CDK9_rep2.0494953ef0ec511912857fa9ef81833c.mugqic.done
)
macs2_callpeak_48_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE &&  $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_48_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_49_JOB_ID: macs2_callpeak_report
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_report
JOB_DEPENDENCIES=$macs2_callpeak_1_JOB_ID:$macs2_callpeak_3_JOB_ID:$macs2_callpeak_5_JOB_ID:$macs2_callpeak_7_JOB_ID:$macs2_callpeak_9_JOB_ID:$macs2_callpeak_11_JOB_ID:$macs2_callpeak_13_JOB_ID:$macs2_callpeak_15_JOB_ID:$macs2_callpeak_17_JOB_ID:$macs2_callpeak_19_JOB_ID:$macs2_callpeak_21_JOB_ID:$macs2_callpeak_23_JOB_ID:$macs2_callpeak_25_JOB_ID:$macs2_callpeak_27_JOB_ID:$macs2_callpeak_29_JOB_ID:$macs2_callpeak_31_JOB_ID:$macs2_callpeak_33_JOB_ID:$macs2_callpeak_35_JOB_ID:$macs2_callpeak_37_JOB_ID:$macs2_callpeak_39_JOB_ID:$macs2_callpeak_41_JOB_ID:$macs2_callpeak_43_JOB_ID:$macs2_callpeak_45_JOB_ID:$macs2_callpeak_47_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_report.b82bf30e5b8592ff02e519aa3dc36e1d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_report.b82bf30e5b8592ff02e519aa3dc36e1d.mugqic.done'
mkdir -p report && \
cp /home/efournie/genpipes/bfx/report/ChipSeq.macs2_callpeak.md report/ && \
for contrast in A549_CTRL_BRD4 A549_CTRL_CDK9 A549_CTRL_MED1 A549_CTRL_NIPBL A549_CTRL_SMC1A A549_DEX_BRD4 A549_DEX_CDK9 A549_DEX_MED1 A549_DEX_NIPBL A549_DEX_SMC1A A549_CTRL_MED1_rep1 A549_CTRL_NIPBL_rep1 A549_CTRL_SMC1A_rep1 A549_CTRL_MED1_rep2 A549_CTRL_NIPBL_rep2 A549_CTRL_SMC1A_rep2 A549_CTRL_BRD4_rep1 A549_CTRL_BRD4_rep2 A549_CTRL_CDK9_rep1 A549_CTRL_CDK9_rep2 A549_DEX_BRD4_rep1 A549_DEX_BRD4_rep2 A549_DEX_CDK9_rep1 A549_DEX_CDK9_rep2
do
  cp -a --parents peak_call/$contrast/ report/ && \
  echo -e "* [Peak Calls File for Design $contrast](peak_call/$contrast/${contrast}_peaks.xls)" \
  >> report/ChipSeq.macs2_callpeak.md
done
macs2_callpeak_report.b82bf30e5b8592ff02e519aa3dc36e1d.mugqic.done
)
macs2_callpeak_49_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE &&  $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE

if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_49_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# Call home with pipeline statistics
#-------------------------------------------------------------------------------
wget "http://mugqic.hpc.mcgill.ca/cgi-bin/pipeline.cgi?hostname=cedar5.cedar.computecanada.ca&ip=206.12.124.6&pipeline=ChipSeq&steps=macs2_callpeak&samples=22&AnonymizedList=cf29f0b79a6fdf85172b3fba4d749551,e827d701807f55d324ab6e04828e0e1b,23380dc1332bbc71a76b52587e6749e9,8a887031ec5ae523ccf57f159ba1cc8e,c943fea81479b8ac1d8802627f59f77b,7624f07d7c46c52b4184bf71629b3bd1,01b37bdc8112c6742f60de53a15612e3,9bb4e4c9a92e657455c6b0465d159bfb,1080ac22c0bbbde9594f8886c2c9c1c8,b4f3f1bd5062aafce8ecc8161749f295,2b9ab9928bce37761cbd5cf78e038e2e,70199dc535d02a760cd016f1288ae793,b86cc53fe7583000c5e891c089286325,e44fbb4967bf7b09a43edd9f40dc6a8c,6e58b87e2df9fee5c171ae930728d26d,151ebda0d13caffcb3815f9dac129675,0b276224b0b1bec51d03b4abd2f1835e,a71351e758a57b1032b63c1faa673fcd,b0d72e85cd7e4b49306a3cf0c9316c20,5113063641ec3d728149cb841aeb453f,858366472bf4ef3fea6d4ccdcd7cdb5c,5511f97f8b7cbf85c702ffced20d11dd" --quiet --output-document=/dev/null

