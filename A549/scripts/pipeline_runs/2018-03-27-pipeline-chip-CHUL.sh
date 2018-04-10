#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

#-------------------------------------------------------------------------------
# ChipSeq SlurmScheduler Job Submission Bash script
# Version: 3.0.1-beta
# Created on: 2018-03-27T13:58:30
# Steps:
#   macs2_callpeak: 19 jobs
#   TOTAL: 19 jobs
#-------------------------------------------------------------------------------

OUTPUT_DIR=/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38
JOB_OUTPUT_DIR=$OUTPUT_DIR/job_output
TIMESTAMP=`date +%FT%H.%M.%S`
JOB_LIST=$JOB_OUTPUT_DIR/ChipSeq_job_list_$TIMESTAMP
export CONFIG_FILES="/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini"
mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR


#-------------------------------------------------------------------------------
# STEP: macs2_callpeak
#-------------------------------------------------------------------------------
STEP=macs2_callpeak
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: macs2_callpeak_1_JOB_ID: macs2_callpeak.A549_DEX_MED1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_DEX_MED1
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_DEX_MED1.7d4c12dbc2b5cd75b6155679a2cb0867.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_DEX_MED1.7d4c12dbc2b5cd75b6155679a2cb0867.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/A549_DEX_MED1 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/A549_DEX_MED1_rep1/A549_DEX_MED1_rep1.sorted.dup.bam \
  alignment/A549_DEX_MED1_rep2/A549_DEX_MED1_rep2.sorted.dup.bam \
  --control \
  alignment/A549_DEX_WCE_rep1/A549_DEX_WCE_rep1.sorted.dup.bam \
  alignment/A549_DEX_WCE_rep2/A549_DEX_WCE_rep2.sorted.dup.bam \
  --name peak_call/A549_DEX_MED1/A549_DEX_MED1 \
  >& peak_call/A549_DEX_MED1/A549_DEX_MED1.diag.macs.out
macs2_callpeak.A549_DEX_MED1.7d4c12dbc2b5cd75b6155679a2cb0867.mugqic.done
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
# JOB: macs2_callpeak_2_JOB_ID: macs2_callpeak_bigBed.A549_DEX_MED1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.A549_DEX_MED1
JOB_DEPENDENCIES=$macs2_callpeak_1_JOB_ID
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
# JOB: macs2_callpeak_3_JOB_ID: macs2_callpeak.A549_DEX_NIPBL
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_DEX_NIPBL
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_DEX_NIPBL.1b1ca9c75bbc8042a3b8987da1e4fa7f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_DEX_NIPBL.1b1ca9c75bbc8042a3b8987da1e4fa7f.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/A549_DEX_NIPBL && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/A549_DEX_NIPBL_rep1/A549_DEX_NIPBL_rep1.sorted.dup.bam \
  alignment/A549_DEX_NIPBL_rep2/A549_DEX_NIPBL_rep2.sorted.dup.bam \
  --control \
  alignment/A549_DEX_WCE_rep1/A549_DEX_WCE_rep1.sorted.dup.bam \
  alignment/A549_DEX_WCE_rep2/A549_DEX_WCE_rep2.sorted.dup.bam \
  --name peak_call/A549_DEX_NIPBL/A549_DEX_NIPBL \
  >& peak_call/A549_DEX_NIPBL/A549_DEX_NIPBL.diag.macs.out
macs2_callpeak.A549_DEX_NIPBL.1b1ca9c75bbc8042a3b8987da1e4fa7f.mugqic.done
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
# JOB: macs2_callpeak_4_JOB_ID: macs2_callpeak_bigBed.A549_DEX_NIPBL
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.A549_DEX_NIPBL
JOB_DEPENDENCIES=$macs2_callpeak_3_JOB_ID
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
# JOB: macs2_callpeak_5_JOB_ID: macs2_callpeak.A549_DEX_SMC1A
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_DEX_SMC1A
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_DEX_SMC1A.97267df97925a05e014a6fee63b04f37.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_DEX_SMC1A.97267df97925a05e014a6fee63b04f37.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/A549_DEX_SMC1A && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/A549_DEX_SMC1A_rep1/A549_DEX_SMC1A_rep1.sorted.dup.bam \
  alignment/A549_DEX_SMC1A_rep2/A549_DEX_SMC1A_rep2.sorted.dup.bam \
  --control \
  alignment/A549_DEX_WCE_rep1/A549_DEX_WCE_rep1.sorted.dup.bam \
  alignment/A549_DEX_WCE_rep2/A549_DEX_WCE_rep2.sorted.dup.bam \
  --name peak_call/A549_DEX_SMC1A/A549_DEX_SMC1A \
  >& peak_call/A549_DEX_SMC1A/A549_DEX_SMC1A.diag.macs.out
macs2_callpeak.A549_DEX_SMC1A.97267df97925a05e014a6fee63b04f37.mugqic.done
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
# JOB: macs2_callpeak_6_JOB_ID: macs2_callpeak_bigBed.A549_DEX_SMC1A
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.A549_DEX_SMC1A
JOB_DEPENDENCIES=$macs2_callpeak_5_JOB_ID
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
# JOB: macs2_callpeak_7_JOB_ID: macs2_callpeak.A549_DEX_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_DEX_MED1_rep1
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_DEX_MED1_rep1.b870dfd5ca855c4e7f64f9dfa39fb5e1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_DEX_MED1_rep1.b870dfd5ca855c4e7f64f9dfa39fb5e1.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/A549_DEX_MED1_rep1 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/A549_DEX_MED1_rep1/A549_DEX_MED1_rep1.sorted.dup.bam \
  --control \
  alignment/A549_DEX_WCE_rep1/A549_DEX_WCE_rep1.sorted.dup.bam \
  --name peak_call/A549_DEX_MED1_rep1/A549_DEX_MED1_rep1 \
  >& peak_call/A549_DEX_MED1_rep1/A549_DEX_MED1_rep1.diag.macs.out
macs2_callpeak.A549_DEX_MED1_rep1.b870dfd5ca855c4e7f64f9dfa39fb5e1.mugqic.done
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
# JOB: macs2_callpeak_8_JOB_ID: macs2_callpeak_bigBed.A549_DEX_MED1_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.A549_DEX_MED1_rep1
JOB_DEPENDENCIES=$macs2_callpeak_7_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.A549_DEX_MED1_rep1.44a9f29f872de12fc865ef9dda7643ee.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.A549_DEX_MED1_rep1.44a9f29f872de12fc865ef9dda7643ee.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/A549_DEX_MED1_rep1/A549_DEX_MED1_rep1_peaks.narrowPeak > peak_call/A549_DEX_MED1_rep1/A549_DEX_MED1_rep1_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/A549_DEX_MED1_rep1/A549_DEX_MED1_rep1_peaks.narrowPeak.bed \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/A549_DEX_MED1_rep1/A549_DEX_MED1_rep1_peaks.narrowPeak.bb
macs2_callpeak_bigBed.A549_DEX_MED1_rep1.44a9f29f872de12fc865ef9dda7643ee.mugqic.done
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
# JOB: macs2_callpeak_9_JOB_ID: macs2_callpeak.A549_DEX_MED1_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_DEX_MED1_rep2
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_DEX_MED1_rep2.0e7d07671e5233f9c6ea64146ba0ccce.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_DEX_MED1_rep2.0e7d07671e5233f9c6ea64146ba0ccce.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/A549_DEX_MED1_rep2 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/A549_DEX_MED1_rep2/A549_DEX_MED1_rep2.sorted.dup.bam \
  --control \
  alignment/A549_DEX_WCE_rep2/A549_DEX_WCE_rep2.sorted.dup.bam \
  --name peak_call/A549_DEX_MED1_rep2/A549_DEX_MED1_rep2 \
  >& peak_call/A549_DEX_MED1_rep2/A549_DEX_MED1_rep2.diag.macs.out
macs2_callpeak.A549_DEX_MED1_rep2.0e7d07671e5233f9c6ea64146ba0ccce.mugqic.done
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
# JOB: macs2_callpeak_10_JOB_ID: macs2_callpeak_bigBed.A549_DEX_MED1_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.A549_DEX_MED1_rep2
JOB_DEPENDENCIES=$macs2_callpeak_9_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.A549_DEX_MED1_rep2.9a12b78975ee388c0cfc62e8e820e71a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.A549_DEX_MED1_rep2.9a12b78975ee388c0cfc62e8e820e71a.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/A549_DEX_MED1_rep2/A549_DEX_MED1_rep2_peaks.narrowPeak > peak_call/A549_DEX_MED1_rep2/A549_DEX_MED1_rep2_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/A549_DEX_MED1_rep2/A549_DEX_MED1_rep2_peaks.narrowPeak.bed \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/A549_DEX_MED1_rep2/A549_DEX_MED1_rep2_peaks.narrowPeak.bb
macs2_callpeak_bigBed.A549_DEX_MED1_rep2.9a12b78975ee388c0cfc62e8e820e71a.mugqic.done
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
# JOB: macs2_callpeak_11_JOB_ID: macs2_callpeak.A549_DEX_NIPBL_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_DEX_NIPBL_rep1
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_DEX_NIPBL_rep1.9def8f07baeb330f9c8b54b20591d1ca.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_DEX_NIPBL_rep1.9def8f07baeb330f9c8b54b20591d1ca.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/A549_DEX_NIPBL_rep1 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/A549_DEX_NIPBL_rep1/A549_DEX_NIPBL_rep1.sorted.dup.bam \
  --control \
  alignment/A549_DEX_WCE_rep1/A549_DEX_WCE_rep1.sorted.dup.bam \
  --name peak_call/A549_DEX_NIPBL_rep1/A549_DEX_NIPBL_rep1 \
  >& peak_call/A549_DEX_NIPBL_rep1/A549_DEX_NIPBL_rep1.diag.macs.out
macs2_callpeak.A549_DEX_NIPBL_rep1.9def8f07baeb330f9c8b54b20591d1ca.mugqic.done
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
# JOB: macs2_callpeak_12_JOB_ID: macs2_callpeak_bigBed.A549_DEX_NIPBL_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.A549_DEX_NIPBL_rep1
JOB_DEPENDENCIES=$macs2_callpeak_11_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.A549_DEX_NIPBL_rep1.effb972b10503fb086575d40439f3d83.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.A549_DEX_NIPBL_rep1.effb972b10503fb086575d40439f3d83.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/A549_DEX_NIPBL_rep1/A549_DEX_NIPBL_rep1_peaks.narrowPeak > peak_call/A549_DEX_NIPBL_rep1/A549_DEX_NIPBL_rep1_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/A549_DEX_NIPBL_rep1/A549_DEX_NIPBL_rep1_peaks.narrowPeak.bed \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/A549_DEX_NIPBL_rep1/A549_DEX_NIPBL_rep1_peaks.narrowPeak.bb
macs2_callpeak_bigBed.A549_DEX_NIPBL_rep1.effb972b10503fb086575d40439f3d83.mugqic.done
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
# JOB: macs2_callpeak_13_JOB_ID: macs2_callpeak.A549_DEX_NIPBL_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_DEX_NIPBL_rep2
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_DEX_NIPBL_rep2.504bef5cedb5ca3df2b3851d572ecec1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_DEX_NIPBL_rep2.504bef5cedb5ca3df2b3851d572ecec1.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/A549_DEX_NIPBL_rep2 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/A549_DEX_NIPBL_rep2/A549_DEX_NIPBL_rep2.sorted.dup.bam \
  --control \
  alignment/A549_DEX_WCE_rep2/A549_DEX_WCE_rep2.sorted.dup.bam \
  --name peak_call/A549_DEX_NIPBL_rep2/A549_DEX_NIPBL_rep2 \
  >& peak_call/A549_DEX_NIPBL_rep2/A549_DEX_NIPBL_rep2.diag.macs.out
macs2_callpeak.A549_DEX_NIPBL_rep2.504bef5cedb5ca3df2b3851d572ecec1.mugqic.done
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
# JOB: macs2_callpeak_14_JOB_ID: macs2_callpeak_bigBed.A549_DEX_NIPBL_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.A549_DEX_NIPBL_rep2
JOB_DEPENDENCIES=$macs2_callpeak_13_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.A549_DEX_NIPBL_rep2.008a2e758788dde7e024d1aefccae83f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.A549_DEX_NIPBL_rep2.008a2e758788dde7e024d1aefccae83f.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/A549_DEX_NIPBL_rep2/A549_DEX_NIPBL_rep2_peaks.narrowPeak > peak_call/A549_DEX_NIPBL_rep2/A549_DEX_NIPBL_rep2_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/A549_DEX_NIPBL_rep2/A549_DEX_NIPBL_rep2_peaks.narrowPeak.bed \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/A549_DEX_NIPBL_rep2/A549_DEX_NIPBL_rep2_peaks.narrowPeak.bb
macs2_callpeak_bigBed.A549_DEX_NIPBL_rep2.008a2e758788dde7e024d1aefccae83f.mugqic.done
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
# JOB: macs2_callpeak_15_JOB_ID: macs2_callpeak.A549_DEX_SMC1A_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_DEX_SMC1A_rep1
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_DEX_SMC1A_rep1.19708e6ecf239793eadd91f8b17feb86.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_DEX_SMC1A_rep1.19708e6ecf239793eadd91f8b17feb86.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/A549_DEX_SMC1A_rep1 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/A549_DEX_SMC1A_rep1/A549_DEX_SMC1A_rep1.sorted.dup.bam \
  --control \
  alignment/A549_DEX_WCE_rep1/A549_DEX_WCE_rep1.sorted.dup.bam \
  --name peak_call/A549_DEX_SMC1A_rep1/A549_DEX_SMC1A_rep1 \
  >& peak_call/A549_DEX_SMC1A_rep1/A549_DEX_SMC1A_rep1.diag.macs.out
macs2_callpeak.A549_DEX_SMC1A_rep1.19708e6ecf239793eadd91f8b17feb86.mugqic.done
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
# JOB: macs2_callpeak_16_JOB_ID: macs2_callpeak_bigBed.A549_DEX_SMC1A_rep1
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.A549_DEX_SMC1A_rep1
JOB_DEPENDENCIES=$macs2_callpeak_15_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.A549_DEX_SMC1A_rep1.08dd73d7b1a98f044a12a36836ad7fff.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.A549_DEX_SMC1A_rep1.08dd73d7b1a98f044a12a36836ad7fff.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/A549_DEX_SMC1A_rep1/A549_DEX_SMC1A_rep1_peaks.narrowPeak > peak_call/A549_DEX_SMC1A_rep1/A549_DEX_SMC1A_rep1_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/A549_DEX_SMC1A_rep1/A549_DEX_SMC1A_rep1_peaks.narrowPeak.bed \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/A549_DEX_SMC1A_rep1/A549_DEX_SMC1A_rep1_peaks.narrowPeak.bb
macs2_callpeak_bigBed.A549_DEX_SMC1A_rep1.08dd73d7b1a98f044a12a36836ad7fff.mugqic.done
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
# JOB: macs2_callpeak_17_JOB_ID: macs2_callpeak.A549_DEX_SMC1A_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak.A549_DEX_SMC1A_rep2
JOB_DEPENDENCIES=
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak.A549_DEX_SMC1A_rep2.dd547240d72c0fd10d390ddc3502dd0b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak.A549_DEX_SMC1A_rep2.dd547240d72c0fd10d390ddc3502dd0b.mugqic.done'
module load mugqic/python/2.7.13 mugqic/MACS2/2.1.1.20160309 && \
mkdir -p peak_call/A549_DEX_SMC1A_rep2 && \
macs2 callpeak --format BAM --nomodel \
  --tempdir ${SLURM_TMPDIR} \
  --gsize 2479938032.8 \
  --treatment \
  alignment/A549_DEX_SMC1A_rep2/A549_DEX_SMC1A_rep2.sorted.dup.bam \
  --control \
  alignment/A549_DEX_WCE_rep2/A549_DEX_WCE_rep2.sorted.dup.bam \
  --name peak_call/A549_DEX_SMC1A_rep2/A549_DEX_SMC1A_rep2 \
  >& peak_call/A549_DEX_SMC1A_rep2/A549_DEX_SMC1A_rep2.diag.macs.out
macs2_callpeak.A549_DEX_SMC1A_rep2.dd547240d72c0fd10d390ddc3502dd0b.mugqic.done
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
# JOB: macs2_callpeak_18_JOB_ID: macs2_callpeak_bigBed.A549_DEX_SMC1A_rep2
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_bigBed.A549_DEX_SMC1A_rep2
JOB_DEPENDENCIES=$macs2_callpeak_17_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_bigBed.A549_DEX_SMC1A_rep2.6c633ea6c0420ea18a7ed10bad5bc988.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_bigBed.A549_DEX_SMC1A_rep2.6c633ea6c0420ea18a7ed10bad5bc988.mugqic.done'
module load mugqic/ucsc/v346 && \
awk ' {if ($9 > 1000) {$9 = 1000} ; printf( "%s\t%s\t%s\t%s\t%0.f\n", $1,$2,$3,$4,$9)} ' peak_call/A549_DEX_SMC1A_rep2/A549_DEX_SMC1A_rep2_peaks.narrowPeak > peak_call/A549_DEX_SMC1A_rep2/A549_DEX_SMC1A_rep2_peaks.narrowPeak.bed && \
bedToBigBed \
  peak_call/A549_DEX_SMC1A_rep2/A549_DEX_SMC1A_rep2_peaks.narrowPeak.bed \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  peak_call/A549_DEX_SMC1A_rep2/A549_DEX_SMC1A_rep2_peaks.narrowPeak.bb
macs2_callpeak_bigBed.A549_DEX_SMC1A_rep2.6c633ea6c0420ea18a7ed10bad5bc988.mugqic.done
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
# JOB: macs2_callpeak_19_JOB_ID: macs2_callpeak_report
#-------------------------------------------------------------------------------
JOB_NAME=macs2_callpeak_report
JOB_DEPENDENCIES=$macs2_callpeak_1_JOB_ID:$macs2_callpeak_3_JOB_ID:$macs2_callpeak_5_JOB_ID:$macs2_callpeak_7_JOB_ID:$macs2_callpeak_9_JOB_ID:$macs2_callpeak_11_JOB_ID:$macs2_callpeak_13_JOB_ID:$macs2_callpeak_15_JOB_ID:$macs2_callpeak_17_JOB_ID
JOB_DONE=job_output/macs2_callpeak/macs2_callpeak_report.bad3813388e42a5b9e6316c2d30cb724.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'macs2_callpeak_report.bad3813388e42a5b9e6316c2d30cb724.mugqic.done'
mkdir -p report && \
cp /home/efournie/genpipes/bfx/report/ChipSeq.macs2_callpeak.md report/ && \
for contrast in A549_CTRL_BRD4 A549_CTRL_CDK9 A549_CTRL_MED1 A549_CTRL_NIPBL A549_CTRL_SMC1A A549_DEX_BRD4 A549_DEX_CDK9 A549_DEX_MED1 A549_DEX_NIPBL A549_DEX_SMC1A A549_CTRL_MED1_rep1 A549_CTRL_NIPBL_rep1 A549_CTRL_SMC1A_rep1 A549_CTRL_MED1_rep2 A549_CTRL_NIPBL_rep2 A549_CTRL_SMC1A_rep2 A549_CTRL_BRD4_rep1 A549_CTRL_BRD4_rep2 A549_CTRL_CDK9_rep1 A549_CTRL_CDK9_rep2 A549_DEX_BRD4_rep1 A549_DEX_BRD4_rep2 A549_DEX_CDK9_rep1 A549_DEX_CDK9_rep2 A549_DEX_MED1_rep1 A549_DEX_MED1_rep2 A549_DEX_NIPBL_rep1 A549_DEX_NIPBL_rep2 A549_DEX_SMC1A_rep1 A549_DEX_SMC1A_rep2
do
  cp -a --parents peak_call/$contrast/ report/ && \
  echo -e "* [Peak Calls File for Design $contrast](peak_call/$contrast/${contrast}_peaks.xls)" \
  >> report/ChipSeq.macs2_callpeak.md
done
macs2_callpeak_report.bad3813388e42a5b9e6316c2d30cb724.mugqic.done
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$macs2_callpeak_19_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# Call home with pipeline statistics
#-------------------------------------------------------------------------------
wget "http://mugqic.hpc.mcgill.ca/cgi-bin/pipeline.cgi?hostname=cedar5.cedar.computecanada.ca&ip=206.12.124.6&pipeline=ChipSeq&steps=macs2_callpeak&samples=25&AnonymizedList=6e58b87e2df9fee5c171ae930728d26d,e827d701807f55d324ab6e04828e0e1b,23380dc1332bbc71a76b52587e6749e9,8a887031ec5ae523ccf57f159ba1cc8e,11728c95eb3cea1f33a00588ac772178,d8a59c316384745ed7db9dbeb5f1e489,c943fea81479b8ac1d8802627f59f77b,7624f07d7c46c52b4184bf71629b3bd1,01b37bdc8112c6742f60de53a15612e3,9bb4e4c9a92e657455c6b0465d159bfb,1080ac22c0bbbde9594f8886c2c9c1c8,b4f3f1bd5062aafce8ecc8161749f295,2b9ab9928bce37761cbd5cf78e038e2e,70199dc535d02a760cd016f1288ae793,b86cc53fe7583000c5e891c089286325,e44fbb4967bf7b09a43edd9f40dc6a8c,ce51d19e1027f42536c485b28dc74e8a,cf29f0b79a6fdf85172b3fba4d749551,151ebda0d13caffcb3815f9dac129675,0b276224b0b1bec51d03b4abd2f1835e,a71351e758a57b1032b63c1faa673fcd,b0d72e85cd7e4b49306a3cf0c9316c20,5113063641ec3d728149cb841aeb453f,858366472bf4ef3fea6d4ccdcd7cdb5c,5511f97f8b7cbf85c702ffced20d11dd" --quiet --output-document=/dev/null

