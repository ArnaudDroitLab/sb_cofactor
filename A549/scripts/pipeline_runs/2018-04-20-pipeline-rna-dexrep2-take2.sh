#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

#-------------------------------------------------------------------------------
# RnaSeq SlurmScheduler Job Submission Bash script
# Version: 3.0.1-beta
# Created on: 2018-04-20T19:20:35
# Steps:
#   picard_sam_to_fastq: 0 job... skipping
#   trimmomatic: 0 job... skipping
#   merge_trimmomatic_stats: 0 job... skipping
#   star: 5 jobs
#   picard_merge_sam_files: 1 job
#   picard_sort_sam: 5 jobs
#   picard_mark_duplicates: 5 jobs
#   picard_rna_metrics: 5 jobs
#   estimate_ribosomal_rna: 6 jobs
#   bam_hard_clip: 5 jobs
#   rnaseqc: 2 jobs
#   wiggle: 30 jobs
#   raw_counts: 5 jobs
#   raw_counts_metrics: 4 jobs
#   cufflinks: 5 jobs
#   cuffmerge: 1 job
#   cuffquant: 5 jobs
#   cuffdiff: 16 jobs
#   cuffnorm: 1 job
#   TOTAL: 101 jobs
#-------------------------------------------------------------------------------

OUTPUT_DIR=/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38
JOB_OUTPUT_DIR=$OUTPUT_DIR/job_output
TIMESTAMP=`date +%FT%H.%M.%S`
JOB_LIST=$JOB_OUTPUT_DIR/RnaSeq_job_list_$TIMESTAMP
export CONFIG_FILES="/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini"
mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR


#-------------------------------------------------------------------------------
# STEP: star
#-------------------------------------------------------------------------------
STEP=star
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: star_1_JOB_ID: star_align.2.HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2
#-------------------------------------------------------------------------------
JOB_NAME=star_align.2.HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2
JOB_DEPENDENCIES=
JOB_DONE=job_output/star/star_align.2.HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2.8ecff57d5be6a2c0eef193cfb1f6a8be.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.2.HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2.8ecff57d5be6a2c0eef193cfb1f6a8be.mugqic.done'
module load mugqic/star/2.5.3a && \
mkdir -p alignment/DEX_shNIPBL_2/HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2 && \
STAR --runMode alignReads \
  --genomeDir reference.Merged \
  --readFilesIn \
    trim/DEX_shNIPBL_2/HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2.trim.pair1.fastq.gz \
    trim/DEX_shNIPBL_2/HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix alignment/DEX_shNIPBL_2/HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2/ \
  --outSAMattrRGline ID:"HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2" 	PL:"ILLUMINA" 			SM:"DEX_shNIPBL_2" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 100000000000 \
  --limitBAMsortRAM 100000000000 \
  --limitIObufferSize 4000000000 \
  --outWigType wiggle read1_5p --outWigStrand Stranded --outWigReferencesPrefix chr \
  --chimSegmentMin 21 && \
ln -s -f HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2/Aligned.sortedByCoord.out.bam alignment/DEX_shNIPBL_2/DEX_shNIPBL_2.sorted.bam
star_align.2.HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2.8ecff57d5be6a2c0eef193cfb1f6a8be.mugqic.done
)
star_1_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"star\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"star\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=128G -N 1 -n 16 | grep "[0-9]" | cut -d\  -f4)
echo "$star_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: star_2_JOB_ID: star_align.2.HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2
#-------------------------------------------------------------------------------
JOB_NAME=star_align.2.HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2
JOB_DEPENDENCIES=
JOB_DONE=job_output/star/star_align.2.HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2.6594585bda9a33609823ce56616d7fb9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.2.HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2.6594585bda9a33609823ce56616d7fb9.mugqic.done'
module load mugqic/star/2.5.3a && \
mkdir -p alignment/DEX_shMED1_2/HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2 && \
STAR --runMode alignReads \
  --genomeDir reference.Merged \
  --readFilesIn \
    trim/DEX_shMED1_2/HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2.trim.pair1.fastq.gz \
    trim/DEX_shMED1_2/HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix alignment/DEX_shMED1_2/HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2/ \
  --outSAMattrRGline ID:"HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2" 	PL:"ILLUMINA" 			SM:"DEX_shMED1_2" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 100000000000 \
  --limitBAMsortRAM 100000000000 \
  --limitIObufferSize 4000000000 \
  --outWigType wiggle read1_5p --outWigStrand Stranded --outWigReferencesPrefix chr \
  --chimSegmentMin 21 && \
ln -s -f HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2/Aligned.sortedByCoord.out.bam alignment/DEX_shMED1_2/DEX_shMED1_2.sorted.bam
star_align.2.HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2.6594585bda9a33609823ce56616d7fb9.mugqic.done
)
star_2_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"star\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"star\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=128G -N 1 -n 16 | grep "[0-9]" | cut -d\  -f4)
echo "$star_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: star_3_JOB_ID: star_align.2.HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2
#-------------------------------------------------------------------------------
JOB_NAME=star_align.2.HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2
JOB_DEPENDENCIES=
JOB_DONE=job_output/star/star_align.2.HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2.679f18c9ea26a0352c8a4c202b14c469.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.2.HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2.679f18c9ea26a0352c8a4c202b14c469.mugqic.done'
module load mugqic/star/2.5.3a && \
mkdir -p alignment/DEX_shSMC1A_4/HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2 && \
STAR --runMode alignReads \
  --genomeDir reference.Merged \
  --readFilesIn \
    trim/DEX_shSMC1A_4/HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2.trim.pair1.fastq.gz \
    trim/DEX_shSMC1A_4/HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix alignment/DEX_shSMC1A_4/HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2/ \
  --outSAMattrRGline ID:"HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2" 	PL:"ILLUMINA" 			SM:"DEX_shSMC1A_4" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 100000000000 \
  --limitBAMsortRAM 100000000000 \
  --limitIObufferSize 4000000000 \
  --outWigType wiggle read1_5p --outWigStrand Stranded --outWigReferencesPrefix chr \
  --chimSegmentMin 21 && \
ln -s -f HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2/Aligned.sortedByCoord.out.bam alignment/DEX_shSMC1A_4/DEX_shSMC1A_4.sorted.bam
star_align.2.HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2.679f18c9ea26a0352c8a4c202b14c469.mugqic.done
)
star_3_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"star\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"star\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=128G -N 1 -n 16 | grep "[0-9]" | cut -d\  -f4)
echo "$star_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: star_4_JOB_ID: star_align.2.HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2
#-------------------------------------------------------------------------------
JOB_NAME=star_align.2.HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2
JOB_DEPENDENCIES=
JOB_DONE=job_output/star/star_align.2.HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2.ec7dfcc78ec02487d7169d1202753a3f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_align.2.HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2.ec7dfcc78ec02487d7169d1202753a3f.mugqic.done'
module load mugqic/star/2.5.3a && \
mkdir -p alignment/DEX_nonMamm_2/HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2 && \
STAR --runMode alignReads \
  --genomeDir reference.Merged \
  --readFilesIn \
    trim/DEX_nonMamm_2/HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2.trim.pair1.fastq.gz \
    trim/DEX_nonMamm_2/HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2.trim.pair2.fastq.gz \
  --runThreadN 16 \
  --readFilesCommand zcat \
  --outStd Log \
  --outSAMunmapped Within \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix alignment/DEX_nonMamm_2/HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2/ \
  --outSAMattrRGline ID:"HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2" 	PL:"ILLUMINA" 			SM:"DEX_nonMamm_2" 	CN:"McGill University and Genome Quebec Innovation Centre"  \
  --limitGenomeGenerateRAM 100000000000 \
  --limitBAMsortRAM 100000000000 \
  --limitIObufferSize 4000000000 \
  --outWigType wiggle read1_5p --outWigStrand Stranded --outWigReferencesPrefix chr \
  --chimSegmentMin 21 && \
ln -s -f HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2/Aligned.sortedByCoord.out.bam alignment/DEX_nonMamm_2/DEX_nonMamm_2.sorted.bam
star_align.2.HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2.ec7dfcc78ec02487d7169d1202753a3f.mugqic.done
)
star_4_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"star\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"star\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=128G -N 1 -n 16 | grep "[0-9]" | cut -d\  -f4)
echo "$star_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: star_5_JOB_ID: star_report
#-------------------------------------------------------------------------------
JOB_NAME=star_report
JOB_DEPENDENCIES=$star_1_JOB_ID:$star_2_JOB_ID:$star_3_JOB_ID:$star_4_JOB_ID
JOB_DONE=job_output/star/star_report.20b949a56558f0d07a6baecd0297fd46.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'star_report.20b949a56558f0d07a6baecd0297fd46.mugqic.done'
module load mugqic/pandoc/1.15.2 && \
mkdir -p report && \
pandoc --to=markdown \
  --template /home/efournie/genpipes/bfx/report/RnaSeq.star.md \
  --variable scientific_name="Homo_sapiens" \
  --variable assembly="GRCh38" \
  /home/efournie/genpipes/bfx/report/RnaSeq.star.md \
  > report/RnaSeq.star.md
star_report.20b949a56558f0d07a6baecd0297fd46.mugqic.done
)
star_5_JOB_ID=$(echo "#! /bin/bash 
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
echo "$star_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# STEP: picard_merge_sam_files
#-------------------------------------------------------------------------------
STEP=picard_merge_sam_files
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_1_JOB_ID: picard_merge_sam_files.DEX_nonMamm_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_merge_sam_files.DEX_nonMamm_1
JOB_DEPENDENCIES=
JOB_DONE=job_output/picard_merge_sam_files/picard_merge_sam_files.DEX_nonMamm_1.3ce73cdb65578ab1834e4f3846d32f90.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_merge_sam_files.DEX_nonMamm_1.3ce73cdb65578ab1834e4f3846d32f90.mugqic.done'
module load java/1.8.0_121 mugqic/picard/2.9.0 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=4 -Dsamjdk.buffer_size=1048576 -Xmx40G -jar $PICARD_HOME/picard.jar MergeSamFiles \
 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/DEX_nonMamm_1/DEX_nonMamm_1_1/Aligned.sortedByCoord.out.bam \
 INPUT=alignment/DEX_nonMamm_1/DEX_nonMamm_1_2/Aligned.sortedByCoord.out.bam \
 OUTPUT=alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.bam \
 MAX_RECORDS_IN_RAM=5750000
picard_merge_sam_files.DEX_nonMamm_1.3ce73cdb65578ab1834e4f3846d32f90.mugqic.done
)
picard_merge_sam_files_1_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_merge_sam_files\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_1.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_merge_sam_files\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_1.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=48G -N 1 -n 12 | grep "[0-9]" | cut -d\  -f4)
echo "$picard_merge_sam_files_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# STEP: picard_sort_sam
#-------------------------------------------------------------------------------
STEP=picard_sort_sam
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: picard_sort_sam_1_JOB_ID: picard_sort_sam.DEX_shNIPBL_2
#-------------------------------------------------------------------------------
JOB_NAME=picard_sort_sam.DEX_shNIPBL_2
JOB_DEPENDENCIES=$star_1_JOB_ID
JOB_DONE=job_output/picard_sort_sam/picard_sort_sam.DEX_shNIPBL_2.b926075a49ce1786cc5b9b67137facce.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sort_sam.DEX_shNIPBL_2.b926075a49ce1786cc5b9b67137facce.mugqic.done'
module load java/1.8.0_121 mugqic/picard/2.9.0 && \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=4 -Dsamjdk.buffer_size=1048576 -Xmx40G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/DEX_shNIPBL_2/DEX_shNIPBL_2.sorted.bam \
 OUTPUT=alignment/DEX_shNIPBL_2/DEX_shNIPBL_2.QueryNameSorted.bam \
 SORT_ORDER=queryname \
 MAX_RECORDS_IN_RAM=5750000
picard_sort_sam.DEX_shNIPBL_2.b926075a49ce1786cc5b9b67137facce.mugqic.done
)
picard_sort_sam_1_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_sort_sam\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_sort_sam\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=48G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$picard_sort_sam_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: picard_sort_sam_2_JOB_ID: picard_sort_sam.DEX_shMED1_2
#-------------------------------------------------------------------------------
JOB_NAME=picard_sort_sam.DEX_shMED1_2
JOB_DEPENDENCIES=$star_2_JOB_ID
JOB_DONE=job_output/picard_sort_sam/picard_sort_sam.DEX_shMED1_2.a3aff1793dc437fb22660ac33ede1bd5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sort_sam.DEX_shMED1_2.a3aff1793dc437fb22660ac33ede1bd5.mugqic.done'
module load java/1.8.0_121 mugqic/picard/2.9.0 && \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=4 -Dsamjdk.buffer_size=1048576 -Xmx40G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/DEX_shMED1_2/DEX_shMED1_2.sorted.bam \
 OUTPUT=alignment/DEX_shMED1_2/DEX_shMED1_2.QueryNameSorted.bam \
 SORT_ORDER=queryname \
 MAX_RECORDS_IN_RAM=5750000
picard_sort_sam.DEX_shMED1_2.a3aff1793dc437fb22660ac33ede1bd5.mugqic.done
)
picard_sort_sam_2_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_sort_sam\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_sort_sam\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=48G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$picard_sort_sam_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: picard_sort_sam_3_JOB_ID: picard_sort_sam.DEX_shSMC1A_4
#-------------------------------------------------------------------------------
JOB_NAME=picard_sort_sam.DEX_shSMC1A_4
JOB_DEPENDENCIES=$star_3_JOB_ID
JOB_DONE=job_output/picard_sort_sam/picard_sort_sam.DEX_shSMC1A_4.021d4c59c9533a361aa9322c7769da04.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sort_sam.DEX_shSMC1A_4.021d4c59c9533a361aa9322c7769da04.mugqic.done'
module load java/1.8.0_121 mugqic/picard/2.9.0 && \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=4 -Dsamjdk.buffer_size=1048576 -Xmx40G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/DEX_shSMC1A_4/DEX_shSMC1A_4.sorted.bam \
 OUTPUT=alignment/DEX_shSMC1A_4/DEX_shSMC1A_4.QueryNameSorted.bam \
 SORT_ORDER=queryname \
 MAX_RECORDS_IN_RAM=5750000
picard_sort_sam.DEX_shSMC1A_4.021d4c59c9533a361aa9322c7769da04.mugqic.done
)
picard_sort_sam_3_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_sort_sam\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_sort_sam\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=48G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$picard_sort_sam_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: picard_sort_sam_4_JOB_ID: picard_sort_sam.DEX_nonMamm_2
#-------------------------------------------------------------------------------
JOB_NAME=picard_sort_sam.DEX_nonMamm_2
JOB_DEPENDENCIES=$star_4_JOB_ID
JOB_DONE=job_output/picard_sort_sam/picard_sort_sam.DEX_nonMamm_2.92fbf653003f66b8f917ca37fbdf4c68.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sort_sam.DEX_nonMamm_2.92fbf653003f66b8f917ca37fbdf4c68.mugqic.done'
module load java/1.8.0_121 mugqic/picard/2.9.0 && \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=4 -Dsamjdk.buffer_size=1048576 -Xmx40G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/DEX_nonMamm_2/DEX_nonMamm_2.sorted.bam \
 OUTPUT=alignment/DEX_nonMamm_2/DEX_nonMamm_2.QueryNameSorted.bam \
 SORT_ORDER=queryname \
 MAX_RECORDS_IN_RAM=5750000
picard_sort_sam.DEX_nonMamm_2.92fbf653003f66b8f917ca37fbdf4c68.mugqic.done
)
picard_sort_sam_4_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_sort_sam\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_sort_sam\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=48G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$picard_sort_sam_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: picard_sort_sam_5_JOB_ID: picard_sort_sam.DEX_nonMamm_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_sort_sam.DEX_nonMamm_1
JOB_DEPENDENCIES=$picard_merge_sam_files_1_JOB_ID
JOB_DONE=job_output/picard_sort_sam/picard_sort_sam.DEX_nonMamm_1.3120f823e40044e5ea933dbff9a319fa.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_sort_sam.DEX_nonMamm_1.3120f823e40044e5ea933dbff9a319fa.mugqic.done'
module load java/1.8.0_121 mugqic/picard/2.9.0 && \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=4 -Dsamjdk.buffer_size=1048576 -Xmx40G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.bam \
 OUTPUT=alignment/DEX_nonMamm_1/DEX_nonMamm_1.QueryNameSorted.bam \
 SORT_ORDER=queryname \
 MAX_RECORDS_IN_RAM=5750000
picard_sort_sam.DEX_nonMamm_1.3120f823e40044e5ea933dbff9a319fa.mugqic.done
)
picard_sort_sam_5_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_sort_sam\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_1.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_sort_sam\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_1.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=48G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$picard_sort_sam_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# STEP: picard_mark_duplicates
#-------------------------------------------------------------------------------
STEP=picard_mark_duplicates
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_1_JOB_ID: picard_mark_duplicates.DEX_shNIPBL_2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.DEX_shNIPBL_2
JOB_DEPENDENCIES=$star_1_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.DEX_shNIPBL_2.f34787f1d3cdeda4a77108b09a2247c9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.DEX_shNIPBL_2.f34787f1d3cdeda4a77108b09a2247c9.mugqic.done'
module load java/1.8.0_121 mugqic/picard/2.9.0 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=4 -Dsamjdk.buffer_size=1048576 -Xmx14G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/DEX_shNIPBL_2/DEX_shNIPBL_2.sorted.bam \
 OUTPUT=alignment/DEX_shNIPBL_2/DEX_shNIPBL_2.sorted.mdup.bam \
 METRICS_FILE=alignment/DEX_shNIPBL_2/DEX_shNIPBL_2.sorted.mdup.metrics \
 MAX_RECORDS_IN_RAM=3500000
picard_mark_duplicates.DEX_shNIPBL_2.f34787f1d3cdeda4a77108b09a2247c9.mugqic.done
)
picard_mark_duplicates_1_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_mark_duplicates\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_mark_duplicates\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=48:00:0 --mem=20G -N 1 -n 5 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$picard_mark_duplicates_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_2_JOB_ID: picard_mark_duplicates.DEX_shMED1_2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.DEX_shMED1_2
JOB_DEPENDENCIES=$star_2_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.DEX_shMED1_2.33ac3b331204ec14c3095208b6a8222f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.DEX_shMED1_2.33ac3b331204ec14c3095208b6a8222f.mugqic.done'
module load java/1.8.0_121 mugqic/picard/2.9.0 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=4 -Dsamjdk.buffer_size=1048576 -Xmx14G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/DEX_shMED1_2/DEX_shMED1_2.sorted.bam \
 OUTPUT=alignment/DEX_shMED1_2/DEX_shMED1_2.sorted.mdup.bam \
 METRICS_FILE=alignment/DEX_shMED1_2/DEX_shMED1_2.sorted.mdup.metrics \
 MAX_RECORDS_IN_RAM=3500000
picard_mark_duplicates.DEX_shMED1_2.33ac3b331204ec14c3095208b6a8222f.mugqic.done
)
picard_mark_duplicates_2_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_mark_duplicates\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_mark_duplicates\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=48:00:0 --mem=20G -N 1 -n 5 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$picard_mark_duplicates_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_3_JOB_ID: picard_mark_duplicates.DEX_shSMC1A_4
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.DEX_shSMC1A_4
JOB_DEPENDENCIES=$star_3_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.DEX_shSMC1A_4.b17ffd1ddda3856df37db5d10fbbca4d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.DEX_shSMC1A_4.b17ffd1ddda3856df37db5d10fbbca4d.mugqic.done'
module load java/1.8.0_121 mugqic/picard/2.9.0 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=4 -Dsamjdk.buffer_size=1048576 -Xmx14G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/DEX_shSMC1A_4/DEX_shSMC1A_4.sorted.bam \
 OUTPUT=alignment/DEX_shSMC1A_4/DEX_shSMC1A_4.sorted.mdup.bam \
 METRICS_FILE=alignment/DEX_shSMC1A_4/DEX_shSMC1A_4.sorted.mdup.metrics \
 MAX_RECORDS_IN_RAM=3500000
picard_mark_duplicates.DEX_shSMC1A_4.b17ffd1ddda3856df37db5d10fbbca4d.mugqic.done
)
picard_mark_duplicates_3_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_mark_duplicates\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_mark_duplicates\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=48:00:0 --mem=20G -N 1 -n 5 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$picard_mark_duplicates_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_4_JOB_ID: picard_mark_duplicates.DEX_nonMamm_2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.DEX_nonMamm_2
JOB_DEPENDENCIES=$star_4_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.DEX_nonMamm_2.7293865370216ee1634bd5fe0b85ffc6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.DEX_nonMamm_2.7293865370216ee1634bd5fe0b85ffc6.mugqic.done'
module load java/1.8.0_121 mugqic/picard/2.9.0 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=4 -Dsamjdk.buffer_size=1048576 -Xmx14G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/DEX_nonMamm_2/DEX_nonMamm_2.sorted.bam \
 OUTPUT=alignment/DEX_nonMamm_2/DEX_nonMamm_2.sorted.mdup.bam \
 METRICS_FILE=alignment/DEX_nonMamm_2/DEX_nonMamm_2.sorted.mdup.metrics \
 MAX_RECORDS_IN_RAM=3500000
picard_mark_duplicates.DEX_nonMamm_2.7293865370216ee1634bd5fe0b85ffc6.mugqic.done
)
picard_mark_duplicates_4_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_mark_duplicates\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_mark_duplicates\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=48:00:0 --mem=20G -N 1 -n 5 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$picard_mark_duplicates_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_5_JOB_ID: picard_mark_duplicates.DEX_nonMamm_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.DEX_nonMamm_1
JOB_DEPENDENCIES=$picard_merge_sam_files_1_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.DEX_nonMamm_1.0dbf63f12f06a87386607f6ffa49bb11.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.DEX_nonMamm_1.0dbf63f12f06a87386607f6ffa49bb11.mugqic.done'
module load java/1.8.0_121 mugqic/picard/2.9.0 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=4 -Dsamjdk.buffer_size=1048576 -Xmx14G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.bam \
 OUTPUT=alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.bam \
 METRICS_FILE=alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.metrics \
 MAX_RECORDS_IN_RAM=3500000
picard_mark_duplicates.DEX_nonMamm_1.0dbf63f12f06a87386607f6ffa49bb11.mugqic.done
)
picard_mark_duplicates_5_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_mark_duplicates\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_1.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_mark_duplicates\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_1.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=48:00:0 --mem=20G -N 1 -n 5 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$picard_mark_duplicates_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# STEP: picard_rna_metrics
#-------------------------------------------------------------------------------
STEP=picard_rna_metrics
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: picard_rna_metrics_1_JOB_ID: picard_rna_metrics.DEX_shNIPBL_2
#-------------------------------------------------------------------------------
JOB_NAME=picard_rna_metrics.DEX_shNIPBL_2
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID
JOB_DONE=job_output/picard_rna_metrics/picard_rna_metrics.DEX_shNIPBL_2.8213251a9eacace8b67ae1878aa82fc9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_rna_metrics.DEX_shNIPBL_2.8213251a9eacace8b67ae1878aa82fc9.mugqic.done'
module load java/1.8.0_121 mugqic/picard/2.9.0 mugqic/R_Bioconductor/3.4.2_3.6 && \
mkdir -p metrics/DEX_shNIPBL_2 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=4 -Dsamjdk.buffer_size=1048576 -Xmx40G -jar $PICARD_HOME/picard.jar CollectMultipleMetrics \
 PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT \
 TMP_DIR=${SLURM_TMPDIR} \
 REFERENCE_SEQUENCE=/project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
 INPUT=alignment/DEX_shNIPBL_2/DEX_shNIPBL_2.sorted.mdup.bam \
 OUTPUT=metrics/DEX_shNIPBL_2/DEX_shNIPBL_2 \
 MAX_RECORDS_IN_RAM=5750000 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=4 -Dsamjdk.buffer_size=1048576 -Xmx40G -jar $PICARD_HOME/picard.jar CollectRnaSeqMetrics \
 VALIDATION_STRINGENCY=SILENT  \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/DEX_shNIPBL_2/DEX_shNIPBL_2.sorted.mdup.bam \
 OUTPUT=metrics/DEX_shNIPBL_2/DEX_shNIPBL_2.picard_rna_metrics \
 REF_FLAT=$MUGQIC_INSTALL_HOME_DEV/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.ref_flat.tsv \
 STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
 MINIMUM_LENGTH=200 \
 REFERENCE_SEQUENCE=$MUGQIC_INSTALL_HOME_DEV/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
 MAX_RECORDS_IN_RAM=5750000
picard_rna_metrics.DEX_shNIPBL_2.8213251a9eacace8b67ae1878aa82fc9.mugqic.done
)
picard_rna_metrics_1_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_rna_metrics\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_rna_metrics\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=48G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$picard_rna_metrics_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: picard_rna_metrics_2_JOB_ID: picard_rna_metrics.DEX_shMED1_2
#-------------------------------------------------------------------------------
JOB_NAME=picard_rna_metrics.DEX_shMED1_2
JOB_DEPENDENCIES=$picard_mark_duplicates_2_JOB_ID
JOB_DONE=job_output/picard_rna_metrics/picard_rna_metrics.DEX_shMED1_2.06ff5d74cadd6075459b31760d045ae9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_rna_metrics.DEX_shMED1_2.06ff5d74cadd6075459b31760d045ae9.mugqic.done'
module load java/1.8.0_121 mugqic/picard/2.9.0 mugqic/R_Bioconductor/3.4.2_3.6 && \
mkdir -p metrics/DEX_shMED1_2 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=4 -Dsamjdk.buffer_size=1048576 -Xmx40G -jar $PICARD_HOME/picard.jar CollectMultipleMetrics \
 PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT \
 TMP_DIR=${SLURM_TMPDIR} \
 REFERENCE_SEQUENCE=/project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
 INPUT=alignment/DEX_shMED1_2/DEX_shMED1_2.sorted.mdup.bam \
 OUTPUT=metrics/DEX_shMED1_2/DEX_shMED1_2 \
 MAX_RECORDS_IN_RAM=5750000 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=4 -Dsamjdk.buffer_size=1048576 -Xmx40G -jar $PICARD_HOME/picard.jar CollectRnaSeqMetrics \
 VALIDATION_STRINGENCY=SILENT  \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/DEX_shMED1_2/DEX_shMED1_2.sorted.mdup.bam \
 OUTPUT=metrics/DEX_shMED1_2/DEX_shMED1_2.picard_rna_metrics \
 REF_FLAT=$MUGQIC_INSTALL_HOME_DEV/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.ref_flat.tsv \
 STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
 MINIMUM_LENGTH=200 \
 REFERENCE_SEQUENCE=$MUGQIC_INSTALL_HOME_DEV/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
 MAX_RECORDS_IN_RAM=5750000
picard_rna_metrics.DEX_shMED1_2.06ff5d74cadd6075459b31760d045ae9.mugqic.done
)
picard_rna_metrics_2_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_rna_metrics\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_rna_metrics\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=48G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$picard_rna_metrics_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: picard_rna_metrics_3_JOB_ID: picard_rna_metrics.DEX_shSMC1A_4
#-------------------------------------------------------------------------------
JOB_NAME=picard_rna_metrics.DEX_shSMC1A_4
JOB_DEPENDENCIES=$picard_mark_duplicates_3_JOB_ID
JOB_DONE=job_output/picard_rna_metrics/picard_rna_metrics.DEX_shSMC1A_4.3de1547eb5aea8149a929d2556775b9c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_rna_metrics.DEX_shSMC1A_4.3de1547eb5aea8149a929d2556775b9c.mugqic.done'
module load java/1.8.0_121 mugqic/picard/2.9.0 mugqic/R_Bioconductor/3.4.2_3.6 && \
mkdir -p metrics/DEX_shSMC1A_4 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=4 -Dsamjdk.buffer_size=1048576 -Xmx40G -jar $PICARD_HOME/picard.jar CollectMultipleMetrics \
 PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT \
 TMP_DIR=${SLURM_TMPDIR} \
 REFERENCE_SEQUENCE=/project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
 INPUT=alignment/DEX_shSMC1A_4/DEX_shSMC1A_4.sorted.mdup.bam \
 OUTPUT=metrics/DEX_shSMC1A_4/DEX_shSMC1A_4 \
 MAX_RECORDS_IN_RAM=5750000 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=4 -Dsamjdk.buffer_size=1048576 -Xmx40G -jar $PICARD_HOME/picard.jar CollectRnaSeqMetrics \
 VALIDATION_STRINGENCY=SILENT  \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/DEX_shSMC1A_4/DEX_shSMC1A_4.sorted.mdup.bam \
 OUTPUT=metrics/DEX_shSMC1A_4/DEX_shSMC1A_4.picard_rna_metrics \
 REF_FLAT=$MUGQIC_INSTALL_HOME_DEV/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.ref_flat.tsv \
 STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
 MINIMUM_LENGTH=200 \
 REFERENCE_SEQUENCE=$MUGQIC_INSTALL_HOME_DEV/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
 MAX_RECORDS_IN_RAM=5750000
picard_rna_metrics.DEX_shSMC1A_4.3de1547eb5aea8149a929d2556775b9c.mugqic.done
)
picard_rna_metrics_3_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_rna_metrics\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_rna_metrics\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=48G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$picard_rna_metrics_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: picard_rna_metrics_4_JOB_ID: picard_rna_metrics.DEX_nonMamm_2
#-------------------------------------------------------------------------------
JOB_NAME=picard_rna_metrics.DEX_nonMamm_2
JOB_DEPENDENCIES=$picard_mark_duplicates_4_JOB_ID
JOB_DONE=job_output/picard_rna_metrics/picard_rna_metrics.DEX_nonMamm_2.8f985d389029d6bffd96928b542a97b4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_rna_metrics.DEX_nonMamm_2.8f985d389029d6bffd96928b542a97b4.mugqic.done'
module load java/1.8.0_121 mugqic/picard/2.9.0 mugqic/R_Bioconductor/3.4.2_3.6 && \
mkdir -p metrics/DEX_nonMamm_2 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=4 -Dsamjdk.buffer_size=1048576 -Xmx40G -jar $PICARD_HOME/picard.jar CollectMultipleMetrics \
 PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT \
 TMP_DIR=${SLURM_TMPDIR} \
 REFERENCE_SEQUENCE=/project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
 INPUT=alignment/DEX_nonMamm_2/DEX_nonMamm_2.sorted.mdup.bam \
 OUTPUT=metrics/DEX_nonMamm_2/DEX_nonMamm_2 \
 MAX_RECORDS_IN_RAM=5750000 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=4 -Dsamjdk.buffer_size=1048576 -Xmx40G -jar $PICARD_HOME/picard.jar CollectRnaSeqMetrics \
 VALIDATION_STRINGENCY=SILENT  \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/DEX_nonMamm_2/DEX_nonMamm_2.sorted.mdup.bam \
 OUTPUT=metrics/DEX_nonMamm_2/DEX_nonMamm_2.picard_rna_metrics \
 REF_FLAT=$MUGQIC_INSTALL_HOME_DEV/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.ref_flat.tsv \
 STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
 MINIMUM_LENGTH=200 \
 REFERENCE_SEQUENCE=$MUGQIC_INSTALL_HOME_DEV/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
 MAX_RECORDS_IN_RAM=5750000
picard_rna_metrics.DEX_nonMamm_2.8f985d389029d6bffd96928b542a97b4.mugqic.done
)
picard_rna_metrics_4_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_rna_metrics\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_rna_metrics\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=48G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$picard_rna_metrics_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: picard_rna_metrics_5_JOB_ID: picard_rna_metrics.DEX_nonMamm_1
#-------------------------------------------------------------------------------
JOB_NAME=picard_rna_metrics.DEX_nonMamm_1
JOB_DEPENDENCIES=$picard_mark_duplicates_5_JOB_ID
JOB_DONE=job_output/picard_rna_metrics/picard_rna_metrics.DEX_nonMamm_1.8b6cccf3e2525af6e7e0d3a1bc4c2b45.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_rna_metrics.DEX_nonMamm_1.8b6cccf3e2525af6e7e0d3a1bc4c2b45.mugqic.done'
module load java/1.8.0_121 mugqic/picard/2.9.0 mugqic/R_Bioconductor/3.4.2_3.6 && \
mkdir -p metrics/DEX_nonMamm_1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=4 -Dsamjdk.buffer_size=1048576 -Xmx40G -jar $PICARD_HOME/picard.jar CollectMultipleMetrics \
 PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics VALIDATION_STRINGENCY=SILENT \
 TMP_DIR=${SLURM_TMPDIR} \
 REFERENCE_SEQUENCE=/project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
 INPUT=alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.bam \
 OUTPUT=metrics/DEX_nonMamm_1/DEX_nonMamm_1 \
 MAX_RECORDS_IN_RAM=5750000 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=4 -Dsamjdk.buffer_size=1048576 -Xmx40G -jar $PICARD_HOME/picard.jar CollectRnaSeqMetrics \
 VALIDATION_STRINGENCY=SILENT  \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.bam \
 OUTPUT=metrics/DEX_nonMamm_1/DEX_nonMamm_1.picard_rna_metrics \
 REF_FLAT=$MUGQIC_INSTALL_HOME_DEV/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.ref_flat.tsv \
 STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \
 MINIMUM_LENGTH=200 \
 REFERENCE_SEQUENCE=$MUGQIC_INSTALL_HOME_DEV/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
 MAX_RECORDS_IN_RAM=5750000
picard_rna_metrics.DEX_nonMamm_1.8b6cccf3e2525af6e7e0d3a1bc4c2b45.mugqic.done
)
picard_rna_metrics_5_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_rna_metrics\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_1.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"picard_rna_metrics\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_1.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=48G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$picard_rna_metrics_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# STEP: estimate_ribosomal_rna
#-------------------------------------------------------------------------------
STEP=estimate_ribosomal_rna
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: estimate_ribosomal_rna_1_JOB_ID: bwa_mem_rRNA.HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_rRNA.HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2
JOB_DEPENDENCIES=$star_1_JOB_ID
JOB_DONE=job_output/estimate_ribosomal_rna/bwa_mem_rRNA.HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2.52bdcc4c8771fe5d76124d57592cc56e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_rRNA.HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2.52bdcc4c8771fe5d76124d57592cc56e.mugqic.done'
module load java/1.8.0_121 mugqic/bvatools/1.6 mugqic/bwa/0.7.15 mugqic/picard/2.9.0 mugqic/mugqic_tools/2.1.9 mugqic/python/2.7.13 && \
mkdir -p alignment/DEX_shNIPBL_2/HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2 metrics/DEX_shNIPBL_2/HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2 && \
java -XX:ParallelGCThreads=4 -Dsamjdk.buffer_size=1048576 -Xmx10G -jar $BVATOOLS_JAR \
  bam2fq --mapped ONLY \
  --bam alignment/DEX_shNIPBL_2/HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2/Aligned.sortedByCoord.out.bam    | \
bwa mem  \
  -M -t 10 \
  -R '@RG	ID:HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2	SM:DEX_shNIPBL_2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  $MUGQIC_INSTALL_HOME_DEV/genomes/species/Homo_sapiens.GRCh38/annotations/rrna_bwa_index/Homo_sapiens.GRCh38.Ensembl87.rrna.fa \
  /dev/stdin | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1  -Dsamjdk.buffer_size=1048576 -Xmx7G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=metrics/DEX_shNIPBL_2/HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2/HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2rRNA.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=1750000 && \
python $PYTHON_TOOLS/rrnaBAMcounter.py \
  -i metrics/DEX_shNIPBL_2/HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2/HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2rRNA.bam \
  -g $MUGQIC_INSTALL_HOME_DEV/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  -o metrics/DEX_shNIPBL_2/HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2/HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2rRNA.stats.tsv \
  -t transcript
bwa_mem_rRNA.HI.4620.001.NEBNext_Index_12.A549_RNA_NA_shNIPBL-3_Dex_Rep2.52bdcc4c8771fe5d76124d57592cc56e.mugqic.done
)
estimate_ribosomal_rna_1_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"estimate_ribosomal_rna\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"estimate_ribosomal_rna\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=64G -N 1 -n 16 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$estimate_ribosomal_rna_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: estimate_ribosomal_rna_2_JOB_ID: bwa_mem_rRNA.HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_rRNA.HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2
JOB_DEPENDENCIES=$star_2_JOB_ID
JOB_DONE=job_output/estimate_ribosomal_rna/bwa_mem_rRNA.HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2.df7e3b912b16aa61b954f6996d66179b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_rRNA.HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2.df7e3b912b16aa61b954f6996d66179b.mugqic.done'
module load java/1.8.0_121 mugqic/bvatools/1.6 mugqic/bwa/0.7.15 mugqic/picard/2.9.0 mugqic/mugqic_tools/2.1.9 mugqic/python/2.7.13 && \
mkdir -p alignment/DEX_shMED1_2/HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2 metrics/DEX_shMED1_2/HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2 && \
java -XX:ParallelGCThreads=4 -Dsamjdk.buffer_size=1048576 -Xmx10G -jar $BVATOOLS_JAR \
  bam2fq --mapped ONLY \
  --bam alignment/DEX_shMED1_2/HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2/Aligned.sortedByCoord.out.bam    | \
bwa mem  \
  -M -t 10 \
  -R '@RG	ID:HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2	SM:DEX_shMED1_2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  $MUGQIC_INSTALL_HOME_DEV/genomes/species/Homo_sapiens.GRCh38/annotations/rrna_bwa_index/Homo_sapiens.GRCh38.Ensembl87.rrna.fa \
  /dev/stdin | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1  -Dsamjdk.buffer_size=1048576 -Xmx7G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=metrics/DEX_shMED1_2/HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2/HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2rRNA.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=1750000 && \
python $PYTHON_TOOLS/rrnaBAMcounter.py \
  -i metrics/DEX_shMED1_2/HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2/HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2rRNA.bam \
  -g $MUGQIC_INSTALL_HOME_DEV/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  -o metrics/DEX_shMED1_2/HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2/HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2rRNA.stats.tsv \
  -t transcript
bwa_mem_rRNA.HI.4620.001.NEBNext_Index_19.A549_RNA_NA_shMED1-2_Dex_Rep2.df7e3b912b16aa61b954f6996d66179b.mugqic.done
)
estimate_ribosomal_rna_2_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"estimate_ribosomal_rna\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"estimate_ribosomal_rna\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=64G -N 1 -n 16 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$estimate_ribosomal_rna_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: estimate_ribosomal_rna_3_JOB_ID: bwa_mem_rRNA.HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_rRNA.HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2
JOB_DEPENDENCIES=$star_3_JOB_ID
JOB_DONE=job_output/estimate_ribosomal_rna/bwa_mem_rRNA.HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2.274cad3d56f6888235ebfdbdc5c12cf8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_rRNA.HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2.274cad3d56f6888235ebfdbdc5c12cf8.mugqic.done'
module load java/1.8.0_121 mugqic/bvatools/1.6 mugqic/bwa/0.7.15 mugqic/picard/2.9.0 mugqic/mugqic_tools/2.1.9 mugqic/python/2.7.13 && \
mkdir -p alignment/DEX_shSMC1A_4/HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2 metrics/DEX_shSMC1A_4/HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2 && \
java -XX:ParallelGCThreads=4 -Dsamjdk.buffer_size=1048576 -Xmx10G -jar $BVATOOLS_JAR \
  bam2fq --mapped ONLY \
  --bam alignment/DEX_shSMC1A_4/HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2/Aligned.sortedByCoord.out.bam    | \
bwa mem  \
  -M -t 10 \
  -R '@RG	ID:HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2	SM:DEX_shSMC1A_4	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  $MUGQIC_INSTALL_HOME_DEV/genomes/species/Homo_sapiens.GRCh38/annotations/rrna_bwa_index/Homo_sapiens.GRCh38.Ensembl87.rrna.fa \
  /dev/stdin | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1  -Dsamjdk.buffer_size=1048576 -Xmx7G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=metrics/DEX_shSMC1A_4/HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2/HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2rRNA.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=1750000 && \
python $PYTHON_TOOLS/rrnaBAMcounter.py \
  -i metrics/DEX_shSMC1A_4/HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2/HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2rRNA.bam \
  -g $MUGQIC_INSTALL_HOME_DEV/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  -o metrics/DEX_shSMC1A_4/HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2/HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2rRNA.stats.tsv \
  -t transcript
bwa_mem_rRNA.HI.4620.001.NEBNext_Index_1.A549_RNA_NA_shSMC1A-2_Dex_Rep2.274cad3d56f6888235ebfdbdc5c12cf8.mugqic.done
)
estimate_ribosomal_rna_3_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"estimate_ribosomal_rna\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"estimate_ribosomal_rna\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=64G -N 1 -n 16 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$estimate_ribosomal_rna_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: estimate_ribosomal_rna_4_JOB_ID: bwa_mem_rRNA.HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_rRNA.HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2
JOB_DEPENDENCIES=$star_4_JOB_ID
JOB_DONE=job_output/estimate_ribosomal_rna/bwa_mem_rRNA.HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2.fd73b6b0b4dd8a78ef360a9af7f3ebdf.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_rRNA.HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2.fd73b6b0b4dd8a78ef360a9af7f3ebdf.mugqic.done'
module load java/1.8.0_121 mugqic/bvatools/1.6 mugqic/bwa/0.7.15 mugqic/picard/2.9.0 mugqic/mugqic_tools/2.1.9 mugqic/python/2.7.13 && \
mkdir -p alignment/DEX_nonMamm_2/HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2 metrics/DEX_nonMamm_2/HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2 && \
java -XX:ParallelGCThreads=4 -Dsamjdk.buffer_size=1048576 -Xmx10G -jar $BVATOOLS_JAR \
  bam2fq --mapped ONLY \
  --bam alignment/DEX_nonMamm_2/HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2/Aligned.sortedByCoord.out.bam    | \
bwa mem  \
  -M -t 10 \
  -R '@RG	ID:HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2	SM:DEX_nonMamm_2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  $MUGQIC_INSTALL_HOME_DEV/genomes/species/Homo_sapiens.GRCh38/annotations/rrna_bwa_index/Homo_sapiens.GRCh38.Ensembl87.rrna.fa \
  /dev/stdin | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1  -Dsamjdk.buffer_size=1048576 -Xmx7G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=metrics/DEX_nonMamm_2/HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2/HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2rRNA.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=1750000 && \
python $PYTHON_TOOLS/rrnaBAMcounter.py \
  -i metrics/DEX_nonMamm_2/HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2/HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2rRNA.bam \
  -g $MUGQIC_INSTALL_HOME_DEV/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  -o metrics/DEX_nonMamm_2/HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2/HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2rRNA.stats.tsv \
  -t transcript
bwa_mem_rRNA.HI.4620.001.NEBNext_Index_5.A549_RNA_NA_shCRTL-2_Dex_Rep2.fd73b6b0b4dd8a78ef360a9af7f3ebdf.mugqic.done
)
estimate_ribosomal_rna_4_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"estimate_ribosomal_rna\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"estimate_ribosomal_rna\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=64G -N 1 -n 16 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$estimate_ribosomal_rna_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: estimate_ribosomal_rna_5_JOB_ID: bwa_mem_rRNA.DEX_nonMamm_1_1
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_rRNA.DEX_nonMamm_1_1
JOB_DEPENDENCIES=
JOB_DONE=job_output/estimate_ribosomal_rna/bwa_mem_rRNA.DEX_nonMamm_1_1.121004d448cbc5b43ea84f2a792dd27b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_rRNA.DEX_nonMamm_1_1.121004d448cbc5b43ea84f2a792dd27b.mugqic.done'
module load java/1.8.0_121 mugqic/bvatools/1.6 mugqic/bwa/0.7.15 mugqic/picard/2.9.0 mugqic/mugqic_tools/2.1.9 mugqic/python/2.7.13 && \
mkdir -p alignment/DEX_nonMamm_1/DEX_nonMamm_1_1 metrics/DEX_nonMamm_1/DEX_nonMamm_1_1 && \
java -XX:ParallelGCThreads=4 -Dsamjdk.buffer_size=1048576 -Xmx10G -jar $BVATOOLS_JAR \
  bam2fq --mapped ONLY \
  --bam alignment/DEX_nonMamm_1/DEX_nonMamm_1_1/Aligned.sortedByCoord.out.bam    | \
bwa mem  \
  -M -t 10 \
  -R '@RG	ID:DEX_nonMamm_1_1	SM:DEX_nonMamm_1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  $MUGQIC_INSTALL_HOME_DEV/genomes/species/Homo_sapiens.GRCh38/annotations/rrna_bwa_index/Homo_sapiens.GRCh38.Ensembl87.rrna.fa \
  /dev/stdin | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1  -Dsamjdk.buffer_size=1048576 -Xmx7G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=metrics/DEX_nonMamm_1/DEX_nonMamm_1_1/DEX_nonMamm_1_1rRNA.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=1750000 && \
python $PYTHON_TOOLS/rrnaBAMcounter.py \
  -i metrics/DEX_nonMamm_1/DEX_nonMamm_1_1/DEX_nonMamm_1_1rRNA.bam \
  -g $MUGQIC_INSTALL_HOME_DEV/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  -o metrics/DEX_nonMamm_1/DEX_nonMamm_1_1/DEX_nonMamm_1_1rRNA.stats.tsv \
  -t transcript
bwa_mem_rRNA.DEX_nonMamm_1_1.121004d448cbc5b43ea84f2a792dd27b.mugqic.done
)
estimate_ribosomal_rna_5_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"estimate_ribosomal_rna\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_1.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"estimate_ribosomal_rna\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_1.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=64G -N 1 -n 16 | grep "[0-9]" | cut -d\  -f4)
echo "$estimate_ribosomal_rna_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: estimate_ribosomal_rna_6_JOB_ID: bwa_mem_rRNA.DEX_nonMamm_1_2
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_rRNA.DEX_nonMamm_1_2
JOB_DEPENDENCIES=
JOB_DONE=job_output/estimate_ribosomal_rna/bwa_mem_rRNA.DEX_nonMamm_1_2.5bed0b1a8d0cbc367bd2a7931c98ee4a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_rRNA.DEX_nonMamm_1_2.5bed0b1a8d0cbc367bd2a7931c98ee4a.mugqic.done'
module load java/1.8.0_121 mugqic/bvatools/1.6 mugqic/bwa/0.7.15 mugqic/picard/2.9.0 mugqic/mugqic_tools/2.1.9 mugqic/python/2.7.13 && \
mkdir -p alignment/DEX_nonMamm_1/DEX_nonMamm_1_2 metrics/DEX_nonMamm_1/DEX_nonMamm_1_2 && \
java -XX:ParallelGCThreads=4 -Dsamjdk.buffer_size=1048576 -Xmx10G -jar $BVATOOLS_JAR \
  bam2fq --mapped ONLY \
  --bam alignment/DEX_nonMamm_1/DEX_nonMamm_1_2/Aligned.sortedByCoord.out.bam    | \
bwa mem  \
  -M -t 10 \
  -R '@RG	ID:DEX_nonMamm_1_2	SM:DEX_nonMamm_1	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  $MUGQIC_INSTALL_HOME_DEV/genomes/species/Homo_sapiens.GRCh38/annotations/rrna_bwa_index/Homo_sapiens.GRCh38.Ensembl87.rrna.fa \
  /dev/stdin | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1  -Dsamjdk.buffer_size=1048576 -Xmx7G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=metrics/DEX_nonMamm_1/DEX_nonMamm_1_2/DEX_nonMamm_1_2rRNA.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=1750000 && \
python $PYTHON_TOOLS/rrnaBAMcounter.py \
  -i metrics/DEX_nonMamm_1/DEX_nonMamm_1_2/DEX_nonMamm_1_2rRNA.bam \
  -g $MUGQIC_INSTALL_HOME_DEV/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  -o metrics/DEX_nonMamm_1/DEX_nonMamm_1_2/DEX_nonMamm_1_2rRNA.stats.tsv \
  -t transcript
bwa_mem_rRNA.DEX_nonMamm_1_2.5bed0b1a8d0cbc367bd2a7931c98ee4a.mugqic.done
)
estimate_ribosomal_rna_6_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"estimate_ribosomal_rna\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_1.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"estimate_ribosomal_rna\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_1.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=64G -N 1 -n 16 | grep "[0-9]" | cut -d\  -f4)
echo "$estimate_ribosomal_rna_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# STEP: bam_hard_clip
#-------------------------------------------------------------------------------
STEP=bam_hard_clip
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: bam_hard_clip_1_JOB_ID: tuxedo_hard_clip.DEX_shNIPBL_2
#-------------------------------------------------------------------------------
JOB_NAME=tuxedo_hard_clip.DEX_shNIPBL_2
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID
JOB_DONE=job_output/bam_hard_clip/tuxedo_hard_clip.DEX_shNIPBL_2.7a845305f786d1622bcb41e4ab987b85.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'tuxedo_hard_clip.DEX_shNIPBL_2.7a845305f786d1622bcb41e4ab987b85.mugqic.done'
module load mugqic/samtools/1.4.1 && \
samtools view -h \
  alignment/DEX_shNIPBL_2/DEX_shNIPBL_2.sorted.mdup.bam | \
awk 'BEGIN {OFS="\t"} {if (substr($1,1,1)=="@") {print;next}; split($6,C,/[0-9]*/); split($6,L,/[SMDIN]/); if (C[2]=="S") {$10=substr($10,L[1]+1); $11=substr($11,L[1]+1)}; if (C[length(C)]=="S") {L1=length($10)-L[length(L)-1]; $10=substr($10,1,L1); $11=substr($11,1,L1); }; gsub(/[0-9]*S/,"",$6); print}'  | \
samtools view -hbS \
  - \
  > alignment/DEX_shNIPBL_2/DEX_shNIPBL_2.sorted.mdup.hardClip.bam
tuxedo_hard_clip.DEX_shNIPBL_2.7a845305f786d1622bcb41e4ab987b85.mugqic.done
)
bam_hard_clip_1_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"bam_hard_clip\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"bam_hard_clip\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=32G -N 1 -n 8 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$bam_hard_clip_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: bam_hard_clip_2_JOB_ID: tuxedo_hard_clip.DEX_shMED1_2
#-------------------------------------------------------------------------------
JOB_NAME=tuxedo_hard_clip.DEX_shMED1_2
JOB_DEPENDENCIES=$picard_mark_duplicates_2_JOB_ID
JOB_DONE=job_output/bam_hard_clip/tuxedo_hard_clip.DEX_shMED1_2.038dd80769c0a6711dfb75d1968cbef5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'tuxedo_hard_clip.DEX_shMED1_2.038dd80769c0a6711dfb75d1968cbef5.mugqic.done'
module load mugqic/samtools/1.4.1 && \
samtools view -h \
  alignment/DEX_shMED1_2/DEX_shMED1_2.sorted.mdup.bam | \
awk 'BEGIN {OFS="\t"} {if (substr($1,1,1)=="@") {print;next}; split($6,C,/[0-9]*/); split($6,L,/[SMDIN]/); if (C[2]=="S") {$10=substr($10,L[1]+1); $11=substr($11,L[1]+1)}; if (C[length(C)]=="S") {L1=length($10)-L[length(L)-1]; $10=substr($10,1,L1); $11=substr($11,1,L1); }; gsub(/[0-9]*S/,"",$6); print}'  | \
samtools view -hbS \
  - \
  > alignment/DEX_shMED1_2/DEX_shMED1_2.sorted.mdup.hardClip.bam
tuxedo_hard_clip.DEX_shMED1_2.038dd80769c0a6711dfb75d1968cbef5.mugqic.done
)
bam_hard_clip_2_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"bam_hard_clip\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"bam_hard_clip\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=32G -N 1 -n 8 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$bam_hard_clip_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: bam_hard_clip_3_JOB_ID: tuxedo_hard_clip.DEX_shSMC1A_4
#-------------------------------------------------------------------------------
JOB_NAME=tuxedo_hard_clip.DEX_shSMC1A_4
JOB_DEPENDENCIES=$picard_mark_duplicates_3_JOB_ID
JOB_DONE=job_output/bam_hard_clip/tuxedo_hard_clip.DEX_shSMC1A_4.1ac9f3f2d4f87e095e6f1c9c2c9f3d7d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'tuxedo_hard_clip.DEX_shSMC1A_4.1ac9f3f2d4f87e095e6f1c9c2c9f3d7d.mugqic.done'
module load mugqic/samtools/1.4.1 && \
samtools view -h \
  alignment/DEX_shSMC1A_4/DEX_shSMC1A_4.sorted.mdup.bam | \
awk 'BEGIN {OFS="\t"} {if (substr($1,1,1)=="@") {print;next}; split($6,C,/[0-9]*/); split($6,L,/[SMDIN]/); if (C[2]=="S") {$10=substr($10,L[1]+1); $11=substr($11,L[1]+1)}; if (C[length(C)]=="S") {L1=length($10)-L[length(L)-1]; $10=substr($10,1,L1); $11=substr($11,1,L1); }; gsub(/[0-9]*S/,"",$6); print}'  | \
samtools view -hbS \
  - \
  > alignment/DEX_shSMC1A_4/DEX_shSMC1A_4.sorted.mdup.hardClip.bam
tuxedo_hard_clip.DEX_shSMC1A_4.1ac9f3f2d4f87e095e6f1c9c2c9f3d7d.mugqic.done
)
bam_hard_clip_3_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"bam_hard_clip\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"bam_hard_clip\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=32G -N 1 -n 8 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$bam_hard_clip_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: bam_hard_clip_4_JOB_ID: tuxedo_hard_clip.DEX_nonMamm_2
#-------------------------------------------------------------------------------
JOB_NAME=tuxedo_hard_clip.DEX_nonMamm_2
JOB_DEPENDENCIES=$picard_mark_duplicates_4_JOB_ID
JOB_DONE=job_output/bam_hard_clip/tuxedo_hard_clip.DEX_nonMamm_2.5660eb08e8244fcff5867784426e86b4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'tuxedo_hard_clip.DEX_nonMamm_2.5660eb08e8244fcff5867784426e86b4.mugqic.done'
module load mugqic/samtools/1.4.1 && \
samtools view -h \
  alignment/DEX_nonMamm_2/DEX_nonMamm_2.sorted.mdup.bam | \
awk 'BEGIN {OFS="\t"} {if (substr($1,1,1)=="@") {print;next}; split($6,C,/[0-9]*/); split($6,L,/[SMDIN]/); if (C[2]=="S") {$10=substr($10,L[1]+1); $11=substr($11,L[1]+1)}; if (C[length(C)]=="S") {L1=length($10)-L[length(L)-1]; $10=substr($10,1,L1); $11=substr($11,1,L1); }; gsub(/[0-9]*S/,"",$6); print}'  | \
samtools view -hbS \
  - \
  > alignment/DEX_nonMamm_2/DEX_nonMamm_2.sorted.mdup.hardClip.bam
tuxedo_hard_clip.DEX_nonMamm_2.5660eb08e8244fcff5867784426e86b4.mugqic.done
)
bam_hard_clip_4_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"bam_hard_clip\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"bam_hard_clip\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=32G -N 1 -n 8 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$bam_hard_clip_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: bam_hard_clip_5_JOB_ID: tuxedo_hard_clip.DEX_nonMamm_1
#-------------------------------------------------------------------------------
JOB_NAME=tuxedo_hard_clip.DEX_nonMamm_1
JOB_DEPENDENCIES=$picard_mark_duplicates_5_JOB_ID
JOB_DONE=job_output/bam_hard_clip/tuxedo_hard_clip.DEX_nonMamm_1.327c6f9596282e972fcfa8d6558912c9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'tuxedo_hard_clip.DEX_nonMamm_1.327c6f9596282e972fcfa8d6558912c9.mugqic.done'
module load mugqic/samtools/1.4.1 && \
samtools view -h \
  alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.bam | \
awk 'BEGIN {OFS="\t"} {if (substr($1,1,1)=="@") {print;next}; split($6,C,/[0-9]*/); split($6,L,/[SMDIN]/); if (C[2]=="S") {$10=substr($10,L[1]+1); $11=substr($11,L[1]+1)}; if (C[length(C)]=="S") {L1=length($10)-L[length(L)-1]; $10=substr($10,1,L1); $11=substr($11,1,L1); }; gsub(/[0-9]*S/,"",$6); print}'  | \
samtools view -hbS \
  - \
  > alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.hardClip.bam
tuxedo_hard_clip.DEX_nonMamm_1.327c6f9596282e972fcfa8d6558912c9.mugqic.done
)
bam_hard_clip_5_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"bam_hard_clip\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_1.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"bam_hard_clip\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_1.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=32G -N 1 -n 8 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$bam_hard_clip_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# STEP: rnaseqc
#-------------------------------------------------------------------------------
STEP=rnaseqc
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: rnaseqc_1_JOB_ID: rnaseqc
#-------------------------------------------------------------------------------
JOB_NAME=rnaseqc
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_2_JOB_ID:$picard_mark_duplicates_3_JOB_ID:$picard_mark_duplicates_4_JOB_ID:$picard_mark_duplicates_5_JOB_ID
JOB_DONE=job_output/rnaseqc/rnaseqc.0bbf9099a1d23b1dad3b3e916996e204.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'rnaseqc.0bbf9099a1d23b1dad3b3e916996e204.mugqic.done'
module load mugqic/java/openjdk-jdk1.7.0_60 mugqic/bwa/0.7.15 mugqic/rnaseqc/1.1.8 && \
mkdir -p metrics/rnaseqRep && \
echo "Sample	BamFile	Note
DEX_shNIPBL_2	alignment/DEX_shNIPBL_2/DEX_shNIPBL_2.sorted.mdup.bam	RNAseq
DEX_shMED1_2	alignment/DEX_shMED1_2/DEX_shMED1_2.sorted.mdup.bam	RNAseq
DEX_shSMC1A_4	alignment/DEX_shSMC1A_4/DEX_shSMC1A_4.sorted.mdup.bam	RNAseq
DEX_nonMamm_2	alignment/DEX_nonMamm_2/DEX_nonMamm_2.sorted.mdup.bam	RNAseq
DEX_nonMamm_1	alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.bam	RNAseq" \
  > alignment/rnaseqc.samples.txt && \
touch dummy_rRNA.fa && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=4 -Dsamjdk.buffer_size=1048576 -Xmx40G -jar $RNASEQC_JAR \
  -n 1000 \
  -o metrics/rnaseqRep \
  -r /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  -s alignment/rnaseqc.samples.txt \
  -t /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.transcript_id.gtf \
  -ttype 2\
  -BWArRNA dummy_rRNA.fa && \
zip -r metrics/rnaseqRep.zip metrics/rnaseqRep
rnaseqc.0bbf9099a1d23b1dad3b3e916996e204.mugqic.done
)
rnaseqc_1_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"rnaseqc\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_1.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"rnaseqc\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_1.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=72:00:0 --mem=48G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$rnaseqc_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: rnaseqc_2_JOB_ID: rnaseqc_report
#-------------------------------------------------------------------------------
JOB_NAME=rnaseqc_report
JOB_DEPENDENCIES=$rnaseqc_1_JOB_ID
JOB_DONE=job_output/rnaseqc/rnaseqc_report.67821e702049231d38161c5db01ef648.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'rnaseqc_report.67821e702049231d38161c5db01ef648.mugqic.done'
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
  /home/efournie/genpipes/bfx/report/RnaSeq.rnaseqc.md \
  --template /home/efournie/genpipes/bfx/report/RnaSeq.rnaseqc.md \
  --variable trim_alignment_table="$trim_alignment_table_md" \
  --to markdown \
  > report/RnaSeq.rnaseqc.md
rnaseqc_report.67821e702049231d38161c5db01ef648.mugqic.done
)
rnaseqc_2_JOB_ID=$(echo "#! /bin/bash 
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
echo "$rnaseqc_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# STEP: wiggle
#-------------------------------------------------------------------------------
STEP=wiggle
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: wiggle_1_JOB_ID: wiggle.DEX_shNIPBL_2.forward_strandspec
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.DEX_shNIPBL_2.forward_strandspec
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.DEX_shNIPBL_2.forward_strandspec.241386d6d04eb0173f1f6f330ffe3b28.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.DEX_shNIPBL_2.forward_strandspec.241386d6d04eb0173f1f6f330ffe3b28.mugqic.done'
module load mugqic/samtools/1.4.1 java/1.8.0_121 mugqic/picard/2.9.0 && \
samtools view -bh -F 256 -f 81 \
  alignment/DEX_shNIPBL_2/DEX_shNIPBL_2.sorted.mdup.bam \
  > alignment/DEX_shNIPBL_2/DEX_shNIPBL_2.sorted.mdup.tmp1.forward.bam && \
samtools view -bh -F 256 -f 161 \
  alignment/DEX_shNIPBL_2/DEX_shNIPBL_2.sorted.mdup.bam \
  > alignment/DEX_shNIPBL_2/DEX_shNIPBL_2.sorted.mdup.tmp2.forward.bam && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=4 -Dsamjdk.buffer_size=1048576 -Xmx40G -jar $PICARD_HOME/picard.jar MergeSamFiles \
 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/DEX_shNIPBL_2/DEX_shNIPBL_2.sorted.mdup.tmp1.forward.bam \
 INPUT=alignment/DEX_shNIPBL_2/DEX_shNIPBL_2.sorted.mdup.tmp2.forward.bam \
 OUTPUT=alignment/DEX_shNIPBL_2/DEX_shNIPBL_2.sorted.mdup.forward.bam \
 MAX_RECORDS_IN_RAM=5750000 && \
rm alignment/DEX_shNIPBL_2/DEX_shNIPBL_2.sorted.mdup.tmp1.forward.bam alignment/DEX_shNIPBL_2/DEX_shNIPBL_2.sorted.mdup.tmp2.forward.bam
wiggle.DEX_shNIPBL_2.forward_strandspec.241386d6d04eb0173f1f6f330ffe3b28.mugqic.done
)
wiggle_1_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:0 --mem=48G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$wiggle_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: wiggle_2_JOB_ID: wiggle.DEX_shNIPBL_2.reverse_strandspec
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.DEX_shNIPBL_2.reverse_strandspec
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.DEX_shNIPBL_2.reverse_strandspec.0f0970c4f8ac970772915f2491d04f4f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.DEX_shNIPBL_2.reverse_strandspec.0f0970c4f8ac970772915f2491d04f4f.mugqic.done'
module load mugqic/samtools/1.4.1 java/1.8.0_121 mugqic/picard/2.9.0 && \
mkdir -p tracks/DEX_shNIPBL_2 tracks/bigWig && \
samtools view -bh -F 256 -f 97 \
  alignment/DEX_shNIPBL_2/DEX_shNIPBL_2.sorted.mdup.bam \
  > alignment/DEX_shNIPBL_2/DEX_shNIPBL_2.sorted.mdup.tmp1.reverse.bam && \
samtools view -bh -F 256 -f 145 \
  alignment/DEX_shNIPBL_2/DEX_shNIPBL_2.sorted.mdup.bam \
  > alignment/DEX_shNIPBL_2/DEX_shNIPBL_2.sorted.mdup.tmp2.reverse.bam && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=4 -Dsamjdk.buffer_size=1048576 -Xmx40G -jar $PICARD_HOME/picard.jar MergeSamFiles \
 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/DEX_shNIPBL_2/DEX_shNIPBL_2.sorted.mdup.tmp1.reverse.bam \
 INPUT=alignment/DEX_shNIPBL_2/DEX_shNIPBL_2.sorted.mdup.tmp2.reverse.bam \
 OUTPUT=alignment/DEX_shNIPBL_2/DEX_shNIPBL_2.sorted.mdup.reverse.bam \
 MAX_RECORDS_IN_RAM=5750000 && \
rm alignment/DEX_shNIPBL_2/DEX_shNIPBL_2.sorted.mdup.tmp1.reverse.bam alignment/DEX_shNIPBL_2/DEX_shNIPBL_2.sorted.mdup.tmp2.reverse.bam
wiggle.DEX_shNIPBL_2.reverse_strandspec.0f0970c4f8ac970772915f2491d04f4f.mugqic.done
)
wiggle_2_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:0 --mem=48G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$wiggle_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: wiggle_3_JOB_ID: bed_graph.DEX_shNIPBL_2.forward
#-------------------------------------------------------------------------------
JOB_NAME=bed_graph.DEX_shNIPBL_2.forward
JOB_DEPENDENCIES=$wiggle_1_JOB_ID
JOB_DONE=job_output/wiggle/bed_graph.DEX_shNIPBL_2.forward.8ab6190ba3ee94acdf0935abd2be301f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bed_graph.DEX_shNIPBL_2.forward.8ab6190ba3ee94acdf0935abd2be301f.mugqic.done'
module load mugqic/samtools/1.4.1 mugqic/bedtools/2.26.0 && \
mkdir -p tracks/DEX_shNIPBL_2  && \
nmblines=$(samtools view -F 256 -f 81  alignment/DEX_shNIPBL_2/DEX_shNIPBL_2.sorted.mdup.forward.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed  -bg -split -scale $scalefactor \
  -ibam alignment/DEX_shNIPBL_2/DEX_shNIPBL_2.sorted.mdup.forward.bam \
  -g /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  > tracks/DEX_shNIPBL_2/DEX_shNIPBL_2.forward.bedGraph
bed_graph.DEX_shNIPBL_2.forward.8ab6190ba3ee94acdf0935abd2be301f.mugqic.done
)
wiggle_3_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:0 --mem=38G -N 1 -n 8 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$wiggle_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: wiggle_4_JOB_ID: wiggle.DEX_shNIPBL_2.forward
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.DEX_shNIPBL_2.forward
JOB_DEPENDENCIES=$wiggle_3_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.DEX_shNIPBL_2.forward.09f141961f774f550c0f7d116b9d85b2.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.DEX_shNIPBL_2.forward.09f141961f774f550c0f7d116b9d85b2.mugqic.done'
module load mugqic/ucsc/v346 && \
mkdir -p tracks/bigWig && \
cat tracks/DEX_shNIPBL_2/DEX_shNIPBL_2.forward.bedGraph | sort --temporary-directory=${SLURM_TMPDIR} -k1,1 -k2,2n | \
awk '{if($0 !~ /^[A-W]/) print ""$0; else print $0}' | grep -vP "GL|lambda|pUC19" | sed 's/MT/chrM/' | \
awk '{printf "%s\t%d\t%d\t%4.4g\n", $1,$2,$3,$4}' > tracks/DEX_shNIPBL_2/DEX_shNIPBL_2.forward.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/DEX_shNIPBL_2/DEX_shNIPBL_2.forward.bedGraph.sorted \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  tracks/bigWig/DEX_shNIPBL_2.forward.bw
wiggle.DEX_shNIPBL_2.forward.09f141961f774f550c0f7d116b9d85b2.mugqic.done
)
wiggle_4_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:0 --mem=48G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$wiggle_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: wiggle_5_JOB_ID: bed_graph.DEX_shNIPBL_2.reverse
#-------------------------------------------------------------------------------
JOB_NAME=bed_graph.DEX_shNIPBL_2.reverse
JOB_DEPENDENCIES=$wiggle_2_JOB_ID
JOB_DONE=job_output/wiggle/bed_graph.DEX_shNIPBL_2.reverse.8c6a56a220923281cbef780c7f81205a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bed_graph.DEX_shNIPBL_2.reverse.8c6a56a220923281cbef780c7f81205a.mugqic.done'
module load mugqic/samtools/1.4.1 mugqic/bedtools/2.26.0 && \
mkdir -p tracks/DEX_shNIPBL_2  && \
nmblines=$(samtools view -F 256 -f 97  alignment/DEX_shNIPBL_2/DEX_shNIPBL_2.sorted.mdup.reverse.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed  -bg -split -scale $scalefactor \
  -ibam alignment/DEX_shNIPBL_2/DEX_shNIPBL_2.sorted.mdup.reverse.bam \
  -g /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  > tracks/DEX_shNIPBL_2/DEX_shNIPBL_2.reverse.bedGraph
bed_graph.DEX_shNIPBL_2.reverse.8c6a56a220923281cbef780c7f81205a.mugqic.done
)
wiggle_5_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:0 --mem=38G -N 1 -n 8 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$wiggle_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: wiggle_6_JOB_ID: wiggle.DEX_shNIPBL_2.reverse
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.DEX_shNIPBL_2.reverse
JOB_DEPENDENCIES=$wiggle_5_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.DEX_shNIPBL_2.reverse.8dfc81d6075b45eda464659ff4b961e0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.DEX_shNIPBL_2.reverse.8dfc81d6075b45eda464659ff4b961e0.mugqic.done'
module load mugqic/ucsc/v346 && \
mkdir -p tracks/bigWig && \
cat tracks/DEX_shNIPBL_2/DEX_shNIPBL_2.reverse.bedGraph | sort --temporary-directory=${SLURM_TMPDIR} -k1,1 -k2,2n | \
awk '{if($0 !~ /^[A-W]/) print ""$0; else print $0}' | grep -vP "GL|lambda|pUC19" | sed 's/MT/chrM/' | \
awk '{printf "%s\t%d\t%d\t%4.4g\n", $1,$2,$3,$4}' > tracks/DEX_shNIPBL_2/DEX_shNIPBL_2.reverse.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/DEX_shNIPBL_2/DEX_shNIPBL_2.reverse.bedGraph.sorted \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  tracks/bigWig/DEX_shNIPBL_2.reverse.bw
wiggle.DEX_shNIPBL_2.reverse.8dfc81d6075b45eda464659ff4b961e0.mugqic.done
)
wiggle_6_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:0 --mem=48G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$wiggle_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: wiggle_7_JOB_ID: wiggle.DEX_shMED1_2.forward_strandspec
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.DEX_shMED1_2.forward_strandspec
JOB_DEPENDENCIES=$picard_mark_duplicates_2_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.DEX_shMED1_2.forward_strandspec.0bc41f37b8eeeeabf907abe0207bc948.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.DEX_shMED1_2.forward_strandspec.0bc41f37b8eeeeabf907abe0207bc948.mugqic.done'
module load mugqic/samtools/1.4.1 java/1.8.0_121 mugqic/picard/2.9.0 && \
samtools view -bh -F 256 -f 81 \
  alignment/DEX_shMED1_2/DEX_shMED1_2.sorted.mdup.bam \
  > alignment/DEX_shMED1_2/DEX_shMED1_2.sorted.mdup.tmp1.forward.bam && \
samtools view -bh -F 256 -f 161 \
  alignment/DEX_shMED1_2/DEX_shMED1_2.sorted.mdup.bam \
  > alignment/DEX_shMED1_2/DEX_shMED1_2.sorted.mdup.tmp2.forward.bam && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=4 -Dsamjdk.buffer_size=1048576 -Xmx40G -jar $PICARD_HOME/picard.jar MergeSamFiles \
 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/DEX_shMED1_2/DEX_shMED1_2.sorted.mdup.tmp1.forward.bam \
 INPUT=alignment/DEX_shMED1_2/DEX_shMED1_2.sorted.mdup.tmp2.forward.bam \
 OUTPUT=alignment/DEX_shMED1_2/DEX_shMED1_2.sorted.mdup.forward.bam \
 MAX_RECORDS_IN_RAM=5750000 && \
rm alignment/DEX_shMED1_2/DEX_shMED1_2.sorted.mdup.tmp1.forward.bam alignment/DEX_shMED1_2/DEX_shMED1_2.sorted.mdup.tmp2.forward.bam
wiggle.DEX_shMED1_2.forward_strandspec.0bc41f37b8eeeeabf907abe0207bc948.mugqic.done
)
wiggle_7_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:0 --mem=48G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$wiggle_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: wiggle_8_JOB_ID: wiggle.DEX_shMED1_2.reverse_strandspec
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.DEX_shMED1_2.reverse_strandspec
JOB_DEPENDENCIES=$picard_mark_duplicates_2_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.DEX_shMED1_2.reverse_strandspec.057d56fbe91e0477226c199de4a505d4.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.DEX_shMED1_2.reverse_strandspec.057d56fbe91e0477226c199de4a505d4.mugqic.done'
module load mugqic/samtools/1.4.1 java/1.8.0_121 mugqic/picard/2.9.0 && \
mkdir -p tracks/DEX_shMED1_2 tracks/bigWig && \
samtools view -bh -F 256 -f 97 \
  alignment/DEX_shMED1_2/DEX_shMED1_2.sorted.mdup.bam \
  > alignment/DEX_shMED1_2/DEX_shMED1_2.sorted.mdup.tmp1.reverse.bam && \
samtools view -bh -F 256 -f 145 \
  alignment/DEX_shMED1_2/DEX_shMED1_2.sorted.mdup.bam \
  > alignment/DEX_shMED1_2/DEX_shMED1_2.sorted.mdup.tmp2.reverse.bam && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=4 -Dsamjdk.buffer_size=1048576 -Xmx40G -jar $PICARD_HOME/picard.jar MergeSamFiles \
 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/DEX_shMED1_2/DEX_shMED1_2.sorted.mdup.tmp1.reverse.bam \
 INPUT=alignment/DEX_shMED1_2/DEX_shMED1_2.sorted.mdup.tmp2.reverse.bam \
 OUTPUT=alignment/DEX_shMED1_2/DEX_shMED1_2.sorted.mdup.reverse.bam \
 MAX_RECORDS_IN_RAM=5750000 && \
rm alignment/DEX_shMED1_2/DEX_shMED1_2.sorted.mdup.tmp1.reverse.bam alignment/DEX_shMED1_2/DEX_shMED1_2.sorted.mdup.tmp2.reverse.bam
wiggle.DEX_shMED1_2.reverse_strandspec.057d56fbe91e0477226c199de4a505d4.mugqic.done
)
wiggle_8_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:0 --mem=48G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$wiggle_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: wiggle_9_JOB_ID: bed_graph.DEX_shMED1_2.forward
#-------------------------------------------------------------------------------
JOB_NAME=bed_graph.DEX_shMED1_2.forward
JOB_DEPENDENCIES=$wiggle_7_JOB_ID
JOB_DONE=job_output/wiggle/bed_graph.DEX_shMED1_2.forward.892ba75c57a875bc2a354f17aa5584fd.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bed_graph.DEX_shMED1_2.forward.892ba75c57a875bc2a354f17aa5584fd.mugqic.done'
module load mugqic/samtools/1.4.1 mugqic/bedtools/2.26.0 && \
mkdir -p tracks/DEX_shMED1_2  && \
nmblines=$(samtools view -F 256 -f 81  alignment/DEX_shMED1_2/DEX_shMED1_2.sorted.mdup.forward.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed  -bg -split -scale $scalefactor \
  -ibam alignment/DEX_shMED1_2/DEX_shMED1_2.sorted.mdup.forward.bam \
  -g /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  > tracks/DEX_shMED1_2/DEX_shMED1_2.forward.bedGraph
bed_graph.DEX_shMED1_2.forward.892ba75c57a875bc2a354f17aa5584fd.mugqic.done
)
wiggle_9_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:0 --mem=38G -N 1 -n 8 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$wiggle_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: wiggle_10_JOB_ID: wiggle.DEX_shMED1_2.forward
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.DEX_shMED1_2.forward
JOB_DEPENDENCIES=$wiggle_9_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.DEX_shMED1_2.forward.da9acc2351e0046e5fa1e49b463bd1ba.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.DEX_shMED1_2.forward.da9acc2351e0046e5fa1e49b463bd1ba.mugqic.done'
module load mugqic/ucsc/v346 && \
mkdir -p tracks/bigWig && \
cat tracks/DEX_shMED1_2/DEX_shMED1_2.forward.bedGraph | sort --temporary-directory=${SLURM_TMPDIR} -k1,1 -k2,2n | \
awk '{if($0 !~ /^[A-W]/) print ""$0; else print $0}' | grep -vP "GL|lambda|pUC19" | sed 's/MT/chrM/' | \
awk '{printf "%s\t%d\t%d\t%4.4g\n", $1,$2,$3,$4}' > tracks/DEX_shMED1_2/DEX_shMED1_2.forward.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/DEX_shMED1_2/DEX_shMED1_2.forward.bedGraph.sorted \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  tracks/bigWig/DEX_shMED1_2.forward.bw
wiggle.DEX_shMED1_2.forward.da9acc2351e0046e5fa1e49b463bd1ba.mugqic.done
)
wiggle_10_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:0 --mem=48G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$wiggle_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: wiggle_11_JOB_ID: bed_graph.DEX_shMED1_2.reverse
#-------------------------------------------------------------------------------
JOB_NAME=bed_graph.DEX_shMED1_2.reverse
JOB_DEPENDENCIES=$wiggle_8_JOB_ID
JOB_DONE=job_output/wiggle/bed_graph.DEX_shMED1_2.reverse.5d206520878e19206ac2f73fb10b91c5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bed_graph.DEX_shMED1_2.reverse.5d206520878e19206ac2f73fb10b91c5.mugqic.done'
module load mugqic/samtools/1.4.1 mugqic/bedtools/2.26.0 && \
mkdir -p tracks/DEX_shMED1_2  && \
nmblines=$(samtools view -F 256 -f 97  alignment/DEX_shMED1_2/DEX_shMED1_2.sorted.mdup.reverse.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed  -bg -split -scale $scalefactor \
  -ibam alignment/DEX_shMED1_2/DEX_shMED1_2.sorted.mdup.reverse.bam \
  -g /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  > tracks/DEX_shMED1_2/DEX_shMED1_2.reverse.bedGraph
bed_graph.DEX_shMED1_2.reverse.5d206520878e19206ac2f73fb10b91c5.mugqic.done
)
wiggle_11_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:0 --mem=38G -N 1 -n 8 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$wiggle_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: wiggle_12_JOB_ID: wiggle.DEX_shMED1_2.reverse
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.DEX_shMED1_2.reverse
JOB_DEPENDENCIES=$wiggle_11_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.DEX_shMED1_2.reverse.1154f3480d2e8184fb896b3ef8c87774.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.DEX_shMED1_2.reverse.1154f3480d2e8184fb896b3ef8c87774.mugqic.done'
module load mugqic/ucsc/v346 && \
mkdir -p tracks/bigWig && \
cat tracks/DEX_shMED1_2/DEX_shMED1_2.reverse.bedGraph | sort --temporary-directory=${SLURM_TMPDIR} -k1,1 -k2,2n | \
awk '{if($0 !~ /^[A-W]/) print ""$0; else print $0}' | grep -vP "GL|lambda|pUC19" | sed 's/MT/chrM/' | \
awk '{printf "%s\t%d\t%d\t%4.4g\n", $1,$2,$3,$4}' > tracks/DEX_shMED1_2/DEX_shMED1_2.reverse.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/DEX_shMED1_2/DEX_shMED1_2.reverse.bedGraph.sorted \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  tracks/bigWig/DEX_shMED1_2.reverse.bw
wiggle.DEX_shMED1_2.reverse.1154f3480d2e8184fb896b3ef8c87774.mugqic.done
)
wiggle_12_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:0 --mem=48G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$wiggle_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: wiggle_13_JOB_ID: wiggle.DEX_shSMC1A_4.forward_strandspec
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.DEX_shSMC1A_4.forward_strandspec
JOB_DEPENDENCIES=$picard_mark_duplicates_3_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.DEX_shSMC1A_4.forward_strandspec.6ca97a8ed4b400fa1f906b61adec9831.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.DEX_shSMC1A_4.forward_strandspec.6ca97a8ed4b400fa1f906b61adec9831.mugqic.done'
module load mugqic/samtools/1.4.1 java/1.8.0_121 mugqic/picard/2.9.0 && \
samtools view -bh -F 256 -f 81 \
  alignment/DEX_shSMC1A_4/DEX_shSMC1A_4.sorted.mdup.bam \
  > alignment/DEX_shSMC1A_4/DEX_shSMC1A_4.sorted.mdup.tmp1.forward.bam && \
samtools view -bh -F 256 -f 161 \
  alignment/DEX_shSMC1A_4/DEX_shSMC1A_4.sorted.mdup.bam \
  > alignment/DEX_shSMC1A_4/DEX_shSMC1A_4.sorted.mdup.tmp2.forward.bam && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=4 -Dsamjdk.buffer_size=1048576 -Xmx40G -jar $PICARD_HOME/picard.jar MergeSamFiles \
 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/DEX_shSMC1A_4/DEX_shSMC1A_4.sorted.mdup.tmp1.forward.bam \
 INPUT=alignment/DEX_shSMC1A_4/DEX_shSMC1A_4.sorted.mdup.tmp2.forward.bam \
 OUTPUT=alignment/DEX_shSMC1A_4/DEX_shSMC1A_4.sorted.mdup.forward.bam \
 MAX_RECORDS_IN_RAM=5750000 && \
rm alignment/DEX_shSMC1A_4/DEX_shSMC1A_4.sorted.mdup.tmp1.forward.bam alignment/DEX_shSMC1A_4/DEX_shSMC1A_4.sorted.mdup.tmp2.forward.bam
wiggle.DEX_shSMC1A_4.forward_strandspec.6ca97a8ed4b400fa1f906b61adec9831.mugqic.done
)
wiggle_13_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:0 --mem=48G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$wiggle_13_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: wiggle_14_JOB_ID: wiggle.DEX_shSMC1A_4.reverse_strandspec
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.DEX_shSMC1A_4.reverse_strandspec
JOB_DEPENDENCIES=$picard_mark_duplicates_3_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.DEX_shSMC1A_4.reverse_strandspec.54881c5425b3d2c1eec80f4a5a965d31.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.DEX_shSMC1A_4.reverse_strandspec.54881c5425b3d2c1eec80f4a5a965d31.mugqic.done'
module load mugqic/samtools/1.4.1 java/1.8.0_121 mugqic/picard/2.9.0 && \
mkdir -p tracks/DEX_shSMC1A_4 tracks/bigWig && \
samtools view -bh -F 256 -f 97 \
  alignment/DEX_shSMC1A_4/DEX_shSMC1A_4.sorted.mdup.bam \
  > alignment/DEX_shSMC1A_4/DEX_shSMC1A_4.sorted.mdup.tmp1.reverse.bam && \
samtools view -bh -F 256 -f 145 \
  alignment/DEX_shSMC1A_4/DEX_shSMC1A_4.sorted.mdup.bam \
  > alignment/DEX_shSMC1A_4/DEX_shSMC1A_4.sorted.mdup.tmp2.reverse.bam && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=4 -Dsamjdk.buffer_size=1048576 -Xmx40G -jar $PICARD_HOME/picard.jar MergeSamFiles \
 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/DEX_shSMC1A_4/DEX_shSMC1A_4.sorted.mdup.tmp1.reverse.bam \
 INPUT=alignment/DEX_shSMC1A_4/DEX_shSMC1A_4.sorted.mdup.tmp2.reverse.bam \
 OUTPUT=alignment/DEX_shSMC1A_4/DEX_shSMC1A_4.sorted.mdup.reverse.bam \
 MAX_RECORDS_IN_RAM=5750000 && \
rm alignment/DEX_shSMC1A_4/DEX_shSMC1A_4.sorted.mdup.tmp1.reverse.bam alignment/DEX_shSMC1A_4/DEX_shSMC1A_4.sorted.mdup.tmp2.reverse.bam
wiggle.DEX_shSMC1A_4.reverse_strandspec.54881c5425b3d2c1eec80f4a5a965d31.mugqic.done
)
wiggle_14_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:0 --mem=48G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$wiggle_14_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: wiggle_15_JOB_ID: bed_graph.DEX_shSMC1A_4.forward
#-------------------------------------------------------------------------------
JOB_NAME=bed_graph.DEX_shSMC1A_4.forward
JOB_DEPENDENCIES=$wiggle_13_JOB_ID
JOB_DONE=job_output/wiggle/bed_graph.DEX_shSMC1A_4.forward.edcb55ba0e36fef336d0b5ba3dcc73e8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bed_graph.DEX_shSMC1A_4.forward.edcb55ba0e36fef336d0b5ba3dcc73e8.mugqic.done'
module load mugqic/samtools/1.4.1 mugqic/bedtools/2.26.0 && \
mkdir -p tracks/DEX_shSMC1A_4  && \
nmblines=$(samtools view -F 256 -f 81  alignment/DEX_shSMC1A_4/DEX_shSMC1A_4.sorted.mdup.forward.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed  -bg -split -scale $scalefactor \
  -ibam alignment/DEX_shSMC1A_4/DEX_shSMC1A_4.sorted.mdup.forward.bam \
  -g /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  > tracks/DEX_shSMC1A_4/DEX_shSMC1A_4.forward.bedGraph
bed_graph.DEX_shSMC1A_4.forward.edcb55ba0e36fef336d0b5ba3dcc73e8.mugqic.done
)
wiggle_15_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:0 --mem=38G -N 1 -n 8 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$wiggle_15_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: wiggle_16_JOB_ID: wiggle.DEX_shSMC1A_4.forward
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.DEX_shSMC1A_4.forward
JOB_DEPENDENCIES=$wiggle_15_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.DEX_shSMC1A_4.forward.5c931bb60afb49891f9ec483ae72b9e2.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.DEX_shSMC1A_4.forward.5c931bb60afb49891f9ec483ae72b9e2.mugqic.done'
module load mugqic/ucsc/v346 && \
mkdir -p tracks/bigWig && \
cat tracks/DEX_shSMC1A_4/DEX_shSMC1A_4.forward.bedGraph | sort --temporary-directory=${SLURM_TMPDIR} -k1,1 -k2,2n | \
awk '{if($0 !~ /^[A-W]/) print ""$0; else print $0}' | grep -vP "GL|lambda|pUC19" | sed 's/MT/chrM/' | \
awk '{printf "%s\t%d\t%d\t%4.4g\n", $1,$2,$3,$4}' > tracks/DEX_shSMC1A_4/DEX_shSMC1A_4.forward.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/DEX_shSMC1A_4/DEX_shSMC1A_4.forward.bedGraph.sorted \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  tracks/bigWig/DEX_shSMC1A_4.forward.bw
wiggle.DEX_shSMC1A_4.forward.5c931bb60afb49891f9ec483ae72b9e2.mugqic.done
)
wiggle_16_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:0 --mem=48G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$wiggle_16_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: wiggle_17_JOB_ID: bed_graph.DEX_shSMC1A_4.reverse
#-------------------------------------------------------------------------------
JOB_NAME=bed_graph.DEX_shSMC1A_4.reverse
JOB_DEPENDENCIES=$wiggle_14_JOB_ID
JOB_DONE=job_output/wiggle/bed_graph.DEX_shSMC1A_4.reverse.f3a49783e4c98a21b35a0e9583dbe11b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bed_graph.DEX_shSMC1A_4.reverse.f3a49783e4c98a21b35a0e9583dbe11b.mugqic.done'
module load mugqic/samtools/1.4.1 mugqic/bedtools/2.26.0 && \
mkdir -p tracks/DEX_shSMC1A_4  && \
nmblines=$(samtools view -F 256 -f 97  alignment/DEX_shSMC1A_4/DEX_shSMC1A_4.sorted.mdup.reverse.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed  -bg -split -scale $scalefactor \
  -ibam alignment/DEX_shSMC1A_4/DEX_shSMC1A_4.sorted.mdup.reverse.bam \
  -g /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  > tracks/DEX_shSMC1A_4/DEX_shSMC1A_4.reverse.bedGraph
bed_graph.DEX_shSMC1A_4.reverse.f3a49783e4c98a21b35a0e9583dbe11b.mugqic.done
)
wiggle_17_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:0 --mem=38G -N 1 -n 8 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$wiggle_17_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: wiggle_18_JOB_ID: wiggle.DEX_shSMC1A_4.reverse
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.DEX_shSMC1A_4.reverse
JOB_DEPENDENCIES=$wiggle_17_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.DEX_shSMC1A_4.reverse.953993519f447d0297d210f96d4c3646.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.DEX_shSMC1A_4.reverse.953993519f447d0297d210f96d4c3646.mugqic.done'
module load mugqic/ucsc/v346 && \
mkdir -p tracks/bigWig && \
cat tracks/DEX_shSMC1A_4/DEX_shSMC1A_4.reverse.bedGraph | sort --temporary-directory=${SLURM_TMPDIR} -k1,1 -k2,2n | \
awk '{if($0 !~ /^[A-W]/) print ""$0; else print $0}' | grep -vP "GL|lambda|pUC19" | sed 's/MT/chrM/' | \
awk '{printf "%s\t%d\t%d\t%4.4g\n", $1,$2,$3,$4}' > tracks/DEX_shSMC1A_4/DEX_shSMC1A_4.reverse.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/DEX_shSMC1A_4/DEX_shSMC1A_4.reverse.bedGraph.sorted \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  tracks/bigWig/DEX_shSMC1A_4.reverse.bw
wiggle.DEX_shSMC1A_4.reverse.953993519f447d0297d210f96d4c3646.mugqic.done
)
wiggle_18_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:0 --mem=48G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$wiggle_18_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: wiggle_19_JOB_ID: wiggle.DEX_nonMamm_2.forward_strandspec
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.DEX_nonMamm_2.forward_strandspec
JOB_DEPENDENCIES=$picard_mark_duplicates_4_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.DEX_nonMamm_2.forward_strandspec.98b881d13727208bf4de46a65fdfac3b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.DEX_nonMamm_2.forward_strandspec.98b881d13727208bf4de46a65fdfac3b.mugqic.done'
module load mugqic/samtools/1.4.1 java/1.8.0_121 mugqic/picard/2.9.0 && \
samtools view -bh -F 256 -f 81 \
  alignment/DEX_nonMamm_2/DEX_nonMamm_2.sorted.mdup.bam \
  > alignment/DEX_nonMamm_2/DEX_nonMamm_2.sorted.mdup.tmp1.forward.bam && \
samtools view -bh -F 256 -f 161 \
  alignment/DEX_nonMamm_2/DEX_nonMamm_2.sorted.mdup.bam \
  > alignment/DEX_nonMamm_2/DEX_nonMamm_2.sorted.mdup.tmp2.forward.bam && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=4 -Dsamjdk.buffer_size=1048576 -Xmx40G -jar $PICARD_HOME/picard.jar MergeSamFiles \
 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/DEX_nonMamm_2/DEX_nonMamm_2.sorted.mdup.tmp1.forward.bam \
 INPUT=alignment/DEX_nonMamm_2/DEX_nonMamm_2.sorted.mdup.tmp2.forward.bam \
 OUTPUT=alignment/DEX_nonMamm_2/DEX_nonMamm_2.sorted.mdup.forward.bam \
 MAX_RECORDS_IN_RAM=5750000 && \
rm alignment/DEX_nonMamm_2/DEX_nonMamm_2.sorted.mdup.tmp1.forward.bam alignment/DEX_nonMamm_2/DEX_nonMamm_2.sorted.mdup.tmp2.forward.bam
wiggle.DEX_nonMamm_2.forward_strandspec.98b881d13727208bf4de46a65fdfac3b.mugqic.done
)
wiggle_19_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:0 --mem=48G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$wiggle_19_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: wiggle_20_JOB_ID: wiggle.DEX_nonMamm_2.reverse_strandspec
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.DEX_nonMamm_2.reverse_strandspec
JOB_DEPENDENCIES=$picard_mark_duplicates_4_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.DEX_nonMamm_2.reverse_strandspec.d1c530369a106c976d2b1771ad0a2f15.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.DEX_nonMamm_2.reverse_strandspec.d1c530369a106c976d2b1771ad0a2f15.mugqic.done'
module load mugqic/samtools/1.4.1 java/1.8.0_121 mugqic/picard/2.9.0 && \
mkdir -p tracks/DEX_nonMamm_2 tracks/bigWig && \
samtools view -bh -F 256 -f 97 \
  alignment/DEX_nonMamm_2/DEX_nonMamm_2.sorted.mdup.bam \
  > alignment/DEX_nonMamm_2/DEX_nonMamm_2.sorted.mdup.tmp1.reverse.bam && \
samtools view -bh -F 256 -f 145 \
  alignment/DEX_nonMamm_2/DEX_nonMamm_2.sorted.mdup.bam \
  > alignment/DEX_nonMamm_2/DEX_nonMamm_2.sorted.mdup.tmp2.reverse.bam && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=4 -Dsamjdk.buffer_size=1048576 -Xmx40G -jar $PICARD_HOME/picard.jar MergeSamFiles \
 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/DEX_nonMamm_2/DEX_nonMamm_2.sorted.mdup.tmp1.reverse.bam \
 INPUT=alignment/DEX_nonMamm_2/DEX_nonMamm_2.sorted.mdup.tmp2.reverse.bam \
 OUTPUT=alignment/DEX_nonMamm_2/DEX_nonMamm_2.sorted.mdup.reverse.bam \
 MAX_RECORDS_IN_RAM=5750000 && \
rm alignment/DEX_nonMamm_2/DEX_nonMamm_2.sorted.mdup.tmp1.reverse.bam alignment/DEX_nonMamm_2/DEX_nonMamm_2.sorted.mdup.tmp2.reverse.bam
wiggle.DEX_nonMamm_2.reverse_strandspec.d1c530369a106c976d2b1771ad0a2f15.mugqic.done
)
wiggle_20_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:0 --mem=48G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$wiggle_20_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: wiggle_21_JOB_ID: bed_graph.DEX_nonMamm_2.forward
#-------------------------------------------------------------------------------
JOB_NAME=bed_graph.DEX_nonMamm_2.forward
JOB_DEPENDENCIES=$wiggle_19_JOB_ID
JOB_DONE=job_output/wiggle/bed_graph.DEX_nonMamm_2.forward.33ef2ef579359c5dfdb65082348d95e9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bed_graph.DEX_nonMamm_2.forward.33ef2ef579359c5dfdb65082348d95e9.mugqic.done'
module load mugqic/samtools/1.4.1 mugqic/bedtools/2.26.0 && \
mkdir -p tracks/DEX_nonMamm_2  && \
nmblines=$(samtools view -F 256 -f 81  alignment/DEX_nonMamm_2/DEX_nonMamm_2.sorted.mdup.forward.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed  -bg -split -scale $scalefactor \
  -ibam alignment/DEX_nonMamm_2/DEX_nonMamm_2.sorted.mdup.forward.bam \
  -g /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  > tracks/DEX_nonMamm_2/DEX_nonMamm_2.forward.bedGraph
bed_graph.DEX_nonMamm_2.forward.33ef2ef579359c5dfdb65082348d95e9.mugqic.done
)
wiggle_21_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:0 --mem=38G -N 1 -n 8 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$wiggle_21_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: wiggle_22_JOB_ID: wiggle.DEX_nonMamm_2.forward
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.DEX_nonMamm_2.forward
JOB_DEPENDENCIES=$wiggle_21_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.DEX_nonMamm_2.forward.efe0d56b29a7044ac7351568265d9d9a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.DEX_nonMamm_2.forward.efe0d56b29a7044ac7351568265d9d9a.mugqic.done'
module load mugqic/ucsc/v346 && \
mkdir -p tracks/bigWig && \
cat tracks/DEX_nonMamm_2/DEX_nonMamm_2.forward.bedGraph | sort --temporary-directory=${SLURM_TMPDIR} -k1,1 -k2,2n | \
awk '{if($0 !~ /^[A-W]/) print ""$0; else print $0}' | grep -vP "GL|lambda|pUC19" | sed 's/MT/chrM/' | \
awk '{printf "%s\t%d\t%d\t%4.4g\n", $1,$2,$3,$4}' > tracks/DEX_nonMamm_2/DEX_nonMamm_2.forward.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/DEX_nonMamm_2/DEX_nonMamm_2.forward.bedGraph.sorted \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  tracks/bigWig/DEX_nonMamm_2.forward.bw
wiggle.DEX_nonMamm_2.forward.efe0d56b29a7044ac7351568265d9d9a.mugqic.done
)
wiggle_22_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:0 --mem=48G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$wiggle_22_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: wiggle_23_JOB_ID: bed_graph.DEX_nonMamm_2.reverse
#-------------------------------------------------------------------------------
JOB_NAME=bed_graph.DEX_nonMamm_2.reverse
JOB_DEPENDENCIES=$wiggle_20_JOB_ID
JOB_DONE=job_output/wiggle/bed_graph.DEX_nonMamm_2.reverse.f77a3903d1bd50fe45646949dca2e245.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bed_graph.DEX_nonMamm_2.reverse.f77a3903d1bd50fe45646949dca2e245.mugqic.done'
module load mugqic/samtools/1.4.1 mugqic/bedtools/2.26.0 && \
mkdir -p tracks/DEX_nonMamm_2  && \
nmblines=$(samtools view -F 256 -f 97  alignment/DEX_nonMamm_2/DEX_nonMamm_2.sorted.mdup.reverse.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed  -bg -split -scale $scalefactor \
  -ibam alignment/DEX_nonMamm_2/DEX_nonMamm_2.sorted.mdup.reverse.bam \
  -g /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  > tracks/DEX_nonMamm_2/DEX_nonMamm_2.reverse.bedGraph
bed_graph.DEX_nonMamm_2.reverse.f77a3903d1bd50fe45646949dca2e245.mugqic.done
)
wiggle_23_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:0 --mem=38G -N 1 -n 8 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$wiggle_23_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: wiggle_24_JOB_ID: wiggle.DEX_nonMamm_2.reverse
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.DEX_nonMamm_2.reverse
JOB_DEPENDENCIES=$wiggle_23_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.DEX_nonMamm_2.reverse.268c4efc9715b0e4f934ecdeec5c2b6c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.DEX_nonMamm_2.reverse.268c4efc9715b0e4f934ecdeec5c2b6c.mugqic.done'
module load mugqic/ucsc/v346 && \
mkdir -p tracks/bigWig && \
cat tracks/DEX_nonMamm_2/DEX_nonMamm_2.reverse.bedGraph | sort --temporary-directory=${SLURM_TMPDIR} -k1,1 -k2,2n | \
awk '{if($0 !~ /^[A-W]/) print ""$0; else print $0}' | grep -vP "GL|lambda|pUC19" | sed 's/MT/chrM/' | \
awk '{printf "%s\t%d\t%d\t%4.4g\n", $1,$2,$3,$4}' > tracks/DEX_nonMamm_2/DEX_nonMamm_2.reverse.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/DEX_nonMamm_2/DEX_nonMamm_2.reverse.bedGraph.sorted \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  tracks/bigWig/DEX_nonMamm_2.reverse.bw
wiggle.DEX_nonMamm_2.reverse.268c4efc9715b0e4f934ecdeec5c2b6c.mugqic.done
)
wiggle_24_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:0 --mem=48G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$wiggle_24_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: wiggle_25_JOB_ID: wiggle.DEX_nonMamm_1.forward_strandspec
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.DEX_nonMamm_1.forward_strandspec
JOB_DEPENDENCIES=$picard_mark_duplicates_5_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.DEX_nonMamm_1.forward_strandspec.8cfd0b948012016043beaf0632807b36.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.DEX_nonMamm_1.forward_strandspec.8cfd0b948012016043beaf0632807b36.mugqic.done'
module load mugqic/samtools/1.4.1 java/1.8.0_121 mugqic/picard/2.9.0 && \
samtools view -bh -F 256 -f 81 \
  alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.bam \
  > alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.tmp1.forward.bam && \
samtools view -bh -F 256 -f 161 \
  alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.bam \
  > alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.tmp2.forward.bam && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=4 -Dsamjdk.buffer_size=1048576 -Xmx40G -jar $PICARD_HOME/picard.jar MergeSamFiles \
 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.tmp1.forward.bam \
 INPUT=alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.tmp2.forward.bam \
 OUTPUT=alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.forward.bam \
 MAX_RECORDS_IN_RAM=5750000 && \
rm alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.tmp1.forward.bam alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.tmp2.forward.bam
wiggle.DEX_nonMamm_1.forward_strandspec.8cfd0b948012016043beaf0632807b36.mugqic.done
)
wiggle_25_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_1.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_1.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:0 --mem=48G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$wiggle_25_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: wiggle_26_JOB_ID: wiggle.DEX_nonMamm_1.reverse_strandspec
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.DEX_nonMamm_1.reverse_strandspec
JOB_DEPENDENCIES=$picard_mark_duplicates_5_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.DEX_nonMamm_1.reverse_strandspec.1d78cadbc746ca633592e019fba0a049.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.DEX_nonMamm_1.reverse_strandspec.1d78cadbc746ca633592e019fba0a049.mugqic.done'
module load mugqic/samtools/1.4.1 java/1.8.0_121 mugqic/picard/2.9.0 && \
mkdir -p tracks/DEX_nonMamm_1 tracks/bigWig && \
samtools view -bh -F 256 -f 97 \
  alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.bam \
  > alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.tmp1.reverse.bam && \
samtools view -bh -F 256 -f 145 \
  alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.bam \
  > alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.tmp2.reverse.bam && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=4 -Dsamjdk.buffer_size=1048576 -Xmx40G -jar $PICARD_HOME/picard.jar MergeSamFiles \
 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.tmp1.reverse.bam \
 INPUT=alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.tmp2.reverse.bam \
 OUTPUT=alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.reverse.bam \
 MAX_RECORDS_IN_RAM=5750000 && \
rm alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.tmp1.reverse.bam alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.tmp2.reverse.bam
wiggle.DEX_nonMamm_1.reverse_strandspec.1d78cadbc746ca633592e019fba0a049.mugqic.done
)
wiggle_26_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_1.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_1.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:0 --mem=48G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$wiggle_26_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: wiggle_27_JOB_ID: bed_graph.DEX_nonMamm_1.forward
#-------------------------------------------------------------------------------
JOB_NAME=bed_graph.DEX_nonMamm_1.forward
JOB_DEPENDENCIES=$wiggle_25_JOB_ID
JOB_DONE=job_output/wiggle/bed_graph.DEX_nonMamm_1.forward.c6a6a4047b36b4ce0db6250101c37a3f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bed_graph.DEX_nonMamm_1.forward.c6a6a4047b36b4ce0db6250101c37a3f.mugqic.done'
module load mugqic/samtools/1.4.1 mugqic/bedtools/2.26.0 && \
mkdir -p tracks/DEX_nonMamm_1  && \
nmblines=$(samtools view -F 256 -f 81  alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.forward.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed  -bg -split -scale $scalefactor \
  -ibam alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.forward.bam \
  -g /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  > tracks/DEX_nonMamm_1/DEX_nonMamm_1.forward.bedGraph
bed_graph.DEX_nonMamm_1.forward.c6a6a4047b36b4ce0db6250101c37a3f.mugqic.done
)
wiggle_27_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_1.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_1.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:0 --mem=38G -N 1 -n 8 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$wiggle_27_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: wiggle_28_JOB_ID: wiggle.DEX_nonMamm_1.forward
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.DEX_nonMamm_1.forward
JOB_DEPENDENCIES=$wiggle_27_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.DEX_nonMamm_1.forward.5a7082d0c5fd18bcec121c94ddd0e7c9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.DEX_nonMamm_1.forward.5a7082d0c5fd18bcec121c94ddd0e7c9.mugqic.done'
module load mugqic/ucsc/v346 && \
mkdir -p tracks/bigWig && \
cat tracks/DEX_nonMamm_1/DEX_nonMamm_1.forward.bedGraph | sort --temporary-directory=${SLURM_TMPDIR} -k1,1 -k2,2n | \
awk '{if($0 !~ /^[A-W]/) print ""$0; else print $0}' | grep -vP "GL|lambda|pUC19" | sed 's/MT/chrM/' | \
awk '{printf "%s\t%d\t%d\t%4.4g\n", $1,$2,$3,$4}' > tracks/DEX_nonMamm_1/DEX_nonMamm_1.forward.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/DEX_nonMamm_1/DEX_nonMamm_1.forward.bedGraph.sorted \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  tracks/bigWig/DEX_nonMamm_1.forward.bw
wiggle.DEX_nonMamm_1.forward.5a7082d0c5fd18bcec121c94ddd0e7c9.mugqic.done
)
wiggle_28_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_1.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_1.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:0 --mem=48G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$wiggle_28_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: wiggle_29_JOB_ID: bed_graph.DEX_nonMamm_1.reverse
#-------------------------------------------------------------------------------
JOB_NAME=bed_graph.DEX_nonMamm_1.reverse
JOB_DEPENDENCIES=$wiggle_26_JOB_ID
JOB_DONE=job_output/wiggle/bed_graph.DEX_nonMamm_1.reverse.b0047d7d102086a23327a739b59962bb.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bed_graph.DEX_nonMamm_1.reverse.b0047d7d102086a23327a739b59962bb.mugqic.done'
module load mugqic/samtools/1.4.1 mugqic/bedtools/2.26.0 && \
mkdir -p tracks/DEX_nonMamm_1  && \
nmblines=$(samtools view -F 256 -f 97  alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.reverse.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed  -bg -split -scale $scalefactor \
  -ibam alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.reverse.bam \
  -g /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  > tracks/DEX_nonMamm_1/DEX_nonMamm_1.reverse.bedGraph
bed_graph.DEX_nonMamm_1.reverse.b0047d7d102086a23327a739b59962bb.mugqic.done
)
wiggle_29_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_1.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_1.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:0 --mem=38G -N 1 -n 8 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$wiggle_29_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: wiggle_30_JOB_ID: wiggle.DEX_nonMamm_1.reverse
#-------------------------------------------------------------------------------
JOB_NAME=wiggle.DEX_nonMamm_1.reverse
JOB_DEPENDENCIES=$wiggle_29_JOB_ID
JOB_DONE=job_output/wiggle/wiggle.DEX_nonMamm_1.reverse.b070a87361e1f2cd9aaf130366cb9b2a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'wiggle.DEX_nonMamm_1.reverse.b070a87361e1f2cd9aaf130366cb9b2a.mugqic.done'
module load mugqic/ucsc/v346 && \
mkdir -p tracks/bigWig && \
cat tracks/DEX_nonMamm_1/DEX_nonMamm_1.reverse.bedGraph | sort --temporary-directory=${SLURM_TMPDIR} -k1,1 -k2,2n | \
awk '{if($0 !~ /^[A-W]/) print ""$0; else print $0}' | grep -vP "GL|lambda|pUC19" | sed 's/MT/chrM/' | \
awk '{printf "%s\t%d\t%d\t%4.4g\n", $1,$2,$3,$4}' > tracks/DEX_nonMamm_1/DEX_nonMamm_1.reverse.bedGraph.sorted && \
bedGraphToBigWig \
  tracks/DEX_nonMamm_1/DEX_nonMamm_1.reverse.bedGraph.sorted \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  tracks/bigWig/DEX_nonMamm_1.reverse.bw
wiggle.DEX_nonMamm_1.reverse.b070a87361e1f2cd9aaf130366cb9b2a.mugqic.done
)
wiggle_30_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_1.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"wiggle\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_1.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=12:00:0 --mem=48G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$wiggle_30_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# STEP: raw_counts
#-------------------------------------------------------------------------------
STEP=raw_counts
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: raw_counts_1_JOB_ID: htseq_count.DEX_shNIPBL_2
#-------------------------------------------------------------------------------
JOB_NAME=htseq_count.DEX_shNIPBL_2
JOB_DEPENDENCIES=$picard_sort_sam_1_JOB_ID
JOB_DONE=job_output/raw_counts/htseq_count.DEX_shNIPBL_2.3cd2a88c92e9c82a29d259bf3934226e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'htseq_count.DEX_shNIPBL_2.3cd2a88c92e9c82a29d259bf3934226e.mugqic.done'
module load mugqic/samtools/1.4.1 mugqic/python/2.7.13 && \
mkdir -p raw_counts && \
samtools view -F 4 \
  alignment/DEX_shNIPBL_2/DEX_shNIPBL_2.QueryNameSorted.bam | \
htseq-count - \
  -m intersection-nonempty \
  --stranded=reverse \
  --format=sam \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  > raw_counts/DEX_shNIPBL_2.readcounts.csv
htseq_count.DEX_shNIPBL_2.3cd2a88c92e9c82a29d259bf3934226e.mugqic.done
)
raw_counts_1_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"raw_counts\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"raw_counts\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=24G -N 1 -n 6 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$raw_counts_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: raw_counts_2_JOB_ID: htseq_count.DEX_shMED1_2
#-------------------------------------------------------------------------------
JOB_NAME=htseq_count.DEX_shMED1_2
JOB_DEPENDENCIES=$picard_sort_sam_2_JOB_ID
JOB_DONE=job_output/raw_counts/htseq_count.DEX_shMED1_2.67ec5966d57d4f0850e851a325c53d8b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'htseq_count.DEX_shMED1_2.67ec5966d57d4f0850e851a325c53d8b.mugqic.done'
module load mugqic/samtools/1.4.1 mugqic/python/2.7.13 && \
mkdir -p raw_counts && \
samtools view -F 4 \
  alignment/DEX_shMED1_2/DEX_shMED1_2.QueryNameSorted.bam | \
htseq-count - \
  -m intersection-nonempty \
  --stranded=reverse \
  --format=sam \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  > raw_counts/DEX_shMED1_2.readcounts.csv
htseq_count.DEX_shMED1_2.67ec5966d57d4f0850e851a325c53d8b.mugqic.done
)
raw_counts_2_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"raw_counts\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"raw_counts\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=24G -N 1 -n 6 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$raw_counts_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: raw_counts_3_JOB_ID: htseq_count.DEX_shSMC1A_4
#-------------------------------------------------------------------------------
JOB_NAME=htseq_count.DEX_shSMC1A_4
JOB_DEPENDENCIES=$picard_sort_sam_3_JOB_ID
JOB_DONE=job_output/raw_counts/htseq_count.DEX_shSMC1A_4.3bceb650f6342218da3a19d0ae3624a0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'htseq_count.DEX_shSMC1A_4.3bceb650f6342218da3a19d0ae3624a0.mugqic.done'
module load mugqic/samtools/1.4.1 mugqic/python/2.7.13 && \
mkdir -p raw_counts && \
samtools view -F 4 \
  alignment/DEX_shSMC1A_4/DEX_shSMC1A_4.QueryNameSorted.bam | \
htseq-count - \
  -m intersection-nonempty \
  --stranded=reverse \
  --format=sam \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  > raw_counts/DEX_shSMC1A_4.readcounts.csv
htseq_count.DEX_shSMC1A_4.3bceb650f6342218da3a19d0ae3624a0.mugqic.done
)
raw_counts_3_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"raw_counts\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"raw_counts\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=24G -N 1 -n 6 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$raw_counts_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: raw_counts_4_JOB_ID: htseq_count.DEX_nonMamm_2
#-------------------------------------------------------------------------------
JOB_NAME=htseq_count.DEX_nonMamm_2
JOB_DEPENDENCIES=$picard_sort_sam_4_JOB_ID
JOB_DONE=job_output/raw_counts/htseq_count.DEX_nonMamm_2.2730d73214d52e4216f47d857b03f58a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'htseq_count.DEX_nonMamm_2.2730d73214d52e4216f47d857b03f58a.mugqic.done'
module load mugqic/samtools/1.4.1 mugqic/python/2.7.13 && \
mkdir -p raw_counts && \
samtools view -F 4 \
  alignment/DEX_nonMamm_2/DEX_nonMamm_2.QueryNameSorted.bam | \
htseq-count - \
  -m intersection-nonempty \
  --stranded=reverse \
  --format=sam \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  > raw_counts/DEX_nonMamm_2.readcounts.csv
htseq_count.DEX_nonMamm_2.2730d73214d52e4216f47d857b03f58a.mugqic.done
)
raw_counts_4_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"raw_counts\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"raw_counts\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=24G -N 1 -n 6 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$raw_counts_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: raw_counts_5_JOB_ID: htseq_count.DEX_nonMamm_1
#-------------------------------------------------------------------------------
JOB_NAME=htseq_count.DEX_nonMamm_1
JOB_DEPENDENCIES=$picard_sort_sam_5_JOB_ID
JOB_DONE=job_output/raw_counts/htseq_count.DEX_nonMamm_1.ecab7dbe77dd2318e5f5d269412f714e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'htseq_count.DEX_nonMamm_1.ecab7dbe77dd2318e5f5d269412f714e.mugqic.done'
module load mugqic/samtools/1.4.1 mugqic/python/2.7.13 && \
mkdir -p raw_counts && \
samtools view -F 4 \
  alignment/DEX_nonMamm_1/DEX_nonMamm_1.QueryNameSorted.bam | \
htseq-count - \
  -m intersection-nonempty \
  --stranded=reverse \
  --format=sam \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  > raw_counts/DEX_nonMamm_1.readcounts.csv
htseq_count.DEX_nonMamm_1.ecab7dbe77dd2318e5f5d269412f714e.mugqic.done
)
raw_counts_5_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"raw_counts\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_1.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"raw_counts\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_1.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=24G -N 1 -n 6 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$raw_counts_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# STEP: raw_counts_metrics
#-------------------------------------------------------------------------------
STEP=raw_counts_metrics
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: raw_counts_metrics_1_JOB_ID: metrics.matrix
#-------------------------------------------------------------------------------
JOB_NAME=metrics.matrix
JOB_DEPENDENCIES=$raw_counts_1_JOB_ID:$raw_counts_2_JOB_ID:$raw_counts_3_JOB_ID:$raw_counts_4_JOB_ID:$raw_counts_5_JOB_ID
JOB_DONE=job_output/raw_counts_metrics/metrics.matrix.a365de06e8e5b5383f5b46b9428f661e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'metrics.matrix.a365de06e8e5b5383f5b46b9428f661e.mugqic.done'
module load mugqic/mugqic_tools/2.1.9 && \
mkdir -p DGE && \
gtf2tmpMatrix.awk \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  DGE/tmpMatrix.txt && \
HEAD='Gene\tSymbol' && \
for read_count_file in \
  raw_counts/DEX_shNIPBL_2.readcounts.csv \
  raw_counts/DEX_shMED1_2.readcounts.csv \
  raw_counts/DEX_shSMC1A_4.readcounts.csv \
  raw_counts/DEX_nonMamm_2.readcounts.csv \
  raw_counts/DEX_nonMamm_1.readcounts.csv
do
  sort -k1,1 $read_count_file > DGE/tmpSort.txt && \
  join -1 1 -2 1 <(sort -k1,1 DGE/tmpMatrix.txt) DGE/tmpSort.txt > DGE/tmpMatrix.2.txt && \
  mv DGE/tmpMatrix.2.txt DGE/tmpMatrix.txt && \
  na=$(basename $read_count_file | rev | cut -d. -f3- | rev) && \
  HEAD="$HEAD\t$na"
done && \
echo -e $HEAD | cat - DGE/tmpMatrix.txt | tr ' ' '\t' > DGE/rawCountMatrix.csv && \
rm DGE/tmpSort.txt DGE/tmpMatrix.txt
metrics.matrix.a365de06e8e5b5383f5b46b9428f661e.mugqic.done
)
raw_counts_metrics_1_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"raw_counts_metrics\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_1.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"raw_counts_metrics\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_1.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=5:00:0 --mem=4G -N 1 -n 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$raw_counts_metrics_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: raw_counts_metrics_2_JOB_ID: metrics.wigzip
#-------------------------------------------------------------------------------
JOB_NAME=metrics.wigzip
JOB_DEPENDENCIES=$wiggle_4_JOB_ID:$wiggle_6_JOB_ID:$wiggle_10_JOB_ID:$wiggle_12_JOB_ID:$wiggle_16_JOB_ID:$wiggle_18_JOB_ID:$wiggle_22_JOB_ID:$wiggle_24_JOB_ID:$wiggle_28_JOB_ID:$wiggle_30_JOB_ID
JOB_DONE=job_output/raw_counts_metrics/metrics.wigzip.edc4e268c60b94d072db60163e113f9c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'metrics.wigzip.edc4e268c60b94d072db60163e113f9c.mugqic.done'
zip -r tracks.zip tracks/bigWig
metrics.wigzip.edc4e268c60b94d072db60163e113f9c.mugqic.done
)
raw_counts_metrics_2_JOB_ID=$(echo "#! /bin/bash 
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=5:00:0 --mem=4G -N 1 -n 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$raw_counts_metrics_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: raw_counts_metrics_3_JOB_ID: rpkm_saturation
#-------------------------------------------------------------------------------
JOB_NAME=rpkm_saturation
JOB_DEPENDENCIES=$raw_counts_metrics_1_JOB_ID
JOB_DONE=job_output/raw_counts_metrics/rpkm_saturation.9e0e736f78a92ef5802a96526d9f3c9e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'rpkm_saturation.9e0e736f78a92ef5802a96526d9f3c9e.mugqic.done'
module load mugqic/R_Bioconductor/3.4.2_3.6 mugqic/mugqic_tools/2.1.9 && \
mkdir -p metrics/saturation && \
Rscript $R_TOOLS/rpkmSaturation.R \
  DGE/rawCountMatrix.csv \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.genes.length.tsv \
  raw_counts \
  metrics/saturation \
  15 \
  1 && \
zip -r metrics/saturation.zip metrics/saturation
rpkm_saturation.9e0e736f78a92ef5802a96526d9f3c9e.mugqic.done
)
raw_counts_metrics_3_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"raw_counts_metrics\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_1.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"raw_counts_metrics\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_1.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=128G -N 1 -n 16 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$raw_counts_metrics_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: raw_counts_metrics_4_JOB_ID: raw_count_metrics_report
#-------------------------------------------------------------------------------
JOB_NAME=raw_count_metrics_report
JOB_DEPENDENCIES=$rnaseqc_1_JOB_ID:$raw_counts_metrics_2_JOB_ID:$raw_counts_metrics_3_JOB_ID
JOB_DONE=job_output/raw_counts_metrics/raw_count_metrics_report.e4faf9a0ac8b5cf62e4257c28ff6a716.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'raw_count_metrics_report.e4faf9a0ac8b5cf62e4257c28ff6a716.mugqic.done'
module load mugqic/pandoc/1.15.2 && \
mkdir -p report && \
cp metrics/rnaseqRep/corrMatrixSpearman.txt report/corrMatrixSpearman.tsv && \
cp tracks.zip report/ && \
cp metrics/saturation.zip report/ && \
pandoc --to=markdown \
  --template /home/efournie/genpipes/bfx/report/RnaSeq.raw_counts_metrics.md \
  --variable corr_matrix_spearman_table="`head -16 report/corrMatrixSpearman.tsv | cut -f-16| awk -F"	" '{OFS="	"; if (NR==1) {$0="Vs"$0; print; gsub(/[^	]/, "-"); print} else {printf $1; for (i=2; i<=NF; i++) {printf "	"sprintf("%.2f", $i)}; print ""}}' | sed 's/	/|/g'`" \
  /home/efournie/genpipes/bfx/report/RnaSeq.raw_counts_metrics.md \
  > report/RnaSeq.raw_counts_metrics.md
raw_count_metrics_report.e4faf9a0ac8b5cf62e4257c28ff6a716.mugqic.done
)
raw_counts_metrics_4_JOB_ID=$(echo "#! /bin/bash 
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
echo "$raw_counts_metrics_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# STEP: cufflinks
#-------------------------------------------------------------------------------
STEP=cufflinks
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: cufflinks_1_JOB_ID: cufflinks.DEX_shNIPBL_2
#-------------------------------------------------------------------------------
JOB_NAME=cufflinks.DEX_shNIPBL_2
JOB_DEPENDENCIES=$bam_hard_clip_1_JOB_ID
JOB_DONE=job_output/cufflinks/cufflinks.DEX_shNIPBL_2.068f320db7b200af4ef9e9313449b665.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cufflinks.DEX_shNIPBL_2.068f320db7b200af4ef9e9313449b665.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/DEX_shNIPBL_2 && \
cufflinks -q  \
  --GTF-guide /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/DEX_shNIPBL_2 \
  --num-threads 8 \
  alignment/DEX_shNIPBL_2/DEX_shNIPBL_2.sorted.mdup.hardClip.bam
cufflinks.DEX_shNIPBL_2.068f320db7b200af4ef9e9313449b665.mugqic.done
)
cufflinks_1_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"cufflinks\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"cufflinks\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=48:00:0 --mem=32G -N 1 -n 8 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$cufflinks_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: cufflinks_2_JOB_ID: cufflinks.DEX_shMED1_2
#-------------------------------------------------------------------------------
JOB_NAME=cufflinks.DEX_shMED1_2
JOB_DEPENDENCIES=$bam_hard_clip_2_JOB_ID
JOB_DONE=job_output/cufflinks/cufflinks.DEX_shMED1_2.4d6692562a42618fbeee2a2c551e81a1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cufflinks.DEX_shMED1_2.4d6692562a42618fbeee2a2c551e81a1.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/DEX_shMED1_2 && \
cufflinks -q  \
  --GTF-guide /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/DEX_shMED1_2 \
  --num-threads 8 \
  alignment/DEX_shMED1_2/DEX_shMED1_2.sorted.mdup.hardClip.bam
cufflinks.DEX_shMED1_2.4d6692562a42618fbeee2a2c551e81a1.mugqic.done
)
cufflinks_2_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"cufflinks\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"cufflinks\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=48:00:0 --mem=32G -N 1 -n 8 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$cufflinks_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: cufflinks_3_JOB_ID: cufflinks.DEX_shSMC1A_4
#-------------------------------------------------------------------------------
JOB_NAME=cufflinks.DEX_shSMC1A_4
JOB_DEPENDENCIES=$bam_hard_clip_3_JOB_ID
JOB_DONE=job_output/cufflinks/cufflinks.DEX_shSMC1A_4.980a5d7d33add214dde7417ea77eeb58.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cufflinks.DEX_shSMC1A_4.980a5d7d33add214dde7417ea77eeb58.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/DEX_shSMC1A_4 && \
cufflinks -q  \
  --GTF-guide /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/DEX_shSMC1A_4 \
  --num-threads 8 \
  alignment/DEX_shSMC1A_4/DEX_shSMC1A_4.sorted.mdup.hardClip.bam
cufflinks.DEX_shSMC1A_4.980a5d7d33add214dde7417ea77eeb58.mugqic.done
)
cufflinks_3_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"cufflinks\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"cufflinks\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=48:00:0 --mem=32G -N 1 -n 8 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$cufflinks_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: cufflinks_4_JOB_ID: cufflinks.DEX_nonMamm_2
#-------------------------------------------------------------------------------
JOB_NAME=cufflinks.DEX_nonMamm_2
JOB_DEPENDENCIES=$bam_hard_clip_4_JOB_ID
JOB_DONE=job_output/cufflinks/cufflinks.DEX_nonMamm_2.cc0910d53cf7696ae7de8cd2a9f86c2f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cufflinks.DEX_nonMamm_2.cc0910d53cf7696ae7de8cd2a9f86c2f.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/DEX_nonMamm_2 && \
cufflinks -q  \
  --GTF-guide /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/DEX_nonMamm_2 \
  --num-threads 8 \
  alignment/DEX_nonMamm_2/DEX_nonMamm_2.sorted.mdup.hardClip.bam
cufflinks.DEX_nonMamm_2.cc0910d53cf7696ae7de8cd2a9f86c2f.mugqic.done
)
cufflinks_4_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"cufflinks\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"cufflinks\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=48:00:0 --mem=32G -N 1 -n 8 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$cufflinks_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: cufflinks_5_JOB_ID: cufflinks.DEX_nonMamm_1
#-------------------------------------------------------------------------------
JOB_NAME=cufflinks.DEX_nonMamm_1
JOB_DEPENDENCIES=$bam_hard_clip_5_JOB_ID
JOB_DONE=job_output/cufflinks/cufflinks.DEX_nonMamm_1.8b18c4e8c05a63fa380d9a42c057eb25.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cufflinks.DEX_nonMamm_1.8b18c4e8c05a63fa380d9a42c057eb25.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/DEX_nonMamm_1 && \
cufflinks -q  \
  --GTF-guide /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/DEX_nonMamm_1 \
  --num-threads 8 \
  alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.hardClip.bam
cufflinks.DEX_nonMamm_1.8b18c4e8c05a63fa380d9a42c057eb25.mugqic.done
)
cufflinks_5_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"cufflinks\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_1.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"cufflinks\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_1.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=48:00:0 --mem=32G -N 1 -n 8 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$cufflinks_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# STEP: cuffmerge
#-------------------------------------------------------------------------------
STEP=cuffmerge
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: cuffmerge_1_JOB_ID: cuffmerge
#-------------------------------------------------------------------------------
JOB_NAME=cuffmerge
JOB_DEPENDENCIES=$cufflinks_1_JOB_ID:$cufflinks_2_JOB_ID:$cufflinks_3_JOB_ID:$cufflinks_4_JOB_ID:$cufflinks_5_JOB_ID
JOB_DONE=job_output/cuffmerge/cuffmerge.31395245413d7dacbc3c72892eee296a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffmerge.31395245413d7dacbc3c72892eee296a.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/AllSamples && \
`cat > cufflinks/cuffmerge.samples.txt << END
cufflinks/DEX_shNIPBL_2/transcripts.gtf
cufflinks/DEX_shMED1_2/transcripts.gtf
cufflinks/DEX_shSMC1A_4/transcripts.gtf
cufflinks/DEX_nonMamm_2/transcripts.gtf
cufflinks/DEX_nonMamm_1/transcripts.gtf
END

` && \
cuffmerge  \
  --ref-gtf /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl87.gtf \
  --ref-sequence /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  -o cufflinks/AllSamples \
  --num-threads 8 \
  cufflinks/cuffmerge.samples.txt
cuffmerge.31395245413d7dacbc3c72892eee296a.mugqic.done
)
cuffmerge_1_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"cuffmerge\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_1.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"cuffmerge\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_1.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=48:00:0 --mem=32G -N 1 -n 8 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$cuffmerge_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# STEP: cuffquant
#-------------------------------------------------------------------------------
STEP=cuffquant
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: cuffquant_1_JOB_ID: cuffquant.DEX_shNIPBL_2
#-------------------------------------------------------------------------------
JOB_NAME=cuffquant.DEX_shNIPBL_2
JOB_DEPENDENCIES=$bam_hard_clip_1_JOB_ID:$cuffmerge_1_JOB_ID
JOB_DONE=job_output/cuffquant/cuffquant.DEX_shNIPBL_2.ba1ee68408781a03600259b45076df16.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffquant.DEX_shNIPBL_2.ba1ee68408781a03600259b45076df16.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/DEX_shNIPBL_2 && \
cuffquant -q  \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/DEX_shNIPBL_2 \
  --num-threads 8 \
  cufflinks/AllSamples/merged.gtf \
  alignment/DEX_shNIPBL_2/DEX_shNIPBL_2.sorted.mdup.hardClip.bam
cuffquant.DEX_shNIPBL_2.ba1ee68408781a03600259b45076df16.mugqic.done
)
cuffquant_1_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"cuffquant\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"cuffquant\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=48:00:0 --mem=32G -N 1 -n 8 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$cuffquant_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: cuffquant_2_JOB_ID: cuffquant.DEX_shMED1_2
#-------------------------------------------------------------------------------
JOB_NAME=cuffquant.DEX_shMED1_2
JOB_DEPENDENCIES=$bam_hard_clip_2_JOB_ID:$cuffmerge_1_JOB_ID
JOB_DONE=job_output/cuffquant/cuffquant.DEX_shMED1_2.a7c18d73011f7c57f0d63758f6cc61a2.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffquant.DEX_shMED1_2.a7c18d73011f7c57f0d63758f6cc61a2.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/DEX_shMED1_2 && \
cuffquant -q  \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/DEX_shMED1_2 \
  --num-threads 8 \
  cufflinks/AllSamples/merged.gtf \
  alignment/DEX_shMED1_2/DEX_shMED1_2.sorted.mdup.hardClip.bam
cuffquant.DEX_shMED1_2.a7c18d73011f7c57f0d63758f6cc61a2.mugqic.done
)
cuffquant_2_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"cuffquant\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"cuffquant\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=48:00:0 --mem=32G -N 1 -n 8 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$cuffquant_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: cuffquant_3_JOB_ID: cuffquant.DEX_shSMC1A_4
#-------------------------------------------------------------------------------
JOB_NAME=cuffquant.DEX_shSMC1A_4
JOB_DEPENDENCIES=$bam_hard_clip_3_JOB_ID:$cuffmerge_1_JOB_ID
JOB_DONE=job_output/cuffquant/cuffquant.DEX_shSMC1A_4.c7fc3b5243fcf8166b649d634a29d7d8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffquant.DEX_shSMC1A_4.c7fc3b5243fcf8166b649d634a29d7d8.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/DEX_shSMC1A_4 && \
cuffquant -q  \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/DEX_shSMC1A_4 \
  --num-threads 8 \
  cufflinks/AllSamples/merged.gtf \
  alignment/DEX_shSMC1A_4/DEX_shSMC1A_4.sorted.mdup.hardClip.bam
cuffquant.DEX_shSMC1A_4.c7fc3b5243fcf8166b649d634a29d7d8.mugqic.done
)
cuffquant_3_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"cuffquant\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"cuffquant\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=48:00:0 --mem=32G -N 1 -n 8 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$cuffquant_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: cuffquant_4_JOB_ID: cuffquant.DEX_nonMamm_2
#-------------------------------------------------------------------------------
JOB_NAME=cuffquant.DEX_nonMamm_2
JOB_DEPENDENCIES=$bam_hard_clip_4_JOB_ID:$cuffmerge_1_JOB_ID
JOB_DONE=job_output/cuffquant/cuffquant.DEX_nonMamm_2.0c8652e691c2965ea58297abf72a087a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffquant.DEX_nonMamm_2.0c8652e691c2965ea58297abf72a087a.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/DEX_nonMamm_2 && \
cuffquant -q  \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/DEX_nonMamm_2 \
  --num-threads 8 \
  cufflinks/AllSamples/merged.gtf \
  alignment/DEX_nonMamm_2/DEX_nonMamm_2.sorted.mdup.hardClip.bam
cuffquant.DEX_nonMamm_2.0c8652e691c2965ea58297abf72a087a.mugqic.done
)
cuffquant_4_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"cuffquant\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"cuffquant\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=48:00:0 --mem=32G -N 1 -n 8 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$cuffquant_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: cuffquant_5_JOB_ID: cuffquant.DEX_nonMamm_1
#-------------------------------------------------------------------------------
JOB_NAME=cuffquant.DEX_nonMamm_1
JOB_DEPENDENCIES=$bam_hard_clip_5_JOB_ID:$cuffmerge_1_JOB_ID
JOB_DONE=job_output/cuffquant/cuffquant.DEX_nonMamm_1.c461f59178793703a318963b05d922ae.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffquant.DEX_nonMamm_1.c461f59178793703a318963b05d922ae.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/DEX_nonMamm_1 && \
cuffquant -q  \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/DEX_nonMamm_1 \
  --num-threads 8 \
  cufflinks/AllSamples/merged.gtf \
  alignment/DEX_nonMamm_1/DEX_nonMamm_1.sorted.mdup.hardClip.bam
cuffquant.DEX_nonMamm_1.c461f59178793703a318963b05d922ae.mugqic.done
)
cuffquant_5_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"cuffquant\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_1.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"cuffquant\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_1.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=48:00:0 --mem=32G -N 1 -n 8 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$cuffquant_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# STEP: cuffdiff
#-------------------------------------------------------------------------------
STEP=cuffdiff
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: cuffdiff_1_JOB_ID: cuffdiff.Dex_vs_EtOH-All
#-------------------------------------------------------------------------------
JOB_NAME=cuffdiff.Dex_vs_EtOH-All
JOB_DEPENDENCIES=$cuffmerge_1_JOB_ID:$cuffquant_1_JOB_ID:$cuffquant_2_JOB_ID:$cuffquant_3_JOB_ID:$cuffquant_4_JOB_ID
JOB_DONE=job_output/cuffdiff/cuffdiff.Dex_vs_EtOH-All.4e4dd850ebb20c21a861a52445748bf1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffdiff.Dex_vs_EtOH-All.4e4dd850ebb20c21a861a52445748bf1.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cuffdiff/Dex_vs_EtOH-All && \
cuffdiff -u \
  --frag-bias-correct /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --library-type fr-firststrand \
  --output-dir cuffdiff/Dex_vs_EtOH-All \
  --num-threads 8 \
  cufflinks/AllSamples/merged.gtf \
   \
  cufflinks/DEX_shNIPBL_2/abundances.cxb,cufflinks/DEX_shMED1_2/abundances.cxb,cufflinks/DEX_shSMC1A_4/abundances.cxb,cufflinks/DEX_nonMamm_2/abundances.cxb
cuffdiff.Dex_vs_EtOH-All.4e4dd850ebb20c21a861a52445748bf1.mugqic.done
)
cuffdiff_1_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"cuffdiff\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"cuffdiff\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=48:00:0 --mem=32G -N 1 -n 8 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$cuffdiff_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: cuffdiff_2_JOB_ID: cuffdiff.Dex_vs_EtOH-nonMamm
#-------------------------------------------------------------------------------
JOB_NAME=cuffdiff.Dex_vs_EtOH-nonMamm
JOB_DEPENDENCIES=$cuffmerge_1_JOB_ID:$cuffquant_4_JOB_ID
JOB_DONE=job_output/cuffdiff/cuffdiff.Dex_vs_EtOH-nonMamm.5260562f5e3edf6bf7122a86ff68f944.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffdiff.Dex_vs_EtOH-nonMamm.5260562f5e3edf6bf7122a86ff68f944.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cuffdiff/Dex_vs_EtOH-nonMamm && \
cuffdiff -u \
  --frag-bias-correct /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --library-type fr-firststrand \
  --output-dir cuffdiff/Dex_vs_EtOH-nonMamm \
  --num-threads 8 \
  cufflinks/AllSamples/merged.gtf \
   \
  cufflinks/DEX_nonMamm_2/abundances.cxb
cuffdiff.Dex_vs_EtOH-nonMamm.5260562f5e3edf6bf7122a86ff68f944.mugqic.done
)
cuffdiff_2_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"cuffdiff\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"cuffdiff\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=48:00:0 --mem=32G -N 1 -n 8 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$cuffdiff_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: cuffdiff_3_JOB_ID: cuffdiff.Dex_vs_EtOH-shMED1
#-------------------------------------------------------------------------------
JOB_NAME=cuffdiff.Dex_vs_EtOH-shMED1
JOB_DEPENDENCIES=$cuffmerge_1_JOB_ID:$cuffquant_2_JOB_ID
JOB_DONE=job_output/cuffdiff/cuffdiff.Dex_vs_EtOH-shMED1.2c2f203ef31fb70c9fb9155e530304cd.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffdiff.Dex_vs_EtOH-shMED1.2c2f203ef31fb70c9fb9155e530304cd.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cuffdiff/Dex_vs_EtOH-shMED1 && \
cuffdiff -u \
  --frag-bias-correct /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --library-type fr-firststrand \
  --output-dir cuffdiff/Dex_vs_EtOH-shMED1 \
  --num-threads 8 \
  cufflinks/AllSamples/merged.gtf \
   \
  cufflinks/DEX_shMED1_2/abundances.cxb
cuffdiff.Dex_vs_EtOH-shMED1.2c2f203ef31fb70c9fb9155e530304cd.mugqic.done
)
cuffdiff_3_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"cuffdiff\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"cuffdiff\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=48:00:0 --mem=32G -N 1 -n 8 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$cuffdiff_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: cuffdiff_4_JOB_ID: cuffdiff.Dex_vs_EtOH-shNIPBL
#-------------------------------------------------------------------------------
JOB_NAME=cuffdiff.Dex_vs_EtOH-shNIPBL
JOB_DEPENDENCIES=$cuffmerge_1_JOB_ID:$cuffquant_1_JOB_ID
JOB_DONE=job_output/cuffdiff/cuffdiff.Dex_vs_EtOH-shNIPBL.b38c769c01e6b08ce3bedb2a28aa7ee2.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffdiff.Dex_vs_EtOH-shNIPBL.b38c769c01e6b08ce3bedb2a28aa7ee2.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cuffdiff/Dex_vs_EtOH-shNIPBL && \
cuffdiff -u \
  --frag-bias-correct /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --library-type fr-firststrand \
  --output-dir cuffdiff/Dex_vs_EtOH-shNIPBL \
  --num-threads 8 \
  cufflinks/AllSamples/merged.gtf \
   \
  cufflinks/DEX_shNIPBL_2/abundances.cxb
cuffdiff.Dex_vs_EtOH-shNIPBL.b38c769c01e6b08ce3bedb2a28aa7ee2.mugqic.done
)
cuffdiff_4_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"cuffdiff\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"cuffdiff\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=48:00:0 --mem=32G -N 1 -n 8 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$cuffdiff_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: cuffdiff_5_JOB_ID: cuffdiff.Dex_vs_EtOH-shSMC1A
#-------------------------------------------------------------------------------
JOB_NAME=cuffdiff.Dex_vs_EtOH-shSMC1A
JOB_DEPENDENCIES=$cuffmerge_1_JOB_ID:$cuffquant_3_JOB_ID
JOB_DONE=job_output/cuffdiff/cuffdiff.Dex_vs_EtOH-shSMC1A.3d2210f968b6cf88757ef5e269e66091.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffdiff.Dex_vs_EtOH-shSMC1A.3d2210f968b6cf88757ef5e269e66091.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cuffdiff/Dex_vs_EtOH-shSMC1A && \
cuffdiff -u \
  --frag-bias-correct /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --library-type fr-firststrand \
  --output-dir cuffdiff/Dex_vs_EtOH-shSMC1A \
  --num-threads 8 \
  cufflinks/AllSamples/merged.gtf \
   \
  cufflinks/DEX_shSMC1A_4/abundances.cxb
cuffdiff.Dex_vs_EtOH-shSMC1A.3d2210f968b6cf88757ef5e269e66091.mugqic.done
)
cuffdiff_5_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"cuffdiff\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"cuffdiff\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=48:00:0 --mem=32G -N 1 -n 8 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$cuffdiff_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: cuffdiff_6_JOB_ID: cuffdiff.Dex_vs_EtOH-shSMC1A-3
#-------------------------------------------------------------------------------
JOB_NAME=cuffdiff.Dex_vs_EtOH-shSMC1A-3
JOB_DEPENDENCIES=$cuffmerge_1_JOB_ID
JOB_DONE=job_output/cuffdiff/cuffdiff.Dex_vs_EtOH-shSMC1A-3.e83461cf172618b0af23c211cb7e3327.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffdiff.Dex_vs_EtOH-shSMC1A-3.e83461cf172618b0af23c211cb7e3327.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cuffdiff/Dex_vs_EtOH-shSMC1A-3 && \
cuffdiff -u \
  --frag-bias-correct /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --library-type fr-firststrand \
  --output-dir cuffdiff/Dex_vs_EtOH-shSMC1A-3 \
  --num-threads 8 \
  cufflinks/AllSamples/merged.gtf \
   \
  
cuffdiff.Dex_vs_EtOH-shSMC1A-3.e83461cf172618b0af23c211cb7e3327.mugqic.done
)
cuffdiff_6_JOB_ID=$(echo "#! /bin/bash 
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=48:00:0 --mem=32G -N 1 -n 8 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$cuffdiff_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: cuffdiff_7_JOB_ID: cuffdiff.Dex_vs_EtOH-shSMC1A-2
#-------------------------------------------------------------------------------
JOB_NAME=cuffdiff.Dex_vs_EtOH-shSMC1A-2
JOB_DEPENDENCIES=$cuffmerge_1_JOB_ID:$cuffquant_3_JOB_ID
JOB_DONE=job_output/cuffdiff/cuffdiff.Dex_vs_EtOH-shSMC1A-2.f8980ecc28b09837fc666839de507eda.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffdiff.Dex_vs_EtOH-shSMC1A-2.f8980ecc28b09837fc666839de507eda.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cuffdiff/Dex_vs_EtOH-shSMC1A-2 && \
cuffdiff -u \
  --frag-bias-correct /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --library-type fr-firststrand \
  --output-dir cuffdiff/Dex_vs_EtOH-shSMC1A-2 \
  --num-threads 8 \
  cufflinks/AllSamples/merged.gtf \
   \
  cufflinks/DEX_shSMC1A_4/abundances.cxb
cuffdiff.Dex_vs_EtOH-shSMC1A-2.f8980ecc28b09837fc666839de507eda.mugqic.done
)
cuffdiff_7_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"cuffdiff\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"cuffdiff\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=48:00:0 --mem=32G -N 1 -n 8 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$cuffdiff_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: cuffdiff_8_JOB_ID: cuffdiff.shMED1_vs_shMamm-Dex
#-------------------------------------------------------------------------------
JOB_NAME=cuffdiff.shMED1_vs_shMamm-Dex
JOB_DEPENDENCIES=$cuffmerge_1_JOB_ID:$cuffquant_2_JOB_ID:$cuffquant_4_JOB_ID
JOB_DONE=job_output/cuffdiff/cuffdiff.shMED1_vs_shMamm-Dex.9cb5783edbe16f9eb03feadf4c50373f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffdiff.shMED1_vs_shMamm-Dex.9cb5783edbe16f9eb03feadf4c50373f.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cuffdiff/shMED1_vs_shMamm-Dex && \
cuffdiff -u \
  --frag-bias-correct /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --library-type fr-firststrand \
  --output-dir cuffdiff/shMED1_vs_shMamm-Dex \
  --num-threads 8 \
  cufflinks/AllSamples/merged.gtf \
  cufflinks/DEX_nonMamm_2/abundances.cxb \
  cufflinks/DEX_shMED1_2/abundances.cxb
cuffdiff.shMED1_vs_shMamm-Dex.9cb5783edbe16f9eb03feadf4c50373f.mugqic.done
)
cuffdiff_8_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"cuffdiff\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"cuffdiff\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=48:00:0 --mem=32G -N 1 -n 8 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$cuffdiff_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: cuffdiff_9_JOB_ID: cuffdiff.shMED1_vs_shMamm-EtOH
#-------------------------------------------------------------------------------
JOB_NAME=cuffdiff.shMED1_vs_shMamm-EtOH
JOB_DEPENDENCIES=$cuffmerge_1_JOB_ID
JOB_DONE=job_output/cuffdiff/cuffdiff.shMED1_vs_shMamm-EtOH.0d4acb8d5c8f5ae82cc59d3ccc73fa48.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffdiff.shMED1_vs_shMamm-EtOH.0d4acb8d5c8f5ae82cc59d3ccc73fa48.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cuffdiff/shMED1_vs_shMamm-EtOH && \
cuffdiff -u \
  --frag-bias-correct /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --library-type fr-firststrand \
  --output-dir cuffdiff/shMED1_vs_shMamm-EtOH \
  --num-threads 8 \
  cufflinks/AllSamples/merged.gtf \
   \
  
cuffdiff.shMED1_vs_shMamm-EtOH.0d4acb8d5c8f5ae82cc59d3ccc73fa48.mugqic.done
)
cuffdiff_9_JOB_ID=$(echo "#! /bin/bash 
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=48:00:0 --mem=32G -N 1 -n 8 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$cuffdiff_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: cuffdiff_10_JOB_ID: cuffdiff.shNIPBL_vs_shMamm-Dex
#-------------------------------------------------------------------------------
JOB_NAME=cuffdiff.shNIPBL_vs_shMamm-Dex
JOB_DEPENDENCIES=$cuffmerge_1_JOB_ID:$cuffquant_1_JOB_ID:$cuffquant_4_JOB_ID
JOB_DONE=job_output/cuffdiff/cuffdiff.shNIPBL_vs_shMamm-Dex.c95fe085c152cbb883f796b9205d0380.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffdiff.shNIPBL_vs_shMamm-Dex.c95fe085c152cbb883f796b9205d0380.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cuffdiff/shNIPBL_vs_shMamm-Dex && \
cuffdiff -u \
  --frag-bias-correct /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --library-type fr-firststrand \
  --output-dir cuffdiff/shNIPBL_vs_shMamm-Dex \
  --num-threads 8 \
  cufflinks/AllSamples/merged.gtf \
  cufflinks/DEX_nonMamm_2/abundances.cxb \
  cufflinks/DEX_shNIPBL_2/abundances.cxb
cuffdiff.shNIPBL_vs_shMamm-Dex.c95fe085c152cbb883f796b9205d0380.mugqic.done
)
cuffdiff_10_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"cuffdiff\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"cuffdiff\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=48:00:0 --mem=32G -N 1 -n 8 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$cuffdiff_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: cuffdiff_11_JOB_ID: cuffdiff.shNIPBL_vs_shMamm-EtOH
#-------------------------------------------------------------------------------
JOB_NAME=cuffdiff.shNIPBL_vs_shMamm-EtOH
JOB_DEPENDENCIES=$cuffmerge_1_JOB_ID
JOB_DONE=job_output/cuffdiff/cuffdiff.shNIPBL_vs_shMamm-EtOH.76ebadaf4058642a0d794b722dae1dae.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffdiff.shNIPBL_vs_shMamm-EtOH.76ebadaf4058642a0d794b722dae1dae.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cuffdiff/shNIPBL_vs_shMamm-EtOH && \
cuffdiff -u \
  --frag-bias-correct /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --library-type fr-firststrand \
  --output-dir cuffdiff/shNIPBL_vs_shMamm-EtOH \
  --num-threads 8 \
  cufflinks/AllSamples/merged.gtf \
   \
  
cuffdiff.shNIPBL_vs_shMamm-EtOH.76ebadaf4058642a0d794b722dae1dae.mugqic.done
)
cuffdiff_11_JOB_ID=$(echo "#! /bin/bash 
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=48:00:0 --mem=32G -N 1 -n 8 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$cuffdiff_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: cuffdiff_12_JOB_ID: cuffdiff.shSMC1A-3_vs_shMamm-Dex
#-------------------------------------------------------------------------------
JOB_NAME=cuffdiff.shSMC1A-3_vs_shMamm-Dex
JOB_DEPENDENCIES=$cuffmerge_1_JOB_ID:$cuffquant_4_JOB_ID
JOB_DONE=job_output/cuffdiff/cuffdiff.shSMC1A-3_vs_shMamm-Dex.4ecc998f1f4b143da651205963499f47.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffdiff.shSMC1A-3_vs_shMamm-Dex.4ecc998f1f4b143da651205963499f47.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cuffdiff/shSMC1A-3_vs_shMamm-Dex && \
cuffdiff -u \
  --frag-bias-correct /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --library-type fr-firststrand \
  --output-dir cuffdiff/shSMC1A-3_vs_shMamm-Dex \
  --num-threads 8 \
  cufflinks/AllSamples/merged.gtf \
  cufflinks/DEX_nonMamm_2/abundances.cxb \
  
cuffdiff.shSMC1A-3_vs_shMamm-Dex.4ecc998f1f4b143da651205963499f47.mugqic.done
)
cuffdiff_12_JOB_ID=$(echo "#! /bin/bash 
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=48:00:0 --mem=32G -N 1 -n 8 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$cuffdiff_12_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: cuffdiff_13_JOB_ID: cuffdiff.shSMC1A-2_vs_shMamm-Dex
#-------------------------------------------------------------------------------
JOB_NAME=cuffdiff.shSMC1A-2_vs_shMamm-Dex
JOB_DEPENDENCIES=$cuffmerge_1_JOB_ID:$cuffquant_3_JOB_ID:$cuffquant_4_JOB_ID
JOB_DONE=job_output/cuffdiff/cuffdiff.shSMC1A-2_vs_shMamm-Dex.0fa16c6cdbf31b0c409e003140addbfb.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffdiff.shSMC1A-2_vs_shMamm-Dex.0fa16c6cdbf31b0c409e003140addbfb.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cuffdiff/shSMC1A-2_vs_shMamm-Dex && \
cuffdiff -u \
  --frag-bias-correct /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --library-type fr-firststrand \
  --output-dir cuffdiff/shSMC1A-2_vs_shMamm-Dex \
  --num-threads 8 \
  cufflinks/AllSamples/merged.gtf \
  cufflinks/DEX_nonMamm_2/abundances.cxb \
  cufflinks/DEX_shSMC1A_4/abundances.cxb
cuffdiff.shSMC1A-2_vs_shMamm-Dex.0fa16c6cdbf31b0c409e003140addbfb.mugqic.done
)
cuffdiff_13_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"cuffdiff\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"cuffdiff\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=48:00:0 --mem=32G -N 1 -n 8 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$cuffdiff_13_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: cuffdiff_14_JOB_ID: cuffdiff.shSMC1A_vs_shMamm-EtOH
#-------------------------------------------------------------------------------
JOB_NAME=cuffdiff.shSMC1A_vs_shMamm-EtOH
JOB_DEPENDENCIES=$cuffmerge_1_JOB_ID
JOB_DONE=job_output/cuffdiff/cuffdiff.shSMC1A_vs_shMamm-EtOH.460aae1332e11d8b33127195e1129b77.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffdiff.shSMC1A_vs_shMamm-EtOH.460aae1332e11d8b33127195e1129b77.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cuffdiff/shSMC1A_vs_shMamm-EtOH && \
cuffdiff -u \
  --frag-bias-correct /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --library-type fr-firststrand \
  --output-dir cuffdiff/shSMC1A_vs_shMamm-EtOH \
  --num-threads 8 \
  cufflinks/AllSamples/merged.gtf \
   \
  
cuffdiff.shSMC1A_vs_shMamm-EtOH.460aae1332e11d8b33127195e1129b77.mugqic.done
)
cuffdiff_14_JOB_ID=$(echo "#! /bin/bash 
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=48:00:0 --mem=32G -N 1 -n 8 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$cuffdiff_14_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: cuffdiff_15_JOB_ID: cuffdiff.shSMC1A-3_vs_shMamm-EtOH
#-------------------------------------------------------------------------------
JOB_NAME=cuffdiff.shSMC1A-3_vs_shMamm-EtOH
JOB_DEPENDENCIES=$cuffmerge_1_JOB_ID
JOB_DONE=job_output/cuffdiff/cuffdiff.shSMC1A-3_vs_shMamm-EtOH.3fa1c91fc7f8373620c70cd947e68b70.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffdiff.shSMC1A-3_vs_shMamm-EtOH.3fa1c91fc7f8373620c70cd947e68b70.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cuffdiff/shSMC1A-3_vs_shMamm-EtOH && \
cuffdiff -u \
  --frag-bias-correct /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --library-type fr-firststrand \
  --output-dir cuffdiff/shSMC1A-3_vs_shMamm-EtOH \
  --num-threads 8 \
  cufflinks/AllSamples/merged.gtf \
   \
  
cuffdiff.shSMC1A-3_vs_shMamm-EtOH.3fa1c91fc7f8373620c70cd947e68b70.mugqic.done
)
cuffdiff_15_JOB_ID=$(echo "#! /bin/bash 
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=48:00:0 --mem=32G -N 1 -n 8 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$cuffdiff_15_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: cuffdiff_16_JOB_ID: cuffdiff.shSMC1A-2_vs_shMamm-EtOH
#-------------------------------------------------------------------------------
JOB_NAME=cuffdiff.shSMC1A-2_vs_shMamm-EtOH
JOB_DEPENDENCIES=$cuffmerge_1_JOB_ID
JOB_DONE=job_output/cuffdiff/cuffdiff.shSMC1A-2_vs_shMamm-EtOH.fb51fab063a3d667faa856ab020846ee.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffdiff.shSMC1A-2_vs_shMamm-EtOH.fb51fab063a3d667faa856ab020846ee.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cuffdiff/shSMC1A-2_vs_shMamm-EtOH && \
cuffdiff -u \
  --frag-bias-correct /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
  --library-type fr-firststrand \
  --output-dir cuffdiff/shSMC1A-2_vs_shMamm-EtOH \
  --num-threads 8 \
  cufflinks/AllSamples/merged.gtf \
   \
  
cuffdiff.shSMC1A-2_vs_shMamm-EtOH.fb51fab063a3d667faa856ab020846ee.mugqic.done
)
cuffdiff_16_JOB_ID=$(echo "#! /bin/bash 
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=48:00:0 --mem=32G -N 1 -n 8 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$cuffdiff_16_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# STEP: cuffnorm
#-------------------------------------------------------------------------------
STEP=cuffnorm
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: cuffnorm_1_JOB_ID: cuffnorm
#-------------------------------------------------------------------------------
JOB_NAME=cuffnorm
JOB_DEPENDENCIES=$cuffmerge_1_JOB_ID:$cuffquant_1_JOB_ID:$cuffquant_2_JOB_ID:$cuffquant_3_JOB_ID:$cuffquant_4_JOB_ID:$cuffquant_5_JOB_ID
JOB_DONE=job_output/cuffnorm/cuffnorm.68196c88b0ddf747d721114f818c586f.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'cuffnorm.68196c88b0ddf747d721114f818c586f.mugqic.done'
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cuffnorm && \
cuffnorm -q  \
  --library-type fr-firststrand \
  --output-dir cuffnorm \
  --num-threads 8 \
  --labels DEX_shNIPBL_2,DEX_shMED1_2,DEX_shSMC1A_4,DEX_nonMamm_2,DEX_nonMamm_1 \
  cufflinks/AllSamples/merged.gtf \
  cufflinks/DEX_shNIPBL_2/abundances.cxb \
  cufflinks/DEX_shMED1_2/abundances.cxb \
  cufflinks/DEX_shSMC1A_4/abundances.cxb \
  cufflinks/DEX_nonMamm_2/abundances.cxb \
  cufflinks/DEX_nonMamm_1/abundances.cxb
cuffnorm.68196c88b0ddf747d721114f818c586f.mugqic.done
)
cuffnorm_1_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"cuffnorm\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_1.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini,/home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"cuffnorm\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shNIPBL_2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shMED1_2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_shSMC1A_4.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/rna-pipeline-GRCh38/json/DEX_nonMamm_1.json\" \
  -f \$MUGQIC_STATE
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
echo '#######################################'
echo 'SLURM FAKE EPILOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
exit \$MUGQIC_STATE" | \
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=48:00:0 --mem=32G -N 1 -n 8 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$cuffnorm_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# Call home with pipeline statistics
#-------------------------------------------------------------------------------
wget "http://mugqic.hpc.mcgill.ca/cgi-bin/pipeline.cgi?hostname=cdr592.int.cedar.computecanada.ca&ip=172.16.138.79&pipeline=RnaSeq&steps=picard_sam_to_fastq,trimmomatic,merge_trimmomatic_stats,star,picard_merge_sam_files,picard_sort_sam,picard_mark_duplicates,picard_rna_metrics,estimate_ribosomal_rna,bam_hard_clip,rnaseqc,wiggle,raw_counts,raw_counts_metrics,cufflinks,cuffmerge,cuffquant,cuffdiff,cuffnorm&samples=5&AnonymizedList=288bbdd32a25fd2689b78e7a576829ec,8befc711a91bdcd48cc37e62a68b2feb,15b56eae2cdf5203323e847cafca2478,cdb870911b7cbb4b8a8795618f3e80e9,74f036ab16e84e94b1e03a92ba25555d" --quiet --output-document=/dev/null

