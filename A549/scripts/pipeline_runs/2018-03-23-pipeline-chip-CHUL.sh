#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

#-------------------------------------------------------------------------------
# ChipSeq SlurmScheduler Job Submission Bash script
# Version: 3.0.1-beta
# Created on: 2018-03-23T19:55:36
# Steps:
#   picard_sam_to_fastq: 0 job... skipping
#   trimmomatic: 8 jobs
#   merge_trimmomatic_stats: 1 job
#   bwa_mem_picard_sort_sam: 9 jobs
#   samtools_view_filter: 9 jobs
#   picard_merge_sam_files: 3 jobs
#   picard_mark_duplicates: 4 jobs
#   metrics: 2 jobs
#   homer_make_tag_directory: 3 jobs
#   qc_metrics: 1 job
#   homer_make_ucsc_file: 7 jobs
#   TOTAL: 47 jobs
#-------------------------------------------------------------------------------

OUTPUT_DIR=/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38
JOB_OUTPUT_DIR=$OUTPUT_DIR/job_output
TIMESTAMP=`date +%FT%H.%M.%S`
JOB_LIST=$JOB_OUTPUT_DIR/ChipSeq_job_list_$TIMESTAMP
export CONFIG_FILES="/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini"
mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR


#-------------------------------------------------------------------------------
# STEP: trimmomatic
#-------------------------------------------------------------------------------
STEP=trimmomatic
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: trimmomatic_1_JOB_ID: trimmomatic.chipseq-KAPA2_GCCAAT_L002_R1_001
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.chipseq-KAPA2_GCCAAT_L002_R1_001
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.chipseq-KAPA2_GCCAAT_L002_R1_001.2e9145f2d99df874a40067d64055865c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.chipseq-KAPA2_GCCAAT_L002_R1_001.2e9145f2d99df874a40067d64055865c.mugqic.done'
module load java/1.8.0_121 mugqic/trimmomatic/0.36 && \
mkdir -p trim/A549_DEX_MED1_rep2 && \
`cat > trim/A549_DEX_MED1_rep2/chipseq-KAPA2_GCCAAT_L002_R1_001.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /project/6001942/Working_Directory/Eric/CofactorHR/A549/raw/chip-seq/chipseq-KAPA2_GCCAAT_L002_R1_001.fastq.gz \
  trim/A549_DEX_MED1_rep2/chipseq-KAPA2_GCCAAT_L002_R1_001.trim.single.fastq.gz \
  ILLUMINACLIP:trim/A549_DEX_MED1_rep2/chipseq-KAPA2_GCCAAT_L002_R1_001.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/A549_DEX_MED1_rep2/chipseq-KAPA2_GCCAAT_L002_R1_001.trim.log
trimmomatic.chipseq-KAPA2_GCCAAT_L002_R1_001.2e9145f2d99df874a40067d64055865c.mugqic.done
)
trimmomatic_1_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"trimmomatic\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_MED1_rep2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"trimmomatic\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_MED1_rep2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=24G -N 1 -n 6 | grep "[0-9]" | cut -d\  -f4)
echo "$trimmomatic_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: trimmomatic_2_JOB_ID: trimmomatic.chipseq-KAPA3_CTTGTA_L002_R1_001
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.chipseq-KAPA3_CTTGTA_L002_R1_001
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.chipseq-KAPA3_CTTGTA_L002_R1_001.06cfaea0a5a5abd709d7ab87c7cde084.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.chipseq-KAPA3_CTTGTA_L002_R1_001.06cfaea0a5a5abd709d7ab87c7cde084.mugqic.done'
module load java/1.8.0_121 mugqic/trimmomatic/0.36 && \
mkdir -p trim/A549_DEX_SMC1A_rep2 && \
`cat > trim/A549_DEX_SMC1A_rep2/chipseq-KAPA3_CTTGTA_L002_R1_001.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /project/6001942/Working_Directory/Eric/CofactorHR/A549/raw/chip-seq/chipseq-KAPA3_CTTGTA_L002_R1_001.fastq.gz \
  trim/A549_DEX_SMC1A_rep2/chipseq-KAPA3_CTTGTA_L002_R1_001.trim.single.fastq.gz \
  ILLUMINACLIP:trim/A549_DEX_SMC1A_rep2/chipseq-KAPA3_CTTGTA_L002_R1_001.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/A549_DEX_SMC1A_rep2/chipseq-KAPA3_CTTGTA_L002_R1_001.trim.log
trimmomatic.chipseq-KAPA3_CTTGTA_L002_R1_001.06cfaea0a5a5abd709d7ab87c7cde084.mugqic.done
)
trimmomatic_2_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"trimmomatic\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_SMC1A_rep2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"trimmomatic\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_SMC1A_rep2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=24G -N 1 -n 6 | grep "[0-9]" | cut -d\  -f4)
echo "$trimmomatic_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: trimmomatic_3_JOB_ID: trimmomatic.chipseq-KAPA4_GTGAAA_L002_R1_001
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.chipseq-KAPA4_GTGAAA_L002_R1_001
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.chipseq-KAPA4_GTGAAA_L002_R1_001.9d5948d3a026c76f25e40ac4eb493213.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.chipseq-KAPA4_GTGAAA_L002_R1_001.9d5948d3a026c76f25e40ac4eb493213.mugqic.done'
module load java/1.8.0_121 mugqic/trimmomatic/0.36 && \
mkdir -p trim/A549_DEX_SMC1A_rep2 && \
`cat > trim/A549_DEX_SMC1A_rep2/chipseq-KAPA4_GTGAAA_L002_R1_001.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /project/6001942/Working_Directory/Eric/CofactorHR/A549/raw/chip-seq/chipseq-KAPA4_GTGAAA_L002_R1_001.fastq.gz \
  trim/A549_DEX_SMC1A_rep2/chipseq-KAPA4_GTGAAA_L002_R1_001.trim.single.fastq.gz \
  ILLUMINACLIP:trim/A549_DEX_SMC1A_rep2/chipseq-KAPA4_GTGAAA_L002_R1_001.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/A549_DEX_SMC1A_rep2/chipseq-KAPA4_GTGAAA_L002_R1_001.trim.log
trimmomatic.chipseq-KAPA4_GTGAAA_L002_R1_001.9d5948d3a026c76f25e40ac4eb493213.mugqic.done
)
trimmomatic_3_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"trimmomatic\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_SMC1A_rep2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"trimmomatic\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_SMC1A_rep2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=24G -N 1 -n 6 | grep "[0-9]" | cut -d\  -f4)
echo "$trimmomatic_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: trimmomatic_4_JOB_ID: trimmomatic.chipseq-NEB-1_ATCACG_L002_R1_001
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.chipseq-NEB-1_ATCACG_L002_R1_001
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.chipseq-NEB-1_ATCACG_L002_R1_001.82c5bbfc9bd6536d296f1e2b97b436a6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.chipseq-NEB-1_ATCACG_L002_R1_001.82c5bbfc9bd6536d296f1e2b97b436a6.mugqic.done'
module load java/1.8.0_121 mugqic/trimmomatic/0.36 && \
mkdir -p trim/A549_DEX_MED1_rep2 && \
`cat > trim/A549_DEX_MED1_rep2/chipseq-NEB-1_ATCACG_L002_R1_001.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /project/6001942/Working_Directory/Eric/CofactorHR/A549/raw/chip-seq/chipseq-NEB-1_ATCACG_L002_R1_001.fastq.gz \
  trim/A549_DEX_MED1_rep2/chipseq-NEB-1_ATCACG_L002_R1_001.trim.single.fastq.gz \
  ILLUMINACLIP:trim/A549_DEX_MED1_rep2/chipseq-NEB-1_ATCACG_L002_R1_001.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/A549_DEX_MED1_rep2/chipseq-NEB-1_ATCACG_L002_R1_001.trim.log
trimmomatic.chipseq-NEB-1_ATCACG_L002_R1_001.82c5bbfc9bd6536d296f1e2b97b436a6.mugqic.done
)
trimmomatic_4_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"trimmomatic\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_MED1_rep2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"trimmomatic\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_MED1_rep2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=24G -N 1 -n 6 | grep "[0-9]" | cut -d\  -f4)
echo "$trimmomatic_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: trimmomatic_5_JOB_ID: trimmomatic.chipseq-NEB-2_TTAGGC_L002_R1_001
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.chipseq-NEB-2_TTAGGC_L002_R1_001
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.chipseq-NEB-2_TTAGGC_L002_R1_001.d20c380fd0a1cca1ff7e2f84c64d91cc.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.chipseq-NEB-2_TTAGGC_L002_R1_001.d20c380fd0a1cca1ff7e2f84c64d91cc.mugqic.done'
module load java/1.8.0_121 mugqic/trimmomatic/0.36 && \
mkdir -p trim/A549_DEX_MED1_rep2 && \
`cat > trim/A549_DEX_MED1_rep2/chipseq-NEB-2_TTAGGC_L002_R1_001.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /project/6001942/Working_Directory/Eric/CofactorHR/A549/raw/chip-seq/chipseq-NEB-2_TTAGGC_L002_R1_001.fastq.gz \
  trim/A549_DEX_MED1_rep2/chipseq-NEB-2_TTAGGC_L002_R1_001.trim.single.fastq.gz \
  ILLUMINACLIP:trim/A549_DEX_MED1_rep2/chipseq-NEB-2_TTAGGC_L002_R1_001.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/A549_DEX_MED1_rep2/chipseq-NEB-2_TTAGGC_L002_R1_001.trim.log
trimmomatic.chipseq-NEB-2_TTAGGC_L002_R1_001.d20c380fd0a1cca1ff7e2f84c64d91cc.mugqic.done
)
trimmomatic_5_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"trimmomatic\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_MED1_rep2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"trimmomatic\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_MED1_rep2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=24G -N 1 -n 6 | grep "[0-9]" | cut -d\  -f4)
echo "$trimmomatic_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: trimmomatic_6_JOB_ID: trimmomatic.chipseq-NEB-3_TGACCA_L002_R1_001
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.chipseq-NEB-3_TGACCA_L002_R1_001
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.chipseq-NEB-3_TGACCA_L002_R1_001.42697005a29c70cf3d5245800bf71310.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.chipseq-NEB-3_TGACCA_L002_R1_001.42697005a29c70cf3d5245800bf71310.mugqic.done'
module load java/1.8.0_121 mugqic/trimmomatic/0.36 && \
mkdir -p trim/A549_DEX_SMC1A_rep2 && \
`cat > trim/A549_DEX_SMC1A_rep2/chipseq-NEB-3_TGACCA_L002_R1_001.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /project/6001942/Working_Directory/Eric/CofactorHR/A549/raw/chip-seq/chipseq-NEB-3_TGACCA_L002_R1_001.fastq.gz \
  trim/A549_DEX_SMC1A_rep2/chipseq-NEB-3_TGACCA_L002_R1_001.trim.single.fastq.gz \
  ILLUMINACLIP:trim/A549_DEX_SMC1A_rep2/chipseq-NEB-3_TGACCA_L002_R1_001.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/A549_DEX_SMC1A_rep2/chipseq-NEB-3_TGACCA_L002_R1_001.trim.log
trimmomatic.chipseq-NEB-3_TGACCA_L002_R1_001.42697005a29c70cf3d5245800bf71310.mugqic.done
)
trimmomatic_6_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"trimmomatic\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_SMC1A_rep2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"trimmomatic\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_SMC1A_rep2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=24G -N 1 -n 6 | grep "[0-9]" | cut -d\  -f4)
echo "$trimmomatic_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: trimmomatic_7_JOB_ID: trimmomatic.chipseq-NEB-4_CAGATC_L002_R1_001
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.chipseq-NEB-4_CAGATC_L002_R1_001
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.chipseq-NEB-4_CAGATC_L002_R1_001.09910caa8075017ba4f9349594b2932c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.chipseq-NEB-4_CAGATC_L002_R1_001.09910caa8075017ba4f9349594b2932c.mugqic.done'
module load java/1.8.0_121 mugqic/trimmomatic/0.36 && \
mkdir -p trim/A549_DEX_SMC1A_rep2 && \
`cat > trim/A549_DEX_SMC1A_rep2/chipseq-NEB-4_CAGATC_L002_R1_001.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /project/6001942/Working_Directory/Eric/CofactorHR/A549/raw/chip-seq/chipseq-NEB-4_CAGATC_L002_R1_001.fastq.gz \
  trim/A549_DEX_SMC1A_rep2/chipseq-NEB-4_CAGATC_L002_R1_001.trim.single.fastq.gz \
  ILLUMINACLIP:trim/A549_DEX_SMC1A_rep2/chipseq-NEB-4_CAGATC_L002_R1_001.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/A549_DEX_SMC1A_rep2/chipseq-NEB-4_CAGATC_L002_R1_001.trim.log
trimmomatic.chipseq-NEB-4_CAGATC_L002_R1_001.09910caa8075017ba4f9349594b2932c.mugqic.done
)
trimmomatic_7_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"trimmomatic\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_SMC1A_rep2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"trimmomatic\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_SMC1A_rep2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=24G -N 1 -n 6 | grep "[0-9]" | cut -d\  -f4)
echo "$trimmomatic_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: trimmomatic_8_JOB_ID: trimmomatic.chipseq-NEB-5_GATCAG_L002_R1_001
#-------------------------------------------------------------------------------
JOB_NAME=trimmomatic.chipseq-NEB-5_GATCAG_L002_R1_001
JOB_DEPENDENCIES=
JOB_DONE=job_output/trimmomatic/trimmomatic.chipseq-NEB-5_GATCAG_L002_R1_001.ec4927419f3e4853747ba4da34b8985d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'trimmomatic.chipseq-NEB-5_GATCAG_L002_R1_001.ec4927419f3e4853747ba4da34b8985d.mugqic.done'
module load java/1.8.0_121 mugqic/trimmomatic/0.36 && \
mkdir -p trim/A549_DEX_NIPBL_rep2 && \
`cat > trim/A549_DEX_NIPBL_rep2/chipseq-NEB-5_GATCAG_L002_R1_001.trim.adapters.fa << END
>Single
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
END
` && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR SE \
  -threads 1 \
  -phred33 \
  /project/6001942/Working_Directory/Eric/CofactorHR/A549/raw/chip-seq/chipseq-NEB-5_GATCAG_L002_R1_001.fastq.gz \
  trim/A549_DEX_NIPBL_rep2/chipseq-NEB-5_GATCAG_L002_R1_001.trim.single.fastq.gz \
  ILLUMINACLIP:trim/A549_DEX_NIPBL_rep2/chipseq-NEB-5_GATCAG_L002_R1_001.trim.adapters.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:50 \
  2> trim/A549_DEX_NIPBL_rep2/chipseq-NEB-5_GATCAG_L002_R1_001.trim.log
trimmomatic.chipseq-NEB-5_GATCAG_L002_R1_001.ec4927419f3e4853747ba4da34b8985d.mugqic.done
)
trimmomatic_8_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"trimmomatic\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_NIPBL_rep2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"trimmomatic\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_NIPBL_rep2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:0 --mem=24G -N 1 -n 6 | grep "[0-9]" | cut -d\  -f4)
echo "$trimmomatic_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# STEP: merge_trimmomatic_stats
#-------------------------------------------------------------------------------
STEP=merge_trimmomatic_stats
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: merge_trimmomatic_stats_1_JOB_ID: merge_trimmomatic_stats
#-------------------------------------------------------------------------------
JOB_NAME=merge_trimmomatic_stats
JOB_DEPENDENCIES=$trimmomatic_1_JOB_ID:$trimmomatic_2_JOB_ID:$trimmomatic_3_JOB_ID:$trimmomatic_4_JOB_ID:$trimmomatic_5_JOB_ID:$trimmomatic_6_JOB_ID:$trimmomatic_7_JOB_ID:$trimmomatic_8_JOB_ID
JOB_DONE=job_output/merge_trimmomatic_stats/merge_trimmomatic_stats.8f738780d1cf5677197f04cef832251b.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'merge_trimmomatic_stats.8f738780d1cf5677197f04cef832251b.mugqic.done'
module load mugqic/pandoc/1.15.2 && \
mkdir -p metrics && \
echo 'Sample	Readset	Raw Single Reads #	Surviving Single Reads #	Surviving Single Reads %' > metrics/trimReadsetTable.tsv && \
grep ^Input trim/A549_DEX_MED1_rep2/chipseq-KAPA2_GCCAAT_L002_R1_001.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/A549_DEX_MED1_rep2	chipseq-KAPA2_GCCAAT_L002_R1_001	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/A549_DEX_SMC1A_rep2/chipseq-KAPA3_CTTGTA_L002_R1_001.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/A549_DEX_SMC1A_rep2	chipseq-KAPA3_CTTGTA_L002_R1_001	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/A549_DEX_SMC1A_rep2/chipseq-KAPA4_GTGAAA_L002_R1_001.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/A549_DEX_SMC1A_rep2	chipseq-KAPA4_GTGAAA_L002_R1_001	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/A549_DEX_MED1_rep2/chipseq-NEB-1_ATCACG_L002_R1_001.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/A549_DEX_MED1_rep2	chipseq-NEB-1_ATCACG_L002_R1_001	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/A549_DEX_MED1_rep2/chipseq-NEB-2_TTAGGC_L002_R1_001.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/A549_DEX_MED1_rep2	chipseq-NEB-2_TTAGGC_L002_R1_001	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/A549_DEX_SMC1A_rep2/chipseq-NEB-3_TGACCA_L002_R1_001.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/A549_DEX_SMC1A_rep2	chipseq-NEB-3_TGACCA_L002_R1_001	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/A549_DEX_SMC1A_rep2/chipseq-NEB-4_CAGATC_L002_R1_001.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/A549_DEX_SMC1A_rep2	chipseq-NEB-4_CAGATC_L002_R1_001	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
grep ^Input trim/A549_DEX_NIPBL_rep2/chipseq-NEB-5_GATCAG_L002_R1_001.trim.log | \
perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/A549_DEX_NIPBL_rep2	chipseq-NEB-5_GATCAG_L002_R1_001	\1	\2/' | \
awk '{OFS="	"; print $0, $4 / $3 * 100}' \
  >> metrics/trimReadsetTable.tsv && \
cut -f1,3- metrics/trimReadsetTable.tsv | awk -F"	" '{OFS="	"; if (NR==1) {if ($2=="Raw Paired Reads #") {paired=1};print "Sample", "Raw Reads #", "Surviving Reads #", "Surviving %"} else {if (paired) {$2=$2*2; $3=$3*2}; raw[$1]+=$2; surviving[$1]+=$3}}END{for (sample in raw){print sample, raw[sample], surviving[sample], surviving[sample] / raw[sample] * 100}}' \
  > metrics/trimSampleTable.tsv && \
mkdir -p report && \
cp metrics/trimReadsetTable.tsv metrics/trimSampleTable.tsv report/ && \
trim_readset_table_md=`LC_NUMERIC=en_CA awk -F "	" '{OFS="|"; if (NR == 1) {$1 = $1; print $0; print "-----|-----|-----:|-----:|-----:"} else {print $1, $2, sprintf("%\47d", $3), sprintf("%\47d", $4), sprintf("%.1f", $5)}}' metrics/trimReadsetTable.tsv` && \
pandoc \
  /home/efournie/genpipes/bfx/report/Illumina.merge_trimmomatic_stats.md \
  --template /home/efournie/genpipes/bfx/report/Illumina.merge_trimmomatic_stats.md \
  --variable trailing_min_quality=30 \
  --variable min_length=50 \
  --variable read_type=Single \
  --variable trim_readset_table="$trim_readset_table_md" \
  --to markdown \
  > report/Illumina.merge_trimmomatic_stats.md
merge_trimmomatic_stats.8f738780d1cf5677197f04cef832251b.mugqic.done
)
merge_trimmomatic_stats_1_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"merge_trimmomatic_stats\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_MED1_rep2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_SMC1A_rep2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_NIPBL_rep2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"merge_trimmomatic_stats\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_MED1_rep2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_SMC1A_rep2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_NIPBL_rep2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$merge_trimmomatic_stats_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# STEP: bwa_mem_picard_sort_sam
#-------------------------------------------------------------------------------
STEP=bwa_mem_picard_sort_sam
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_1_JOB_ID: bwa_mem_picard_sort_sam.chipseq-KAPA2_GCCAAT_L002_R1_001
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.chipseq-KAPA2_GCCAAT_L002_R1_001
JOB_DEPENDENCIES=$trimmomatic_1_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.chipseq-KAPA2_GCCAAT_L002_R1_001.0bcb9f9cfc678ff96af3fe01fa6c31c1.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.chipseq-KAPA2_GCCAAT_L002_R1_001.0bcb9f9cfc678ff96af3fe01fa6c31c1.mugqic.done'
module load mugqic/bwa/0.7.12 java/1.8.0_121 mugqic/picard/2.0.1 && \
mkdir -p alignment/A549_DEX_MED1_rep2/chipseq-KAPA2_GCCAAT_L002_R1_001 && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:chipseq-KAPA2_GCCAAT_L002_R1_001	SM:A549_DEX_MED1_rep2	LB:A549_DEX_MED1_rep2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/A549_DEX_MED1_rep2/chipseq-KAPA2_GCCAAT_L002_R1_001.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/A549_DEX_MED1_rep2/chipseq-KAPA2_GCCAAT_L002_R1_001/chipseq-KAPA2_GCCAAT_L002_R1_001.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.chipseq-KAPA2_GCCAAT_L002_R1_001.0bcb9f9cfc678ff96af3fe01fa6c31c1.mugqic.done
)
bwa_mem_picard_sort_sam_1_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"bwa_mem_picard_sort_sam\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_MED1_rep2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"bwa_mem_picard_sort_sam\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_MED1_rep2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=48G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$bwa_mem_picard_sort_sam_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_2_JOB_ID: bwa_mem_picard_sort_sam.chipseq-KAPA3_CTTGTA_L002_R1_001
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.chipseq-KAPA3_CTTGTA_L002_R1_001
JOB_DEPENDENCIES=$trimmomatic_2_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.chipseq-KAPA3_CTTGTA_L002_R1_001.d9abae60e62c183d04acec8e9aed7e70.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.chipseq-KAPA3_CTTGTA_L002_R1_001.d9abae60e62c183d04acec8e9aed7e70.mugqic.done'
module load mugqic/bwa/0.7.12 java/1.8.0_121 mugqic/picard/2.0.1 && \
mkdir -p alignment/A549_DEX_SMC1A_rep2/chipseq-KAPA3_CTTGTA_L002_R1_001 && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:chipseq-KAPA3_CTTGTA_L002_R1_001	SM:A549_DEX_SMC1A_rep2	LB:A549_DEX_SMC1A_rep2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/A549_DEX_SMC1A_rep2/chipseq-KAPA3_CTTGTA_L002_R1_001.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/A549_DEX_SMC1A_rep2/chipseq-KAPA3_CTTGTA_L002_R1_001/chipseq-KAPA3_CTTGTA_L002_R1_001.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.chipseq-KAPA3_CTTGTA_L002_R1_001.d9abae60e62c183d04acec8e9aed7e70.mugqic.done
)
bwa_mem_picard_sort_sam_2_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"bwa_mem_picard_sort_sam\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_SMC1A_rep2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"bwa_mem_picard_sort_sam\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_SMC1A_rep2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=48G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$bwa_mem_picard_sort_sam_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_3_JOB_ID: bwa_mem_picard_sort_sam.chipseq-KAPA4_GTGAAA_L002_R1_001
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.chipseq-KAPA4_GTGAAA_L002_R1_001
JOB_DEPENDENCIES=$trimmomatic_3_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.chipseq-KAPA4_GTGAAA_L002_R1_001.50d808628330b187d3fe032823609a33.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.chipseq-KAPA4_GTGAAA_L002_R1_001.50d808628330b187d3fe032823609a33.mugqic.done'
module load mugqic/bwa/0.7.12 java/1.8.0_121 mugqic/picard/2.0.1 && \
mkdir -p alignment/A549_DEX_SMC1A_rep2/chipseq-KAPA4_GTGAAA_L002_R1_001 && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:chipseq-KAPA4_GTGAAA_L002_R1_001	SM:A549_DEX_SMC1A_rep2	LB:A549_DEX_SMC1A_rep2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/A549_DEX_SMC1A_rep2/chipseq-KAPA4_GTGAAA_L002_R1_001.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/A549_DEX_SMC1A_rep2/chipseq-KAPA4_GTGAAA_L002_R1_001/chipseq-KAPA4_GTGAAA_L002_R1_001.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.chipseq-KAPA4_GTGAAA_L002_R1_001.50d808628330b187d3fe032823609a33.mugqic.done
)
bwa_mem_picard_sort_sam_3_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"bwa_mem_picard_sort_sam\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_SMC1A_rep2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"bwa_mem_picard_sort_sam\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_SMC1A_rep2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=48G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$bwa_mem_picard_sort_sam_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_4_JOB_ID: bwa_mem_picard_sort_sam.chipseq-NEB-1_ATCACG_L002_R1_001
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.chipseq-NEB-1_ATCACG_L002_R1_001
JOB_DEPENDENCIES=$trimmomatic_4_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.chipseq-NEB-1_ATCACG_L002_R1_001.22d60d85eb72216115abc94f5a1e7409.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.chipseq-NEB-1_ATCACG_L002_R1_001.22d60d85eb72216115abc94f5a1e7409.mugqic.done'
module load mugqic/bwa/0.7.12 java/1.8.0_121 mugqic/picard/2.0.1 && \
mkdir -p alignment/A549_DEX_MED1_rep2/chipseq-NEB-1_ATCACG_L002_R1_001 && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:chipseq-NEB-1_ATCACG_L002_R1_001	SM:A549_DEX_MED1_rep2	LB:A549_DEX_MED1_rep2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/A549_DEX_MED1_rep2/chipseq-NEB-1_ATCACG_L002_R1_001.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/A549_DEX_MED1_rep2/chipseq-NEB-1_ATCACG_L002_R1_001/chipseq-NEB-1_ATCACG_L002_R1_001.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.chipseq-NEB-1_ATCACG_L002_R1_001.22d60d85eb72216115abc94f5a1e7409.mugqic.done
)
bwa_mem_picard_sort_sam_4_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"bwa_mem_picard_sort_sam\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_MED1_rep2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"bwa_mem_picard_sort_sam\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_MED1_rep2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=48G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$bwa_mem_picard_sort_sam_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_5_JOB_ID: bwa_mem_picard_sort_sam.chipseq-NEB-2_TTAGGC_L002_R1_001
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.chipseq-NEB-2_TTAGGC_L002_R1_001
JOB_DEPENDENCIES=$trimmomatic_5_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.chipseq-NEB-2_TTAGGC_L002_R1_001.4c64983f747ce1270bf9576b7dd9ec7c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.chipseq-NEB-2_TTAGGC_L002_R1_001.4c64983f747ce1270bf9576b7dd9ec7c.mugqic.done'
module load mugqic/bwa/0.7.12 java/1.8.0_121 mugqic/picard/2.0.1 && \
mkdir -p alignment/A549_DEX_MED1_rep2/chipseq-NEB-2_TTAGGC_L002_R1_001 && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:chipseq-NEB-2_TTAGGC_L002_R1_001	SM:A549_DEX_MED1_rep2	LB:A549_DEX_MED1_rep2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/A549_DEX_MED1_rep2/chipseq-NEB-2_TTAGGC_L002_R1_001.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/A549_DEX_MED1_rep2/chipseq-NEB-2_TTAGGC_L002_R1_001/chipseq-NEB-2_TTAGGC_L002_R1_001.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.chipseq-NEB-2_TTAGGC_L002_R1_001.4c64983f747ce1270bf9576b7dd9ec7c.mugqic.done
)
bwa_mem_picard_sort_sam_5_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"bwa_mem_picard_sort_sam\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_MED1_rep2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"bwa_mem_picard_sort_sam\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_MED1_rep2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=48G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$bwa_mem_picard_sort_sam_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_6_JOB_ID: bwa_mem_picard_sort_sam.chipseq-NEB-3_TGACCA_L002_R1_001
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.chipseq-NEB-3_TGACCA_L002_R1_001
JOB_DEPENDENCIES=$trimmomatic_6_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.chipseq-NEB-3_TGACCA_L002_R1_001.214a1bebeac91b8331efc03cebe0e858.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.chipseq-NEB-3_TGACCA_L002_R1_001.214a1bebeac91b8331efc03cebe0e858.mugqic.done'
module load mugqic/bwa/0.7.12 java/1.8.0_121 mugqic/picard/2.0.1 && \
mkdir -p alignment/A549_DEX_SMC1A_rep2/chipseq-NEB-3_TGACCA_L002_R1_001 && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:chipseq-NEB-3_TGACCA_L002_R1_001	SM:A549_DEX_SMC1A_rep2	LB:A549_DEX_SMC1A_rep2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/A549_DEX_SMC1A_rep2/chipseq-NEB-3_TGACCA_L002_R1_001.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/A549_DEX_SMC1A_rep2/chipseq-NEB-3_TGACCA_L002_R1_001/chipseq-NEB-3_TGACCA_L002_R1_001.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.chipseq-NEB-3_TGACCA_L002_R1_001.214a1bebeac91b8331efc03cebe0e858.mugqic.done
)
bwa_mem_picard_sort_sam_6_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"bwa_mem_picard_sort_sam\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_SMC1A_rep2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"bwa_mem_picard_sort_sam\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_SMC1A_rep2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=48G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$bwa_mem_picard_sort_sam_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_7_JOB_ID: bwa_mem_picard_sort_sam.chipseq-NEB-4_CAGATC_L002_R1_001
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.chipseq-NEB-4_CAGATC_L002_R1_001
JOB_DEPENDENCIES=$trimmomatic_7_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.chipseq-NEB-4_CAGATC_L002_R1_001.e73c19497d4cb2754155ead565accab6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.chipseq-NEB-4_CAGATC_L002_R1_001.e73c19497d4cb2754155ead565accab6.mugqic.done'
module load mugqic/bwa/0.7.12 java/1.8.0_121 mugqic/picard/2.0.1 && \
mkdir -p alignment/A549_DEX_SMC1A_rep2/chipseq-NEB-4_CAGATC_L002_R1_001 && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:chipseq-NEB-4_CAGATC_L002_R1_001	SM:A549_DEX_SMC1A_rep2	LB:A549_DEX_SMC1A_rep2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/A549_DEX_SMC1A_rep2/chipseq-NEB-4_CAGATC_L002_R1_001.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/A549_DEX_SMC1A_rep2/chipseq-NEB-4_CAGATC_L002_R1_001/chipseq-NEB-4_CAGATC_L002_R1_001.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.chipseq-NEB-4_CAGATC_L002_R1_001.e73c19497d4cb2754155ead565accab6.mugqic.done
)
bwa_mem_picard_sort_sam_7_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"bwa_mem_picard_sort_sam\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_SMC1A_rep2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"bwa_mem_picard_sort_sam\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_SMC1A_rep2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=48G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$bwa_mem_picard_sort_sam_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_8_JOB_ID: bwa_mem_picard_sort_sam.chipseq-NEB-5_GATCAG_L002_R1_001
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam.chipseq-NEB-5_GATCAG_L002_R1_001
JOB_DEPENDENCIES=$trimmomatic_8_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam.chipseq-NEB-5_GATCAG_L002_R1_001.420da78fa14d62dc47841cc4665597c3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam.chipseq-NEB-5_GATCAG_L002_R1_001.420da78fa14d62dc47841cc4665597c3.mugqic.done'
module load mugqic/bwa/0.7.12 java/1.8.0_121 mugqic/picard/2.0.1 && \
mkdir -p alignment/A549_DEX_NIPBL_rep2/chipseq-NEB-5_GATCAG_L002_R1_001 && \
bwa mem  \
  -M -t 11 \
  -R '@RG	ID:chipseq-NEB-5_GATCAG_L002_R1_001	SM:A549_DEX_NIPBL_rep2	LB:A549_DEX_NIPBL_rep2	CN:McGill University and Genome Quebec Innovation Centre	PL:Illumina' \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/bwa_index/Homo_sapiens.GRCh38.fa \
  trim/A549_DEX_NIPBL_rep2/chipseq-NEB-5_GATCAG_L002_R1_001.trim.single.fastq.gz | \
 java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx54G -jar $PICARD_HOME/picard.jar SortSam \
 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=/dev/stdin \
 OUTPUT=alignment/A549_DEX_NIPBL_rep2/chipseq-NEB-5_GATCAG_L002_R1_001/chipseq-NEB-5_GATCAG_L002_R1_001.sorted.bam \
 SORT_ORDER=coordinate \
 MAX_RECORDS_IN_RAM=13500000
bwa_mem_picard_sort_sam.chipseq-NEB-5_GATCAG_L002_R1_001.420da78fa14d62dc47841cc4665597c3.mugqic.done
)
bwa_mem_picard_sort_sam_8_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"bwa_mem_picard_sort_sam\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_NIPBL_rep2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"bwa_mem_picard_sort_sam\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_NIPBL_rep2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=48G -N 1 -n 12 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$bwa_mem_picard_sort_sam_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: bwa_mem_picard_sort_sam_9_JOB_ID: bwa_mem_picard_sort_sam_report
#-------------------------------------------------------------------------------
JOB_NAME=bwa_mem_picard_sort_sam_report
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_1_JOB_ID:$bwa_mem_picard_sort_sam_2_JOB_ID:$bwa_mem_picard_sort_sam_3_JOB_ID:$bwa_mem_picard_sort_sam_4_JOB_ID:$bwa_mem_picard_sort_sam_5_JOB_ID:$bwa_mem_picard_sort_sam_6_JOB_ID:$bwa_mem_picard_sort_sam_7_JOB_ID:$bwa_mem_picard_sort_sam_8_JOB_ID
JOB_DONE=job_output/bwa_mem_picard_sort_sam/bwa_mem_picard_sort_sam_report.22d7ff6198a068bd5763589f2515f7f6.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'bwa_mem_picard_sort_sam_report.22d7ff6198a068bd5763589f2515f7f6.mugqic.done'
module load mugqic/pandoc/1.15.2 && \
mkdir -p report && \
pandoc --to=markdown \
  --template /home/efournie/genpipes/bfx/report/DnaSeq.bwa_mem_picard_sort_sam.md \
  --variable scientific_name="Homo_sapiens" \
  --variable assembly="GRCh38" \
  /home/efournie/genpipes/bfx/report/DnaSeq.bwa_mem_picard_sort_sam.md \
  > report/DnaSeq.bwa_mem_picard_sort_sam.md
bwa_mem_picard_sort_sam_report.22d7ff6198a068bd5763589f2515f7f6.mugqic.done
)
bwa_mem_picard_sort_sam_9_JOB_ID=$(echo "#! /bin/bash 
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
echo "$bwa_mem_picard_sort_sam_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# STEP: samtools_view_filter
#-------------------------------------------------------------------------------
STEP=samtools_view_filter
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_1_JOB_ID: samtools_view_filter.chipseq-KAPA2_GCCAAT_L002_R1_001
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.chipseq-KAPA2_GCCAAT_L002_R1_001
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_1_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.chipseq-KAPA2_GCCAAT_L002_R1_001.2935c239f85202375560e2c2b1669b15.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.chipseq-KAPA2_GCCAAT_L002_R1_001.2935c239f85202375560e2c2b1669b15.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/A549_DEX_MED1_rep2/chipseq-KAPA2_GCCAAT_L002_R1_001/chipseq-KAPA2_GCCAAT_L002_R1_001.sorted.bam \
  > alignment/A549_DEX_MED1_rep2/chipseq-KAPA2_GCCAAT_L002_R1_001/chipseq-KAPA2_GCCAAT_L002_R1_001.sorted.filtered.bam
samtools_view_filter.chipseq-KAPA2_GCCAAT_L002_R1_001.2935c239f85202375560e2c2b1669b15.mugqic.done
)
samtools_view_filter_1_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"samtools_view_filter\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_MED1_rep2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"samtools_view_filter\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_MED1_rep2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$samtools_view_filter_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_2_JOB_ID: samtools_view_filter.chipseq-KAPA3_CTTGTA_L002_R1_001
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.chipseq-KAPA3_CTTGTA_L002_R1_001
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_2_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.chipseq-KAPA3_CTTGTA_L002_R1_001.e46f362dedceaf9803699ffa06539eb7.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.chipseq-KAPA3_CTTGTA_L002_R1_001.e46f362dedceaf9803699ffa06539eb7.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/A549_DEX_SMC1A_rep2/chipseq-KAPA3_CTTGTA_L002_R1_001/chipseq-KAPA3_CTTGTA_L002_R1_001.sorted.bam \
  > alignment/A549_DEX_SMC1A_rep2/chipseq-KAPA3_CTTGTA_L002_R1_001/chipseq-KAPA3_CTTGTA_L002_R1_001.sorted.filtered.bam
samtools_view_filter.chipseq-KAPA3_CTTGTA_L002_R1_001.e46f362dedceaf9803699ffa06539eb7.mugqic.done
)
samtools_view_filter_2_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"samtools_view_filter\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_SMC1A_rep2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"samtools_view_filter\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_SMC1A_rep2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$samtools_view_filter_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_3_JOB_ID: samtools_view_filter.chipseq-KAPA4_GTGAAA_L002_R1_001
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.chipseq-KAPA4_GTGAAA_L002_R1_001
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_3_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.chipseq-KAPA4_GTGAAA_L002_R1_001.f595ae51fb247003e5f3cacec2f39e38.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.chipseq-KAPA4_GTGAAA_L002_R1_001.f595ae51fb247003e5f3cacec2f39e38.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/A549_DEX_SMC1A_rep2/chipseq-KAPA4_GTGAAA_L002_R1_001/chipseq-KAPA4_GTGAAA_L002_R1_001.sorted.bam \
  > alignment/A549_DEX_SMC1A_rep2/chipseq-KAPA4_GTGAAA_L002_R1_001/chipseq-KAPA4_GTGAAA_L002_R1_001.sorted.filtered.bam
samtools_view_filter.chipseq-KAPA4_GTGAAA_L002_R1_001.f595ae51fb247003e5f3cacec2f39e38.mugqic.done
)
samtools_view_filter_3_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"samtools_view_filter\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_SMC1A_rep2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"samtools_view_filter\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_SMC1A_rep2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$samtools_view_filter_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_4_JOB_ID: samtools_view_filter.chipseq-NEB-1_ATCACG_L002_R1_001
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.chipseq-NEB-1_ATCACG_L002_R1_001
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_4_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.chipseq-NEB-1_ATCACG_L002_R1_001.7e592c7a2f43bce880dc0ec298a5a00e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.chipseq-NEB-1_ATCACG_L002_R1_001.7e592c7a2f43bce880dc0ec298a5a00e.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/A549_DEX_MED1_rep2/chipseq-NEB-1_ATCACG_L002_R1_001/chipseq-NEB-1_ATCACG_L002_R1_001.sorted.bam \
  > alignment/A549_DEX_MED1_rep2/chipseq-NEB-1_ATCACG_L002_R1_001/chipseq-NEB-1_ATCACG_L002_R1_001.sorted.filtered.bam
samtools_view_filter.chipseq-NEB-1_ATCACG_L002_R1_001.7e592c7a2f43bce880dc0ec298a5a00e.mugqic.done
)
samtools_view_filter_4_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"samtools_view_filter\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_MED1_rep2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"samtools_view_filter\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_MED1_rep2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$samtools_view_filter_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_5_JOB_ID: samtools_view_filter.chipseq-NEB-2_TTAGGC_L002_R1_001
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.chipseq-NEB-2_TTAGGC_L002_R1_001
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_5_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.chipseq-NEB-2_TTAGGC_L002_R1_001.7d122d29440302a92399ac9e07353866.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.chipseq-NEB-2_TTAGGC_L002_R1_001.7d122d29440302a92399ac9e07353866.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/A549_DEX_MED1_rep2/chipseq-NEB-2_TTAGGC_L002_R1_001/chipseq-NEB-2_TTAGGC_L002_R1_001.sorted.bam \
  > alignment/A549_DEX_MED1_rep2/chipseq-NEB-2_TTAGGC_L002_R1_001/chipseq-NEB-2_TTAGGC_L002_R1_001.sorted.filtered.bam
samtools_view_filter.chipseq-NEB-2_TTAGGC_L002_R1_001.7d122d29440302a92399ac9e07353866.mugqic.done
)
samtools_view_filter_5_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"samtools_view_filter\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_MED1_rep2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"samtools_view_filter\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_MED1_rep2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$samtools_view_filter_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_6_JOB_ID: samtools_view_filter.chipseq-NEB-3_TGACCA_L002_R1_001
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.chipseq-NEB-3_TGACCA_L002_R1_001
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_6_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.chipseq-NEB-3_TGACCA_L002_R1_001.a617063655e620e09c4d1fb57d360abb.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.chipseq-NEB-3_TGACCA_L002_R1_001.a617063655e620e09c4d1fb57d360abb.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/A549_DEX_SMC1A_rep2/chipseq-NEB-3_TGACCA_L002_R1_001/chipseq-NEB-3_TGACCA_L002_R1_001.sorted.bam \
  > alignment/A549_DEX_SMC1A_rep2/chipseq-NEB-3_TGACCA_L002_R1_001/chipseq-NEB-3_TGACCA_L002_R1_001.sorted.filtered.bam
samtools_view_filter.chipseq-NEB-3_TGACCA_L002_R1_001.a617063655e620e09c4d1fb57d360abb.mugqic.done
)
samtools_view_filter_6_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"samtools_view_filter\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_SMC1A_rep2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"samtools_view_filter\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_SMC1A_rep2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$samtools_view_filter_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_7_JOB_ID: samtools_view_filter.chipseq-NEB-4_CAGATC_L002_R1_001
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.chipseq-NEB-4_CAGATC_L002_R1_001
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_7_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.chipseq-NEB-4_CAGATC_L002_R1_001.656b3ce28ec4f68315c0203d4fa386cf.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.chipseq-NEB-4_CAGATC_L002_R1_001.656b3ce28ec4f68315c0203d4fa386cf.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/A549_DEX_SMC1A_rep2/chipseq-NEB-4_CAGATC_L002_R1_001/chipseq-NEB-4_CAGATC_L002_R1_001.sorted.bam \
  > alignment/A549_DEX_SMC1A_rep2/chipseq-NEB-4_CAGATC_L002_R1_001/chipseq-NEB-4_CAGATC_L002_R1_001.sorted.filtered.bam
samtools_view_filter.chipseq-NEB-4_CAGATC_L002_R1_001.656b3ce28ec4f68315c0203d4fa386cf.mugqic.done
)
samtools_view_filter_7_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"samtools_view_filter\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_SMC1A_rep2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"samtools_view_filter\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_SMC1A_rep2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$samtools_view_filter_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_8_JOB_ID: samtools_view_filter.chipseq-NEB-5_GATCAG_L002_R1_001
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter.chipseq-NEB-5_GATCAG_L002_R1_001
JOB_DEPENDENCIES=$bwa_mem_picard_sort_sam_8_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter.chipseq-NEB-5_GATCAG_L002_R1_001.0f3c60ef09c7944ea1b8f6821525a212.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter.chipseq-NEB-5_GATCAG_L002_R1_001.0f3c60ef09c7944ea1b8f6821525a212.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools view -b -F4 -q 20 \
  alignment/A549_DEX_NIPBL_rep2/chipseq-NEB-5_GATCAG_L002_R1_001/chipseq-NEB-5_GATCAG_L002_R1_001.sorted.bam \
  > alignment/A549_DEX_NIPBL_rep2/chipseq-NEB-5_GATCAG_L002_R1_001/chipseq-NEB-5_GATCAG_L002_R1_001.sorted.filtered.bam
samtools_view_filter.chipseq-NEB-5_GATCAG_L002_R1_001.0f3c60ef09c7944ea1b8f6821525a212.mugqic.done
)
samtools_view_filter_8_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"samtools_view_filter\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_NIPBL_rep2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"samtools_view_filter\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_NIPBL_rep2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$samtools_view_filter_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: samtools_view_filter_9_JOB_ID: samtools_view_filter_report
#-------------------------------------------------------------------------------
JOB_NAME=samtools_view_filter_report
JOB_DEPENDENCIES=$samtools_view_filter_1_JOB_ID:$samtools_view_filter_2_JOB_ID:$samtools_view_filter_3_JOB_ID:$samtools_view_filter_4_JOB_ID:$samtools_view_filter_5_JOB_ID:$samtools_view_filter_6_JOB_ID:$samtools_view_filter_7_JOB_ID:$samtools_view_filter_8_JOB_ID
JOB_DONE=job_output/samtools_view_filter/samtools_view_filter_report.adda516c06f797351611d7331668a33d.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'samtools_view_filter_report.adda516c06f797351611d7331668a33d.mugqic.done'
module load mugqic/pandoc/1.15.2 && \
mkdir -p report && \
pandoc --to=markdown \
  --template /home/efournie/genpipes/bfx/report/ChipSeq.samtools_view_filter.md \
  --variable min_mapq="20" \
  /home/efournie/genpipes/bfx/report/ChipSeq.samtools_view_filter.md \
  > report/ChipSeq.samtools_view_filter.md
samtools_view_filter_report.adda516c06f797351611d7331668a33d.mugqic.done
)
samtools_view_filter_9_JOB_ID=$(echo "#! /bin/bash 
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
echo "$samtools_view_filter_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# STEP: picard_merge_sam_files
#-------------------------------------------------------------------------------
STEP=picard_merge_sam_files
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_1_JOB_ID: picard_merge_sam_files.A549_DEX_MED1_rep2
#-------------------------------------------------------------------------------
JOB_NAME=picard_merge_sam_files.A549_DEX_MED1_rep2
JOB_DEPENDENCIES=$samtools_view_filter_1_JOB_ID:$samtools_view_filter_4_JOB_ID:$samtools_view_filter_5_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/picard_merge_sam_files.A549_DEX_MED1_rep2.1c310cf6f60e5384c343468a7ec12a99.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_merge_sam_files.A549_DEX_MED1_rep2.1c310cf6f60e5384c343468a7ec12a99.mugqic.done'
module load java/1.8.0_121 mugqic/picard/2.0.1 && \
mkdir -p alignment/A549_DEX_MED1_rep2 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx1700M -jar $PICARD_HOME/picard.jar MergeSamFiles \
 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/A549_DEX_MED1_rep2/chipseq-KAPA2_GCCAAT_L002_R1_001/chipseq-KAPA2_GCCAAT_L002_R1_001.sorted.filtered.bam \
 INPUT=alignment/A549_DEX_MED1_rep2/chipseq-NEB-1_ATCACG_L002_R1_001/chipseq-NEB-1_ATCACG_L002_R1_001.sorted.filtered.bam \
 INPUT=alignment/A549_DEX_MED1_rep2/chipseq-NEB-2_TTAGGC_L002_R1_001/chipseq-NEB-2_TTAGGC_L002_R1_001.sorted.filtered.bam \
 OUTPUT=alignment/A549_DEX_MED1_rep2/A549_DEX_MED1_rep2.merged.bam \
 MAX_RECORDS_IN_RAM=250000
picard_merge_sam_files.A549_DEX_MED1_rep2.1c310cf6f60e5384c343468a7ec12a99.mugqic.done
)
picard_merge_sam_files_1_JOB_ID=$(echo "#! /bin/bash 
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:0 --mem=8G -N 1 -n 2 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$picard_merge_sam_files_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_2_JOB_ID: picard_merge_sam_files.A549_DEX_SMC1A_rep2
#-------------------------------------------------------------------------------
JOB_NAME=picard_merge_sam_files.A549_DEX_SMC1A_rep2
JOB_DEPENDENCIES=$samtools_view_filter_2_JOB_ID:$samtools_view_filter_3_JOB_ID:$samtools_view_filter_6_JOB_ID:$samtools_view_filter_7_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/picard_merge_sam_files.A549_DEX_SMC1A_rep2.ff5844bd65b97bcd6cbea0e39359a2f9.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_merge_sam_files.A549_DEX_SMC1A_rep2.ff5844bd65b97bcd6cbea0e39359a2f9.mugqic.done'
module load java/1.8.0_121 mugqic/picard/2.0.1 && \
mkdir -p alignment/A549_DEX_SMC1A_rep2 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx1700M -jar $PICARD_HOME/picard.jar MergeSamFiles \
 VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/A549_DEX_SMC1A_rep2/chipseq-KAPA3_CTTGTA_L002_R1_001/chipseq-KAPA3_CTTGTA_L002_R1_001.sorted.filtered.bam \
 INPUT=alignment/A549_DEX_SMC1A_rep2/chipseq-KAPA4_GTGAAA_L002_R1_001/chipseq-KAPA4_GTGAAA_L002_R1_001.sorted.filtered.bam \
 INPUT=alignment/A549_DEX_SMC1A_rep2/chipseq-NEB-3_TGACCA_L002_R1_001/chipseq-NEB-3_TGACCA_L002_R1_001.sorted.filtered.bam \
 INPUT=alignment/A549_DEX_SMC1A_rep2/chipseq-NEB-4_CAGATC_L002_R1_001/chipseq-NEB-4_CAGATC_L002_R1_001.sorted.filtered.bam \
 OUTPUT=alignment/A549_DEX_SMC1A_rep2/A549_DEX_SMC1A_rep2.merged.bam \
 MAX_RECORDS_IN_RAM=250000
picard_merge_sam_files.A549_DEX_SMC1A_rep2.ff5844bd65b97bcd6cbea0e39359a2f9.mugqic.done
)
picard_merge_sam_files_2_JOB_ID=$(echo "#! /bin/bash 
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=35:00:0 --mem=8G -N 1 -n 2 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$picard_merge_sam_files_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: picard_merge_sam_files_3_JOB_ID: symlink_readset_sample_bam.A549_DEX_NIPBL_rep2
#-------------------------------------------------------------------------------
JOB_NAME=symlink_readset_sample_bam.A549_DEX_NIPBL_rep2
JOB_DEPENDENCIES=$samtools_view_filter_8_JOB_ID
JOB_DONE=job_output/picard_merge_sam_files/symlink_readset_sample_bam.A549_DEX_NIPBL_rep2.96609edb1858e61a3f18890e05240d82.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'symlink_readset_sample_bam.A549_DEX_NIPBL_rep2.96609edb1858e61a3f18890e05240d82.mugqic.done'
mkdir -p alignment/A549_DEX_NIPBL_rep2 && \
ln -s -f chipseq-NEB-5_GATCAG_L002_R1_001/chipseq-NEB-5_GATCAG_L002_R1_001.sorted.filtered.bam alignment/A549_DEX_NIPBL_rep2/A549_DEX_NIPBL_rep2.merged.bam
symlink_readset_sample_bam.A549_DEX_NIPBL_rep2.96609edb1858e61a3f18890e05240d82.mugqic.done
)
picard_merge_sam_files_3_JOB_ID=$(echo "#! /bin/bash 
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
echo "$picard_merge_sam_files_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# STEP: picard_mark_duplicates
#-------------------------------------------------------------------------------
STEP=picard_mark_duplicates
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_1_JOB_ID: picard_mark_duplicates.A549_DEX_MED1_rep2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.A549_DEX_MED1_rep2
JOB_DEPENDENCIES=$picard_merge_sam_files_1_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.A549_DEX_MED1_rep2.1fed2620d9143610c9562e6442fadd00.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.A549_DEX_MED1_rep2.1fed2620d9143610c9562e6442fadd00.mugqic.done'
module load java/1.8.0_121 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/A549_DEX_MED1_rep2/A549_DEX_MED1_rep2.merged.bam \
 OUTPUT=alignment/A549_DEX_MED1_rep2/A549_DEX_MED1_rep2.sorted.dup.bam \
 METRICS_FILE=alignment/A549_DEX_MED1_rep2/A549_DEX_MED1_rep2.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.A549_DEX_MED1_rep2.1fed2620d9143610c9562e6442fadd00.mugqic.done
)
picard_mark_duplicates_1_JOB_ID=$(echo "#! /bin/bash 
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=48:00:0 --mem=8G -N 1 -n 2 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$picard_mark_duplicates_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_2_JOB_ID: picard_mark_duplicates.A549_DEX_SMC1A_rep2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.A549_DEX_SMC1A_rep2
JOB_DEPENDENCIES=$picard_merge_sam_files_2_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.A549_DEX_SMC1A_rep2.9141aae28d8bbc3ed1a9a5dd5112f352.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.A549_DEX_SMC1A_rep2.9141aae28d8bbc3ed1a9a5dd5112f352.mugqic.done'
module load java/1.8.0_121 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/A549_DEX_SMC1A_rep2/A549_DEX_SMC1A_rep2.merged.bam \
 OUTPUT=alignment/A549_DEX_SMC1A_rep2/A549_DEX_SMC1A_rep2.sorted.dup.bam \
 METRICS_FILE=alignment/A549_DEX_SMC1A_rep2/A549_DEX_SMC1A_rep2.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.A549_DEX_SMC1A_rep2.9141aae28d8bbc3ed1a9a5dd5112f352.mugqic.done
)
picard_mark_duplicates_2_JOB_ID=$(echo "#! /bin/bash 
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=48:00:0 --mem=8G -N 1 -n 2 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$picard_mark_duplicates_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_3_JOB_ID: picard_mark_duplicates.A549_DEX_NIPBL_rep2
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates.A549_DEX_NIPBL_rep2
JOB_DEPENDENCIES=$picard_merge_sam_files_3_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates.A549_DEX_NIPBL_rep2.81429933aa303fce625574c2449af224.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates.A549_DEX_NIPBL_rep2.81429933aa303fce625574c2449af224.mugqic.done'
module load java/1.8.0_121 mugqic/picard/2.0.1 && \
java -Djava.io.tmpdir=${SLURM_TMPDIR} -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
 REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
 TMP_DIR=${SLURM_TMPDIR} \
 INPUT=alignment/A549_DEX_NIPBL_rep2/A549_DEX_NIPBL_rep2.merged.bam \
 OUTPUT=alignment/A549_DEX_NIPBL_rep2/A549_DEX_NIPBL_rep2.sorted.dup.bam \
 METRICS_FILE=alignment/A549_DEX_NIPBL_rep2/A549_DEX_NIPBL_rep2.sorted.dup.metrics \
 MAX_RECORDS_IN_RAM=1000000
picard_mark_duplicates.A549_DEX_NIPBL_rep2.81429933aa303fce625574c2449af224.mugqic.done
)
picard_mark_duplicates_3_JOB_ID=$(echo "#! /bin/bash 
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=48:00:0 --mem=8G -N 1 -n 2 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$picard_mark_duplicates_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates_4_JOB_ID: picard_mark_duplicates_report
#-------------------------------------------------------------------------------
JOB_NAME=picard_mark_duplicates_report
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_2_JOB_ID:$picard_mark_duplicates_3_JOB_ID
JOB_DONE=job_output/picard_mark_duplicates/picard_mark_duplicates_report.8bfd4b03ed7ee119ffcd84a77e108045.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'picard_mark_duplicates_report.8bfd4b03ed7ee119ffcd84a77e108045.mugqic.done'
mkdir -p report && \
cp \
  /home/efournie/genpipes/bfx/report/ChipSeq.picard_mark_duplicates.md \
  report/ChipSeq.picard_mark_duplicates.md
picard_mark_duplicates_report.8bfd4b03ed7ee119ffcd84a77e108045.mugqic.done
)
picard_mark_duplicates_4_JOB_ID=$(echo "#! /bin/bash 
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
echo "$picard_mark_duplicates_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# STEP: metrics
#-------------------------------------------------------------------------------
STEP=metrics
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: metrics_1_JOB_ID: metrics.flagstat
#-------------------------------------------------------------------------------
JOB_NAME=metrics.flagstat
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID:$picard_mark_duplicates_2_JOB_ID:$picard_mark_duplicates_3_JOB_ID
JOB_DONE=job_output/metrics/metrics.flagstat.4a8c2927f5741d458261ec37aaf0cf44.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'metrics.flagstat.4a8c2927f5741d458261ec37aaf0cf44.mugqic.done'
module load mugqic/samtools/1.3.1 && \
samtools flagstat \
  alignment/A549_DEX_MED1_rep2/A549_DEX_MED1_rep2.sorted.dup.bam \
  > alignment/A549_DEX_MED1_rep2/A549_DEX_MED1_rep2.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/A549_DEX_SMC1A_rep2/A549_DEX_SMC1A_rep2.sorted.dup.bam \
  > alignment/A549_DEX_SMC1A_rep2/A549_DEX_SMC1A_rep2.sorted.dup.bam.flagstat && \
samtools flagstat \
  alignment/A549_DEX_NIPBL_rep2/A549_DEX_NIPBL_rep2.sorted.dup.bam \
  > alignment/A549_DEX_NIPBL_rep2/A549_DEX_NIPBL_rep2.sorted.dup.bam.flagstat
metrics.flagstat.4a8c2927f5741d458261ec37aaf0cf44.mugqic.done
)
metrics_1_JOB_ID=$(echo "#! /bin/bash 
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
echo "$metrics_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: metrics_2_JOB_ID: metrics_report
#-------------------------------------------------------------------------------
JOB_NAME=metrics_report
JOB_DEPENDENCIES=$metrics_1_JOB_ID
JOB_DONE=job_output/metrics/metrics_report.93c783539badb5118acf6fc6375c616c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'metrics_report.93c783539badb5118acf6fc6375c616c.mugqic.done'
module load mugqic/pandoc/1.15.2 && \
for sample in A549_DEX_MED1_rep2 A549_DEX_SMC1A_rep2 A549_DEX_NIPBL_rep2
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
  --template /home/efournie/genpipes/bfx/report/ChipSeq.metrics.md \
  --variable trim_mem_sample_table="$trim_mem_sample_table" \
  /home/efournie/genpipes/bfx/report/ChipSeq.metrics.md \
  > report/ChipSeq.metrics.md

metrics_report.93c783539badb5118acf6fc6375c616c.mugqic.done
)
metrics_2_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"metrics\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_MED1_rep2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_SMC1A_rep2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_NIPBL_rep2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"metrics\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_MED1_rep2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_SMC1A_rep2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_NIPBL_rep2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$metrics_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# STEP: homer_make_tag_directory
#-------------------------------------------------------------------------------
STEP=homer_make_tag_directory
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_1_JOB_ID: homer_make_tag_directory.A549_DEX_MED1_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.A549_DEX_MED1_rep2
JOB_DEPENDENCIES=$picard_mark_duplicates_1_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.A549_DEX_MED1_rep2.05311f9da02591c574b2bf637a1be7b0.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.A549_DEX_MED1_rep2.05311f9da02591c574b2bf637a1be7b0.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/A549_DEX_MED1_rep2 \
            alignment/A549_DEX_MED1_rep2/A549_DEX_MED1_rep2.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.A549_DEX_MED1_rep2.05311f9da02591c574b2bf637a1be7b0.mugqic.done
)
homer_make_tag_directory_1_JOB_ID=$(echo "#! /bin/bash 
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
echo "$homer_make_tag_directory_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_2_JOB_ID: homer_make_tag_directory.A549_DEX_SMC1A_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.A549_DEX_SMC1A_rep2
JOB_DEPENDENCIES=$picard_mark_duplicates_2_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.A549_DEX_SMC1A_rep2.f9c22d757345a90afc63bb1374b9d64c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.A549_DEX_SMC1A_rep2.f9c22d757345a90afc63bb1374b9d64c.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/A549_DEX_SMC1A_rep2 \
            alignment/A549_DEX_SMC1A_rep2/A549_DEX_SMC1A_rep2.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.A549_DEX_SMC1A_rep2.f9c22d757345a90afc63bb1374b9d64c.mugqic.done
)
homer_make_tag_directory_2_JOB_ID=$(echo "#! /bin/bash 
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
echo "$homer_make_tag_directory_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: homer_make_tag_directory_3_JOB_ID: homer_make_tag_directory.A549_DEX_NIPBL_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_tag_directory.A549_DEX_NIPBL_rep2
JOB_DEPENDENCIES=$picard_mark_duplicates_3_JOB_ID
JOB_DONE=job_output/homer_make_tag_directory/homer_make_tag_directory.A549_DEX_NIPBL_rep2.2e95feb80984cd7bc8865d02de483aae.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_tag_directory.A549_DEX_NIPBL_rep2.2e95feb80984cd7bc8865d02de483aae.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 mugqic/samtools/1.3.1 && \
makeTagDirectory tags/A549_DEX_NIPBL_rep2 \
            alignment/A549_DEX_NIPBL_rep2/A549_DEX_NIPBL_rep2.sorted.dup.bam \
            -genome hg38 \
            -checkGC \
 
homer_make_tag_directory.A549_DEX_NIPBL_rep2.2e95feb80984cd7bc8865d02de483aae.mugqic.done
)
homer_make_tag_directory_3_JOB_ID=$(echo "#! /bin/bash 
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
echo "$homer_make_tag_directory_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# STEP: qc_metrics
#-------------------------------------------------------------------------------
STEP=qc_metrics
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: qc_metrics_1_JOB_ID: qc_plots_R
#-------------------------------------------------------------------------------
JOB_NAME=qc_plots_R
JOB_DEPENDENCIES=$homer_make_tag_directory_1_JOB_ID:$homer_make_tag_directory_2_JOB_ID:$homer_make_tag_directory_3_JOB_ID
JOB_DONE=job_output/qc_metrics/qc_plots_R.6454129ca20a7ad00c28aee8e77889c3.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'qc_plots_R.6454129ca20a7ad00c28aee8e77889c3.mugqic.done'
module load mugqic/mugqic_tools/2.1.9 mugqic/R_Bioconductor/3.2.3_3.2 && \
mkdir -p graphs && \
Rscript $R_TOOLS/chipSeqGenerateQCMetrics.R \
  ../../raw/chip-seq/design_partial_CHUL.txt \
  /project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38 && \
cp /home/efournie/genpipes/bfx/report/ChipSeq.qc_metrics.md report/ChipSeq.qc_metrics.md && \
for sample in A549_DEX_MED1_rep2 A549_DEX_SMC1A_rep2 A549_DEX_NIPBL_rep2
do
  cp --parents graphs/${sample}_QC_Metrics.ps report/
  convert -rotate 90 graphs/${sample}_QC_Metrics.ps report/graphs/${sample}_QC_Metrics.png
  echo -e "----

![QC Metrics for Sample $sample ([download high-res image](graphs/${sample}_QC_Metrics.ps))](graphs/${sample}_QC_Metrics.png)
" \
  >> report/ChipSeq.qc_metrics.md
done
qc_plots_R.6454129ca20a7ad00c28aee8e77889c3.mugqic.done
)
qc_metrics_1_JOB_ID=$(echo "#! /bin/bash 
echo '#######################################'
echo 'SLURM FAKE PROLOGUE (MUGQIC)'
date 
scontrol show job \$SLURM_JOBID
sstat -j \$SLURM_JOBID.batch 
echo '#######################################'
rm -f $JOB_DONE && module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"qc_metrics\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_MED1_rep2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_SMC1A_rep2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_NIPBL_rep2.json\" \
  -f \"running\"
module unload mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9 && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
module load mugqic/python/2.7.13 mugqic/mugqic_tools/2.1.9
/home/efournie/genpipes/utils/job2json.py \
  -u \"efournie\" \
  -c \"/home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini,/home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini,/cvmfs/soft.mugqic/CentOS6/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini\" \
  -s \"qc_metrics\" \
  -j \"$JOB_NAME\" \
  -d \"$JOB_DONE\" \
  -l \"$JOB_OUTPUT\" \
  -o \"/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_MED1_rep2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_SMC1A_rep2.json,/project/6001942/Working_Directory/Eric/CofactorHR/A549/output/chip-pipeline-GRCh38/json/A549_DEX_NIPBL_rep2.json\" \
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
sbatch --mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID -D $OUTPUT_DIR -o $JOB_OUTPUT -J $JOB_NAME --time=24:00:00 --mem=4G -n 1 -N 1 --depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]" | cut -d\  -f4)
echo "$qc_metrics_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# STEP: homer_make_ucsc_file
#-------------------------------------------------------------------------------
STEP=homer_make_ucsc_file
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_1_JOB_ID: homer_make_ucsc_file.A549_DEX_MED1_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.A549_DEX_MED1_rep2
JOB_DEPENDENCIES=$homer_make_tag_directory_1_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.A549_DEX_MED1_rep2.7e44441a9b28bd92b739b33de9cb5763.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.A549_DEX_MED1_rep2.7e44441a9b28bd92b739b33de9cb5763.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 && \
mkdir -p tracks/A549_DEX_MED1_rep2 && \
makeUCSCfile \
        tags/A549_DEX_MED1_rep2 > tracks/A549_DEX_MED1_rep2/A549_DEX_MED1_rep2.ucsc.bedGraph && \
        gzip -c tracks/A549_DEX_MED1_rep2/A549_DEX_MED1_rep2.ucsc.bedGraph > tracks/A549_DEX_MED1_rep2/A549_DEX_MED1_rep2.ucsc.bedGraph.gz
homer_make_ucsc_file.A549_DEX_MED1_rep2.7e44441a9b28bd92b739b33de9cb5763.mugqic.done
)
homer_make_ucsc_file_1_JOB_ID=$(echo "#! /bin/bash 
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
echo "$homer_make_ucsc_file_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_2_JOB_ID: homer_make_ucsc_file_bigWig.A549_DEX_MED1_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file_bigWig.A549_DEX_MED1_rep2
JOB_DEPENDENCIES=$homer_make_ucsc_file_1_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file_bigWig.A549_DEX_MED1_rep2.04e5de4b13b341c4a12defa83d13870e.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file_bigWig.A549_DEX_MED1_rep2.04e5de4b13b341c4a12defa83d13870e.mugqic.done'
module load mugqic/ucsc/v346 && \
mkdir -p tracks/A549_DEX_MED1_rep2/bigWig && \
export TMPDIR=${SLURM_TMPDIR} && \
cat tracks/A549_DEX_MED1_rep2/A549_DEX_MED1_rep2.ucsc.bedGraph | head -n 1 > tracks/A549_DEX_MED1_rep2/A549_DEX_MED1_rep2.ucsc.bedGraph.head.tmp && \
cat tracks/A549_DEX_MED1_rep2/A549_DEX_MED1_rep2.ucsc.bedGraph | awk ' NR > 1 ' | sort  --temporary-directory=${SLURM_TMPDIR} -k1,1 -k2,2n | \
awk '{if($0 !~ /^[A-W]/) print ""$0; else print $0}' | grep -vP "GL|lambda|pUC19" | sed 's/MT/chrM/' | \
awk '{printf "%s\t%d\t%d\t%4.4g\n", $1,$2,$3,$4}' > tracks/A549_DEX_MED1_rep2/A549_DEX_MED1_rep2.ucsc.bedGraph.body.tmp && \
cat tracks/A549_DEX_MED1_rep2/A549_DEX_MED1_rep2.ucsc.bedGraph.head.tmp tracks/A549_DEX_MED1_rep2/A549_DEX_MED1_rep2.ucsc.bedGraph.body.tmp > tracks/A549_DEX_MED1_rep2/A549_DEX_MED1_rep2.ucsc.bedGraph.sorted && \
rm tracks/A549_DEX_MED1_rep2/A549_DEX_MED1_rep2.ucsc.bedGraph.head.tmp tracks/A549_DEX_MED1_rep2/A549_DEX_MED1_rep2.ucsc.bedGraph.body.tmp && \
bedGraphToBigWig \
  tracks/A549_DEX_MED1_rep2/A549_DEX_MED1_rep2.ucsc.bedGraph.sorted \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  tracks/A549_DEX_MED1_rep2/bigWig/A549_DEX_MED1_rep2.bw
homer_make_ucsc_file_bigWig.A549_DEX_MED1_rep2.04e5de4b13b341c4a12defa83d13870e.mugqic.done
)
homer_make_ucsc_file_2_JOB_ID=$(echo "#! /bin/bash 
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
echo "$homer_make_ucsc_file_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_3_JOB_ID: homer_make_ucsc_file.A549_DEX_SMC1A_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.A549_DEX_SMC1A_rep2
JOB_DEPENDENCIES=$homer_make_tag_directory_2_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.A549_DEX_SMC1A_rep2.139bc21531540ee3766866406c647b14.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.A549_DEX_SMC1A_rep2.139bc21531540ee3766866406c647b14.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 && \
mkdir -p tracks/A549_DEX_SMC1A_rep2 && \
makeUCSCfile \
        tags/A549_DEX_SMC1A_rep2 > tracks/A549_DEX_SMC1A_rep2/A549_DEX_SMC1A_rep2.ucsc.bedGraph && \
        gzip -c tracks/A549_DEX_SMC1A_rep2/A549_DEX_SMC1A_rep2.ucsc.bedGraph > tracks/A549_DEX_SMC1A_rep2/A549_DEX_SMC1A_rep2.ucsc.bedGraph.gz
homer_make_ucsc_file.A549_DEX_SMC1A_rep2.139bc21531540ee3766866406c647b14.mugqic.done
)
homer_make_ucsc_file_3_JOB_ID=$(echo "#! /bin/bash 
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
echo "$homer_make_ucsc_file_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_4_JOB_ID: homer_make_ucsc_file_bigWig.A549_DEX_SMC1A_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file_bigWig.A549_DEX_SMC1A_rep2
JOB_DEPENDENCIES=$homer_make_ucsc_file_3_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file_bigWig.A549_DEX_SMC1A_rep2.b2bf8241ac054928eb0878ac99e7494a.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file_bigWig.A549_DEX_SMC1A_rep2.b2bf8241ac054928eb0878ac99e7494a.mugqic.done'
module load mugqic/ucsc/v346 && \
mkdir -p tracks/A549_DEX_SMC1A_rep2/bigWig && \
export TMPDIR=${SLURM_TMPDIR} && \
cat tracks/A549_DEX_SMC1A_rep2/A549_DEX_SMC1A_rep2.ucsc.bedGraph | head -n 1 > tracks/A549_DEX_SMC1A_rep2/A549_DEX_SMC1A_rep2.ucsc.bedGraph.head.tmp && \
cat tracks/A549_DEX_SMC1A_rep2/A549_DEX_SMC1A_rep2.ucsc.bedGraph | awk ' NR > 1 ' | sort  --temporary-directory=${SLURM_TMPDIR} -k1,1 -k2,2n | \
awk '{if($0 !~ /^[A-W]/) print ""$0; else print $0}' | grep -vP "GL|lambda|pUC19" | sed 's/MT/chrM/' | \
awk '{printf "%s\t%d\t%d\t%4.4g\n", $1,$2,$3,$4}' > tracks/A549_DEX_SMC1A_rep2/A549_DEX_SMC1A_rep2.ucsc.bedGraph.body.tmp && \
cat tracks/A549_DEX_SMC1A_rep2/A549_DEX_SMC1A_rep2.ucsc.bedGraph.head.tmp tracks/A549_DEX_SMC1A_rep2/A549_DEX_SMC1A_rep2.ucsc.bedGraph.body.tmp > tracks/A549_DEX_SMC1A_rep2/A549_DEX_SMC1A_rep2.ucsc.bedGraph.sorted && \
rm tracks/A549_DEX_SMC1A_rep2/A549_DEX_SMC1A_rep2.ucsc.bedGraph.head.tmp tracks/A549_DEX_SMC1A_rep2/A549_DEX_SMC1A_rep2.ucsc.bedGraph.body.tmp && \
bedGraphToBigWig \
  tracks/A549_DEX_SMC1A_rep2/A549_DEX_SMC1A_rep2.ucsc.bedGraph.sorted \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  tracks/A549_DEX_SMC1A_rep2/bigWig/A549_DEX_SMC1A_rep2.bw
homer_make_ucsc_file_bigWig.A549_DEX_SMC1A_rep2.b2bf8241ac054928eb0878ac99e7494a.mugqic.done
)
homer_make_ucsc_file_4_JOB_ID=$(echo "#! /bin/bash 
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
echo "$homer_make_ucsc_file_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_5_JOB_ID: homer_make_ucsc_file.A549_DEX_NIPBL_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file.A549_DEX_NIPBL_rep2
JOB_DEPENDENCIES=$homer_make_tag_directory_3_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file.A549_DEX_NIPBL_rep2.b9f901931424a3bbf445e6531ebca056.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file.A549_DEX_NIPBL_rep2.b9f901931424a3bbf445e6531ebca056.mugqic.done'
module load mugqic/perl/5.22.1 mugqic/homer/4.9.1 && \
mkdir -p tracks/A549_DEX_NIPBL_rep2 && \
makeUCSCfile \
        tags/A549_DEX_NIPBL_rep2 > tracks/A549_DEX_NIPBL_rep2/A549_DEX_NIPBL_rep2.ucsc.bedGraph && \
        gzip -c tracks/A549_DEX_NIPBL_rep2/A549_DEX_NIPBL_rep2.ucsc.bedGraph > tracks/A549_DEX_NIPBL_rep2/A549_DEX_NIPBL_rep2.ucsc.bedGraph.gz
homer_make_ucsc_file.A549_DEX_NIPBL_rep2.b9f901931424a3bbf445e6531ebca056.mugqic.done
)
homer_make_ucsc_file_5_JOB_ID=$(echo "#! /bin/bash 
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
echo "$homer_make_ucsc_file_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_6_JOB_ID: homer_make_ucsc_file_bigWig.A549_DEX_NIPBL_rep2
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file_bigWig.A549_DEX_NIPBL_rep2
JOB_DEPENDENCIES=$homer_make_ucsc_file_5_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file_bigWig.A549_DEX_NIPBL_rep2.f4df00be3e031ab54ee37dbb0cd2dc94.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file_bigWig.A549_DEX_NIPBL_rep2.f4df00be3e031ab54ee37dbb0cd2dc94.mugqic.done'
module load mugqic/ucsc/v346 && \
mkdir -p tracks/A549_DEX_NIPBL_rep2/bigWig && \
export TMPDIR=${SLURM_TMPDIR} && \
cat tracks/A549_DEX_NIPBL_rep2/A549_DEX_NIPBL_rep2.ucsc.bedGraph | head -n 1 > tracks/A549_DEX_NIPBL_rep2/A549_DEX_NIPBL_rep2.ucsc.bedGraph.head.tmp && \
cat tracks/A549_DEX_NIPBL_rep2/A549_DEX_NIPBL_rep2.ucsc.bedGraph | awk ' NR > 1 ' | sort  --temporary-directory=${SLURM_TMPDIR} -k1,1 -k2,2n | \
awk '{if($0 !~ /^[A-W]/) print ""$0; else print $0}' | grep -vP "GL|lambda|pUC19" | sed 's/MT/chrM/' | \
awk '{printf "%s\t%d\t%d\t%4.4g\n", $1,$2,$3,$4}' > tracks/A549_DEX_NIPBL_rep2/A549_DEX_NIPBL_rep2.ucsc.bedGraph.body.tmp && \
cat tracks/A549_DEX_NIPBL_rep2/A549_DEX_NIPBL_rep2.ucsc.bedGraph.head.tmp tracks/A549_DEX_NIPBL_rep2/A549_DEX_NIPBL_rep2.ucsc.bedGraph.body.tmp > tracks/A549_DEX_NIPBL_rep2/A549_DEX_NIPBL_rep2.ucsc.bedGraph.sorted && \
rm tracks/A549_DEX_NIPBL_rep2/A549_DEX_NIPBL_rep2.ucsc.bedGraph.head.tmp tracks/A549_DEX_NIPBL_rep2/A549_DEX_NIPBL_rep2.ucsc.bedGraph.body.tmp && \
bedGraphToBigWig \
  tracks/A549_DEX_NIPBL_rep2/A549_DEX_NIPBL_rep2.ucsc.bedGraph.sorted \
  /project/6007512/C3G/analyste_dev/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai \
  tracks/A549_DEX_NIPBL_rep2/bigWig/A549_DEX_NIPBL_rep2.bw
homer_make_ucsc_file_bigWig.A549_DEX_NIPBL_rep2.f4df00be3e031ab54ee37dbb0cd2dc94.mugqic.done
)
homer_make_ucsc_file_6_JOB_ID=$(echo "#! /bin/bash 
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
echo "$homer_make_ucsc_file_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# JOB: homer_make_ucsc_file_7_JOB_ID: homer_make_ucsc_file_report
#-------------------------------------------------------------------------------
JOB_NAME=homer_make_ucsc_file_report
JOB_DEPENDENCIES=$homer_make_ucsc_file_5_JOB_ID
JOB_DONE=job_output/homer_make_ucsc_file/homer_make_ucsc_file_report.90b1cb89146cadaf805ea66244074321.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'homer_make_ucsc_file_report.90b1cb89146cadaf805ea66244074321.mugqic.done'
mkdir -p report && \
zip -r report/tracks.zip tracks/*/*.ucsc.bedGraph.gz && \
cp /home/efournie/genpipes/bfx/report/ChipSeq.homer_make_ucsc_file.md report/
homer_make_ucsc_file_report.90b1cb89146cadaf805ea66244074321.mugqic.done
)
homer_make_ucsc_file_7_JOB_ID=$(echo "#! /bin/bash 
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
echo "$homer_make_ucsc_file_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST

sleep 0.5


#-------------------------------------------------------------------------------
# Call home with pipeline statistics
#-------------------------------------------------------------------------------
wget "http://mugqic.hpc.mcgill.ca/cgi-bin/pipeline.cgi?hostname=cedar5.cedar.computecanada.ca&ip=206.12.124.6&pipeline=ChipSeq&steps=picard_sam_to_fastq,trimmomatic,merge_trimmomatic_stats,bwa_mem_picard_sort_sam,samtools_view_filter,picard_merge_sam_files,picard_mark_duplicates,metrics,homer_make_tag_directory,qc_metrics,homer_make_ucsc_file&samples=3&AnonymizedList=d8a59c316384745ed7db9dbeb5f1e489,11728c95eb3cea1f33a00588ac772178,ce51d19e1027f42536c485b28dc74e8a" --quiet --output-document=/dev/null

