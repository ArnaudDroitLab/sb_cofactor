
OUTPUT_DIR=/gs/scratch/efournier/CofactorHR/A549/output/rna-pipeline
JOB_OUTPUT_DIR=$OUTPUT_DIR/job_output
TIMESTAMP=`date +%FT%H.%M.%S`
JOB_LIST=$JOB_OUTPUT_DIR/RnaSeq_job_list_$TIMESTAMP
mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR


#-------------------------------------------------------------------------------
# STEP: raw_counts
#-------------------------------------------------------------------------------
STEP=raw_counts_ENSEMBL
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: raw_counts_1_JOB_ID: htseq_count.DEX_nonMamm_1
#-------------------------------------------------------------------------------
JOB_NAME=htseq_count.DEX_nonMamm_1
JOB_DEPENDENCIES=$picard_sort_sam_1_JOB_ID
JOB_DONE=job_output/raw_counts_ENSEMBL/htseq_count.DEX_nonMamm_1.cda09854fa5c0036c14132a63939f746.mugqic.done
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
  ../../ENSEMBL_UCUSC.gtf \
  > raw_counts_ENSEMBL/DEX_nonMamm_1.readcounts.csv
htseq_count.DEX_nonMamm_1.cda09854fa5c0036c14132a63939f746.mugqic.done
)
raw_counts_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=1700m  | grep "[0-9]")
echo "$raw_counts_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: raw_counts_2_JOB_ID: htseq_count.DEX_shMED1_1
#-------------------------------------------------------------------------------
JOB_NAME=htseq_count.DEX_shMED1_1
JOB_DEPENDENCIES=$picard_sort_sam_2_JOB_ID
JOB_DONE=job_output/raw_counts_ENSEMBL/htseq_count.DEX_shMED1_1.20649cac2eaf3c487e3806d360495b6c.mugqic.done
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
  ../../ENSEMBL_UCUSC.gtf \
  > raw_counts_ENSEMBL/DEX_shMED1_1.readcounts.csv
htseq_count.DEX_shMED1_1.20649cac2eaf3c487e3806d360495b6c.mugqic.done
)
raw_counts_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=1700m  | grep "[0-9]")
echo "$raw_counts_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: raw_counts_3_JOB_ID: htseq_count.DEX_shNIPBL_1
#-------------------------------------------------------------------------------
JOB_NAME=htseq_count.DEX_shNIPBL_1
JOB_DEPENDENCIES=$picard_sort_sam_3_JOB_ID
JOB_DONE=job_output/raw_counts_ENSEMBL/htseq_count.DEX_shNIPBL_1.b4c69e6af7975293fb13169f1f651f5d.mugqic.done
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
  ../../ENSEMBL_UCUSC.gtf \
  > raw_counts_ENSEMBL/DEX_shNIPBL_1.readcounts.csv
htseq_count.DEX_shNIPBL_1.b4c69e6af7975293fb13169f1f651f5d.mugqic.done
)
raw_counts_3_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=1700m  | grep "[0-9]")
echo "$raw_counts_3_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: raw_counts_4_JOB_ID: htseq_count.DEX_shSMC1A_3
#-------------------------------------------------------------------------------
JOB_NAME=htseq_count.DEX_shSMC1A_3
JOB_DEPENDENCIES=$picard_sort_sam_4_JOB_ID
JOB_DONE=job_output/raw_counts_ENSEMBL/htseq_count.DEX_shSMC1A_3.d788f23a20e7cb8aaed60833ebd91556.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'htseq_count.DEX_shSMC1A_3.d788f23a20e7cb8aaed60833ebd91556.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/python/2.7.13 && \
mkdir -p raw_counts && \
samtools view -F 4 \
  alignment/DEX_shSMC1A_3/DEX_shSMC1A_3.QueryNameSorted.bam | \
htseq-count - \
  -m intersection-nonempty \
  --stranded=reverse \
  --format=sam \
  ../../ENSEMBL_UCUSC.gtf \
  > raw_counts_ENSEMBL/DEX_shSMC1A_3.readcounts.csv
htseq_count.DEX_shSMC1A_3.d788f23a20e7cb8aaed60833ebd91556.mugqic.done
)
raw_counts_4_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=1700m  | grep "[0-9]")
echo "$raw_counts_4_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: raw_counts_5_JOB_ID: htseq_count.ETOH_nonMamm_1
#-------------------------------------------------------------------------------
JOB_NAME=htseq_count.ETOH_nonMamm_1
JOB_DEPENDENCIES=$picard_sort_sam_5_JOB_ID
JOB_DONE=job_output/raw_counts_ENSEMBL/htseq_count.ETOH_nonMamm_1.c002606d033d98f723d9720fcb62c27d.mugqic.done
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
  ../../ENSEMBL_UCUSC.gtf \
  > raw_counts_ENSEMBL/ETOH_nonMamm_1.readcounts.csv
htseq_count.ETOH_nonMamm_1.c002606d033d98f723d9720fcb62c27d.mugqic.done
)
raw_counts_5_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=1700m  | grep "[0-9]")
echo "$raw_counts_5_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: raw_counts_6_JOB_ID: htseq_count.ETOH_nonMamm_2
#-------------------------------------------------------------------------------
JOB_NAME=htseq_count.ETOH_nonMamm_2
JOB_DEPENDENCIES=$picard_sort_sam_6_JOB_ID
JOB_DONE=job_output/raw_counts_ENSEMBL/htseq_count.ETOH_nonMamm_2.2f4b039f189ead7f2c504bc5b3c82267.mugqic.done
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
  ../../ENSEMBL_UCUSC.gtf \
  > raw_counts_ENSEMBL/ETOH_nonMamm_2.readcounts.csv
htseq_count.ETOH_nonMamm_2.2f4b039f189ead7f2c504bc5b3c82267.mugqic.done
)
raw_counts_6_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=1700m  | grep "[0-9]")
echo "$raw_counts_6_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: raw_counts_7_JOB_ID: htseq_count.ETOH_shMED1_1
#-------------------------------------------------------------------------------
JOB_NAME=htseq_count.ETOH_shMED1_1
JOB_DEPENDENCIES=$picard_sort_sam_7_JOB_ID
JOB_DONE=job_output/raw_counts_ENSEMBL/htseq_count.ETOH_shMED1_1.1805706205cc6eefcf397f33fc39af06.mugqic.done
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
  ../../ENSEMBL_UCUSC.gtf \
  > raw_counts_ENSEMBL/ETOH_shMED1_1.readcounts.csv
htseq_count.ETOH_shMED1_1.1805706205cc6eefcf397f33fc39af06.mugqic.done
)
raw_counts_7_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=1700m  | grep "[0-9]")
echo "$raw_counts_7_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: raw_counts_8_JOB_ID: htseq_count.ETOH_shNIPBL_2
#-------------------------------------------------------------------------------
JOB_NAME=htseq_count.ETOH_shNIPBL_2
JOB_DEPENDENCIES=$picard_sort_sam_8_JOB_ID
JOB_DONE=job_output/raw_counts_ENSEMBL/htseq_count.ETOH_shNIPBL_2.64d8757ddc896dac291d14f05c36fb4f.mugqic.done
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
  ../../ENSEMBL_UCUSC.gtf \
  > raw_counts_ENSEMBL/ETOH_shNIPBL_2.readcounts.csv
htseq_count.ETOH_shNIPBL_2.64d8757ddc896dac291d14f05c36fb4f.mugqic.done
)
raw_counts_8_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=1700m  | grep "[0-9]")
echo "$raw_counts_8_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: raw_counts_9_JOB_ID: htseq_count.ETOH_shNIPBL_1
#-------------------------------------------------------------------------------
JOB_NAME=htseq_count.ETOH_shNIPBL_1
JOB_DEPENDENCIES=$picard_sort_sam_9_JOB_ID
JOB_DONE=job_output/raw_counts_ENSEMBL/htseq_count.ETOH_shNIPBL_1.ec5717e8bfb63aab5a8b0a3885573b27.mugqic.done
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
  ../../ENSEMBL_UCUSC.gtf \
  > raw_counts_ENSEMBL/ETOH_shNIPBL_1.readcounts.csv
htseq_count.ETOH_shNIPBL_1.ec5717e8bfb63aab5a8b0a3885573b27.mugqic.done
)
raw_counts_9_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=1700m  | grep "[0-9]")
echo "$raw_counts_9_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: raw_counts_10_JOB_ID: htseq_count.ETOH_shMED1_2
#-------------------------------------------------------------------------------
JOB_NAME=htseq_count.ETOH_shMED1_2
JOB_DEPENDENCIES=$picard_sort_sam_10_JOB_ID
JOB_DONE=job_output/raw_counts_ENSEMBL/htseq_count.ETOH_shMED1_2.728a1efbfc6b094f685e20206001f2f8.mugqic.done
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
  ../../ENSEMBL_UCUSC.gtf \
  > raw_counts_ENSEMBL/ETOH_shMED1_2.readcounts.csv
htseq_count.ETOH_shMED1_2.728a1efbfc6b094f685e20206001f2f8.mugqic.done
)
raw_counts_10_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=1700m  | grep "[0-9]")
echo "$raw_counts_10_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: raw_counts_11_JOB_ID: htseq_count.ETOH_shSMC1A_3
#-------------------------------------------------------------------------------
JOB_NAME=htseq_count.ETOH_shSMC1A_3
JOB_DEPENDENCIES=$picard_sort_sam_11_JOB_ID
JOB_DONE=job_output/raw_counts_ENSEMBL/htseq_count.ETOH_shSMC1A_3.148fb66ebed515441801e53b9f0330e7.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'htseq_count.ETOH_shSMC1A_3.148fb66ebed515441801e53b9f0330e7.mugqic.done'
module load mugqic/samtools/1.3.1 mugqic/python/2.7.13 && \
mkdir -p raw_counts && \
samtools view -F 4 \
  alignment/ETOH_shSMC1A_3/ETOH_shSMC1A_3.QueryNameSorted.bam | \
htseq-count - \
  -m intersection-nonempty \
  --stranded=reverse \
  --format=sam \
  ../../ENSEMBL_UCUSC.gtf \
  > raw_counts_ENSEMBL/ETOH_shSMC1A_3.readcounts.csv
htseq_count.ETOH_shSMC1A_3.148fb66ebed515441801e53b9f0330e7.mugqic.done
)
raw_counts_11_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=1700m  | grep "[0-9]")
echo "$raw_counts_11_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: raw_counts_12_JOB_ID: htseq_count.ETOH_shSMC1A_2
#-------------------------------------------------------------------------------
JOB_NAME=htseq_count.ETOH_shSMC1A_2
JOB_DEPENDENCIES=$picard_sort_sam_12_JOB_ID
JOB_DONE=job_output/raw_counts_ENSEMBL/htseq_count.ETOH_shSMC1A_2.bd31be1e8d841acb792d73fc68d648e3.mugqic.done
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
  ../../ENSEMBL_UCUSC.gtf \
  > raw_counts_ENSEMBL/ETOH_shSMC1A_2.readcounts.csv
htseq_count.ETOH_shSMC1A_2.bd31be1e8d841acb792d73fc68d648e3.mugqic.done
)
raw_counts_12_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=24:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=1700m  | grep "[0-9]")
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
JOB_DONE=job_output/raw_counts_metrics/metrics.matrix.73a19052bb7cbcc74617e6452a1f64d7.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'metrics.matrix.73a19052bb7cbcc74617e6452a1f64d7.mugqic.done'
module load mugqic/mugqic_tools/2.1.7 && \
mkdir -p DGE && \
gtf2tmpMatrix.awk \
  ../../ENSEMBL_UCUSC.gtf \
  DGE/tmpMatrix.txt && \
HEAD='Gene\tSymbol' && \
for read_count_file in \
  raw_counts_ENSEMBL/DEX_nonMamm_1.readcounts.csv \
  raw_counts_ENSEMBL/DEX_shMED1_1.readcounts.csv \
  raw_counts_ENSEMBL/DEX_shNIPBL_1.readcounts.csv \
  raw_counts_ENSEMBL/DEX_shSMC1A_3.readcounts.csv \
  raw_counts_ENSEMBL/ETOH_nonMamm_1.readcounts.csv \
  raw_counts_ENSEMBL/ETOH_nonMamm_2.readcounts.csv \
  raw_counts_ENSEMBL/ETOH_shMED1_1.readcounts.csv \
  raw_counts_ENSEMBL/ETOH_shNIPBL_2.readcounts.csv \
  raw_counts_ENSEMBL/ETOH_shNIPBL_1.readcounts.csv \
  raw_counts_ENSEMBL/ETOH_shMED1_2.readcounts.csv \
  raw_counts_ENSEMBL/ETOH_shSMC1A_3.readcounts.csv \
  raw_counts_ENSEMBL/ETOH_shSMC1A_2.readcounts.csv
do
  sort -k1,1 $read_count_file > DGE/tmpSort.txt && \
  join -1 1 -2 1 <(sort -k1,1 DGE/tmpMatrix.txt) DGE/tmpSort.txt > DGE/tmpMatrix.2.txt && \
  mv DGE/tmpMatrix.2.txt DGE/tmpMatrix.txt && \
  na=$(basename $read_count_file | rev | cut -d. -f3- | rev) && \
  HEAD="$HEAD\t$na"
done && \
echo -e $HEAD | cat - DGE/tmpMatrix.txt | tr ' ' '\t' > DGE/rawCountMatrix.csv && \
rm DGE/tmpSort.txt DGE/tmpMatrix.txt
metrics.matrix.73a19052bb7cbcc74617e6452a1f64d7.mugqic.done
)
raw_counts_metrics_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=5:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=1700m -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
echo "$raw_counts_metrics_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST
