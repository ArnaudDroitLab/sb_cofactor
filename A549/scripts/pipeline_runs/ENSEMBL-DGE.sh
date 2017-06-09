OUTPUT_DIR=/gs/scratch/efournier/CofactorHR/A549/output/rna-pipeline
JOB_OUTPUT_DIR=$OUTPUT_DIR/job_output
TIMESTAMP=`date +%FT%H.%M.%S`
JOB_LIST=$JOB_OUTPUT_DIR/RnaSeq_job_list_$TIMESTAMP
mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR


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
JOB_DONE=job_output/differential_expression/differential_expression.471a56e71dafee5bb143cd0a2091d23c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'differential_expression.471a56e71dafee5bb143cd0a2091d23c.mugqic.done'
module load mugqic/mugqic_tools/2.1.7 mugqic/R_Bioconductor/3.1.2_3.0 && \
mkdir -p DGE_ENSEMBL && \
Rscript $R_TOOLS/edger.R \
  -d ../../raw/rna-seq/design.txt \
  -c DGE_ENSEMBL/rawCountMatrix.csv \
  -o DGE_ENSEMBL && \
Rscript $R_TOOLS/deseq.R \
  -d ../../raw/rna-seq/design.txt \
  -c DGE_ENSEMBL/rawCountMatrix.csv \
  -o DGE_ENSEMBL \
  -l
differential_expression.471a56e71dafee5bb143cd0a2091d23c.mugqic.done
)
differential_expression_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -M $JOB_MAIL -W umask=0002 -A $RAP_ID -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=10:00:0 -q metaq -l nodes=1:ppn=1 -l pmem=1700m  | grep "[0-9]")
