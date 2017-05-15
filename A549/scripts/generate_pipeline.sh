export RAP_ID=fhq-091-aa
module load mugqic/mugqic_pipelines/2.3.0

PROJECT_BASE=/gs/scratch/efournier/CofactorHR/A549

mkdir -p $PROJECT_BASE/output/rna-pipeline
$MUGQIC_PIPELINES_HOME/pipelines/rnaseq/rnaseq.py -s '1-22' -l debug \
    -r $PROJECT_BASE/raw/rna-seq/readset.txt \
    -d $PROJECT_BASE/raw/rna-seq/design.txt \
    -o $PROJECT_BASE/output/rna-pipeline \
    --config $MUGQIC_PIPELINES_HOME/pipelines/rnaseq/rnaseq.base.ini \
        $MUGQIC_PIPELINES_HOME/pipelines/rnaseq/rnaseq.guillimin.ini \
        $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/Homo_sapiens.hg19.ini \
        $PROJECT_BASE/input/this_run.ini