mkdir -p output/chip-pipeline

$MUGQIC_PIPELINES_HOME/pipelines/chipseq/chipseq.py -s '1-15' \
    -l debug \
    -r raw/chip-seq/readset.txt \
    -d raw/chip-seq/design.txt \
    -o output/chip-pipeline \
    --config $MUGQIC_PIPELINES_HOME/pipelines/chipseq/chipseq.base.ini \
        $MUGQIC_PIPELINES_HOME/pipelines/chipseq/chipseq.guillimin.ini \
        $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.hg19/Homo_sapiens.hg19.ini \
        input/chipseq.numpy.bug.ini
