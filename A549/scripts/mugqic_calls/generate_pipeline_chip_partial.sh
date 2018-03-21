mkdir -p output/chip-pipeline-GRCh38

$MUGQIC_PIPELINES_HOME/pipelines/chipseq/chipseq.py -s '1-15' \
    -l debug \
    -r raw/chip-seq/readset_partial.txt \
    -d raw/chip-seq/design_partial.txt \
    -o output/chip-pipeline-GRCh38 \
    --config $MUGQIC_PIPELINES_HOME/pipelines/chipseq/chipseq.base.ini \
        $MUGQIC_PIPELINES_HOME/pipelines/chipseq/chipseq.guillimin.ini \
        $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini \
        input/chipseq.numpy.bug.ini
