mkdir -p output/chip-pipeline-ENCODE

$MUGQIC_PIPELINES_HOME/pipelines/chipseq/chipseq.py -s '1-15' \
    -l debug \
    -r raw/ENCODE-chip/readset.txt \
    -d raw/ENCODE-chip/design.txt \
    -o output/ENCODE-chip \
    --config $MUGQIC_PIPELINES_HOME/pipelines/chipseq/chipseq.base.ini \
        $MUGQIC_PIPELINES_HOME/pipelines/chipseq/chipseq.guillimin.ini \
        $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini \
        input/chipseq.numpy.bug.ini