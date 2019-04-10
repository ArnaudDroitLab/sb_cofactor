mkdir -p output/chip-pipeline-GRCm38_GSE94300_B

chipseq.py -j slurm -s '1-10,12' \
    -l debug \
    -r raw/chip-seq/readset_GSE94300_B.txt \
    -d raw/chip-seq/design_GSE94300_B.txt \
    -o output/chip-pipeline-GRCm38_GSE94300_B \
    --config $MUGQIC_PIPELINES_HOME/pipelines/chipseq/chipseq.base.ini \
        $MUGQIC_PIPELINES_HOME/pipelines/chipseq/chipseq.cedar.ini \
        $MUGQIC_INSTALL_HOME/genomes/species/Mus_musculus.GRCm38/Mus_musculus.GRCm38.ini \
        input/cedar_assembly_bugfix.txt \
        input/cedar_T30M30.txt \
        input/cedar_BWA128G.txt
