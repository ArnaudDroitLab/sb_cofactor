mkdir -p output/chip-pipeline-GRCm38_GSE20530

chipseq.py -j slurm -s '1-10,12' \
    -l debug \
    -r raw/chip-seq/readset_GSE20530.txt \
    -d raw/chip-seq/design_GSE20530.txt \
    -o output/chip-pipeline-GRCm38_GSE20530 \
    --config $MUGQIC_PIPELINES_HOME/pipelines/chipseq/chipseq.base.ini \
        $MUGQIC_PIPELINES_HOME/pipelines/chipseq/chipseq.cedar.ini \
        $MUGQIC_INSTALL_HOME/genomes/species/Mus_musculus.GRCm38/Mus_musculus.GRCm38.ini \
        input/cedar_assembly_bugfix.txt \
        input/cedar_T25M25.txt