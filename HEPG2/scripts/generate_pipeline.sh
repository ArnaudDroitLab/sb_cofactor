mkdir -p output/chip-pipeline-GRCh38

chipseq.py -j slurm -s '1-10,12' \
    -l debug \
    -r raw/chip-seq/readset.txt \
    -d raw/chip-seq/design.txt \
    -o output/chip-pipeline-GRCh38 \
    --config $MUGQIC_PIPELINES_HOME/pipelines/chipseq/chipseq.base.ini \
        $MUGQIC_PIPELINES_HOME/pipelines/chipseq/chipseq.cedar.ini \
        $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini \
        input/cedar_assembly_bugfix.txt

