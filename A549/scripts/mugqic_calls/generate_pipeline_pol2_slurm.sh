mkdir -p output/chip-pipeline-PolII-GRCh38

/home/efournie/genpipes/pipelines/chipseq/chipseq.py -j slurm -s '12' \
    -l debug \
    -r raw/readset.txt \
    -d raw/design.txt \
    -o output/chip-pipeline-PolII-GRCh38 \
    --config /home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini \
        /home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini \
        $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini \
        input/chipseq.numpy.bug.ini
