mkdir -p output/chip-pipeline-GRCh38

/home/efournie/genpipes/pipelines/chipseq/chipseq.py -j slurm -s '1-7' \
    -l debug \
    -r raw/chip-seq/readset_SMC1A_Rep1.txt \
    -d raw/chip-seq/design_SMC1A_Rep1.txt \
    -o output/chip-pipeline-GRCh38 \
    --config /home/efournie/genpipes/pipelines/chipseq/chipseq.base.ini \
        /home/efournie/genpipes/pipelines/chipseq/chipseq.cedar.ini \
        $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini
