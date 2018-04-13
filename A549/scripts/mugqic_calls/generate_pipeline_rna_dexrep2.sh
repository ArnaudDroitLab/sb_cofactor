export RAP_ID=def-stbil30

mkdir -p output/rna-pipeline-GRCh38
/home/efournie/genpipes/pipelines/rnaseq/rnaseq.py -j slurm -s '1-7' -l debug \
    -r raw/rna-seq/readset_dexrep2.txt \
    -d raw/rna-seq/design_dexrep2.txt \
    -o output/rna-pipeline-GRCh38 \
    --config /home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini \
        /home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini \
        $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini 