export RAP_ID=def-stbil30
module load mugqic/mugqic_pipelines


mkdir -p output/rna-pipeline-GRCh38
/home/efournie/genpipes/pipelines/rnaseq/rnaseq.py -j slurm -s '10' -l debug \
    -r raw/rna-seq/readset.txt \
    -d raw/rna-seq/design.txt \
    -o output/rna-pipeline-GRCh38 \
    --config /home/efournie/genpipes/pipelines/rnaseq/rnaseq.base.ini \
        /home/efournie/genpipes/pipelines/rnaseq/rnaseq.cedar.ini \
        $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini