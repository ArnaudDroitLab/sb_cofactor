mkdir -p $SCRATCH/sb_cofactor/MCF7/output/chip-pipeline-GRCh38

chipseq.py -j slurm -s '1-10,12' \
    -l debug \
    -r $SCRATCH/sb_cofactor/MCF7/raw/chip-seq/readset_2019-02-23.txt \
    -d $SCRATCH/sb_cofactor/MCF7/raw/chip-seq/design_2019-02-23.txt \
    -o $SCRATCH/sb_cofactor/MCF7/output/chip-pipeline-GRCh38 \
    --config $MUGQIC_PIPELINES_HOME/pipelines/chipseq/chipseq.base.ini \
        $MUGQIC_PIPELINES_HOME/pipelines/chipseq/chipseq.cedar.ini \
        $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.ini \
        input/cedar_assembly_bugfix.txt

