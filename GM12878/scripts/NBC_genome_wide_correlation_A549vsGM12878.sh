#!/bin/bash

OUTPUT_DIR=output/chip-pipeline-GRCh38/analyses/NBC_genome_wide_correlation_A549vsGM12878
mkdir -p $OUTPUT_DIR

BAM_A549_DIR=/home/chris/Bureau/sb_cofactor_hr/A549/output/chip-pipeline-GRCh38/alignment
BAM_GM_DIR=output/chip-pipeline-GRCh38/alignment

COUNT_NAME=NBC_genome_wide_correlation_A549vsGM12878_counts.txt
MATRIX_NAME=NBC_genome_wide_correlation_A549vsGM12878_matrix.npz
HEATMAP_NAME=20181122_A549vsGM12878_CDK9_heatmap.png
OUTPUT_BASENAME=NBC_genome_wide_correlation_A549vsGM12878

# multiBamSummary bins
multiBamSummary bins \
	--bamfiles \
		$BAM_A549_DIR/A549_CTRL_NIPBL_rep1/A549_CTRL_NIPBL_rep1.sorted.dup.bam \
		$BAM_A549_DIR/A549_CTRL_BRD4_rep1/A549_CTRL_BRD4_rep1.sorted.dup.bam \
		$BAM_A549_DIR/A549_CTRL_CDK9_rep1/A549_CTRL_CDK9_rep1.sorted.dup.bam \
		$BAM_A549_DIR/A549_DEX_NIPBL_rep1/A549_DEX_NIPBL_rep1.sorted.dup.bam \
		$BAM_A549_DIR/A549_DEX_BRD4_rep1/A549_DEX_BRD4_rep1.sorted.dup.bam \
		$BAM_A549_DIR/A549_DEX_CDK9_rep1/A549_DEX_CDK9_rep1.sorted.dup.bam \
		$BAM_GM_DIR/GM12878_NIPBL_rep1/GM12878_NIPBL_rep1.sorted.dup.bam \
		$BAM_GM_DIR/SRR1636861_BRD4/SRR1636861_BRD4.sorted.dup.bam \
		$BAM_GM_DIR/GM12878_CDK9_rep1/GM12878_CDK9_rep1.sorted.dup.bam \
	--labels \
		A549_NIPBL_CTRL A549_BRD4_CTRL A549_CDK9_CTRL \
		A549_NIPBL_DEX A549_BRD4_DEX A549_CDK9_DEX \
		GM12878_NIPBL_CTRL GM12878_BRD4_CTRL GM12878_CDK9_CTRL \
	--numberOfProcessors 8 \
	--outRawCounts $OUTPUT_DIR/$COUNT_NAME \
	--outFileName $OUTPUT_DIR/$MATRIX_NAME

# plotCorrelation
plotCorrelation --corData $OUTPUT_DIR/$MATRIX_NAME \
		--corMethod pearson \
		--whatToPlot scatterplot \
		--plotFile $OUTPUT_DIR/NBC_genome_wide_correlation_A549vsGM12878_scatterplot_pearson.png

plotCorrelation --corData $OUTPUT_DIR/$MATRIX_NAME \
		--corMethod spearman \
		--whatToPlot scatterplot \
		--plotFile $OUTPUT_DIR/NBC_genome_wide_correlation_A549vsGM12878_scatterplot_spearman.png

plotCorrelation --corData $OUTPUT_DIR/$MATRIX_NAME \
		--corMethod pearson \
		--whatToPlot heatmap \
		--colorMap RdBu \
		--plotNumbers \
		--plotFile $OUTPUT_DIR/NBC_genome_wide_correlation_A549vsGM12878_heatmap_pearson.png

plotCorrelation --corData $OUTPUT_DIR/$MATRIX_NAME \
		--corMethod spearman \
		--whatToPlot heatmap \
		--colorMap RdBu \
		--plotNumbers \
		--plotFile $OUTPUT_DIR/NBC_genome_wide_correlation_A549vsGM12878_heatmap_spearman.png
