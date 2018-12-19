#!/bin/bash

OUTPUT_DIR=output/analyses/NBC_genome_wide_correlation
mkdir -p $OUTPUT_DIR

BAM_DIR=/home/chris/Bureau/sb_cofactor_hr/A549/output/chip-pipeline-GRCh38/alignment

# multiBamSummary bins
multiBamSummary bins \
	--bamfiles \
		$BAM_DIR/A549_CTRL_NIPBL_rep1/A549_CTRL_NIPBL_rep1.sorted.dup.bam \
		$BAM_DIR/A549_CTRL_BRD4_rep1/A549_CTRL_BRD4_rep1.sorted.dup.bam \
		$BAM_DIR/A549_CTRL_CDK9_rep1/A549_CTRL_CDK9_rep1.sorted.dup.bam \
		$BAM_DIR/A549_DEX_NIPBL_rep1/A549_DEX_NIPBL_rep1.sorted.dup.bam \
		$BAM_DIR/A549_DEX_BRD4_rep1/A549_DEX_BRD4_rep1.sorted.dup.bam \
		$BAM_DIR/A549_DEX_CDK9_rep1/A549_DEX_CDK9_rep1.sorted.dup.bam \
	--labels NIPBL_CTRL BRD4_CTRL CDK9_CTRL NIPBL_DEX BRD4_DEX CDK9_DEX \
	--numberOfProcessors 8 \
	--outRawCounts $OUTPUT_DIR/NBC_genome_wide_correlation_counts.txt \
	--outFileName $OUTPUT_DIR/NBC_genome_wide_correlation_corrmatrix.npz

# plotCorrelation
plotCorrelation --corData $OUTPUT_DIR/NBC_genome_wide_correlation_corrmatrix.npz \
		--corMethod pearson \
		--whatToPlot scatterplot \
		--plotFile $OUTPUT_DIR/NBC_genome_wide_correlation_scatterplot_pearson.png

plotCorrelation --corData $OUTPUT_DIR/NBC_genome_wide_correlation_corrmatrix.npz \
		--corMethod spearman \
		--whatToPlot scatterplot \
		--plotFile $OUTPUT_DIR/NBC_genome_wide_correlation_scatterplot_spearman.png

plotCorrelation --corData $OUTPUT_DIR/NBC_genome_wide_correlation_corrmatrix.npz \
		--corMethod pearson \
		--whatToPlot heatmap \
		--colorMap RdBu \
		--plotNumbers \
		--plotFile $OUTPUT_DIR/NBC_genome_wide_correlation_heatmap_pearson.png

plotCorrelation --corData $OUTPUT_DIR/NBC_genome_wide_correlation_corrmatrix.npz \
		--corMethod spearman \
		--whatToPlot heatmap \
		--colorMap RdBu \
		--plotNumbers \
		--plotFile $OUTPUT_DIR/NBC_genome_wide_correlation_heatmap_spearman.png
