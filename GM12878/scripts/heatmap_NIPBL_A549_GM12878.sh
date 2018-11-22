#!/bin/bash

OUTPUT_DIR=output/chip-pipeline-GRCh38/analysis/heatmap_NIPBL_A549_GM12878
mkdir -p $OUTPUT_DIR

PEAKS_DIR=output/chip-pipeline-GRCh38/peak_call/A549vsGM12878_NIPBL

BW_DIR_A549=/home/chris/Bureau/sb_cofactor_hr/A549/output/chip-pipeline-GRCh38/tracks
BW_DIR_GM12878=output/chip-pipeline-GRCh38/tracks

MATRIX_NAME=20181121_A549vsGM12878_NIPBL_matrix.gzip
HEATMAP_NAME=20181121_A549vsGM12878_NIPBL_heatmap.png

### computeMatrix
computeMatrix reference-point --referencePoint center \
	--regionsFileName \
		$PEAKS_DIR/A549vsGM12878_NIPBL_A_specific.bed \
		$PEAKS_DIR/A549vsGM12878_NIPBL_common.bed \
		$PEAKS_DIR/A549vsGM12878_NIPBL_GM_specific.bed \
	--scoreFileName \
		$BW_DIR_A549/A549_CTRL_NIPBL_rep1.bw \
		$BW_DIR_A549/A549_CTRL_BRD4_rep1.bw \
		$BW_DIR_A549/A549_CTRL_CDK9_rep1.bw \
		$BW_DIR_A549/A549_CTRL_MED1_rep1.bw \
		$BW_DIR_A549/A549_CTRL_SMC1A_rep1.bw \
		$BW_DIR_GM12878/GM12878_NIPBL_rep1.bw \
		$BW_DIR_GM12878/SRR1636861_BRD4.bw \
		$BW_DIR_GM12878/GM12878_CDK9_rep1.bw \
		$BW_DIR_GM12878/GM12878_MED1_rep1.bw \
		$BW_DIR_GM12878/GM12878_SMC1_rep1.bw \
	--upstream 1000 --downstream 1000 -p 8 \
	--outFileName $OUTPUT_DIR/$MATRIX_NAME

### plotHeatmap
plotHeatmap \
	--matrixFile $OUTPUT_DIR/$MATRIX_NAME \
	--colorMap rainbow \
	--regionsLabel A549_NIPBL COMMON_NIPBL GM12878_NIPBL \
	--samplesLabel \
		A549_NIPBL A549_BRD4 A549_CDK9 A549_MED1 A549_SMC1A \
		GM12878_NIPBL GM12878_BRD4 GM12878_CDK9 GM12878_MED1 GM12878_SMC1_rep1 \ 		--yMax 1000 \
	--outFileName $OUTPUT_DIR/$HEATMAP_NAME
