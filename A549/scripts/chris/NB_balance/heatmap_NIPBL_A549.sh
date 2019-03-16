#!/bin/bash

OUTPUT_DIR=output/chip-pipeline-GRCh38/analysis/heatmap_NIPBL_A549
mkdir -p $OUTPUT_DIR

PEAKS_DIR=output/chip-pipeline-GRCh38/peak_call

BW_DIR_A549=/home/chris/Bureau/sb_cofactor_hr/A549/output/chip-pipeline-GRCh38/tracks

MATRIX_NAME=20190416_NIPBL_A549.gzip
HEATMAP_NAME=20190416_NIPBL_A549_heatmap_greens.png

### computeMatrix
time computeMatrix reference-point --referencePoint center \
	--regionsFileName \
		$PEAKS_DIR/A549_CTRL_NIPBL_rep1/A549_CTRL_NIPBL_rep1_peaks.narrowPeak.stdchr.bed \
		$PEAKS_DIR/A549_DEX_NIPBL_rep1/A549_DEX_NIPBL_rep1_peaks.narrowPeak.stdchr.bed \
	--scoreFileName \
		$BW_DIR_A549/A549_CTRL_NIPBL_rep1.bw \
		$BW_DIR_A549/A549_CTRL_BRD4_rep1.bw \
		$BW_DIR_A549/A549_CTRL_CDK9_rep1.bw \
		$BW_DIR_A549/A549_CTRL_MED1_rep1.bw \
		$BW_DIR_A549/A549_CTRL_SMC1A_rep1.bw \
		$BW_DIR_A549/A549_DEX_NIPBL_rep1.bw \
		$BW_DIR_A549/A549_DEX_BRD4_rep1.bw \
		$BW_DIR_A549/A549_DEX_CDK9_rep1.bw \
		$BW_DIR_A549/A549_DEX_MED1_rep1.bw \
		$BW_DIR_A549/A549_DEX_SMC1A_rep1.bw \
	--upstream 1000 --downstream 1000 -p 8 \
	--outFileName $OUTPUT_DIR/$MATRIX_NAME

### plotHeatmap
time plotHeatmap \
	--matrixFile $OUTPUT_DIR/$MATRIX_NAME \
	--colorMap Greens \
	--regionsLabel A549_NIPBL_CTRL A549_NIPBL_DEX \
	--samplesLabel \
		NIPBL_CTRL BRD4_CTRL CDK9_CTRL MED1_CTRL SMC1A_CTRL \
		NIPBL_DEX BRD4_DEX CDK9_DEX MED1_DEX SMC1A_DEX \
	--whatToShow 'heatmap only'\
	--outFileName $OUTPUT_DIR/$HEATMAP_NAME
