#!/bin/bash

OUTPUT_DIR=output/analyses/heatmap_NIPBL_A549
mkdir -p $OUTPUT_DIR

PEAKS_DIR=output/chip-pipeline-GRCh38/peak_call

BW_DIR_A549=/home/chris/Bureau/sb_cofactor_hr/A549/output/chip-pipeline-GRCh38/tracks
BW_DIR_Reddy=/home/chris/Bureau/sb_cofactor_hr/A549/output/chip-Reddy/tracks

MATRIX_NAME=20190425_NB_EP300_A549.gzip
HEATMAP_NAME=20190425_NB_EP300_A549_heatmap.png

### computeMatrix
echo "computeMatrix..."
time computeMatrix reference-point --referencePoint center \
	--regionsFileName \
		$PEAKS_DIR/A549_NB/NB_specific_CTRL.bed \
		$PEAKS_DIR/A549_NB/NB_common.bed \
		$PEAKS_DIR/A549_NB/NB_specific_DEX.bed \
	--scoreFileName \
		$BW_DIR_A549/A549_CTRL_NIPBL_rep1.bw \
		$BW_DIR_A549/A549_CTRL_BRD4_rep1.bw \
		$BW_DIR_Reddy/A549_EP300_CTRL_rep1_ENCFF859CJY.bw \
		$BW_DIR_A549/A549_DEX_NIPBL_rep1.bw \
		$BW_DIR_A549/A549_DEX_BRD4_rep1.bw \
		$BW_DIR_Reddy/A549_EP300_DEX_1h_rep1_ENCFF074CYV.bw \
	--upstream 2500 --downstream 2500 -p 8 \
	--sortRegions keep \
	--outFileName $OUTPUT_DIR/$MATRIX_NAME

### plotHeatmap
echo "plotHeatmap..."
time plotHeatmap \
	--matrixFile $OUTPUT_DIR/$MATRIX_NAME \
	--colorMap Greens Greens Greens Oranges Oranges Oranges \
	--regionsLabel NB_CTRL NB_common NB_DEX \
	--samplesLabel \
		NIPBL_CTRL BRD4_CTRL EP300_CTRL \
		NIPBL_DEX BRD4_DEX EP300_DEX \
	--sortRegions no \
	--outFileName $OUTPUT_DIR/$HEATMAP_NAME