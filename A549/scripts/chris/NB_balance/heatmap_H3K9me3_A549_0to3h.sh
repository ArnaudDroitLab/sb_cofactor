#!/bin/bash

OUTPUT_DIR=output/analyses/heatmap_NIPBL_A549
mkdir -p $OUTPUT_DIR

PEAKS_DIR=output/chip-pipeline-GRCh38/peak_call

BW_DIR_Reddy=/home/chris/Bureau/sb_cofactor_hr/A549/output/chip-Reddy/tracks

MATRIX_NAME=20190422_H3K9me3_A549_0to3h.gzip
HEATMAP_NAME=20190422_H3K9me3_A549_0to3h_heatmap.png

### computeMatrix
echo "computeMatrix..."
time computeMatrix reference-point --referencePoint center \
	--regionsFileName \
		$PEAKS_DIR/A549_NB/NB_specific_CTRL.bed \
		$PEAKS_DIR/A549_NB/NB_common.bed \
		$PEAKS_DIR/A549_NB/NB_specific_DEX.bed \
	--scoreFileName \
		$BW_DIR_Reddy/A549_H3K9me3_CTRL_rep1_ENCFF801BLX.bam.bw \
		$BW_DIR_Reddy/A549_H3K9me3_DEX_30m_rep1_ENCFF104HAJ.bw \
		$BW_DIR_Reddy/A549_H3K9me3_DEX_1h_rep1_ENCFF506RDX.bw \
		$BW_DIR_Reddy/A549_H3K9me3_DEX_2h_rep1_ENCFF255IZP.bw \
		$BW_DIR_Reddy/A549_H3K9me3_DEX_3h_rep1_ENCFF561NAT.bw \
	--upstream 1000 --downstream 1000 -p 8 \
	--sortRegions keep \
	--outFileName $OUTPUT_DIR/$MATRIX_NAME

### plotHeatmap
echo "plotHeatmap..."
time plotHeatmap \
	--matrixFile $OUTPUT_DIR/$MATRIX_NAME \
	--colorMap Blues \
	--regionsLabel NB_CTRL NB_common NB_DEX \
	--samplesLabel \
		H3K9me3_CTRL H3K9me3_DEX_30m H3K9me3_DEX_1h H3K9me3_DEX_2h H3K9me3_DEX_3h \
	--sortRegions no \
	--outFileName $OUTPUT_DIR/$HEATMAP_NAME
