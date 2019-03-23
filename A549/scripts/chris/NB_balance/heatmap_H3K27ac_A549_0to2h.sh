#!/bin/bash

OUTPUT_DIR=output/analyses/heatmap_NIPBL_A549
mkdir -p $OUTPUT_DIR

PEAKS_DIR=output/chip-pipeline-GRCh38/peak_call

BW_DIR_Reddy=/home/chris/Bureau/sb_cofactor_hr/A549/output/chip-Reddy/tracks

MATRIX_NAME=20190423_H3K27acme3_A549_0to2h.gzip
HEATMAP_NAME=20190423_H3K27ac_A549_0to2h_heatmap.png

### computeMatrix
echo "computeMatrix..."
time computeMatrix reference-point --referencePoint center \
	--regionsFileName \
		$PEAKS_DIR/A549_NB/NB_specific_CTRL.bed \
		$PEAKS_DIR/A549_NB/NB_common.bed \
		$PEAKS_DIR/A549_NB/NB_specific_DEX.bed \
	--scoreFileName \
		$BW_DIR_Reddy/A549_H3K27ac_CTRL_rep1_ENCFF808SQL.bw \
		$BW_DIR_Reddy/A549_H3K27ac_DEX_5m_rep1_ENCFF012FFF.bw \
		$BW_DIR_Reddy/A549_H3K27ac_DEX_10m_rep1_ENCFF728NYI.bw \
		$BW_DIR_Reddy/A549_H3K27ac_DEX_15m_rep1_ENCFF189MFM.bw \
		$BW_DIR_Reddy/A549_H3K27ac_DEX_20m_rep1_ENCFF622MPU.bw \
		$BW_DIR_Reddy/A549_H3K27ac_DEX_25m_rep1_ENCFF616ZIX.bw \
		$BW_DIR_Reddy/A549_H3K27ac_DEX_30m_rep1_ENCFF039EBI.bw \
		$BW_DIR_Reddy/A549_H3K27ac_DEX_1h_rep1_ENCFF049RFE.bw \
		$BW_DIR_Reddy/A549_H3K27ac_DEX_2h_rep1_ENCFF255KSG.bw \
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
		H3K27ac_CTRL H3K27ac_DEX_5m H3K27ac_DEX_10m H3K27ac_DEX_15m \
		H3K27ac_DEX_20m H3K27ac_DEX_25m H3K27ac_DEX_30m \
		H3K27ac_DEX_1h H3K27ac_DEX_2h \
	--sortRegions no \
	--outFileName $OUTPUT_DIR/$HEATMAP_NAME
