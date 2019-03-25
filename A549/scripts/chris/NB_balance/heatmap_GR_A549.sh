#!/bin/bash

OUTPUT_DIR=output/analyses/heatmap_NIPBL_A549
mkdir -p $OUTPUT_DIR

PEAKS_DIR=output/chip-pipeline-GRCh38/peak_call

BW_DIR_GR=/home/chris/Bureau/sb_cofactor_hr/A549/output/chip-Reddy/tracks

MATRIX_NAME=20190425_GR_A549.gzip
HEATMAP_NAME=20190425_GR_A549_heatmap.png

### computeMatrix
echo "computeMatrix..."
time computeMatrix reference-point --referencePoint center \
	--regionsFileName \
		$PEAKS_DIR/A549_NB/NB_specific_CTRL.bed \
		$PEAKS_DIR/A549_NB/NB_common.bed \
		$PEAKS_DIR/A549_NB/NB_specific_DEX.bed \
	--scoreFileName \
		$BW_DIR_GR/A549_CTRL_NR3C1_rep2_ENCFF181HLP.bw \
		$BW_DIR_GR/A549_DEX_5m_NR3C1_rep1_ENCFF880BMZ.bw \
		$BW_DIR_GR/A549_DEX_10m_NR3C1_rep1_ENCFF792BBU.bw \
		$BW_DIR_GR/A549_DEX_15m_NR3C1_rep1_ENCFF733DBG.bw \
		$BW_DIR_GR/A549_DEX_20m_NR3C1_rep1_ENCFF645MAL.bw \
		$BW_DIR_GR/A549_DEX_25m_NR3C1_rep1_ENCFF883JID.bw \
		$BW_DIR_GR/A549_DEX_30m_NR3C1_rep1_ENCFF613TLN.bw \
		$BW_DIR_GR/A549_DEX_1h_NR3C1_rep1_ENCFF331QXR.bw \
	--upstream 2500 --downstream 2500 -p 8 \
	--sortRegions keep \
	--outFileName $OUTPUT_DIR/$MATRIX_NAME

### plotHeatmap
echo "plotHeatmap..."
time plotHeatmap \
	--matrixFile $OUTPUT_DIR/$MATRIX_NAME \
	--colorMap Blues Blues Blues Blues Blues Blues Blues Blues \
	--regionsLabel NB_CTRL NB_common NB_DEX \
	--samplesLabel \
		GR_CTRL GR_DEX_5m GR_DEX_10m GR_DEX_15m GR_DEX_20m GR_DEX_25m GR_DEX_30m GR_DEX_1h \
	--sortRegions no \
	--whatToShow 'heatmap only' \
	--outFileName $OUTPUT_DIR/$HEATMAP_NAME
