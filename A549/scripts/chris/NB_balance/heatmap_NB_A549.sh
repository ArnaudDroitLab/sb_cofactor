#!/bin/bash

OUTPUT_DIR=output/heatmap_NIPBL
mkdir -p $OUTPUT_DIR

A549_DIR=/home/chris/Bureau/sb_cofactor_hr/A549/output/chip-pipeline-GRCh38
HEPG2_DIR=/home/chris/Bureau/sb_cofactor_hr/HEPG2/output/chip-pipeline-GRCh38
K562_DIR=/home/chris/Bureau/sb_cofactor_hr/K562/output/chip-pipeline-GRCh38
MCF7_DIR=/home/chris/Bureau/sb_cofactor_hr/MCF7/output/chip-pipeline-GRCh38

##### regionsDir
A549_PEAKS_DIR=$A549_DIR/peak_call
HEPG2_PEAKS_DIR=$HEPG2_DIR/peak_call
K562_PEAKS_DIR=$K562_DIR/peak_call
MCF7_PEAKS_DIR=$MCF7_DIR/peak_call

##### scoreDir
A549_BW_DIR=$A549_DIR/tracks
HEPG2_BW_DIR=$HEPG2_DIR/tracks
K562_BW_DIR=$K562_DIR/tracks
MCF7_BW_DIR=$MCF7_DIR/tracks

##### output filenames
MATRIX_NAME=20190307_heatmap_NIPBL_v3.gzip
HEATMAP_NAME=20190307_heatmap_NIPBL_v3.png

### computeMatrix
echo 'computeMatrix...'
time computeMatrix reference-point --referencePoint center \
	--regionsFileName \
		$A549_PEAKS_DIR/A549_CTRL_NIPBL_rep1/A549_CTRL_NIPBL_rep1_peaks.narrowPeak.bed \
		$HEPG2_PEAKS_DIR/HEPG2_NIPBL_rep1/HEPG2_NIPBL_rep1_peaks.narrowPeak.bed \
		$HEPG2_PEAKS_DIR/HEPG2_NIPBL_rep2/HEPG2_NIPBL_rep2_peaks.narrowPeak.bed \
		$K562_PEAKS_DIR/K562_NIPBL_rep1/K562_NIPBL_rep1_peaks.narrowPeak.bed \
		$MCF7_PEAKS_DIR/MCF7_CTRL_NIPBL_rep1/MCF7_CTRL_NIPBL_rep1_peaks.narrowPeak.bed \
		$MCF7_PEAKS_DIR/MCF7_CTRL_NIPBL_rep2/MCF7_CTRL_NIPBL_rep2_peaks.narrowPeak.bed \
	--scoreFileName \
		$A549_BW_DIR/A549_CTRL_NIPBL_rep1.bw \
		$A549_BW_DIR/A549_CTRL_BRD4_rep1.bw \
		$A549_BW_DIR/A549_CTRL_CDK9_rep1.bw \
		$HEPG2_BW_DIR/HEPG2_NIPBL_rep1.bw \
		$HEPG2_BW_DIR/HEPG2_NIPBL_rep2.bw \
		$HEPG2_BW_DIR/HEPG2_BRD4_rep1.bw \
		$HEPG2_BW_DIR/HEPG2_BRD4_rep2.bw \
		$HEPG2_BW_DIR/HEPG2_CDK9_rep1.bw \
		$K562_BW_DIR/K562_NIPBL_rep1.bw \
		$K562_BW_DIR/K562_BRD4_rep1.bw \
		$K562_BW_DIR/K562_BRD4_rep1.bw \
		$K562_BW_DIR/K562_CDK9_rep1.bw \
		$MCF7_BW_DIR/MCF7_CTRL_NIPBL_rep1.bw \
		$MCF7_BW_DIR/MCF7_CTRL_NIPBL_rep2.bw \
		$MCF7_BW_DIR/MCF7_CTRL_BRD4_rep1.bw \
		$MCF7_BW_DIR/MCF7_CTRL_BRD4_rep2.bw \
		$MCF7_BW_DIR/MCF7_CTRL_BRD4_rep3.bw \
	--upstream 1000 --downstream 1000 -p 8 \
	--outFileName $OUTPUT_DIR/$MATRIX_NAME

### plotHeatmap
echo 'plotHeatmap...'
time plotHeatmap \
	--matrixFile $OUTPUT_DIR/$MATRIX_NAME \
	--colorMap Purples Purples Purples \
						 Oranges Oranges Oranges Oranges Oranges \
						 Blues Blues Blues Blues \
						 Greens Greens Greens Greens Greens \
	--regionsLabel A549 HEPG2_rep1 HEPG2_rep2 K562 MCF7_rep1 MCF7_rep2 \
	--samplesLabel \
		A549_NIPBL A549_BRD4 A549_CDK9 \
		HEPG2_NIPBL_rep1 HEPG2_NIPBL_rep2 HEPG2_BRD4_rep1 HEPG2_BRD4_rep2 HEPG2_CDK9 \
		K562_NIPBL K562_BRD4_rep1 K562_BRD4_rep2 K562_CDK9 \
		MCF7_NIPBL_rep1 MCF7_NIPBL_rep2 MCF7_BRD4_rep1 MCF7_BRD4_rep2 MCF7_BRD4_rep3 \
	--whatToShow 'heatmap only'\
	--outFileName $OUTPUT_DIR/$HEATMAP_NAME

