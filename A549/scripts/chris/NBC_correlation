#!/bin/bash

# 3 sets of regions:
# /home/chris/Bureau/sb_cofactor_hr/A549/output/chip-pipeline-GRCh38/peak_call/A549_NBC/A549_NBC_CTRL_specific.bed  
# /home/chris/Bureau/sb_cofactor_hr/A549/output/chip-pipeline-GRCh38/peak_call/A549_NBC/A549_NBC_common.bed
# /home/chris/Bureau/sb_cofactor_hr/A549/output/chip-pipeline-GRCh38/peak_call/A549_NBC/A549_NBC_DEX_specific.bed

multiBamSummary BED-file \
	--bamfiles /home/chris/Bureau/sb_cofactor_hr/A549/output/chip-pipeline-GRCh38/alignment/A549_CTRL_NIPBL_rep1/A549_CTRL_NIPBL_rep1.sorted.dup.bam \
		  /home/chris/Bureau/sb_cofactor_hr/A549/output/chip-pipeline-GRCh38/alignment/A549_CTRL_BRD4_rep1/A549_CTRL_BRD4_rep1.sorted.dup.bam \
		  /home/chris/Bureau/sb_cofactor_hr/A549/output/chip-pipeline-GRCh38/alignment/A549_CTRL_CDK9_rep1/A549_CTRL_CDK9_rep1.sorted.dup.bam \
		  /home/chris/Bureau/sb_cofactor_hr/A549/output/chip-pipeline-GRCh38/alignment/A549_DEX_NIPBL_rep1/A549_DEX_NIPBL_rep1.sorted.dup.bam \
		  /home/chris/Bureau/sb_cofactor_hr/A549/output/chip-pipeline-GRCh38/alignment/A549_DEX_BRD4_rep1/A549_DEX_BRD4_rep1.sorted.dup.bam \
		  /home/chris/Bureau/sb_cofactor_hr/A549/output/chip-pipeline-GRCh38/alignment/A549_DEX_CDK9_rep1/A549_DEX_CDK9_rep1.sorted.dup.bam \
	--labels NIPBL_CTRL BRD4_CTRL CDK9_CTRL NIPBL_DEX BRD4_DEX CDK9_DEX \
	--BED /home/chris/Bureau/sb_cofactor_hr/A549/output/chip-pipeline-GRCh38/peak_call/A549_NBC/A549_NBC_CTRL_specific.bed \
	--numberOfProcessors 8 \
	--outRawCounts /home/chris/Bureau/sb_cofactor_hr/A549/output/analyses/NBC_correlation/A549_NBC_CTRL_specific_fromBam_counts.txt \
	--outFileName /home/chris/Bureau/sb_cofactor_hr/A549/output/analyses/NBC_correlation/A549_NBC_CTRL_specific_fromBam_corrmatrix.npz

multiBamSummary BED-file \
	--bamfiles /home/chris/Bureau/sb_cofactor_hr/A549/output/chip-pipeline-GRCh38/alignment/A549_CTRL_NIPBL_rep1/A549_CTRL_NIPBL_rep1.sorted.dup.bam \
		  /home/chris/Bureau/sb_cofactor_hr/A549/output/chip-pipeline-GRCh38/alignment/A549_CTRL_BRD4_rep1/A549_CTRL_BRD4_rep1.sorted.dup.bam \
		  /home/chris/Bureau/sb_cofactor_hr/A549/output/chip-pipeline-GRCh38/alignment/A549_CTRL_CDK9_rep1/A549_CTRL_CDK9_rep1.sorted.dup.bam \
		  /home/chris/Bureau/sb_cofactor_hr/A549/output/chip-pipeline-GRCh38/alignment/A549_DEX_NIPBL_rep1/A549_DEX_NIPBL_rep1.sorted.dup.bam \
		  /home/chris/Bureau/sb_cofactor_hr/A549/output/chip-pipeline-GRCh38/alignment/A549_DEX_BRD4_rep1/A549_DEX_BRD4_rep1.sorted.dup.bam \
		  /home/chris/Bureau/sb_cofactor_hr/A549/output/chip-pipeline-GRCh38/alignment/A549_DEX_CDK9_rep1/A549_DEX_CDK9_rep1.sorted.dup.bam \
	--labels NIPBL_CTRL BRD4_CTRL CDK9_CTRL NIPBL_DEX BRD4_DEX CDK9_DEX \
	--BED /home/chris/Bureau/sb_cofactor_hr/A549/output/chip-pipeline-GRCh38/peak_call/A549_NBC/A549_NBC_common.bed \
	--numberOfProcessors 8 \
	--outRawCounts /home/chris/Bureau/sb_cofactor_hr/A549/output/analyses/NBC_correlation/A549_NBC_common_fromBam_counts.txt \
	--outFileName /home/chris/Bureau/sb_cofactor_hr/A549/output/analyses/NBC_correlation/A549_NBC_common_fromBam_corrmatrix.npz

multiBamSummary BED-file \
	--bamfiles /home/chris/Bureau/sb_cofactor_hr/A549/output/chip-pipeline-GRCh38/alignment/A549_CTRL_NIPBL_rep1/A549_CTRL_NIPBL_rep1.sorted.dup.bam \
		  /home/chris/Bureau/sb_cofactor_hr/A549/output/chip-pipeline-GRCh38/alignment/A549_CTRL_BRD4_rep1/A549_CTRL_BRD4_rep1.sorted.dup.bam \
		  /home/chris/Bureau/sb_cofactor_hr/A549/output/chip-pipeline-GRCh38/alignment/A549_CTRL_CDK9_rep1/A549_CTRL_CDK9_rep1.sorted.dup.bam \
		  /home/chris/Bureau/sb_cofactor_hr/A549/output/chip-pipeline-GRCh38/alignment/A549_DEX_NIPBL_rep1/A549_DEX_NIPBL_rep1.sorted.dup.bam \
		  /home/chris/Bureau/sb_cofactor_hr/A549/output/chip-pipeline-GRCh38/alignment/A549_DEX_BRD4_rep1/A549_DEX_BRD4_rep1.sorted.dup.bam \
		  /home/chris/Bureau/sb_cofactor_hr/A549/output/chip-pipeline-GRCh38/alignment/A549_DEX_CDK9_rep1/A549_DEX_CDK9_rep1.sorted.dup.bam \
	--labels NIPBL_CTRL BRD4_CTRL CDK9_CTRL NIPBL_DEX BRD4_DEX CDK9_DEX \
	--BED /home/chris/Bureau/sb_cofactor_hr/A549/output/chip-pipeline-GRCh38/peak_call/A549_NBC/A549_NBC_DEX_specific.bed \
	--numberOfProcessors 8 \
	--outRawCounts /home/chris/Bureau/sb_cofactor_hr/A549/output/analyses/NBC_correlation/A549_NBC_DEX_specific_fromBam_counts.txt \
	--outFileName /home/chris/Bureau/sb_cofactor_hr/A549/output/analyses/NBC_correlation/A549_NBC_DEX_specific_fromBam_corrmatrix.npz

######
## NBC_CTRL_specific
plotCorrelation --corData /home/chris/Bureau/sb_cofactor_hr/A549/output/analyses/NBC_correlation/A549_NBC_CTRL_specific_fromBam_corrmatrix.npz \
		--corMethod pearson \
		--whatToPlot scatterplot \
		--plotFile /home/chris/Bureau/sb_cofactor_hr/A549/output/analyses/NBC_correlation/A549_NBC_CTRL_specific_fromBam_scatterplot_pearson.png

plotCorrelation --corData /home/chris/Bureau/sb_cofactor_hr/A549/output/analyses/NBC_correlation/A549_NBC_CTRL_specific_fromBam_corrmatrix.npz \
		--corMethod spearman \
		--whatToPlot scatterplot \
		--plotFile /home/chris/Bureau/sb_cofactor_hr/A549/output/analyses/NBC_correlation/A549_NBC_CTRL_specific_fromBam_scatterplot_spearman.png

plotCorrelation --corData /home/chris/Bureau/sb_cofactor_hr/A549/output/analyses/NBC_correlation/A549_NBC_CTRL_specific_fromBam_corrmatrix.npz \
		--corMethod spearman \
		--whatToPlot heatmap \
		--colorMap RdBu \
		--plotNumbers \
		--plotFile /home/chris/Bureau/sb_cofactor_hr/A549/output/analyses/NBC_correlation/A549_NBC_CTRL_specific_fromBam_heatmap_spearman.png

plotCorrelation --corData /home/chris/Bureau/sb_cofactor_hr/A549/output/analyses/NBC_correlation/A549_NBC_CTRL_specific_fromBam_corrmatrix.npz \
		--corMethod pearson \
		--whatToPlot heatmap \
		--colorMap RdBu \
		--plotNumbers \
		--plotFile /home/chris/Bureau/sb_cofactor_hr/A549/output/analyses/NBC_correlation/A549_NBC_common_fromBam_heatmap_pearson.png

## NBC_common
plotCorrelation --corData /home/chris/Bureau/sb_cofactor_hr/A549/output/analyses/NBC_correlation/A549_NBC_common_fromBam_corrmatrix.npz \
		--corMethod pearson \
		--whatToPlot scatterplot \
		--plotFile /home/chris/Bureau/sb_cofactor_hr/A549/output/analyses/NBC_correlation/A549_NBC_common_fromBam_scatterplot_pearson.png

plotCorrelation --corData /home/chris/Bureau/sb_cofactor_hr/A549/output/analyses/NBC_correlation/A549_NBC_common_fromBam_corrmatrix.npz \
		--corMethod spearman \
		--whatToPlot scatterplot \
		--plotFile /home/chris/Bureau/sb_cofactor_hr/A549/output/analyses/NBC_correlation/A549_NBC_common_fromBam_scatterplot_spearman.png

plotCorrelation --corData /home/chris/Bureau/sb_cofactor_hr/A549/output/analyses/NBC_correlation/A549_NBC_common_fromBam_corrmatrix.npz \
		--corMethod spearman \
		--whatToPlot heatmap \
		--colorMap RdBu \
		--plotNumbers \
		--plotFile /home/chris/Bureau/sb_cofactor_hr/A549/output/analyses/NBC_correlation/A549_NBC_common_fromBam_heatmap_spearman.png

plotCorrelation --corData /home/chris/Bureau/sb_cofactor_hr/A549/output/analyses/NBC_correlation/A549_NBC_common_fromBam_corrmatrix.npz \
		--corMethod pearson \
		--whatToPlot heatmap \
		--colorMap RdBu \
		--plotNumbers \
		--plotFile /home/chris/Bureau/sb_cofactor_hr/A549/output/analyses/NBC_correlation/A549_NBC_common_fromBam_heatmap_pearson.png

## NBC_DEX_specific
plotCorrelation --corData /home/chris/Bureau/sb_cofactor_hr/A549/output/analyses/NBC_correlation/A549_NBC_DEX_specific_fromBam_corrmatrix.npz \
		--corMethod pearson \
		--whatToPlot scatterplot \
		--plotFile /home/chris/Bureau/sb_cofactor_hr/A549/output/analyses/NBC_correlation/A549_NBC_DEX_specific_fromBam_scatterplot_pearson.png

plotCorrelation --corData /home/chris/Bureau/sb_cofactor_hr/A549/output/analyses/NBC_correlation/A549_NBC_DEX_specific_fromBam_corrmatrix.npz \
		--corMethod spearman \
		--whatToPlot scatterplot \
		--plotFile /home/chris/Bureau/sb_cofactor_hr/A549/output/analyses/NBC_correlation/A549_NBC_DEX_specific_fromBam_scatterplot_spearman.png

plotCorrelation --corData /home/chris/Bureau/sb_cofactor_hr/A549/output/analyses/NBC_correlation/A549_NBC_DEX_specific_fromBam_corrmatrix.npz \
		--corMethod spearman \
		--whatToPlot heatmap \
		--colorMap RdBu \
		--plotNumbers \
		--plotFile /home/chris/Bureau/sb_cofactor_hr/A549/output/analyses/NBC_correlation/A549_NBC_DEX_specific_fromBam_heatmap_spearman.png

plotCorrelation --corData /home/chris/Bureau/sb_cofactor_hr/A549/output/analyses/NBC_correlation/A549_NBC_DEX_specific_fromBam_corrmatrix.npz \
		--corMethod pearson \
		--whatToPlot heatmap \
		--colorMap RdBu \
		--plotNumbers \
		--plotFile /home/chris/Bureau/sb_cofactor_hr/A549/output/analyses/NBC_correlation/A549_NBC_DEX_specific_fromBam_heatmap_pearson.png
