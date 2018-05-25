library(knitr)
library(dplyr)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(rtracklayer)

source("scripts/metagene/function_generate_metagene_object.R")
source("scripts/metagene/function_generate_WCE_metagene_object.R")
source("scripts/metagene/function_resize_peaks.R")
source("scripts/load_reddy.R")

########################################
# Define directory
########################################
bam_dir <- "output/chip-pipeline-GRCh38/alignment"
binding_diff_dir <- "output/chip-pipeline-GRCh38/binding_diff"

###############################################################################
# Loop on all cofactors
###############################################################################
# cofactor_list <- c("BRD4", "CDK9", "NIPBL", "SMC1A", "MED1")

cofactor_list <- c("BRD4", "CDK9")
	
for (cofactor in cofactor_list) {
	message("##########     ", cofactor, "     ##########")
	output_dir <- file.path(binding_diff_dir, paste0("A549_", cofactor), "output_filters")
	
	###############################################################################
	# Define regions over which metagenes will be plotted.
	###############################################################################
	
	CstRegions_path <- file.path(output_dir, paste0("A549_", cofactor, "_unbiased_peaks.bed"))
	BindUpRegions_path <- file.path(output_dir, paste0("A549_DEX_", cofactor, "_rep1_peaks.narrowPeak_M_above_1.0_biased_peaks.bed"))
	BindDownRegions_path <- file.path(output_dir, paste0("A549_CTRL_", cofactor, "_rep1_peaks.narrowPeak_M_below_-1.0_biased_peaks.bed"))
	
	CstRegions <- rtracklayer::import(CstRegions_path)
	BindUpRegions <- rtracklayer::import(BindUpRegions_path)
	BindDownregions <- rtracklayer::import(BindDownRegions_path)
	
	region_list = list(CstRegions = CstRegions,
					   BindUpRegions = BindUpRegions,
					   BindDownregions = BindDownregions)
	
	###############################################################################
	# Generate the metagene object.
	###############################################################################
	generate_WCE_metagene_object(region_list, bin=200)
	generate_metagene_object(cofactor, region_list, bin=200)
}
