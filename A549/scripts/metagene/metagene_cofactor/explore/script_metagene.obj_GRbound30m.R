library(knitr)
library(dplyr)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

source("scripts/metagene/function_generate_metagene_object.R")
source("scripts/metagene/function_resize_peaks.R")
source("scripts/load_reddy.R")

########################################
# Define directory
########################################
bam_dir <- "output/chip-pipeline-GRCh38/alignment"
output_dir <- "output/chip-pipeline-GRCh38/metagene/metagene_cofactor"

###############################################################################
# Define regions over which metagenes will be plotted.
###############################################################################

# Import GR binding regions.
gr_regions = load_reddy_gr_binding_consensus()

#
gr_regions_30 <- gr_regions[["30 minutes"]]
gr_regions_30 <- resize_all_peaks(gr_regions_30, window = 600)

#
region_list = list(GR_Regions_30m = gr_regions_30)

###############################################################################
# Generate the metagene object.
###############################################################################

# cofactor_list <- c("BRD4", "CDK9", "NIPBL", "SMC1A", "MED1")

cofactor_list <- c("MED1")
	
for (cofactor in cofactor_list) {
	message("##########     ", cofactor, "     ##########")
	generate_metagene_object(cofactor, region_list, bin=200)
}


