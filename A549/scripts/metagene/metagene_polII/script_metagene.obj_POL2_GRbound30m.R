library(knitr)
library(dplyr)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

source("scripts/metagene_polII/function_generate_metagene_pol2_object.R")
source("scripts/metagene/function_resize_peaks.R")
source("scripts/load_reddy.R")

########################################
# Define directory
########################################
bam_dir <- "output/chip-pipeline-PolII-GRCh38/alignment"
output_dir <- "output/chip-pipeline-PolII-GRCh38/metagene_polII"

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

target_list <- c("POL2-ser2", "POL2", "WCE")
sh_list <-  c("shCTRL-1", "shCTRL-2", "shNIPBL-3", "shNIPBL-5")

for (target in target_list) {
	for (sh in sh_list) {
	generate_metagene_object(target, sh, region_list, bin=200)
	}
}


