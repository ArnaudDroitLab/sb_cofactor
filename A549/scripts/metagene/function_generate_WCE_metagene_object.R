library(knitr)
library(metagene)

########################################
# Define secondary functions
########################################
generate_WCE_design <- function() {
	Samples <- c("A549_DEX_WCE_rep1", "A549_DEX_WCE_rep2", "A549_CTRL_WCE_rep1", "A549_CTRL_WCE_rep2", "A549_CTRL_WCE_rep3")
	DEX_WCE_rep1 <- c(1, 0, 0, 0, 0)
	DEX_WCE_rep2 <- c(0, 1, 0, 0, 0)
	CTRL_WCE_rep1 <- c(0, 0, 1, 0, 0)
	CTRL_WCE_rep2 <- c(0, 0, 0, 1, 0)
	CTRL_WCE_rep3 <- c(0, 0, 0, 0, 1)
	
	design <- data.frame(Samples, DEX_WCE_rep1, DEX_WCE_rep2, CTRL_WCE_rep1, CTRL_WCE_rep2, CTRL_WCE_rep3)
	
	print(kable(design))
	
	design$Samples = file.path(bam_dir, design$Samples, design$Samples)
	design$Samples = paste0(design$Samples, ".sorted.dup.bam")
	
	return (design)
}

generate_region_names <- function(region_list) {
	region_names <- ""
	for (region in names(region_list)) {
		region_names <- paste0(region_names, "_", region)
	}
	return (region_names)
}

########################################
# Main function
########################################

generate_WCE_metagene_object <- function(cofactor, region_list, bin) {
	mg <- ""
	
	design <- generate_WCE_design()
	
	mg <- metagene$new(regions = region_list, bam_files = design$Samples, assay="chipseq",
					   force_seqlevels = TRUE, verbose = TRUE, cores = 4)
	message("DONE\t", "WCE", " | metagene object")
	
	mg$produce_table(design = design, normalization="RPM", flip_regions=TRUE, bin_count=bin)
	message("DONE\t", "WCE", " | metagene table")
	
	mg$produce_data_frame()
	message("DONE\t", "WCE", " | metagene data fame")

	regions_names <- generate_region_names(region_list)
	
	output_filepath <- file.path(output_dir, paste0("WCE_", cofactor, regions_names, "_metagene_obj.RData"))
	
	save(mg, file = output_filepath)
	message("Saved in ", output_filepath)
}

