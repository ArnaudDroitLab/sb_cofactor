library(knitr)
library(metagene)

########################################
# Define secondary functions
########################################

generate_samples <- function(target, sh) {
	Samples <- vector("character")
	for (condition in c("Dex", "EtOH")) {
			OneSample <- paste0(target, "_", sh, "_", condition, "_Rep1")
			Samples <- c(Samples, OneSample)
	}
	return (Samples)
}

generate_colnames_design <- function(target, sh) {
	colnames_design <- c("Samples")
	for (condition in c("Dex", "EtOh")) {
			OneColname <- paste0(target, "_", sh, "_", condition, "_Rep1")
			colnames_design <- c(colnames_design, OneColname)
	}
	return (colnames_design)
}

generate_design <- function(target, sh) {
	Samples <- generate_samples(target, sh)
	design <- data.frame(Samples, diag(2))
	colnames(design) <- generate_colnames_design(target, sh)
	
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

generate_metagene_object <- function(target, sh, region_list, bin) {
	design <- generate_design(target, sh)
	
	mg <- metagene$new(regions = region_list, bam_files = design$Samples, assay="chipseq",
					   force_seqlevels = TRUE, verbose = TRUE, cores = 4)
	message("DONE\t", target, "_", sh , " | metagene object")
	
	mg$produce_table(design = design, normalization="RPM", flip_regions=TRUE, bin_count=bin)
	message("DONE\t", target, "_", sh , " | metagene table")
	
	mg$produce_data_frame()
	message("DONE\t", target, "_", sh , " | metagene data fame")

	regions_names <- generate_region_names(region_list)
	
	output_filepath <- file.path(output_dir, paste0(target, "_", sh, regions_names, "_metagene_obj.RData"))
	
	save(mg, file = output_filepath)
	message("Saved in ", output_filepath)
}
