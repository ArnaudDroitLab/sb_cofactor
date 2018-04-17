library(knitr)
library(metagene)

########################################
# Define directory
########################################
bam_dir <- "output/chip-pipeline-GRCh38/alignment"
output_dir <- "output/chip-pipeline-GRCh38/metagene_cofactor"

########################################
# Define secondary functions
########################################

generate_samples <- function(cofactor) {
	Samples <- vector("character")
	for (condition in c("DEX", "CTRL")) {
		for (rep in c("rep1","rep2")) {
			OneSample <- paste0("A549", "_", condition, "_", cofactor, "_", rep)
			Samples <- c(Samples, OneSample)
		}
	}
	return (Samples)
}

generate_colnames_design <- function(cofactor) {
	colnames_design <- c("Samples")
	for (condition in c("DEX", "CTRL")) {
		for (rep in c("rep1","rep2")) {
			OneColname <- paste0(condition, "_", cofactor, "_", rep)
			colnames_design <- c(colnames_design, OneColname)
		}
	}
	return (colnames_design)
}

generate_design <- function(cofactor) {
	Samples <- generate_samples(cofactor)
	design <- data.frame(Samples, diag(4))
	colnames(design) <- generate_colnames_design(cofactor)
	
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

generate_metagene_object <- function(cofactor, region_list) {
	design <- generate_design(cofactor)
	
	mg <- metagene$new(regions = region_list, bam_files = design$Samples, assay="chipseq",
					   force_seqlevels = TRUE, verbose = TRUE, cores = 4)
	message("DONE\t", cofactor, " | metagene object")
	
	mg$produce_table(design = design, normalization="RPM", flip_regions=TRUE, bin_count=100)
	message("DONE\t", cofactor, " | metagene table")
	
	mg$produce_data_frame()
	message("DONE\t", cofactor, " | metagene data fame")

	regions_names <- generate_region_names(region_list)
	
	output_filepath <- file.path(output_dir, paste0(cofactor, regions_names, "_metagene_obj.RData"))
	
	save(mg, file = output_filepath)
	message("Saved in ", output_filepath)
}