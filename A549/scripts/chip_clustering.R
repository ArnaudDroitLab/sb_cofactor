source("scripts/load_reddy.R")

# Loads all peaks calls for a given target and build consensus regions
# for each time point.
summarize_GGR_chip <- function(encode_results, diagnostic_dir=NULL) {
    # Make sure the passed in ENCODE subset ahs only one target.
    target = unique(encode_results$target)
    stopifnot(length(target)==1)
    
    # Determine if the target has broad or narrow peaks so we can
    # load the peaks in an appropriate manner.
    if(all(unique(encode_results$file_type)=="bed narrowPeak")) {
        bed_type="narrow"
    } else if(all(unique(encode_results$file_type)=="bed broadPeak")) {
        bed_type="broad"
    } else {
        stop("Not all downloaded files are of teh same type!")
    }
    
    # Download the peak files.
    download_dir = file.path("input/ENCODE/A549/GRCh38/chip-seq", bed_type)
    dir.create(download_dir, recursive=TRUE, showWarnings=FALSE)
    downloaded_files = ENCODExplorer::downloadEncode(encode_results, dir=download_dir, force=FALSE)
    
    # Import the peak files into a GRangesList object.
    names(downloaded_files) <- gsub(".bed.gz", "", basename(downloaded_files))
    all_regions = ef.utils::import_files_into_grl(downloaded_files, file.format=bed_type, 
                                                  file.ext=".bed.gz", discard.metadata=TRUE)
    
    # Determine the time_point for all downloaded files.
    time_points = paste(encode_results$treatment_duration, encode_results$treatment_duration_unit)
    # Files without a treatment (NA) are controls, which we label "0 minute"
    time_points = gsub("NA NA", "0 minute", time_points)
    names(time_points) = encode_results$file_accession
    
    # Group files by time point, then intersect them.
    results = list()
    for(time_point in unique(time_points)) {
        intersect_object = ef.utils::build_intersect(all_regions[names(time_points[time_points==time_point])])
        two_replicates = rowSums(intersect_object$Matrix) >= 2
        results[[time_point]] = intersect_object$Regions[two_replicates]
        
        if(!is.null(diagnostic_dir)) {
            # Turn off VennDiagram logging.
            futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
        
            dir.create(file.path(diagnostic_dir, "Venn diagrams of peak calls"), recursive=TRUE, showWarnings=FALSE)
            ef.utils::intersect_venn_plot(intersect_object, file.path(diagnostic_dir, "Venn diagrams of peak calls", paste0(target, "-", time_point, ".tiff")))
        }
    }
    
    return(GenomicRanges::GRangesList(results))
}

# Load all ENCODE chips.
all_chip = queryEncodeGeneric(biosample_name="A549", file_format="bed", lab="Tim Reddy, Duke", assay="ChIP-seq")

# Keeping the number of replicates at three makes ChIPs more compaable. The latest Reddy
# batch all coems in triplicates. So we'll remove older batches when there are more
# than 3 replicates.
# Remove a couple of extraneous JUN replicates
all_chip = all_chip %>% dplyr::filter(!(target=="JUN" & date_released=="2016-12-20" &is.na(treatment)))

# Do the same thing for some H3K27ac replicates.
all_chip = all_chip %>% dplyr::filter(!(target=="H3K27ac" & date_released=="2016-07-21" &is.na(treatment)))
all_targets = unique(all_chip$target)

all_chip_regions = list()
for(target_name in all_targets) {
    all_chip_regions[[target_name]] = summarize_GGR_chip(all_chip %>% dplyr::filter(target==target_name),
                                                         diagnostic_dir="output/analyses/cofactors")
}

# Load our own ChIP.
lab_chip = load_cofactor_binding()
for(target in names(lab_chip[["DEX"]])) {
    all_chip_regions[[target]] = GRangesList("1 hour"   = lab_chip[["DEX"]][[target]],
                                             "0 minute" = lab_chip[["CTRL"]][[target]])
}


all_time_points = unique(unlist(lapply(all_chip_regions, names)))
for(time_point in all_time_points) {
    target_at_time = list()
    for(target in names(all_chip_regions)) {
        if(time_point %in% names(all_chip_regions[[target]])) {
            target_at_time[[target]] = all_chip_regions[[target]][[time_point]]
        }
    }
    target_at_time = GenomicRanges::GRangesList(target_at_time)
    
    time_intersect = ef.utils::build_intersect(target_at_time)
    pdf(paste0("output/analyses/cofactors/Clustering of cofactors at time ", time_point, ".pdf"))
    plot(hclust(dist(t(time_intersect$Matrix))))
    dev.off()
}
