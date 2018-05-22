source("scripts/load_reddy.R")

summarize_GGR_chip <- function(target, diagnostic_dir=NULL) {
    peak_files = ENCODExplorer::queryEncodeGeneric(biosample_name="A549", treatment="dexamethasone", 
                                                   file_format="bed", lab="Tim Reddy, Duke",
                                                   treatment_amount=100, treatment_amount_unit="nM",
                                                   target=target)
    
    if(all(unique(peak_files$file_type)=="bed narrowPeak")) {
        bed_type="narrow"
    } else if(all(unique(peak_files$file_type)=="bed broadPeak")) {
        bed_type="broad"
    } else {
        stop("Not all downloaded files are of teh same type!")
    }
    
    download_dir = file.path("input/ENCODE/A549/GRCh38/chip-seq", bed_type)
    dir.create(download_dir, recursive=TRUE, showWarnings=FALSE)
    downloaded_files = ENCODExplorer::downloadEncode(peak_files, dir=download_dir, force=FALSE)
    
    names(downloaded_files) <- gsub(".bed.gz", "", basename(downloaded_files))
    all_regions = ef.utils::import_files_into_grl(downloaded_files, file.format=bed_type, 
                                                  file.ext=".bed.gz", discard.metadata=TRUE)
    
    # Group files by time point, then intersect them.
    time_points = paste(peak_files$treatment_duration, peak_files$treatment_duration_unit)
    names(time_points) = peak_files$file_accession
    
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
all_chip = queryEncodeGeneric(biosample_name="A549", treatment="dexamethasone", file_format="bed", treatment_amount=100, treatment_amount_unit="nM", lab="Tim Reddy, Duke", assay="ChIP-seq")
all_targets = unique(all_chip$target)

all_chip_regions = list()
for(target in all_targets) {
    all_chip_regions[[target]] = summarize_GGR_chip(target, diagnostic_dir="output/analyses/cofactors")
}

# Load our own ChIP.
lab_chip = load_cofactor_binding()

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
