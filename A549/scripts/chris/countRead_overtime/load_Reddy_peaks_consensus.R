setwd("/home/chris/Bureau/sb_cofactor_hr/A549")
library(ENCODExplorer)
library(dplyr)
library(knitr)
library(ef.utils)

##### Make ENCODE_Reddy_ChIP_experiments.txt file
make_ENCODE_Reddy_ChIP_experiments_file <- function(target_name) {
  all_chip_bed <- ENCODExplorer::queryEncodeGeneric(biosample_name="A549", file_format = "bed", assay="ChIP-seq")
  target_bed <- all_chip_bed %>% filter(target == target_name, assembly == "GRCh38", lab == "Tim Reddy, Duke")
  report_target_bed <- target_bed %>% dplyr::select(accession, file_accession, submitted_by, file_format, target, treatment, treatment_duration, treatment_duration_unit, biological_replicates, controls)
  report_target_bed$treatment_duration[is.na(report_target_bed$treatment_duration)] <- 0
  report_target_bed[which(report_target_bed$submitted_by == "Alejandro Barrera"), ]$treatment_duration_unit[is.na(report_target_bed[which(report_target_bed$submitted_by == "Alejandro Barrera"), ]$treatment_duration_unit)] <- "minute"
  report_target_bed[which(report_target_bed$submitted_by == "Ian McDowell"), ]$treatment_duration_unit[is.na(report_target_bed[which(report_target_bed$submitted_by == "Ian McDowell"), ]$treatment_duration_unit)] <- "hour"
  report_target_bed <- report_target_bed %>% arrange(desc(treatment_duration_unit), treatment_duration, submitted_by, biological_replicates)
  print(kable(report_target_bed))
  
  exp_tmp <- report_target_bed %>% dplyr::select(accession, treatment_duration, treatment_duration_unit) %>% unique
  exp_target <- data.frame(Experiment = exp_tmp$accession, Time = paste0(exp_tmp$treatment_duration, " ", exp_tmp$treatment_duration_unit), Order = 1:nrow(exp_tmp))
  print(kable(exp_target))
  
  filename <- paste0("ENCODE_Reddy_", target_name, "_ChIP_experiments.txt")
  write.table(exp_target, file = file.path("input", filename), row.names = FALSE, quote = FALSE, sep = "\t")
  message(filename, " saved in ", file.path(getwd(), "input"))
}

##### load_reddy_binding_consensus function
load_reddy_binding_consensus <- function(target_name, diagnostic_dir=NULL) {
  # Get ENCODE accessions for target binding time series
  filename <- paste0("ENCODE_Reddy_", target_name, "_ChIP_experiments.txt")
  target_accession = read.table(file.path("input", filename), header=TRUE, sep="\t")
  target_accession = target_accession[order(target_accession$Order),]
  target_accession$Time = factor(target_accession$Time, levels=target_accession$Time)
  
  chip_dir = "input/ENCODE/A549/GRCh38/chip-seq/narrow"
  dir.create(chip_dir, recursive=TRUE, showWarnings=FALSE)
  
  # Loop over all ENCODE accessions, downloading all replicates and 
  # extracting consensus regions.
  target_regions = list()
  for(i in 1:nrow(target_accession)) {
    encode_accession = target_accession$Experiment[i]
    time_point = as.character(target_accession$Time[i])
    
    # Download and import peak calls for the time point.
    encodeResults = ENCODExplorer::queryEncodeGeneric(accession=encode_accession, file_type="bed narrowPeak", assembly="GRCh38")
    downloaded_files = ENCODExplorer::downloadEncode(encodeResults, dir=chip_dir, force=FALSE)
    names(downloaded_files) = gsub(".*\\/(.*).bed.gz", "\\1", downloaded_files)
    
    extraCols <- c(signalValue = "numeric", pValue = "numeric", qValue = "numeric", peak = "integer")
    binding_sites = GRangesList(lapply(downloaded_files, rtracklayer::import, extraCols=extraCols))
    
    # Only keep peak calls found in at least two replicates.
    intersect_object = build_intersect(binding_sites)
    two_replicates = rowSums(intersect_object$Matrix) >= 2
    target_regions[[time_point]] = intersect_object$Regions[two_replicates]
    if(!is.null(diagnostic_dir)) {
      # Turn off VennDiagram logging.
      futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
      
      dir.create(file.path(diagnostic_dir, "Venn diagrams of peak calls"), recursive=TRUE, showWarnings=FALSE)
      intersect_venn_plot(intersect_object, file.path(diagnostic_dir, "Venn diagrams of peak calls", paste0(time_point, ".tiff")))
    }
    
    # Get rid of incomplete/alternate scaffolds, since they interfere with annotation later.
    target_regions[[time_point]] = target_regions[[time_point]][!grepl("_", seqnames(target_regions[[time_point]]))]
  }
  
  return(target_regions)
}
