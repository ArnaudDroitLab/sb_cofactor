setwd("/home/chris/Bureau/sb_cofactor_hr/A549")
library(ENCODExplorer)
library(dplyr)
library(knitr)

##### Make ENCODE_Reddy_EP300_ChIP_experiments.txt file
all_chip_bed <- ENCODExplorer::queryEncodeGeneric(biosample_name="A549", file_format = "bed", assay="ChIP-seq")
ep300_bed <- all_chip_bed %>% filter(target == "EP300", assembly == "GRCh38", lab == "Tim Reddy, Duke")
report_ep300_bed <- ep300_bed %>% select(accession, file_accession, submitted_by, file_format, target, treatment, treatment_duration, treatment_duration_unit, biological_replicates, controls)
report_ep300_bed$treatment_duration[is.na(report_ep300_bed$treatment_duration)] <- 0
report_ep300_bed[which(report_ep300_bed$submitted_by == "Alejandro Barrera"), ]$treatment_duration_unit[is.na(report_ep300_bed[which(report_ep300_bed$submitted_by == "Alejandro Barrera"), ]$treatment_duration_unit)] <- "minute"
report_ep300_bed[which(report_ep300_bed$submitted_by == "Ian McDowell"), ]$treatment_duration_unit[is.na(report_ep300_bed[which(report_ep300_bed$submitted_by == "Ian McDowell"), ]$treatment_duration_unit)] <- "hour"
report_ep300_bed <- report_ep300_bed %>% arrange(desc(treatment_duration_unit), treatment_duration, submitted_by, biological_replicates)
kable(report_ep300_bed)

exp_tmp <- report_ep300_bed %>% select(accession, treatment_duration, treatment_duration_unit) %>% unique
exp_ep300 <- data.frame(Experiment = exp_tmp$accession, Time = paste0(exp_tmp$treatment_duration, " ", exp_tmp$treatment_duration_unit), Order = 1:nrow(exp_tmp))
kable(exp_ep300)

# write.table(exp_ep300, file = "input/ENCODE_Reddy_EP300_ChIP_experiments.txt", row.names = FALSE, quote = FALSE, sep = "\t")

##### Make load_reddy_ep300_binding_consensus

load_reddy_ep300_binding_consensus <- function(diagnostic_dir=NULL) {
  # Get ENCODE accessions for GR binding time series
  ep300_accession = read.table("input/ENCODE_Reddy_EP300_ChIP_experiments.txt", header=TRUE, sep="\t")
  ep300_accession = ep300_accession[order(ep300_accession$Order),]
  ep300_accession$Time = factor(ep300_accession$Time, levels=ep300_accession$Time)
  
  chip_dir = "input/ENCODE/A549/GRCh38/chip-seq/narrow"
  dir.create(chip_dir, recursive=TRUE, showWarnings=FALSE)
  
  # Loop over all ENCODE accessions, downloading all replicates and 
  # extracting consensus regions.
  ep300_regions = list()
  for(i in 1:nrow(ep300_accession)) {
    encode_accession = ep300_accession$Experiment[i]
    time_point = as.character(ep300_accession$Time[i])
    
    # Download and import peak calls for the time point.
    encodeResults = ENCODExplorer::queryEncodeGeneric(accession=encode_accession, file_type="bed narrowPeak", assembly="GRCh38")
    downloaded_files = ENCODExplorer::downloadEncode(encodeResults, dir=chip_dir, force=FALSE)
    names(downloaded_files) = gsub(".*\\/(.*).bed.gz", "\\1", downloaded_files)
    
    extraCols <- c(signalValue = "numeric", pValue = "numeric", qValue = "numeric", peak = "integer")
    binding_sites = GRangesList(lapply(downloaded_files, rtracklayer::import, extraCols=extraCols))
    
    # Only keep peak calls found in at least two replicates.
    intersect_object = build_intersect(binding_sites)
    two_replicates = rowSums(intersect_object$Matrix) >= 2
    ep300_regions[[time_point]] = intersect_object$Regions[two_replicates]
    if(!is.null(diagnostic_dir)) {
      # Turn off VennDiagram logging.
      futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
      
      dir.create(file.path(diagnostic_dir, "Venn diagrams of peak calls"), recursive=TRUE, showWarnings=FALSE)
      intersect_venn_plot(intersect_object, file.path(diagnostic_dir, "Venn diagrams of peak calls", paste0(time_point, ".tiff")))
    }
    
    # Get rid of incomplete/alternate scaffolds, since they interfere with annotation later.
    ep300_regions[[time_point]] = ep300_regions[[time_point]][!grepl("_", seqnames(ep300_regions[[time_point]]))]
  }
  
  return(ep300_regions)
}