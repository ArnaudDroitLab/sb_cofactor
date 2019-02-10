setwd("/home/chris/Bureau/sb_cofactor_hr/A549")
source("scripts/load_reddy.R")
source("scripts/chris/countRead_overtime/countReads.utils.R")
library(Rsamtools)
library(knitr)
library(dplyr)

######################
# background argument can be: FALSE, background or gaps_background
countInBam <- function(target, indexBam = FALSE, background = FALSE) {
  ### List of bam
  report_bam <- make_report_bam(target_name = target, all_chip_bam)
  
  if (background == "wce") {
    report_wce_bam <- make_report_WCE_bam(report_bam, all_chip_bam)
    bamPath_wce <- generate_bamPath_from_report(report_wce_bam, "/home/chris/Bureau/sb_cofactor_hr/A549/input/ENCODE/A549/GRCh38/chip-seq/bam/")
    if (indexBam == TRUE) {indexBam(bamPath_wce)}
  } else {
    bamPath <- generate_bamPath_from_report(report_bam, "/home/chris/Bureau/sb_cofactor_hr/A549/input/ENCODE/A549/GRCh38/chip-seq/bam/")
    if (indexBam == TRUE) {indexBam(bamPath)}
  }
  
  ### Gather all peaks all over the time frame
  make_ENCODE_Reddy_ChIP_experiments_file(target)
  target_regions <- load_reddy_binding_consensus(target = target)
  all_regions_reduced <- get_reduced_peaks_from_list(target_regions)
  
  if (background == "background") {
    background_regions_reduced <- shift(all_regions_reduced, 10000)
    peaks_background_coordVector <- generate_coordVector(background_regions_reduced)
  } else if (background == "gaps_background") {
    gapbackground_regions_reduced <- gaps(all_regions_reduced)
    peaks_gapbackground_coordVector <- generate_coordVector(gapbackground_regions_reduced)
  } else {
    peaks_coordVector <- generate_coordVector(all_regions_reduced)
  }
  
  ### Count reads
  if (background == "wce") {
    count_total <- countRead(all_regions_reduced, peaks_coordVector, bamPath_wce, report_bam_wce)
  } else if (background == "background") {
    count_total <- countRead(background_regions_reduced, peaks_background_coordVector, bamPath, report_bam)
  } else if (background == "gaps_background") {
    count_total <- countRead(gapbackground_regions_reduced, peaks_gapbackground_coordVector, bamPath, report_bam)
  } else {
    count_total <- countRead(all_regions_reduced, peaks_coordVector, bamPath, report_bam)
  }
  
  ### Save
  filename <- paste0("count_total_", target)
  if (background != FALSE) {
    filename <- paste0(filename, "_", background)
  }
  filename_rdata <- paste0(filename, ".RData")
  filename_txt <- paste0(filename, ".txt")
  
  save(count_total, file = file.path(output_path, filename_rdata))
  write.table(count_total, file = file.path(output_path, filename_txt), sep = "\t", quote = FALSE, row.names = FALSE)
}

######################
all_chip_bam <- ENCODExplorer::queryEncodeGeneric(biosample_name="A549", file_format = "bam", assay="ChIP-seq")
output_path <- "/home/chris/Bureau/sb_cofactor_hr/A549/output/analyses/countTable_overtime"

######################
targets <- c("NR3C1", "EP300", "SMC3", "CTCF", "RAD21", "FOSL2", "BCL3", "JUN", "JUNB", "HES2", "CEBPB")
for (protein in targets) {
  message("#####\t", protein)
  # countInBam(protein, indexBam = FALSE)
  # countInBam(protein, indexBam = FALSE, background = "wce")
  # countInBam(protein, indexBam = FALSE, background = "background")
  # countInBam(protein, indexBam = FALSE, background = "gaps_background")
}
