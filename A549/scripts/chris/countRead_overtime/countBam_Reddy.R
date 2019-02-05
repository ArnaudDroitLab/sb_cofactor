setwd("/home/chris/Bureau/sb_cofactor_hr/A549")
source("scripts/load_reddy.R")
source("scripts/chris/countRead_overtime/load_Reddy_peaks_consensus.R")
source("scripts/chris/countRead_overtime/countReads.utils.R")
library(Rsamtools)
library(knitr)
library(dplyr)

######################
countInBam <- function(target, indexBam = FALSE) {
  ### List of bam
  report_bam <- make_report_bam(target_name = target, all_chip_bam)
  bamPath <- generate_bamPath_from_report(report_bam, "/home/chris/Bureau/sb_cofactor_hr/A549/input/ENCODE/A549/GRCh38/chip-seq/bam/")
  if (indexBam == TRUE) {indexBam(bamPath)}
  
  ### Gather all peaks all over the time frame
  make_ENCODE_Reddy_ChIP_experiments_file(target)
  target_regions <- load_reddy_binding_consensus(target = target)
  all_regions_reduced <- get_reduced_peaks_from_list(target_regions)
  summary(width(all_regions_reduced))
  peaks_coordVector <- generate_coordVector(all_regions_reduced)
  
  ### Count reads
  count_total <- countRead(all_regions_reduced, peaks_coordVector, bamPath, report_bam)
  
  ### Save
  filename <- paste0("count_total_", target)
  filename_rdata <- paste0(filename, ".RData")
  filename_txt <- paste0(filename, ".txt")
  
  save(count_total, file = file.path(output_path, filename_rdata))
  write.table(count_total, file = file.path(output_path, filename_txt), sep = "\t", quote = FALSE, row.names = FALSE)
}

countInBam_gapbackground <- function(target, indexBam = FALSE) {
  report_bam <- make_report_bam(target_name = target, all_chip_bam)
  bamPath <- generate_bamPath_from_report(report_bam, "/home/chris/Bureau/sb_cofactor_hr/A549/input/ENCODE/A549/GRCh38/chip-seq/bam/")
  if (indexBam == TRUE) {indexBam(bamPath)}
  
  ### Gather all peaks all over the time frame
  make_ENCODE_Reddy_ChIP_experiments_file(target)
  target_regions <- load_reddy_binding_consensus(target = target)
  all_regions_reduced <- get_reduced_peaks_from_list(target_regions)

  gapbackground_regions_reduced <- gaps(all_regions_reduced)
  summary(width(gapbackground_regions_reduced))
  length(gapbackground_regions_reduced)
  peaks_gapbackground_coordVector <- generate_coordVector(gapbackground_regions_reduced)
  
  ### Count reads
  count_total <- countRead(gapbackground_regions_reduced, peaks_gapbackground_coordVector, bamPath, report_bam)
  
  ### Save
  filename <- paste0("count_total_", target, "_gaps_background")
  filename_rdata <- paste0(filename, ".RData")
  filename_txt <- paste0(filename, ".txt")
  
  save(count_total, file = file.path(output_path, filename_rdata))
  write.table(count_total, file = file.path(output_path, filename_txt), sep = "\t", quote = FALSE, row.names = FALSE)
}

######################
all_chip_bam <- ENCODExplorer::queryEncodeGeneric(biosample_name="A549", file_format = "bam", assay="ChIP-seq")
output_path <- "/home/chris/Bureau/sb_cofactor_hr/A549/output/analyses/countTable_overtime"

######################
targets <- c("SMC3", "CTCF", "RAD21", "FOSL2", "BCL3", "JUN", "JUNB", "HES2", "CEBPB")

for (protein in targets) {
  message("#####\t", protein)
  # countInBam(protein, indexBam = TRUE)
  countInBam_gapbackground(protein, indexBam = FALSE)
}
