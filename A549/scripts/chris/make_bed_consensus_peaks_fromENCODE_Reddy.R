setwd("/home/chris/Bureau/sb_cofactor_hr/A549")
source("scripts/load_reddy.R")

library(rtracklayer)

##### Aim: Generate consensus peaks from Reddy data (A549 + DEX): 0-25m and 0-12h
# in order to display it on tracks
##### Output: bed files for consensus peaks at each time point for the following targets:

##### Function to format time_point for filename
# Remove spaces
# Transform "minutes" to "m" and "hour" to "h"
format_timepoint_name <- function(time_point) {
  elements <- strsplit(time_point, split = " ")[[1]]
  number <- elements[1]
  unit <- elements[2]
  
  if (unit == "minute") {
    res <- paste0(number, "m")
  } else {
    res <- paste0(number, "h")
  }
  
  return(res)
}

###### Function to make bed file representing consensus peaks from replicates
# From Reddy data
# At several time points
dir_export <- "input/ENCODE/A549/GRCh38/chip-seq/narrow"
make_consensus_bed <- function(target_name) {
  make_ENCODE_Reddy_ChIP_experiments_file(target_name)
  peaks <- load_reddy_binding_consensus(target_name)
  
  for (time_point in names(peaks)) {
    if (length(peaks[[time_point]]) != 0) {
      bed_filename <- paste(target_name, format_timepoint_name(time_point), "consensus_peaks.bed", sep = "_")
      output_path <- file.path(dir_export, bed_filename)
      rtracklayer::export.bed(peaks[[time_point]], con = output_path)
      message("Saved in ", output_path)
    }
  }
}

##### Available at 0-25m and 0-12h
# make_consensus_bed("NR3C1")
# make_consensus_bed("EP300")
# make_consensus_bed("H3K27ac")
# make_consensus_bed("JUN")

##### Available at 0-12h
# make_consensus_bed("CEBPB")
# make_consensus_bed("BCL3")
# make_consensus_bed("FOSL2")
# make_consensus_bed("HES2")
# make_consensus_bed("CTCF") # some timepoint have only   one replicate so does not have a file
# make_consensus_bed("JUNB")
# make_consensus_bed("SMC3")
# make_consensus_bed("RAD21")
# make_consensus_bed("H3K4me1")
# make_consensus_bed("H3K4me2")
# make_consensus_bed("H3K4me3")
# make_consensus_bed("H3K9me3") # not dowloaded yet because it's a broadPeak file
