setwd("/home/chris/Bureau/sb_cofactor_hr/A549")

library(GenomicRanges)

peaks_dir <- "output/chip-pipeline-GRCh38/peak_call"

load_cofactor_peaks <- function(cofactors) {
  cofactors_peaks <- GRangesList()
  name_cofactors_peaks <- c()
  for (cofactor in cofactors) {
    for (condition in c("CTRL", "DEX")) {
      message("####\t", cofactor, " | ", condition)
      basename <- paste0("A549_", condition, "_", cofactor, "_rep1")
      peaks_path <- file.path(peaks_dir, basename, paste0(basename, "_peaks.narrowPeak.bed"))
      message(peaks_path)
      peaks <- rtracklayer::import(peaks_path)
      message("Number of regions : ", length(peaks))
      cofactors_peaks <- append(cofactors_peaks, GRangesList(peaks))
      name_cofactors_peaks <- c(name_cofactors_peaks, paste0(cofactor, "_", condition))
    }
  }
  names(cofactors_peaks) <- name_cofactors_peaks
  return(cofactors_peaks)
}

## test
cofactors <- c("NIPBL", "BRD4", "CDK9", "SMC1A", "MED1")
p <- load_cofactor_peaks(cofactors)
