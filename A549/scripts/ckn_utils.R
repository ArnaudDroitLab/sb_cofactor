# setwd("/home/chris/Bureau/sb_cofactor_hr/A549")

library(GenomicRanges)

load_cofactor_peaks <- function(cofactors = c("NIPBL", "BRD4", "CDK9", "MED1","SMC1A")) {
  peaks_dir <- "output/chip-pipeline-GRCh38/peak_call"
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
  message("#####################################")
  message("Available set of regions: ")
  print(names(cofactors_peaks))
  return(cofactors_peaks)
}

keepStdChr <- function(gr) {
  message("With all chromosomes, including contigs : ", length(gr), " regions")
  stdChr <- paste0("chr", c(seq(1:22), "X", "Y"))
  gr_StdChr <- keepSeqlevels(gr, stdChr[stdChr %in% seqlevels(gr)], pruning.mode = "coarse")
  message("Keeping standard chromosomes : ", length(gr_StdChr), " regions")
  message("\t--> ", length(gr) - length(gr_StdChr), " regions removed")
  return(gr_StdChr)
}
