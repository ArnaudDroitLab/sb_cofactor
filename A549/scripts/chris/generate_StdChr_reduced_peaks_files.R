# setwd("/home/chris/Bureau/sb_cofactor_hr/A549")

source("scripts/ckn_utils.R")

cofactors_stdchr_peaks <- load_cofactor_stdchr_peaks()
peaks_dir <- "output/chip-pipeline-GRCh38/peak_call"

cofactors <- c("NIPBL", "BRD4", "CDK9", "MED1", "SMC1A")

for (cofactor in cofactors) {
  message("#####\t", cofactor)
  ctrl = paste0(cofactor, "_CTRL")
  dex = paste0(cofactor, "_DEX")
  
  cofactor_ctrl_peaks <- cofactors_stdchr_peaks[[ctrl]]
  cat("CTRL :", length(cofactor_ctrl_peaks), "peaks\n")
  cofactor_dex_peaks <- cofactors_stdchr_peaks[[dex]]
  cat("DEX :", length(cofactor_dex_peaks), "peaks\n")
  reduced <- GenomicRanges::reduce(c(cofactor_ctrl_peaks, cofactor_dex_peaks))
  cat("Reduced :", length(reduced), "peaks\n")
  
  foldername <- paste0("A549_", cofactor)
  mkdir_cmd <- paste("mkdir -p", file.path(peaks_dir, foldername), sep = " ")
  cat(mkdir_cmd, "\n")
  system(mkdir_cmd)
  
  bedname <- paste0("A549", "_", cofactor, "_reduced.bed")
  output_path <- file.path(peaks_dir, foldername, bedname)
  rtracklayer::export.bed(reduced, con = output_path, format = "bed")
  cat(output_path, "\n")
}
