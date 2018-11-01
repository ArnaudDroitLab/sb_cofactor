# setwd("/home/chris/Bureau/sb_cofactor_hr/A549")

source("scripts/ckn_utils.R")

cofactors_peaks <- load_cofactor_peaks()
cofactors_stdchr_peaks <- lapply(cofactors_peaks, keepStdChr)

peaks_dir <- "output/chip-pipeline-GRCh38/peak_call"
for (name in names(cofactors_stdchr_peaks)) {
  regions <- cofactors_stdchr_peaks[[name]]
  split_name <- strsplit(name, split = "_")
  target <- split_name[[1]][1]
  condition <- split_name[[1]][2]
  
  basename <- paste0("A549", "_", condition, "_", target, "_rep1")
  output_filename <- paste0(basename, "_peaks.narrowPeak", ".stdchr", ".bed")
  output_path <- file.path(peaks_dir, basename, output_filename)
  message(output_path)
  
  rtracklayer::export.bed(regions, con = output_path, format = "bed")
}
