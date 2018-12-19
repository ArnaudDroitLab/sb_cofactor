# to delete: peaks_nipbl_vs_itself.R
setwd("/home/chris/Bureau/sb_cofactor_hr/A549")

source("scripts/ckn_utils.R")

cofactors <- c("NIPBL", "BRD4", "CDK9", "MED1", "SMC1A")

for (cofactor in cofactors) {
  cofactors_StdChr <- load_cofactor_stdchr_peaks(cofactors = cofactor)
  ctrl <- cofactors_StdChr[[paste0(cofactor, "_CTRL")]]
  dex <- cofactors_StdChr[[paste0(cofactor, "_DEX")]]
  
  ovCTRLvsDEX <- subsetByOverlaps(ctrl, dex)
  ovDEXvsCTRL <- subsetByOverlaps(dex, ctrl)
  
  ov <- sort(reduce(c(ovCTRLvsDEX, ovDEXvsCTRL)))
  speCTRL <- sort(setdiff(ctrl, ovCTRLvsDEX))
  speDEX <- sort(setdiff(dex, ovDEXvsCTRL))
  
  peaks_dir <- file.path("output/chip-pipeline-GRCh38/peak_call", paste0("A549_", cofactor))
  dir.create(peaks_dir, recursive=TRUE, showWarnings=FALSE)
  
  common_path <- file.path(peaks_dir, paste0("A549_", cofactor, "_common.bed"))
  speCTRL_path <- file.path(peaks_dir, paste0("A549_", cofactor, "_CTRL_specific.bed"))
  speDEX_path <- file.path(peaks_dir, paste0("A549_", cofactor, "_DEX_specific.bed"))
  
  message("###### Summary ", cofactor)
  message("#\t", "Specific CTRL\t", length(speCTRL), " peaks")
  print(summary(width(speCTRL)))
  message("#\t", "Common\t", length(ov), " peaks")
  print(summary(width(ov)))
  message("#\t", "Specific DEX\t", length(speDEX), " peaks")
  print(summary(width(speDEX)))
  
  rtracklayer::export.bed(ov, con = common_path, format = "bed")
  rtracklayer::export.bed(speCTRL, con = speCTRL_path, format = "bed")
  rtracklayer::export.bed(speDEX, con = speDEX_path, format = "bed")
  
  message("#####################################")
}



