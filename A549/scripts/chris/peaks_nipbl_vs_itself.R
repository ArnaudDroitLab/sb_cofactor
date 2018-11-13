setwd("/home/chris/Bureau/sb_cofactor_hr/A549")

source("scripts/ckn_utils.R")

cofactors_peaks <- load_cofactor_peaks()
cofactors_StdChr <- lapply(cofactors_peaks, keepStdChr)

NIPBL_CTRL <- cofactors_StdChr[["NIPBL_CTRL"]] # 9470
summary(width(NIPBL_CTRL))
NIPBL_DEX <- cofactors_StdChr[["NIPBL_DEX"]] # 2733
summary(width(NIPBL_DEX))

ovCTRLvsDEX <- subsetByOverlaps(NIPBL_CTRL, NIPBL_DEX) # 1819
summary(width(ovCTRLvsDEX))
ovDEXvsCTRL <- subsetByOverlaps(NIPBL_DEX, NIPBL_CTRL) # 1850
summary(width(ovDEXvsCTRL))

ov <- sort(reduce(c(ovCTRLvsDEX, ovDEXvsCTRL))) # 1802
summary(width(ov))

speCTRL <- sort(setdiff(NIPBL_CTRL, ovCTRLvsDEX)) # 7651
summary(width(speCTRL))
speDEX <- sort(setdiff(NIPBL_DEX, ovDEXvsCTRL)) # 883
summary(width(speDEX))

peaks_NIPBL_dir <- "output/chip-pipeline-GRCh38/peak_call/A549_NIPBL"
dir.create(peaks_NIPBL_dir, recursive=TRUE, showWarnings=FALSE)

commonNIPBL_path <- file.path(peaks_NIPBL_dir, "A549_NIPBL_common.bed")
speNIPBL_CTRL_path <- file.path(peaks_NIPBL_dir, "A549_NIPBL_CTRL_specific.bed")
speNIPBL_DEX_path <- file.path(peaks_NIPBL_dir, "A549_NIPBL_DEX_specific.bed")

rtracklayer::export.bed(ov, con = commonNIPBL_path, format = "bed")
rtracklayer::export.bed(speCTRL, con = speNIPBL_CTRL_path, format = "bed")
rtracklayer::export.bed(speDEX, con = speNIPBL_DEX_path, format = "bed")