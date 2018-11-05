setwd("/home/chris/Bureau/sb_cofactor_hr/A549")

source("scripts/ckn_utils.R")

cofactors_peaks <- load_cofactor_stdchr_peaks()

NIPBL_CTRL <- cofactors_peaks[["NIPBL_CTRL"]] # 9470
BRD4_CTRL <- cofactors_peaks[["BRD4_CTRL"]] # 28084
CDK9_CTRL <- cofactors_peaks[["CDK9_CTRL"]] # 10788

NIPBL_DEX <- cofactors_peaks[["NIPBL_DEX"]] # 2733
BRD4_DEX <- cofactors_peaks[["BRD4_DEX"]] # 21225
CDK9_DEX <- cofactors_peaks[["CDK9_DEX"]] # 6061

NB_CTRL <- subsetByOverlaps(NIPBL_CTRL, BRD4_CTRL) # 8273
NBC_CTRL <- subsetByOverlaps(NB_CTRL, CDK9_CTRL) # 5566

NB_DEX <- subsetByOverlaps(NIPBL_DEX, BRD4_DEX) # 2634
NBC_DEX <- subsetByOverlaps(NB_DEX, CDK9_DEX) # 2225
############

subsetByOverlaps(NBC_CTRL, NBC_DEX)
subsetByOverlaps(NBC_DEX, NBC_CTRL)



peaks_NBC_dir <- "output/chip-pipeline-GRCh38/peak_call/A549_NBC"
dir.create(peaks_NBC_dir, recursive=TRUE, showWarnings=FALSE)

commonNIPBL_path <- file.path(peaks_NIPBL_dir, "A549_NIPBL_common.bed")
speNIPBL_CTRL_path <- file.path(peaks_NIPBL_dir, "A549_NIPBL_CTRL_specific.bed")
speNIPBL_DEX_path <- file.path(peaks_NIPBL_dir, "A549_NIPBL_DEX_specific.bed")

rtracklayer::export.bed(ov, con = commonNIPBL_path, format = "bed")
rtracklayer::export.bed(speCTRL, con = speNIPBL_CTRL_path, format = "bed")
rtracklayer::export.bed(speDEX, con = speNIPBL_DEX_path, format = "bed")