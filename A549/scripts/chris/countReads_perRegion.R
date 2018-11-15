setwd("/home/chris/Bureau/sb_cofactor_hr/A549")
source("scripts/ckn_utils.R")

peaks_NBC_dir <- "output/chip-pipeline-GRCh38/peak_call/A549_NBC"

speNBC_CTRL_path <- file.path(peaks_NBC_dir, "A549_NBC_CTRL_specific.bed")
speNBC_CTRL <- rtracklayer::import(speNBC_CTRL_path)

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
ovCTRLvsDEX <- subsetByOverlaps(NBC_CTRL, NBC_DEX) # 1390
summary(width(ovCTRLvsDEX))
ovDEXvsCTRL <- subsetByOverlaps(NBC_DEX, NBC_CTRL) # 1417
summary(width(ovDEXvsCTRL))

ov <- sort(reduce(c(ovCTRLvsDEX, ovDEXvsCTRL))) # 1376
summary(width(ov))

speCTRL <- sort(setdiff(NBC_CTRL, ovCTRLvsDEX)) # 4176
summary(width(speCTRL))
speDEX <- sort(setdiff(NBC_DEX, ovDEXvsCTRL)) # 808
summary(width(speDEX))


dir.create(peaks_NBC_dir, recursive=TRUE, showWarnings=FALSE)

commonNBC_path <- file.path(peaks_NBC_dir, "A549_NBC_common.bed")

speNBC_DEX_path <- file.path(peaks_NBC_dir, "A549_NBC_DEX_specific.bed")

# rtracklayer::export.bed(ov, con = commonNBC_path, format = "bed")
# rtracklayer::export.bed(speCTRL, con = speNBC_CTRL_path, format = "bed")
# rtracklayer::export.bed(speDEX, con = speNBC_DEX_path, format = "bed")