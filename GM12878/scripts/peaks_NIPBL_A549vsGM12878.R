setwd("/home/chris/Bureau/sb_cofactor_hr/GM12878")
library(GenomicRanges)

###
keepStdChr <- function(gr) {
  message("With all chromosomes, including contigs : ", length(gr), " regions")
  stdChr <- paste0("chr", c(seq(1:22), "X", "Y"))
  gr_StdChr <- keepSeqlevels(gr, stdChr[stdChr %in% seqlevels(gr)], pruning.mode = "coarse")
  message("Keeping standard chromosomes : ", length(gr_StdChr), " regions")
  message("\t--> ", length(gr) - length(gr_StdChr), " regions removed")
  return(gr_StdChr)
}

### All NIPBL peaks are in EtOH condition
# A549
NIPBL_A549_path <- "/home/chris/Bureau/sb_cofactor_hr/A549/output/chip-pipeline-GRCh38/peak_call/A549_CTRL_NIPBL_rep1/A549_CTRL_NIPBL_rep1_peaks.narrowPeak.stdchr.bed"
NIPBL_A549 <- rtracklayer::import(NIPBL_A549_path) # 9470 peaks

# GM12878
NIPBL_GM12878_path <- "output/chip-pipeline-GRCh38/peak_call/NIPBL_rep1/NIPBL_rep1_peaks.narrowPeak.bed"
NIPBL_GM12878_nonStdChr <- rtracklayer::import(NIPBL_GM12878_path) # 14450 peaks
NIPBL_GM12878 <- keepStdChr(NIPBL_GM12878_nonStdChr) # 14392 peaks | 58 removed regions

##### Overlaps
ovAvsGM <- subsetByOverlaps(NIPBL_A549, NIPBL_GM12878) # 1059
summary(width(ovAvsGM)) 
ovGMvsA <- subsetByOverlaps(NIPBL_GM12878, NIPBL_A549) # 1057
summary(width(ovAvsGM)) 

ov <- sort(reduce(c(ovAvsGM, ovGMvsA))) # 1032
summary(width(ov))

speA <- sort(GenomicRanges::setdiff(NIPBL_A549, ovAvsGM)) # 8411
summary(width(speA))
speGM <- sort(GenomicRanges::setdiff(NIPBL_GM12878, ovGMvsA)) # 13335
summary(width(speGM))

peaks_NIPBL_dir <- "output/chip-pipeline-GRCh38/peak_call/A549vsGM12878_NIPBL"
dir.create(peaks_NIPBL_dir, recursive=TRUE, showWarnings=FALSE)

commonNIPBL_path <- file.path(peaks_NIPBL_dir, "A549vsGM12878_NIPBL_common.bed")
speA_path <- file.path(peaks_NIPBL_dir, "A549vsGM12878_NIPBL_A_specific.bed")
speGM_path <- file.path(peaks_NIPBL_dir, "A549vsGM12878_NIPBL_GM_specific.bed")

rtracklayer::export.bed(ov, con = commonNIPBL_path, format = "bed")
rtracklayer::export.bed(speA, con = speA_path, format = "bed")
rtracklayer::export.bed(speGM, con = speGM_path, format = "bed")