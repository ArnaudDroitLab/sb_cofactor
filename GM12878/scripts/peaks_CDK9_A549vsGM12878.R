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

### All CDK9 peaks are in EtOH condition
# A549
CDK9_A549_path <- "/home/chris/Bureau/sb_cofactor_hr/A549/output/chip-pipeline-GRCh38/peak_call/A549_CTRL_CDK9_rep1/A549_CTRL_CDK9_rep1_peaks.narrowPeak.stdchr.bed"
CDK9_A549 <- rtracklayer::import(CDK9_A549_path) # 10788 peaks

# GM12878
CDK9_GM12878_path <- "output/chip-pipeline-GRCh38/peak_call/CDK9_rep1/CDK9_rep1_peaks.narrowPeak.bed"
CDK9_GM12878_nonStdChr <- rtracklayer::import(CDK9_GM12878_path) # 22660 peaks
CDK9_GM12878 <- keepStdChr(CDK9_GM12878_nonStdChr) # 22640 peaks | 20 removed regions

##### Overlaps
ovAvsGM <- subsetByOverlaps(CDK9_A549, CDK9_GM12878) # 3205
summary(width(ovAvsGM)) 
ovGMvsA <- subsetByOverlaps(CDK9_GM12878, CDK9_A549) # 3208
summary(width(ovAvsGM)) 

ov <- sort(reduce(c(ovAvsGM, ovGMvsA))) # 3112
summary(width(ov))

speA <- sort(GenomicRanges::setdiff(CDK9_A549, ovAvsGM)) # 7583
summary(width(speA))
speGM <- sort(GenomicRanges::setdiff(CDK9_GM12878, ovGMvsA)) # 19432
summary(width(speGM))

peaks_CDK9_dir <- "output/chip-pipeline-GRCh38/peak_call/A549vsGM12878_CDK9"
dir.create(peaks_CDK9_dir, recursive=TRUE, showWarnings=FALSE)

commonCDK9_path <- file.path(peaks_CDK9_dir, "A549vsGM12878_CDK9_common.bed")
speA_path <- file.path(peaks_CDK9_dir, "A549vsGM12878_CDK9_A_specific.bed")
speGM_path <- file.path(peaks_CDK9_dir, "A549vsGM12878_CDK9_GM_specific.bed")

rtracklayer::export.bed(ov, con = commonCDK9_path, format = "bed")
rtracklayer::export.bed(speA, con = speA_path, format = "bed")
rtracklayer::export.bed(speGM, con = speGM_path, format = "bed")