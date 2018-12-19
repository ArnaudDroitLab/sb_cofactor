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

### All MED1 peaks are in EtOH condition
# A549
MED1_A549_path <- "/home/chris/Bureau/sb_cofactor_hr/A549/output/chip-pipeline-GRCh38/peak_call/A549_CTRL_MED1_rep1/A549_CTRL_MED1_rep1_peaks.narrowPeak.stdchr.bed"
MED1_A549 <- rtracklayer::import(MED1_A549_path) # 12668 peaks

# GM12878
MED1_GM12878_path <- "output/chip-pipeline-GRCh38/peak_call/MED1_rep1/MED1_rep1_peaks.narrowPeak"
MED1_GM12878_nonStdChr <- rtracklayer::import(MED1_GM12878_path) # 12041 peaks
MED1_GM12878 <- keepStdChr(MED1_GM12878_nonStdChr) # 12030 peaks | 11 removed regions

##### Overlaps
ovAvsGM <- subsetByOverlaps(MED1_A549, MED1_GM12878) # 1598
summary(width(ovAvsGM)) 
ovGMvsA <- subsetByOverlaps(MED1_GM12878, MED1_A549) # 1578
summary(width(ovAvsGM)) 

ov <- sort(reduce(c(ovAvsGM, ovGMvsA))) # 1523
summary(width(ov))

speA <- sort(GenomicRanges::setdiff(MED1_A549, ovAvsGM)) # 11070
summary(width(speA))
speGM <- sort(GenomicRanges::setdiff(MED1_GM12878, ovGMvsA)) # 10452
summary(width(speGM))

peaks_MED1_dir <- "output/chip-pipeline-GRCh38/peak_call/A549vsGM12878_MED1"
dir.create(peaks_MED1_dir, recursive=TRUE, showWarnings=FALSE)

commonMED1_path <- file.path(peaks_MED1_dir, "A549vsGM12878_MED1_common.bed")
speA_path <- file.path(peaks_MED1_dir, "A549vsGM12878_MED1_A_specific.bed")
speGM_path <- file.path(peaks_MED1_dir, "A549vsGM12878_MED1_GM_specific.bed")

rtracklayer::export.bed(ov, con = commonMED1_path, format = "bed")
rtracklayer::export.bed(speA, con = speA_path, format = "bed")
rtracklayer::export.bed(speGM, con = speGM_path, format = "bed")
