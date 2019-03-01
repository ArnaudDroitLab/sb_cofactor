# setwd("/Users/chris/Desktop/sb_cofactor_hr/A549")
setwd("/home/chris/Bureau/sb_cofactor_hr/A549")

source("scripts/ckn_utils.R")
library(ChIPseeker)

# Loading peaks
peaks_dir <- "output/chip-pipeline-GRCh38/peak_call/A549_NBC"
bedname <- "NBC_DEX_to_None_CTRL_ovGR.bed"
peaks <- rtracklayer::import(con = file.path(peaks_dir, bedname)) # 294

# Width
summary(width(peaks))
hist(width(peaks), breaks = 60)

# Annotation
anno_df <- annotatePeaks(peaks, output = "df")
anno <- annotatePeaks(peaks, output = "anno")

# Annotation visualization
plotAnnotation(anno_df)
covplot(peaks)
plotAnnoPie(anno)
plotAnnoBar(anno)
plotDistToTSS(anno)
vennpie(anno)
upsetplot(anno)
upsetplot(anno, vennpie = TRUE)

# input to GREAT: Genomic Regions Enrichment of Annotations Tool
# GREAT predicts functions of cis-regulatory regions
  # liftover to hg39
  # 


# tenter le 3D
# liste des TSS


# une fois qu'il y a une assignation NBC to gene, comment evaluer son effet:
  # Reddy POL2 Reddy: import track Reddy to do, just to see
  # Reddy RNA-Seq: +/- DEX
  # ChIP-Seq POL2 shCTRL shNIPBL
  # RNA-Seq POL2 shCTRL shNIPBL