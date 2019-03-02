# setwd("/Users/chris/Desktop/sb_cofactor_hr/A549")
setwd("/home/chris/Bureau/sb_cofactor_hr/A549")

source("scripts/ckn_utils.R")
library(ChIPseeker)

# Loading peaks
peaks_dir <- "output/chip-pipeline-GRCh38/peak_call/A549_NBC"
bedname <- "NBC_CTRL_to_None_DEX_ovGR_hg38.bed"
peaks <- rtracklayer::import(con = file.path(peaks_dir, bedname)) # 294

# Width
summary(width(peaks))
hist(width(peaks), breaks = 60)

# Annotation
# anno_df <- annotatePeaks(peaks, output = "df")
anno_df2 <- annotatePeaks2(peaks, output = "df")
# anno <- annotatePeaks2(peaks, output = "anno")
anno2 <- annotatePeaks2(peaks, output = "anno")


# Annotation visualization
plotAnnotation(anno_df2)
covplot(peaks)
plotAnnoPie(anno2)
plotAnnoBar(anno2)
plotDistToTSS(anno2)
vennpie(anno2)
upsetplot(anno2)
upsetplot(anno2, vennpie = TRUE)

#
# anno_df %>% filter(abs(distanceToTSS) <= 500)
anno_df2 %>% filter(Annot == "Promoter")

# Action répressive de GR par binding direct