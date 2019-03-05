# setwd("/Users/chris/Desktop/sb_cofactor_hr/A549")
setwd("/home/chris/Bureau/sb_cofactor_hr/A549")

source("scripts/ckn_utils.R")
library(ChIPseeker)

# Loading peaks
peaks_dir <- "output/chip-pipeline-GRCh38/peak_call/A549_NBC"
gainNBC <- rtracklayer::import(con = file.path(peaks_dir, "NBC_DEX_to_None_CTRL_ovGR_hg38.bed")); print(length(gainNBC)) # 294
lossNBC_ovGR <- rtracklayer::import(con = file.path(peaks_dir, "NBC_CTRL_to_None_DEX_ovGR_hg38.bed")); print(length(lossNBC_ovGR)) # 205
lossNBC_notovGR <- rtracklayer::import(con = file.path(peaks_dir, "NBC_CTRL_to_None_DEX_notovGR_hg38.bed")); print(length(lossNBC_notovGR)) # 216

# Width
summary(width(gainNBC))
hist(width(gainNBC), breaks = 60)

summary(width(lossNBC_ovGR))
hist(width(lossNBC_ovGR), breaks = 60)

summary(width(lossNBC_notovGR))
hist(width(lossNBC_notovGR), breaks = 60)

# Annotation
gainNBC_annodf2 <- annotatePeaks2(gainNBC, output = "df")
lossNBC_ovGR_annodf2 <- annotatePeaks2(lossNBC_ovGR, output = "df")
lossNBC_notovGR_annodf2 <- annotatePeaks2(lossNBC_notovGR, output = "df")

# Retrieve genes which gain or lose NBC at the promoters
geneGainNBC <- gainNBC_annodf2 %>% filter(Annot == "Promoter") %>% select(ENSEMBL) %>% unique
geneLossNBC_ovGR <- lossNBC_ovGR_annodf2 %>% filter(Annot == "Promoter") %>% select(ENSEMBL) %>% unique
geneLossNBC_notovGR <- lossNBC_notovGR_annodf2 %>% filter(Annot == "Promoter") %>% select(ENSEMBL) %>% unique

######################
# Draw FC time series
######################
source("scripts/reddy_time_series/draw_graph_log2FC_0-12h.R")

geneGroupList <- list("GainNBC" = geneGainNBC$ENSEMBL,
                      "LossNBC_withGR" = geneLossNBC_ovGR$ENSEMBL,
                      "LossNBC_withoutGR" = geneLossNBC_notovGR$ENSEMBL)

draw_time_course_FC(geneGainNBC$ENSEMBL)
draw_time_course_FC(geneLossNBC_ovGR$ENSEMBL)
draw_time_course_FC(geneLossNBC_notovGR$ENSEMBL)
draw_time_course_pergroup_FC(geneGroupList)

# geneLossNBC_ovGR: Action répressive de GR par binding direct
# geneLossNBC_notovGR: Les premières observations ne montre pas de grand changements dans le niveau de fold change de gene expression, réservoir de cofacteurs?
