# setwd("/Users/chris/Desktop/sb_cofactor_hr/A549")
setwd("/home/chris/Bureau/sb_cofactor_hr/A549")

source("scripts/ckn_utils.R")
library(ChIPseeker)

# Loading peaks
peaks_dir <- "output/chip-pipeline-GRCh38/peak_call/A549_NB"
gainNB_ovGR <- rtracklayer::import(con = file.path(peaks_dir, "NB_DEX_to_None_CTRL_ovGR_hg38.bed")); print(length(gainNB_ovGR)) # 399
gainNB_notovGR <- rtracklayer::import(con = file.path(peaks_dir, "NB_DEX_to_None_CTRL_notovGR_hg38.bed")); print(length(gainNB_notovGR)) # 3
lossNB_ovGR <- rtracklayer::import(con = file.path(peaks_dir, "NB_CTRL_to_None_DEX_ovGR_hg38.bed")); print(length(lossNB_ovGR)) # 561
lossNB_notovGR <- rtracklayer::import(con = file.path(peaks_dir, "NB_CTRL_to_None_DEX_notovGR_hg38.bed")); print(length(lossNB_notovGR)) # 904

# Width
summary(width(gainNB_ovGR))
hist(width(gainNB_ovGR), breaks = 60)

summary(width(gainNB_notovGR))
hist(width(gainNB_notovGR), breaks = 60)

summary(width(lossNB_ovGR))
hist(width(lossNB_ovGR), breaks = 60)

summary(width(lossNB_notovGR))
hist(width(lossNB_notovGR), breaks = 60)

# Annotation
gainNB_ovGR_annodf <- annotatePeaks(gainNB_ovGR, output = "df")
gainNB_notovGR_annodf <- annotatePeaks(gainNB_notovGR, output = "df")
lossNB_ovGR_annodf <- annotatePeaks(lossNB_ovGR, output = "df")
lossNB_notovGR_annodf <- annotatePeaks(lossNB_notovGR, output = "df")

# Retrieve genes which gain or lose NBC at the promoters
geneGainNB_ovGR <- gainNB_ovGR_annodf %>% filter(Annot == "Promoter") %>% pull(geneId) %>% unique
geneGainNB_notovGR <- gainNB_notovGR_annodf %>% filter(Annot == "Promoter") %>% pull(geneId) %>% unique
geneLossNB_ovGR <- lossNB_ovGR_annodf %>% filter(Annot == "Promoter") %>% pull(geneId) %>% unique
geneLossNB_notovGR <- lossNB_notovGR_annodf %>% filter(Annot == "Promoter") %>% pull(geneId) %>% unique

symbol_all_geneGainNB_ovGR <- gainNB_ovGR_annodf %>% pull(SYMBOL) %>% unique
symbol_all_geneGainNB_notovGR <- gainNB_notovGR_annodf %>% pull(SYMBOL) %>% unique
symbol_all_geneLossNB_ovGR <- lossNB_ovGR_annodf %>% pull(SYMBOL) %>% unique
symbol_all_geneLossNB_notovGR <- lossNB_notovGR_annodf %>% pull(SYMBOL) %>% unique

symbol_prom_geneGainNB_ovGR <- gainNB_ovGR_annodf %>% filter(Annot == "Promoter") %>% pull(SYMBOL) %>% unique
symbol_prom_geneGainNB_notovGR <- gainNB_notovGR_annodf %>% filter(Annot == "Promoter") %>% pull(SYMBOL) %>% unique
symbol_prom_geneLossNB_ovGR <- lossNB_ovGR_annodf %>% filter(Annot == "Promoter") %>% pull(SYMBOL) %>% unique
symbol_prom_geneLossNB_notovGR <- lossNB_notovGR_annodf %>% filter(Annot == "Promoter") %>% pull(SYMBOL) %>% unique

#########
upDEX <- c("PER1", "ZFP36", "ERRFI1", "ANGPTL4", "NR1D2", "CRY2")
upDEX_in_gainNB <- upDEX %in% symbol_all_geneGainNB_ovGR; names(upDEX_in_gainNB) <- upDEX
upDEX_in_gainNB

downDEX <- c("IL11")

######################
# Draw FC time series
######################
source("scripts/reddy_time_series/draw_graph_log2FC_0-12h.R")

geneGroupList <- list("GainNB_ovGR_withGR" = geneGainNB_ovGR,
                      "GainNB_notovGR_withGR" = geneGainNB_notovGR,
                      "LossNB_withGR" = geneLossNB_ovGR,
                      "LossNB_withoutGR" = geneLossNB_notovGR)

draw_time_course_FC(geneGainNB_ovGR)
draw_time_course_FC(gainNB_ovGR_annodf %>% pull(geneId) %>% unique)
draw_time_course_FC(geneGainNB_notovGR)
draw_time_course_FC(geneLossNB_ovGR)
draw_time_course_FC(geneLossNB_notovGR)
draw_time_course_pergroup_FC(geneGroupList)

# geneLossNBC_ovGR: Action répressive de GR par binding direct
# geneLossNBC_notovGR: Les premières observations ne montre pas de grand changements dans le niveau de fold change de gene expression, réservoir de cofacteurs?

gainNB_ovGR_annodf %>% filter(distanceToTSS > 500000)
