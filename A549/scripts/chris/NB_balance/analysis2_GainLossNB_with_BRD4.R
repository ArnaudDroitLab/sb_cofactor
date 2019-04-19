# setwd("/Users/chris/Desktop/sb_cofactor_hr/A549")
setwd("/home/chris/Bureau/sb_cofactor_hr/A549")

source("scripts/ckn_utils.R")
source("scripts/load_reddy.R")
library(ChIPseeker)

# Loading peaks
peaks_dir <- "output/chip-pipeline-GRCh38/peak_call/A549_NB"
gainNB <- rtracklayer::import(con = file.path(peaks_dir, "NB_specific_DEX.bed")); print(length(gainNB)) # 803
lossNB <- rtracklayer::import(con = file.path(peaks_dir, "NB_specific_CTRL.bed")); print(length(lossNB)) # 5952
commonNB <- rtracklayer::import(con = file.path(peaks_dir, "NB_common.bed")); print(length(commonNB)) # 1680

### Overlaps with BRD4
BRD4 <- load_cofactor_peaks(cofactors = c("BRD4"))
BRD4_CTRL <- BRD4[["BRD4_CTRL"]]
BRD4_DEX <- BRD4[["BRD4_DEX"]]

gainNB_ovBRD4 <- subsetByOverlaps(gainNB, BRD4_DEX); print(length(gainNB_ovBRD4)) #  ; 803 / 803 = .%
gainNB_notovBRD4 <- gainNB[!(gainNB %in% gainNB_ovBRD4)]; print(length(gainNB_notovBRD4)) # ; 0 / 803 = .%
lossNB_ovBRD4 <- subsetByOverlaps(lossNB, BRD4_CTRL); print(length(lossNB_ovBRD4)) # ; / 5952 = .%
lossNB_notovBRD4 <- lossNB[!(lossNB %in% lossNB_ovBRD4)]; print(length(lossNB_notovBRD4)) #  ; / 5952 = .%






# Annotation
gainNB_ovGR_annodf <- annotatePeaks(gainNB_ovGR, output = "df")
gainNB_notovGR_annodf <- annotatePeaks(gainNB_notovGR, output = "df")
lossNB_ovGR_annodf <- annotatePeaks(lossNB_ovGR, output = "df")
lossNB_notovGR_annodf <- annotatePeaks(lossNB_notovGR, output = "df")
commonNB_ovGR_annodf <- annotatePeaks(commonNB_ovGR, output = "df")
commonNB_notovGR_annodf <- annotatePeaks(commonNB_notovGR, output = "df")

# Retrieve genes which gain or lose NBC at the promoters
geneGainNB_ovGR <- gainNB_ovGR_annodf %>% filter(abs(distanceToTSS) <= 3000) %>% pull(geneId) %>% unique # 20
geneGainNB_notovGR <- gainNB_notovGR_annodf %>% filter(abs(distanceToTSS) <= 3000) %>% pull(geneId) %>% unique # 3
geneLossNB_ovGR <- lossNB_ovGR_annodf %>% filter(abs(distanceToTSS) <= 3000) %>% pull(geneId) %>% unique # 368
geneLossNB_notovGR <- lossNB_notovGR_annodf %>% filter(abs(distanceToTSS) <= 3000) %>% pull(geneId) %>% unique # 1516
geneCommonNB_ovGR <- commonNB_ovGR_annodf %>% filter(abs(distanceToTSS) <= 3000) %>% pull(geneId) %>% unique # 86
geneCommonNB_notovGR <- commonNB_notovGR_annodf %>% filter(abs(distanceToTSS) <= 3000) %>% pull(geneId) %>% unique # 18

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
downDEX_in_gainNB <- downDEX %in% symbol_all_geneGainNB_ovGR; names(downDEX_in_gainNB) <- downDEX
downDEX_in_gainNB

######################
# Draw FC time series
######################
source("scripts/reddy_time_series/draw_graph_log2FC_0-12h.R")

geneGroupList <- list("GainNB_ovGR_withGR" = geneGainNB_ovGR,
                      "GainNB_notovGR_withGR" = geneGainNB_notovGR,
                      "LossNB_withGR" = geneLossNB_ovGR,
                      "LossNB_withoutGR" = geneLossNB_notovGR,
                      "CommonNB_withGR" = geneCommonNB_ovGR,
                      "CommonNB_withoutGR" = geneCommonNB_notovGR)

draw_time_course_FC(geneGainNB_ovGR)
draw_time_course_FC(geneGainNB_notovGR)
draw_time_course_FC(geneLossNB_ovGR)
draw_time_course_FC(geneLossNB_notovGR)
draw_time_course_FC(geneCommonNB_ovGR)
draw_time_course_FC(geneCommonNB_notovGR)
draw_time_course_pergroup_FC(geneGroupList)

# geneLossNBC_ovGR: Action répressive de GR par binding direct
# geneLossNBC_notovGR: Les premières observations ne montre pas de grand changements dans le niveau de fold change de gene expression, réservoir de cofacteurs?

gainNB_ovGR_annodf %>% filter(distanceToTSS > 500000)
