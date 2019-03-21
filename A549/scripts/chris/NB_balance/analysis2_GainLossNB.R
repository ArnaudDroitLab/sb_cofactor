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

# Overlaps with GR
gr_regions <- load_reddy_binding_consensus("NR3C1")
gr_5m_1h <- GRanges()
for (time in names(gr_regions)[2:8]) {
  gr_time <- gr_regions[[time]]
  gr_5m_1h <- append(gr_5m_1h, gr_time)
}

gainNB_ovGR <- subsetByOverlaps(gainNB, gr_5m_1h); print(length(gainNB_ovGR)) # 791 ; 791/803 = 98.51%
gainNB_notovGR <- gainNB[!(gainNB %in% gainNB_ovGR)]; print(length(gainNB_notovGR)) # 12 ; 12/803 = 1.49%

lossNB_ovGR <- subsetByOverlaps(lossNB, gr_5m_1h); print(length(lossNB_ovGR)) # 3600 ; 3600/5952 = 60.48%
lossNB_notovGR <- lossNB[!(lossNB %in% lossNB_ovGR)]; print(length(lossNB_notovGR)) # 2352 ; 2352/803 = 2.93%

commonNB_ovGR <- subsetByOverlaps(commonNB, gr_5m_1h); print(length(commonNB_ovGR)) # 1632 ; 1632/1680 = 97.14%
commonNB_notovGR <- commonNB[!(commonNB %in% commonNB_ovGR)]; print(length(commonNB_notovGR)) # 48 ; 48/1680 = 2.85%

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
geneGainNB_ovGR <- gainNB_ovGR_annodf %>% filter(Annot %in% c("Promoter")) %>% pull(geneId) %>% unique
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
downDEX_in_gainNB <- downDEX %in% symbol_all_geneGainNB_ovGR; names(downDEX_in_gainNB) <- downDEX
downDEX_in_gainNB

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
