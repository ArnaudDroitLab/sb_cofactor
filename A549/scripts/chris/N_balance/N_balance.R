# setwd("/Users/chris/Desktop/sb_cofactor_hr/A549")
setwd("/home/chris/Bureau/sb_cofactor_hr/A549")

library(eulerr)
source("scripts/ckn_utils.R")

###
getVennRegionsNIPBL <- function(grl) {
  inter <- build_intersect(grl)
  regions <- inter$Regions
  M <- as.data.frame(inter$Matrix)
  
  N_CTRL <- (M$NIPBL_CTRL != 0 & M$NIPBL_DEX == 0); print(sum(N_CTRL)) #
  N_CTRL_regions <- regions[N_CTRL,]
  
  N_DEX <- (M$NIPBL_CTRL == 0 & M$NIPBL_DEX != 0); print(sum(N_DEX)) #
  N_DEX_regions <- regions[N_DEX,]
  
  N_common <- (M$NIPBL_CTRL != 0 & M$NIPBL_DEX != 0); print(sum(N_common)) #
  N_common_regions <- regions[N_common,]
  
  res <- list("N_CTRL" = N_CTRL_regions,
              "N_DEX" = N_DEX_regions,
              "N_common" = N_common_regions)
  return(res)
}

### Split NIPBL
NIPBL <- load_cofactor_peaks(cofactors = c("NIPBL"))
NIPBL_CTRL <- NIPBL[["NIPBL_CTRL"]]
NIPBL_DEX <- NIPBL[["NIPBL_DEX"]]

Nlist <- GenomicRangesList("NIPBL_CTRL" = NIPBL_CTRL, "NIPBL_DEX" = NIPBL_DEX)
plotVenn(Nlist)

regions_NIPBL <- getVennRegionsNIPBL(Nlist)
gainN <- regions_NIPBL[["N_DEX"]] # 890
lossN <- regions_NIPBL[["N_CTRL"]] # 7662
commonN <- regions_NIPBL[["N_common"]] # 1804

### Overlaps with GR
gr_regions <- load_reddy_binding_consensus("NR3C1")

gr_5m_1h <- GRanges()
for (time in names(gr_regions)[2:8]) {
  gr_time <- gr_regions[[time]]
  gr_5m_1h <- append(gr_5m_1h, gr_time)
}

gainN_ovGR <- subsetByOverlaps(gainN, gr_5m_1h); print(length(gainN_ovGR)) # 841 ; 841/890 = 94.49%
gainN_notovGR <- gainN[!(gainN %in% gainN_ovGR)]; print(length(gainN_notovGR)) # 49 ; 49/890 = 5.50%

lossN_ovGR <- subsetByOverlaps(lossN, gr_5m_1h); print(length(lossN_ovGR)) # 4271 ; 4271/7662 = 55.74%
lossN_notovGR <- lossN[!(lossN %in% lossN_ovGR)]; print(length(lossN_notovGR)) # 3391 ; 3391/7662 = 44.26%

commonN_ovGR <- subsetByOverlaps(commonN, gr_5m_1h); print(length(commonN_ovGR)) # 1703 ; 1703/1804 = 94.40%
commonN_notovGR <- commonN[!(commonN %in% commonN_ovGR)]; print(length(commonN_notovGR)) # 101 ; 101/1804 = 5.60%

### Annotation
gainN_ovGR_annodf <- annotatePeaks(gainN_ovGR, output = "df")
gainN_notovGR_annodf <- annotatePeaks(gainN_notovGR, output = "df")
lossN_ovGR_annodf <- annotatePeaks(lossN_ovGR, output = "df")
lossN_notovGR_annodf <- annotatePeaks(lossN_notovGR, output = "df")
commonN_ovGR_annodf <- annotatePeaks(commonN_ovGR, output = "df")
commonN_notovGR_annodf <- annotatePeaks(commonN_notovGR, output = "df")

###
gene_gainN_ovGR <- gainN_ovGR_annodf %>% filter(abs(distanceToTSS) <= 5000) %>% pull(geneId) %>% unique # 23
gene_gainN_notovGR <- gainN_notovGR_annodf %>% filter(abs(distanceToTSS) <= 3000) %>% pull(geneId) %>% unique # 6
gene_lossN_ovGR <- lossN_ovGR_annodf %>% filter(abs(distanceToTSS) <= 3000) %>% pull(geneId) %>% unique # 375
gene_lossN_notovGR <- lossN_notovGR_annodf %>% filter(abs(distanceToTSS) <= 3000) %>% pull(geneId) %>% unique # 1829
gene_commonN_ovGR <- commonN_ovGR_annodf %>% filter(abs(distanceToTSS) <= 3000) %>% pull(geneId) %>% unique # 76
gene_commonN_notovGR <- commonN_notovGR_annodf %>% filter(abs(distanceToTSS) <= 3000) %>% pull(geneId) %>% unique # 25

###

######################
# Draw FC time series
######################
source("scripts/reddy_time_series/draw_graph_log2FC_0-12h.R")

geneGroupList <- list("gene_gainN_ovGR" = gene_gainN_ovGR,
                      "gene_gainN_notovGR" = gene_gainN_notovGR,
                      "gene_lossN_ovGR" = gene_lossN_ovGR,
                      "gene_lossN_notovGR" = gene_lossN_notovGR,
                      "gene_commonN_ovGR" = gene_commonN_ovGR,
                      "gene_commonN_notovGR" = gene_commonN_notovGR)

draw_time_course_FC(gene_gainN_ovGR)
draw_time_course_FC(gene_gainN_notovGR)
draw_time_course_FC(gene_lossN_ovGR)
draw_time_course_FC(gene_lossN_notovGR)
draw_time_course_FC(gene_commonN_ovGR)
draw_time_course_FC(gene_commonN_notovGR)
draw_time_course_pergroup_FC(geneGroupList)

### With BRD4
BRD4 <- load_cofactor_peaks(cofactors = c("BRD4"))
BRD4_CTRL <- BRD4[["BRD4_CTRL"]]
BRD4_DEX <- BRD4[["BRD4_DEX"]]

Blist <- GenomicRangesList("BRD4_CTRL" = BRD4_CTRL, "BRD4_DEX" = BRD4_DEX)
plotVenn(Blist)

### Overlaps with BRD4
gainN_ovBRD4 <- subsetByOverlaps(gainN, BRD4_DEX); print(length(gainN_ovBRD4)) # 852 ; 852/890 = 95.73%
gainN_notovBRD4 <- gainN[!(gainN %in% gainN_ovBRD4)]; print(length(gainN_notovBRD4)) # 38 ; 38/890 = 4.26%

lossN_ovBRD4 <- subsetByOverlaps(lossN, BRD4_CTRL); print(length(lossN_ovBRD4)) # 6523 ; 6523/7662 = 85.13%
lossN_notovBRD4 <- lossN[!(lossN %in% lossN_ovBRD4)]; print(length(lossN_notovBRD4)) # 1136 ; 1136/7662 = 14.82%

# commonN_ovBRD4 <- subsetByOverlaps(commonN, c(BRD4_CTRL, BRD4_DEX)); print(length(commonN_ovBRD4)) # 1772 ; 1772/1804 = 98.22%
# commonN_notovBRD4 <- commonN[!(commonN %in% c(BRD4_CTRL, BRD4_DEX))]; print(length(commonN_notovBRD4)) # 101 ; 101/1804 = 5.60%

### Annotation
gainN_ovBRD4_annodf <- annotatePeaks(gainN_ovBRD4, output = "df")
gainN_notovBRD4_annodf <- annotatePeaks(gainN_notovBRD4, output = "df")
lossN_ovBRD4_annodf <- annotatePeaks(lossN_ovBRD4, output = "df")
lossN_notovBRD4_annodf <- annotatePeaks(lossN_notovBRD4, output = "df")

###
gene_gainN_ovBRD4 <- gainN_ovBRD4_annodf %>% filter(abs(distanceToTSS) <= 3000) %>% pull(geneId) %>% unique # 27
gene_gainN_notovBRD4 <- gainN_notovBRD4_annodf %>% filter(abs(distanceToTSS) <= 3000) %>% pull(geneId) %>% unique # 2
gene_lossN_ovBRD4 <- lossN_ovBRD4_annodf %>% filter(abs(distanceToTSS) <= 3000) %>% pull(geneId) %>% unique # 1855
gene_lossN_notovBRD4 <- lossN_notovBRD4_annodf %>% filter(abs(distanceToTSS) <= 3000) %>% pull(geneId) %>% unique # 305

###
draw_time_course_FC(gene_gainN_ovBRD4)
draw_time_course_FC(gene_gainN_notovBRD4)
draw_time_course_FC(gene_lossN_ovBRD4)
draw_time_course_FC(gene_lossN_notovBRD4)
