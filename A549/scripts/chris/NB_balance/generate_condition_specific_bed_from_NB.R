# setwd("/Users/chris/Desktop/sb_cofactor_hr/A549")
setwd("/home/chris/Bureau/sb_cofactor_hr/A549")

library(eulerr)
source("scripts/ckn_utils.R")

######### getVennRegionsDEX
# input: GenomicRangesList
# output: List of GRanges corresponding to each regions of the Venn Diagramm, construted from the GenomicRangesList
getVennRegionsDEX <- function(grl) {
  inter <- build_intersect(grl)
  regions <- inter$Regions
  M <- as.data.frame(inter$Matrix)
  
  NB <- (M$NIPBL_DEX != 0 & M$BRD4_DEX != 0); print(sum(NB)) #
  NB_DEX_regions <- regions[NB,]
  
  N <- (M$NIPBL_DEX != 0 & M$BRD4_DEX == 0); print(sum(N)) #
  N_DEX_regions <- regions[N,]
  
  B <- (M$NIPBL_DEX == 0 & M$BRD4_DEX != 0); print(sum(B)) #
  B_DEX_regions <- regions[B,]
  
  res <- list("NB_DEX" = NB_DEX_regions,
              "N_DEX" = N_DEX_regions,
              "B_DEX" = B_DEX_regions)
  return(res)
}

######### getVennRegionsCTRL
# input: GenomicRangesList
# output: List of GRanges corresponding to each regions of the Venn Diagramm, construted from the GenomicRangesList
getVennRegionsCTRL <- function(grl) {
  inter <- build_intersect(grl)
  regions <- inter$Regions
  M <- as.data.frame(inter$Matrix)
  
  NB <- (M$NIPBL_CTRL != 0 & M$BRD4_CTRL != 0); print(sum(NB)) #
  NB_CTRL_regions <- regions[NB,]
  
  N <- (M$NIPBL_CTRL != 0 & M$BRD4_CTRL == 0); print(sum(N)) #
  N_CTRL_regions <- regions[N,]
  
  B <- (M$NIPBL_CTRL == 0 & M$BRD4_CTRL != 0); print(sum(B)) #
  B_CTRL_regions <- regions[B,]
  
  res <- list("NB_CTRL" = NB_CTRL_regions,
              "N_CTRL" = N_CTRL_regions,
              "B_CTRL" = B_CTRL_regions)
  return(res)
}

######### getVennRegionsCTRLDEX
# input: GenomicRangesList
# output: List of GRanges corresponding to each regions of the Venn Diagramm, construted from the GenomicRangesList
getVennRegionsCTRLDEX <- function(grl) {
  inter <- build_intersect(grl)
  regions <- inter$Regions
  M <- as.data.frame(inter$Matrix)
  
  NB_common <- (M$"NIPBL-BRD4_CTRL" != 0 & M$"NIPBL-BRD4_DEX" != 0); print(sum(NB_common)) #
  NB_common_regions <- regions[NB_common,]
  
  NB_specific_CTRL <- (M$"NIPBL-BRD4_CTRL" != 0 & M$"NIPBL-BRD4_DEX" == 0); print(sum(NB_specific_CTRL)) #
  NB_specific_CTRL_regions <- regions[NB_specific_CTRL,]
  
  NB_specific_DEX <- (M$"NIPBL-BRD4_CTRL" == 0 & M$"NIPBL-BRD4_DEX" != 0); print(sum(NB_specific_DEX)) #
  NB_specific_DEX_regions <- regions[NB_specific_DEX,]
  
  res <- list("NB_common" = NB_common_regions,
              "NB_specific_CTRL" = NB_specific_CTRL_regions,
              "NB_specific_DEX" = NB_specific_DEX_regions)
  return(res)
}

##########
cofactors_peaks <- load_cofactor_stdchr_peaks(cofactors = c("NIPBL", "BRD4"))

NIPBL_CTRL <- cofactors_peaks[["NIPBL_CTRL"]]; print(length(NIPBL_CTRL)) # 9470
BRD4_CTRL <- cofactors_peaks[["BRD4_CTRL"]]; print(length(BRD4_CTRL)) # 28084
NBlist_CTRL <- GenomicRangesList("NIPBL_CTRL" = NIPBL_CTRL, "BRD4_CTRL" = BRD4_CTRL)
plotVenn(NBlist_CTRL) # NB_CTRL 7652

NIPBL_DEX <- cofactors_peaks[["NIPBL_DEX"]]; print(length(NIPBL_DEX)) # 2733
BRD4_DEX <- cofactors_peaks[["BRD4_DEX"]]; print(length(BRD4_DEX)) # 21225
NBlist_DEX <- GenomicRangesList("NIPBL_DEX" = NIPBL_DEX, "BRD4_DEX" = BRD4_DEX)
plotVenn(NBlist_DEX) # NB_DEX 2512

##########
vennRegions_NBlist_CTRL <- getVennRegionsCTRL(NBlist_CTRL)
NB_CTRL <- vennRegions_NBlist_CTRL[["NB_CTRL"]]

vennRegions_NBlist_DEX <- getVennRegionsDEX(NBlist_DEX)
NB_DEX <- vennRegions_NBlist_DEX[["NB_DEX"]]

NBlist_CTRL_DEX <- GenomicRangesList("NIPBL-BRD4_CTRL" = NB_CTRL, "NIPBL-BRD4_DEX" = NB_DEX)
plotVenn(NBlist_CTRL_DEX)

vennRegions_NBlist_CTRL_DEX <- getVennRegionsCTRLDEX(NBlist_CTRL_DEX)
NB_specific_ctrl <- vennRegions_NBlist_CTRL_DEX[["NB_specific_CTRL"]]; print(length(NB_specific_ctrl))
NB_common <- vennRegions_NBlist_CTRL_DEX[["NB_common"]]; print(length(NB_common))
NB_specific_dex <- vennRegions_NBlist_CTRL_DEX[["NB_specific_DEX"]]; print(length(NB_specific_dex))
summary(width(NB_specific_ctrl))
summary(width(NB_common))
summary(width(NB_specific_dex))

########## Create bed files for further analysis
output_dir <- "output/chip-pipeline-GRCh38/peak_call/A549_NB"
command <- paste("mkdir -p", output_dir, sep = " ")
system(command)

rtracklayer::export(NB_specific_ctrl, con = file.path(output_dir, "NB_specific_CTRL.bed"))
rtracklayer::export(NB_common, con = file.path(output_dir, "NB_common.bed"))
rtracklayer::export(NB_specific_dex, con = file.path(output_dir, "NB_specific_DEX.bed"))