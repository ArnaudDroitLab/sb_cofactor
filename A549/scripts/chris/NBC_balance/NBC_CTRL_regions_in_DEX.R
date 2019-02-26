# get OS:http://conjugateprior.org/2015/06/identifying-the-os-from-r/
setwd("/home/chris/Bureau/sb_cofactor_hr/A549")
setwd("/Users/chris/Desktop/sb_cofactor_hr/A549")
source("scripts/ckn_utils.R")
library(ef.utils)
library(GenomicRanges)

##########
cofactors_peaks <- load_cofactor_stdchr_peaks(cofactors = c("NIPBL", "BRD4", "CDK9"))

NIPBL_CTRL <- cofactors_peaks[["NIPBL_CTRL"]]; print(length(NIPBL_CTRL)) # 9470
BRD4_CTRL <- cofactors_peaks[["BRD4_CTRL"]]; print(length(BRD4_CTRL)) # 28084
CDK9_CTRL <- cofactors_peaks[["CDK9_CTRL"]]; print(length(CDK9_CTRL)) # 10788

NBC_CTRL <- GenomicRangesList("NIPBL_CTRL" = NIPBL_CTRL, "BRD4_CTRL" = BRD4_CTRL, "CDK9_CTRL" = CDK9_CTRL)
plotVenn(NBC_CTRL)

#########
# input: GenomicRangesList
# output: regions of each overlaps
noname <- function(grl) {
  inter <- build_intersect(grl)
  regions <- inter$Regions
  M <- as.data.frame(inter$Matrix)
  
  NBC <- (M$NIPBL_CTRL != 0 & M$BRD4_CTRL != 0 & M$CDK9_CTRL != 0); print(sum(NBC)) # 5200
  NBC_regions <- regions[NBC,]
  
  NB <- (M$NIPBL_CTRL != 0 & M$BRD4_CTRL != 0 & M$CDK9_CTRL == 0); print(sum(NB)) # 2439
  NB_regions <- regions[NB,]
  
  NC <- (M$NIPBL_CTRL != 0 & M$BRD4_CTRL == 0 & M$CDK9_CTRL != 0); print(sum(NC)) # 104
  NC_regions <- regions[NC,]
  
  BC <- (M$NIPBL_CTRL == 0 & M$BRD4_CTRL != 0 & M$CDK9_CTRL != 0); print(sum(BC)) # 3683
  BC_regions <- regions[BC,]
  
  N <- (M$NIPBL_CTRL != 0 & M$BRD4_CTRL == 0 & M$CDK9_CTRL == 0); print(sum(N)) # 1080
  N_regions <- regions[N,]
  
  B <- (M$NIPBL_CTRL == 0 & M$BRD4_CTRL != 0 & M$CDK9_CTRL == 0); print(sum(B)) # 16540
  B_regions <- regions[B,]
  
  C <- (M$NIPBL_CTRL == 0 & M$BRD4_CTRL == 0 & M$CDK9_CTRL != 0); print(sum(C)) # 883
  C_regions <- regions[C,]
  
  res <- list("NBC" = NBC_regions,
              "NB" = NB_regions,
              "NC" = NC_regions,
              "BC" = BC_regions,
              "N" = N_regions,
              "B" = B_regions,
              "C" = C_regions)
  return(res)
}

test <- noname(NBC_CTRL)
names(test)

NIPBL_DEX <- cofactors_peaks[["NIPBL_DEX"]] # 2733
BRD4_DEX <- cofactors_peaks[["BRD4_DEX"]] # 21225
CDK9_DEX <- cofactors_peaks[["CDK9_DEX"]] # 6061
NBC_DEX <- list("NIPBL_DEX" = NIPBL_DEX, "BRD4_DEX" = BRD4_DEX, "CDK9_DEX" = CDK9_DEX)

count_vs <- function(regions, region_list) {
  grl <- region_list
  all.regions = regions
  overlap.matrix <- matrix(0, nrow = length(all.regions), ncol = length(grl))
  overlap.list = list()
  
  for (i in 1:length(grl)) {
    overlap.matrix[, i] <- GenomicRanges::countOverlaps(all.regions, 
                                                        grl[[i]], type = "any")
    overlap.list[[names(grl)[i]]] <- which(overlap.matrix[, i] != 0)
  }
  colnames(overlap.matrix) <- names(grl)
  res <- overlap.matrix
  
  resM <- as.data.frame(res); nrow(resM) # 5200
  NBC <- sum(resM$NIPBL_DEX != 0 & resM$BRD4_DEX != 0 & resM$CDK9_DEX != 0); print(NBC) # 1343
  NB <- sum(resM$NIPBL_DEX != 0 & resM$BRD4_DEX != 0 & resM$CDK9_DEX == 0); print(NB) # 136
  NC <- sum(resM$NIPBL_DEX != 0 & resM$BRD4_DEX == 0 & resM$CDK9_DEX != 0); print(NC) # 3
  BC <- sum(resM$NIPBL_DEX == 0 & resM$BRD4_DEX != 0 & resM$CDK9_DEX != 0); print(BC) # 1404
  N <- sum(resM$NIPBL_DEX != 0 & resM$BRD4_DEX == 0 & resM$CDK9_DEX == 0); print(N) # 4
  B <- sum(resM$NIPBL_DEX == 0 & resM$BRD4_DEX != 0 & resM$CDK9_DEX == 0); print(B) # 1833
  C <- sum(resM$NIPBL_DEX == 0 & resM$BRD4_DEX == 0 & resM$CDK9_DEX != 0); print(C) # 56
  zero <- sum(resM$NIPBL_DEX == 0 & resM$BRD4_DEX == 0 & resM$CDK9_DEX == 0); print(zero) # 421
  total <- NBC + NB + NC + BC + N + B + C + zero; print(total) # 5200
  
  vres <- c(total, NBC, NB, NC, BC, N, B, C, zero)
  return(vres)
}

df_res <- data.frame()
for (region_name in names(test)) {
  vres <- count_vs(test[[region_name]], NBC_DEX)
  df_res <- rbind(df_res, vres)
}

rownames(df_res) <- names(test)
colnames(df_res) <- c("total", "NBC", "NB", "NC", "BC", "N", "B", "C", "None")
df_res

NBC_DEX2 <- GenomicRangesList("NIPBL_DEX" = NIPBL_DEX, "BRD4_DEX" = BRD4_DEX, "CDK9_DEX" = CDK9_DEX)
plotVenn(NBC_DEX2)
