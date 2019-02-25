setwd("/home/chris/Bureau/sb_cofactor_hr/A549")
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
# output: 
function(grl)



intersect_NBC_CTRL <- build_intersect(NBC_CTRL)




NIPBL_DEX <- cofactors_peaks[["NIPBL_DEX"]] # 2733
BRD4_DEX <- cofactors_peaks[["BRD4_DEX"]] # 21225
CDK9_DEX <- cofactors_peaks[["CDK9_DEX"]] # 6061



intersect_NBC_CTRL <- build_intersect2(NBC_CTRL)
names(intersect_NBC_CTRL)
intersect_NBC_CTRL$Regions
intersect_NBC_CTRL$Matrix
intersect_NBC_CTRL$List
intersect_NBC_CTRL$Names
intersect_NBC_CTRL$Length

length(intersect_NBC_CTRL$Regions)
nrow(intersect_NBC_CTRL$Matrix)
9470 + 28084 + 10788
M <- as.data.frame(intersect_NBC_CTRL$Matrix); nrow(M)
M_NIPBL <- M[1:9470,]; nrow(M_NIPBL)
M_BRD4 <- M[9471:(9471+28084-1),]; nrow(M_BRD4)
M_CDK9 <- M[(9471+28084):(9471+28084+10788-1),]; nrow(M_CDK9)

# NIPBL point of view
NBC <- sum(M_CDK9$NIPBL_CTRL != 0 & M_CDK9$BRD4_CTRL != 0 & M_CDK9$CDK9_CTRL != 0); print(NBC) # 5476
NB <- sum(M_CDK9$NIPBL_CTRL != 0 & M_CDK9$BRD4_CTRL != 0 & M_CDK9$CDK9_CTRL == 0); print(NB) # 0
NC <- sum(M_CDK9$NIPBL_CTRL != 0 & M_CDK9$BRD4_CTRL == 0 & M_CDK9$CDK9_CTRL != 0); print(NC) # 110
BC <- sum(M_CDK9$NIPBL_CTRL == 0 & M_CDK9$BRD4_CTRL != 0 & M_CDK9$CDK9_CTRL != 0); print(BC) # 4318
N <- sum(M_CDK9$NIPBL_CTRL != 0 & M_CDK9$BRD4_CTRL == 0 & M_CDK9$CDK9_CTRL == 0); print(N) # 0
B <- sum(M_CDK9$NIPBL_CTRL == 0 & M_CDK9$BRD4_CTRL != 0 & M_CDK9$CDK9_CTRL == 0); print(B) # 0
C <- sum(M_CDK9$NIPBL_CTRL == 0 & M_CDK9$BRD4_CTRL == 0 & M_CDK9$CDK9_CTRL != 0); print(C) # 884
total <- NBC + NB + NC + BC + N + B + C; print(total) # 9470

# BRD4 point of view
NBC <- sum(M_BRD4$NIPBL_CTRL != 0 & M_BRD4$BRD4_CTRL != 0 & M_BRD4$CDK9_CTRL != 0); print(NBC) # 5231
NB <- sum(M_BRD4$NIPBL_CTRL != 0 & M_BRD4$BRD4_CTRL != 0 & M_BRD4$CDK9_CTRL == 0); print(NB) # 2507
NC <- sum(M_BRD4$NIPBL_CTRL != 0 & M_BRD4$BRD4_CTRL == 0 & M_BRD4$CDK9_CTRL != 0); print(NC) # 0
BC <- sum(M_BRD4$NIPBL_CTRL == 0 & M_BRD4$BRD4_CTRL != 0 & M_BRD4$CDK9_CTRL != 0); print(BC) # 3805
N <- sum(M_BRD4$NIPBL_CTRL != 0 & M_BRD4$BRD4_CTRL == 0 & M_BRD4$CDK9_CTRL == 0); print(N) # 0
B <- sum(M_BRD4$NIPBL_CTRL == 0 & M_BRD4$BRD4_CTRL != 0 & M_BRD4$CDK9_CTRL == 0); print(B) # 16541
C <- sum(M_BRD4$NIPBL_CTRL == 0 & M_BRD4$BRD4_CTRL == 0 & M_BRD4$CDK9_CTRL != 0); print(C) # 0
total <- NBC + NB + NC + BC + N + B + C; print(total) # 28084

# CDK9 point of view
NBC <- sum(M_BRD4$NIPBL_CTRL != 0 & M_BRD4$BRD4_CTRL != 0 & M_BRD4$CDK9_CTRL != 0); print(NBC) # 5231
NB <- sum(M_BRD4$NIPBL_CTRL != 0 & M_BRD4$BRD4_CTRL != 0 & M_BRD4$CDK9_CTRL == 0); print(NB) # 2507
NC <- sum(M_BRD4$NIPBL_CTRL != 0 & M_BRD4$BRD4_CTRL == 0 & M_BRD4$CDK9_CTRL != 0); print(NC) # 0
BC <- sum(M_BRD4$NIPBL_CTRL == 0 & M_BRD4$BRD4_CTRL != 0 & M_BRD4$CDK9_CTRL != 0); print(BC) # 3805
N <- sum(M_BRD4$NIPBL_CTRL != 0 & M_BRD4$BRD4_CTRL == 0 & M_BRD4$CDK9_CTRL == 0); print(N) # 0
B <- sum(M_BRD4$NIPBL_CTRL == 0 & M_BRD4$BRD4_CTRL != 0 & M_BRD4$CDK9_CTRL == 0); print(B) # 16541
C <- sum(M_BRD4$NIPBL_CTRL == 0 & M_BRD4$BRD4_CTRL == 0 & M_BRD4$CDK9_CTRL != 0); print(C) # 0
total <- NBC + NB + NC + BC + N + B + C; print(total) # 10788

#


###
B <- build_intersect(NBC_CTRL)
regions <- B$Regions
M <- as.data.frame(B$Matrix)

NBC <- (M$NIPBL_CTRL != 0 & M$BRD4_CTRL != 0 & M$CDK9_CTRL != 0); print(sum(NBC)) # 5200
NBC_regions <- regions[NBC,]

NB <- (M$NIPBL_CTRL != 0 & M$BRD4_CTRL != 0 & M$CDK9_CTRL == 0); print(sum(NB)) # 2439
NB_regions <- regions[NB,]

NC <- (M$NIPBL_CTRL != 0 & M$BRD4_CTRL == 0 & M$CDK9_CTRL != 0); print(sum(NC)) # 104
NC_regions <- regions[NC,]

BC <- (M$NIPBL_CTRL == 0 & M$BRD4_CTRL != 0 & M$CDK9_CTRL != 0); print(sum(BC)) # 3683
BC_regions <- regions[BC,]

N <- (M$NIPBL_CTRL != 0 & M$BRD4_CTRL == 0 & M$CDK9_CTRL == 0); print(sum(N)) # 1080
N_regions <- regions[NBC,]

B <- (M$NIPBL_CTRL == 0 & M$BRD4_CTRL != 0 & M$CDK9_CTRL == 0); print(sum(B)) # 16540
B_regions <- regions[NBC,]

C <- (M$NIPBL_CTRL == 0 & M$BRD4_CTRL == 0 & M$CDK9_CTRL != 0); print(sum(C)) # 883
C_regions <- regions[NBC,]

### NBC_regions
NBC_regions
dex_regions <- list("NIPBL_DEX" = NIPBL_DEX, "BRD4_DEX" = BRD4_DEX, "CDK9_DEX" = CDK9_DEX)
# output: matrice, ligne: NBC_regions, colonnes: dex_regions

grl <- dex_regions
all.regions = NBC_regions
overlap.matrix <- matrix(0, nrow = length(all.regions), ncol = length(grl))
overlap.list = list()

for (i in 1:length(grl)) {
  overlap.matrix[, i] <- GenomicRanges::countOverlaps(all.regions, 
                                                      grl[[i]], type = "any")
  overlap.list[[names(grl)[i]]] <- which(overlap.matrix[, i] != 0)
}
colnames(overlap.matrix) <- names(grl)
res <- overlap.matrix

# for (range_names in names(dex_regions)) {
#   data <- dex_regions[[range_names]]
#   message(length(data))
#   subset <- subsetByOverlaps(NBC_regions, data)
#   message(length(subset), " / ",length(NBC_regions))
# }

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
