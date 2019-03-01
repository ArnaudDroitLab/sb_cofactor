# get OS:http://conjugateprior.org/2015/06/identifying-the-os-from-r/
# setwd("/Users/chris/Desktop/sb_cofactor_hr/A549")
setwd("/home/chris/Bureau/sb_cofactor_hr/A549")

output_dir <- "output/analyses/NBC_balance"

source("scripts/ckn_utils.R")
library(ef.utils)
library(GenomicRanges)
library(eulerr)
library(cowplot)

##########
cofactors_peaks <- load_cofactor_stdchr_peaks(cofactors = c("NIPBL", "BRD4", "CDK9"))

NIPBL_CTRL <- cofactors_peaks[["NIPBL_CTRL"]]; print(length(NIPBL_CTRL)) # 9470
BRD4_CTRL <- cofactors_peaks[["BRD4_CTRL"]]; print(length(BRD4_CTRL)) # 28084
CDK9_CTRL <- cofactors_peaks[["CDK9_CTRL"]]; print(length(CDK9_CTRL)) # 10788

NBC_CTRL <- GenomicRangesList("NIPBL_CTRL" = NIPBL_CTRL, "BRD4_CTRL" = BRD4_CTRL, "CDK9_CTRL" = CDK9_CTRL)
plotVenn(NBC_CTRL)

######### getVennRegions
# input: GenomicRangesList
# output: List of GRanges corresponding to each regions of the Venn Diagramm, construted from the GenomicRangesList
getVennRegions <- function(grl) {
  inter <- build_intersect(grl)
  regions <- inter$Regions
  M <- as.data.frame(inter$Matrix)
  
  NBC <- (M$NIPBL_CTRL != 0 & M$BRD4_CTRL != 0 & M$CDK9_CTRL != 0); print(sum(NBC)) # 5200
  NBC_CTRL_regions <- regions[NBC,]
  
  NB <- (M$NIPBL_CTRL != 0 & M$BRD4_CTRL != 0 & M$CDK9_CTRL == 0); print(sum(NB)) # 2439
  NB_CTRL_regions <- regions[NB,]
  
  NC <- (M$NIPBL_CTRL != 0 & M$BRD4_CTRL == 0 & M$CDK9_CTRL != 0); print(sum(NC)) # 104
  NC_CTRL_regions <- regions[NC,]
  
  BC <- (M$NIPBL_CTRL == 0 & M$BRD4_CTRL != 0 & M$CDK9_CTRL != 0); print(sum(BC)) # 3683
  BC_CTRL_regions <- regions[BC,]
  
  N <- (M$NIPBL_CTRL != 0 & M$BRD4_CTRL == 0 & M$CDK9_CTRL == 0); print(sum(N)) # 1080
  N_CTRL_regions <- regions[N,]
  
  B <- (M$NIPBL_CTRL == 0 & M$BRD4_CTRL != 0 & M$CDK9_CTRL == 0); print(sum(B)) # 16540
  B_CTRL_regions <- regions[B,]
  
  C <- (M$NIPBL_CTRL == 0 & M$BRD4_CTRL == 0 & M$CDK9_CTRL != 0); print(sum(C)) # 883
  C_CTRL_regions <- regions[C,]
  
  res <- list("NBC_CTRL" = NBC_CTRL_regions,
              "NB_CTRL" = NB_CTRL_regions,
              "NC_CTRL" = NC_CTRL_regions,
              "BC_CTRL" = BC_CTRL_regions,
              "N_CTRL" = N_CTRL_regions,
              "B_CTRL" = B_CTRL_regions,
              "C_CTRL" = C_CTRL_regions)
  return(res)
}
#########
vennRegions_NBC_CTRL <- getVennRegions(NBC_CTRL)
names(vennRegions_NBC_CTRL)

######### Construct list of GRanges to compare with
NIPBL_DEX <- cofactors_peaks[["NIPBL_DEX"]]; print(length(NIPBL_DEX)) # 2733
BRD4_DEX <- cofactors_peaks[["BRD4_DEX"]];  print(length(BRD4_DEX)) # 21225
CDK9_DEX <- cofactors_peaks[["CDK9_DEX"]];  print(length(CDK9_DEX)) # 6061
NBC_DEX <- list("NIPBL_DEX" = NIPBL_DEX, "BRD4_DEX" = BRD4_DEX, "CDK9_DEX" = CDK9_DEX)

######### getTableCount
# input:  1. GRanges
#         2. List of GRanges to compare against input1
# output: vector of integers, corresponding to the count of genomic ranges in input1 overlapping with each GRanges of input2
getTableCount <- function(regions, region_list) {
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

######### for loop to compare each GRanges of 
df_res <- data.frame()
for (region_name in names(vennRegions_NBC_CTRL)) {
  vres <- getTableCount(vennRegions_NBC_CTRL[[region_name]], NBC_DEX)
  df_res <- rbind(df_res, vres)
}

rownames(df_res) <- names(vennRegions_NBC_CTRL)
colnames(df_res) <- c("total", "NBC_DEX", "NB_DEX", "NC_DEX", "BC_DEX", "N_DEX", "B_DEX", "C_DEX", "None_DEX")
df_res
write.table(df_res, file = file.path(output_dir, "table_NBC_CTRL_regions_in_DEX.txt"),
            sep = "\t", quote = FALSE)

df_res_percent <- round(df_res / df_res$total * 100, 2)
df_res_percent
write.table(df_res_percent, file = file.path(output_dir, "table_NBC_CTRL_regions_in_DEX_percent.txt"),
            sep = "\t", quote = FALSE)

############################
getRegions <- function(regions, region_list) {
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
  NBC_regions <- all.regions[resM$NIPBL_DEX != 0 & resM$BRD4_DEX != 0 & resM$CDK9_DEX != 0]; print(length(NBC_regions))
  NB_regions <- all.regions[resM$NIPBL_DEX != 0 & resM$BRD4_DEX != 0 & resM$CDK9_DEX == 0]; print(length(NB_regions))
  NC_regions <- all.regions[resM$NIPBL_DEX != 0 & resM$BRD4_DEX == 0 & resM$CDK9_DEX != 0]; print(length(NC_regions))
  BC_regions <- all.regions[resM$NIPBL_DEX == 0 & resM$BRD4_DEX != 0 & resM$CDK9_DEX != 0]; print(length(BC_regions))
  N_regions <- all.regions[resM$NIPBL_DEX != 0 & resM$BRD4_DEX == 0 & resM$CDK9_DEX == 0]; print(length(N_regions))
  B_regions <- all.regions[resM$NIPBL_DEX == 0 & resM$BRD4_DEX != 0 & resM$CDK9_DEX == 0]; print(length(B_regions))
  C_regions <- all.regions[resM$NIPBL_DEX == 0 & resM$BRD4_DEX == 0 & resM$CDK9_DEX != 0]; print(length(C_regions))
  zero_regions <- all.regions[resM$NIPBL_DEX == 0 & resM$BRD4_DEX == 0 & resM$CDK9_DEX == 0]; print(length(zero_regions))
  total_regions <- all.regions; print(length(all.regions)) # 5200
  
  individual_list <- list("total_DEX" = total_regions,
                          "NBC_DEX" = NBC_regions,
                          "NB_DEX" = NB_regions,
                          "NC_DEX" = NC_regions,
                          "BC_DEX" =  BC_regions,
                          "N_DEX" =  N_regions,
                          "B_DEX" = B_regions,
                          "C_DEX" = C_regions,
                          "None_DEX" = zero_regions)
  return(individual_list)
}

#########
list_res <- list()
for (region_name in names(vennRegions_NBC_CTRL)) {
  individual_list <- getRegions(vennRegions_NBC_CTRL[[region_name]], NBC_DEX)
  aslist_individual_list <- list(individual_list)
  list_res <- c(list_res, aslist_individual_list)
}

names(list_res) <- names(vennRegions_NBC_CTRL)

##########################################
# Overlaps_with_GR
##########################################

### Now, we want to know, what is the binding of these regions with GR
# At 1h:
gr_regions <- load_reddy_binding_consensus("NR3C1")
gr_1h <- gr_regions[["1 hour"]]

df_res_overlapsGR <- data.frame()
for (ctrl_regions in names(list_res)) {
  message("\n##### ", ctrl_regions)
  
  vres <- c()
  for (dex_regions in names(list_res[[ctrl_regions]])) {
    message("    ### ", ctrl_regions)
    initial_regions <- list_res[[ctrl_regions]][[dex_regions]]
    cat("\t", length(initial_regions), "regions \n")
    ov_with_GR <- subsetByOverlaps(initial_regions, gr_1h)
    cat("\t", length(ov_with_GR), "regions overlapping with GR\n")
    
    nb_regions <- length(initial_regions)
    nb_regions_withGR <- length(ov_with_GR)
    ratio <- round(nb_regions_withGR/nb_regions*100, 2)
    cat("\t", ratio, " %\n")
    
    vres <- c(vres, nb_regions_withGR)
  }
  df_res_overlapsGR <- rbind(df_res_overlapsGR, vres)
}

df_res_overlapsGR
rownames(df_res_overlapsGR) <- names(list_res)
colnames(df_res_overlapsGR) <- names(list_res[["NBC_CTRL"]])
df_res_overlapsGR

write.table(df_res_overlapsGR, file = file.path(output_dir, "table_NBC_CTRL_withGR_regions_in_DEX.txt"),
            sep = "\t", quote = FALSE)

#### ratio table of GR overlapping
df_res_overlapsGR_percent <- round(df_res_overlapsGR/df_res*100, 2)
df_res_overlapsGR_percent

write.table(df_res_overlapsGR_percent, file = file.path(output_dir, "table_NBC_CTRL_withGR_regions_in_DEX_percent.txt"),
            sep = "\t", quote = FALSE)

########################################
#
########################################
