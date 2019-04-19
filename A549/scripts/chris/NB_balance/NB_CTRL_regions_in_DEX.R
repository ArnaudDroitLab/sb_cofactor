# get OS:http://conjugateprior.org/2015/06/identifying-the-os-from-r/
# setwd("/Users/chris/Desktop/sb_cofactor_hr/A549")
setwd("/home/chris/Bureau/sb_cofactor_hr/A549")

output_dir <- "output/analyses/NB_balance"

source("scripts/ckn_utils.R")
library(ef.utils)
library(GenomicRanges)
library(eulerr)
library(cowplot)

##########
cofactors_peaks <- load_cofactor_stdchr_peaks(cofactors = c("NIPBL", "BRD4"))

NIPBL_CTRL <- cofactors_peaks[["NIPBL_CTRL"]]; print(length(NIPBL_CTRL)) # 9470
BRD4_CTRL <- cofactors_peaks[["BRD4_CTRL"]]; print(length(BRD4_CTRL)) # 28084

NB_CTRL <- GenomicRangesList("NIPBL_CTRL" = NIPBL_CTRL, "BRD4_CTRL" = BRD4_CTRL)
plotVenn(NB_CTRL)

######### getVennRegions
# input: GenomicRangesList
# output: List of GRanges corresponding to each regions of the Venn Diagramm, construted from the GenomicRangesList
getVennRegions <- function(grl) {
  inter <- build_intersect(grl)
  regions <- inter$Regions
  M <- as.data.frame(inter$Matrix)
  
  NB <- (M$NIPBL_CTRL != 0 & M$BRD4_CTRL != 0); print(sum(NB)) # 2439
  NB_CTRL_regions <- regions[NB,]
  
  N <- (M$NIPBL_CTRL != 0 & M$BRD4_CTRL == 0); print(sum(N)) # 1080
  N_CTRL_regions <- regions[N,]
  
  B <- (M$NIPBL_CTRL == 0 & M$BRD4_CTRL != 0); print(sum(B)) # 16540
  B_CTRL_regions <- regions[B,]
  
  
  res <- list("NB_CTRL" = NB_CTRL_regions,
              "N_CTRL" = N_CTRL_regions,
              "B_CTRL" = B_CTRL_regions)
  return(res)
}
#########
vennRegions_NB_CTRL <- getVennRegions(NB_CTRL)
names(vennRegions_NB_CTRL)

######### Construct list of GRanges to compare with
NIPBL_DEX <- cofactors_peaks[["NIPBL_DEX"]]; print(length(NIPBL_DEX)) # 2733
BRD4_DEX <- cofactors_peaks[["BRD4_DEX"]];  print(length(BRD4_DEX)) # 21225
NB_DEX <- list("NIPBL_DEX" = NIPBL_DEX, "BRD4_DEX" = BRD4_DEX)

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
  NB <- sum(resM$NIPBL_DEX != 0 & resM$BRD4_DEX != 0); print(NB) # 136
  N <- sum(resM$NIPBL_DEX != 0 & resM$BRD4_DEX == 0); print(N) # 4
  B <- sum(resM$NIPBL_DEX == 0 & resM$BRD4_DEX != 0); print(B) # 1833
  zero <- sum(resM$NIPBL_DEX == 0 & resM$BRD4_DEX == 0); print(zero) # 421
  total <- NB + N + B + zero; print(total) # 5200
  
  vres <- c(total, NB, N, B, zero)
  return(vres)
}

######### for loop to compare each GRanges of 
df_res <- data.frame()
for (region_name in names(vennRegions_NB_CTRL)) {
  vres <- getTableCount(vennRegions_NB_CTRL[[region_name]], NB_DEX)
  df_res <- rbind(df_res, vres)
}

rownames(df_res) <- names(vennRegions_NB_CTRL)
colnames(df_res) <- c("total", "NB_DEX", "N_DEX", "B_DEX", "None_DEX")
df_res
write.table(df_res, file = file.path(output_dir, "table_NB_CTRL_regions_in_DEX.txt"),
            sep = "\t", quote = FALSE)

df_res_percent <- round(df_res / df_res$total * 100, 2)
df_res_percent
write.table(df_res_percent, file = file.path(output_dir, "table_NB_CTRL_regions_in_DEX_percent.txt"),
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
  NB_regions <- all.regions[resM$NIPBL_DEX != 0 & resM$BRD4_DEX != 0]; print(length(NB_regions))
  N_regions <- all.regions[resM$NIPBL_DEX != 0 & resM$BRD4_DEX == 0]; print(length(N_regions))
  B_regions <- all.regions[resM$NIPBL_DEX == 0 & resM$BRD4_DEX != 0]; print(length(B_regions))
  zero_regions <- all.regions[resM$NIPBL_DEX == 0 & resM$BRD4_DEX == 0]; print(length(zero_regions))
  total_regions <- all.regions; print(length(all.regions)) # 5200
  
  individual_list <- list("total" = total_regions,
                          "NB_DEX" = NB_regions,
                          "N_DEX" =  N_regions,
                          "B_DEX" = B_regions,
                          "None_DEX" = zero_regions)
  return(individual_list)
}

#########
list_res <- list()
for (region_name in names(vennRegions_NB_CTRL)) {
  individual_list <- getRegions(vennRegions_NB_CTRL[[region_name]], NB_DEX)
  aslist_individual_list <- list(individual_list)
  list_res <- c(list_res, aslist_individual_list)
}

names(list_res) <- names(vennRegions_NB_CTRL)

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
    message("    ### ", dex_regions)
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
colnames(df_res_overlapsGR) <- names(list_res[["NB_CTRL"]])
df_res_overlapsGR

# write.table(df_res_overlapsGR, file = file.path(output_dir, "table_NBC_CTRL_withGR1h_regions_in_DEX.txt"),
#             sep = "\t", quote = FALSE)

#### ratio table of GR overlapping
df_res_overlapsGR_percent <- round(df_res_overlapsGR/df_res*100, 2)
df_res_overlapsGR_percent

# write.table(df_res_overlapsGR_percent, file = file.path(output_dir, "table_NBC_CTRL_withGR1h_regions_in_DEX_percent.txt"),
#             sep = "\t", quote = FALSE)

########################################
# Perte de NBC:
# Dans NBC_CTRL qui sont None_DEX, voici le binding de GR:
# A 5 min: 29.83 %
# A 10 min: 36.04 %
# A 15 min: 32.35 %
# A 20 min: 27.51 %
# A 25 min: 23.48 %
# A 30 m: 6.08 %
# A 1h: 2.12 %
########################################
########################################
# Il existe donc deux groupes NBC_CTRL>None_DEX: ceux bindés par GR et ceux non bindés par GR
# L'idée est donc de déterminer ces deux groupes
# step1: retrieve the genomic coordinates of NBC_CTRL>None_DEX
genomic_regions <- list_res[["NB_CTRL"]][["None_DEX"]] # 421

# step2: gather all gr regions from 5min to 1h
gr_5m_1h <- GRanges()
for (time in names(gr_regions)[2:8]) {
  gr_time <- gr_regions[[time]]
  gr_5m_1h <- append(gr_5m_1h, gr_time)
}

# step3a: retrieve the genomic coordinates of NBC_CTRL>None_DEX that are binded by GR at least one time between 5m and 1 hour
ov_gr_regions <- subsetByOverlaps(genomic_regions, gr_5m_1h); print(length(ov_gr_regions)) # 561 ; 561/1467 = 38.24%

# step3b
not_ov_gr_regions <- genomic_regions[!(genomic_regions %in% ov_gr_regions)]; print(length(not_ov_gr_regions)) # 904 ; 904/1467 = 61.62%

# step4: export to bed file for further analysis (hg38)
output_dir <- "output/chip-pipeline-GRCh38/peak_call/A549_NB"
command <- paste("mkdir -p", output_dir, sep = " ")
system(command)

bedname_gr_hg38 <- "NB_CTRL_to_None_DEX_ovGR_hg38.bed"
rtracklayer::export(ov_gr_regions, con = file.path(output_dir, bedname_gr_hg38))

bedname_notgr_hg38 <- "NB_CTRL_to_None_DEX_notovGR_hg38.bed"
rtracklayer::export(not_ov_gr_regions, con = file.path(output_dir, bedname_notgr_hg38))

# step4: convert to hg19 in order to be submitted to GREAT
# chain.hg38tog19 <- import.chain("input/hg38ToHg19.over.chain")
# ov_gr_regions_hg19 <- liftOver(ov_gr_regions, chain.hg38tog19)
# not_ov_gr_regions_hg19 <- liftOver(not_ov_gr_regions, chain.hg38tog19)
# 
# bedname_gr_hg19 <- "NBC_CTRL_to_None_DEX_ovGR_hg19.bed"
# rtracklayer::export(ov_gr_regions_hg19, con = file.path(output_dir, bedname_gr_hg19))
# 
# bedname_notgr_hg19 <- "NBC_CTRL_to_None_DEX_notovGR_hg19.bed"
# rtracklayer::export(not_ov_gr_regions_hg19, con = file.path(output_dir, bedname_notgr_hg19))
