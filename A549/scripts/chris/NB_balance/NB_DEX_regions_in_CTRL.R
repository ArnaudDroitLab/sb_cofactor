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

NIPBL_DEX <- cofactors_peaks[["NIPBL_DEX"]]; print(length(NIPBL_DEX)) # 2733
BRD4_DEX <- cofactors_peaks[["BRD4_DEX"]]; print(length(BRD4_DEX)) # 21225

NB_DEX <- GenomicRangesList("NIPBL_DEX" = NIPBL_DEX, "BRD4_DEX" = BRD4_DEX)
plot_grid(plotVenn(NB_DEX), labels = c("A549_DEX"))

######### getVennRegions
# input: GenomicRangesList
# output: List of GRanges corresponding to each regions of the Venn Diagramm, construted from the GenomicRangesList
getVennRegions <- function(grl) {
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

#########
vennRegions_NB_DEX <- getVennRegions(NB_DEX)
names(vennRegions_NB_DEX)

######### Construct list of GRanges to compare with
NIPBL_CTRL <- cofactors_peaks[["NIPBL_CTRL"]]; print(length(NIPBL_CTRL)) # 9470
BRD4_CTRL <- cofactors_peaks[["BRD4_CTRL"]];  print(length(BRD4_CTRL)) # 28084
NB_CTRL <- list("NIPBL_CTRL" = NIPBL_CTRL, "BRD4_CTRL" = BRD4_CTRL)

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

  resM <- as.data.frame(res); nrow(resM) # 
  NB <- sum(resM$NIPBL_CTRL != 0 & resM$BRD4_CTRL != 0); print(NB) # 
  N <- sum(resM$NIPBL_CTRL != 0 & resM$BRD4_CTRL == 0); print(N) # 
  B <- sum(resM$NIPBL_CTRL == 0 & resM$BRD4_CTRL != 0); print(B) # 
  zero <- sum(resM$NIPBL_CTRL == 0 & resM$BRD4_CTRL == 0); print(zero) # 
  total <- NB + N + B + zero; print(total) # 

  vres <- c(total, NB, N, B, zero)
  return(vres)
}

######### for loop to compare each GRanges of
df_res <- data.frame()
for (region_name in names(vennRegions_NB_DEX)) {
  vres <- getTableCount(vennRegions_NB_DEX[[region_name]], NB_CTRL)
  df_res <- rbind(df_res, vres)
}

rownames(df_res) <- names(vennRegions_NB_DEX)
colnames(df_res) <- c("total", "NB_CTRL", "N_CTRL", "B_CTRL", "None_CTRL")
df_res
write.table(df_res, file = file.path(output_dir, "table_NB_DEX_regions_in_CTRL.txt"),
            sep = "\t", quote = FALSE)

df_res_percent <- round(df_res / df_res$total * 100, 2)
df_res_percent
write.table(df_res_percent, file = file.path(output_dir, "table_NB_DEX_regions_in_CTRL_percent.txt"),
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
  NB_regions <- all.regions[resM$NIPBL_CTRL != 0 & resM$BRD4_CTRL != 0]; print(length(NB_regions))
  N_regions <- all.regions[resM$NIPBL_CTRL != 0 & resM$BRD4_CTRL == 0]; print(length(N_regions))
  B_regions <- all.regions[resM$NIPBL_CTRL == 0 & resM$BRD4_CTRL != 0]; print(length(B_regions))
  zero_regions <- all.regions[resM$NIPBL_CTRL == 0 & resM$BRD4_CTRL == 0]; print(length(zero_regions))
  total_regions <- all.regions; print(length(all.regions)) # 5200

  individual_list <- list("total_CTRL" = total_regions,
                          "NB_CTRL" = NB_regions,
                          "N_CTRL" =  N_regions,
                          "B_CTRL" = B_regions,
                          "None_CTRL" = zero_regions)
  return(individual_list)
}

#########
list_res <- list()
for (region_name in names(vennRegions_NB_DEX)) {
  individual_list <- getRegions(vennRegions_NB_DEX[[region_name]], NB_CTRL)
  aslist_individual_list <- list(individual_list)
  list_res <- c(list_res, aslist_individual_list)
}

names(list_res) <- names(vennRegions_NB_DEX)

# list_res is a list of list containing all GRanges extracted from comparisons: NBC_DEX venn regions and what were their binding of cofactors status in control condition

##########################################
# Overlaps_with_GR
##########################################

### Now, we want to know, what is the binding of these regions with GR
# At 1h:
gr_regions <- load_reddy_binding_consensus("NR3C1")
gr_1h <- gr_regions[["1 hour"]]

df_res_overlapsGR <- data.frame()
for (dex_regions in names(list_res)) {
  message("\n##### ", dex_regions)

  vres <- c()
  for (ctrl_regions in names(list_res[[dex_regions]])) {
    message("    ### ", ctrl_regions)
    initial_regions <- list_res[[dex_regions]][[ctrl_regions]]
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
colnames(df_res_overlapsGR) <- names(list_res[["NB_DEX"]])
df_res_overlapsGR

# write.table(df_res_overlapsGR, file = file.path(output_dir, "table_NB_DEX_withGR1h_regions_in_CTRL.txt"),
#             sep = "\t", quote = FALSE)

#### ratio table of GR overlapping
df_res_overlapsGR_percent <- round(df_res_overlapsGR/df_res*100, 2)
df_res_overlapsGR_percent

# write.table(df_res_overlapsGR_percent, file = file.path(output_dir, "table_NB_DEX_withGR1h_regions_in_CTRL_percent.txt"),
#             sep = "\t", quote = FALSE)

########################################
# Gain de NBC:
# Dans NBC_DEX qui étaient None_CTRL, voici le binding de GR:
# A 5 min: 97.51 %
# A 10 min: 99.25 %
# A 15 min: 99 %
# A 20 min: 98.76 %
# A 25 min:  98.76 %
# A 30 m:  96.27 %
# A 1h:  94.53 %
########################################
########################################
# Retrieve sequences NBC_DEX>None_CTRL that overlaps with GR
# Take the set of sequences that maximizes the overlaps with GR between 0 and 1 hour
########################################
# step1: retrieve the genomic coordinates of NBC_DEX>None_CTRL
genomic_regions <- list_res[["NB_DEX"]][["None_CTRL"]]

# step2: gather all gr regions from 5min to 1h
gr_5m_1h <- GRanges()
for (time in names(gr_regions)[2:8]) {
  gr_time <- gr_regions[[time]]
  gr_5m_1h <- append(gr_5m_1h, gr_time)
}

# step3a: retrieve the genomic coordinades of NBC_DEX>None_CTRL that are binded by GR at least one time between 5m and 1 hour
ov_gr_regions <- subsetByOverlaps(genomic_regions, gr_5m_1h); print(length(ov_gr_regions)) # 399 (as expected)
summary(width(ov_gr_regions))

# step3b
not_ov_gr_regions <- genomic_regions[!(genomic_regions %in% ov_gr_regions)]; print(length(not_ov_gr_regions)) # 3 ; 3/402 = 0.74%

# step4: export to bed file for further analysis (hg38)
output_dir <- "output/chip-pipeline-GRCh38/peak_call/A549_NB"
command <- paste("mkdir -p", output_dir, sep = " ")
system(command)

bednamehg38 <- "NB_DEX_to_None_CTRL_ovGR_hg38.bed"
rtracklayer::export(ov_gr_regions, con = file.path(output_dir,bednamehg38))

bedname_notgr_hg38 <- "NB_DEX_to_None_CTRL_notovGR_hg38.bed"
rtracklayer::export(not_ov_gr_regions, con = file.path(output_dir, bedname_notgr_hg38))

# step4: convert to hg19 in order to be submitted to GREAT
# chain.hg38tog19 <- rtracklayer::import.chain("input/hg38ToHg19.over.chain")
# ov_gr_regions_hg19 <- unlist(rtracklayer::liftOver(ov_gr_regions, chain.hg38tog19)); print(length(ov_gr_regions_hg19)) # 293
# bednamehg19 <- "NBC_DEX_to_None_CTRL_ovGR_hg19.bed"
# rtracklayer::export(ov_gr_regions_hg19, con = file.path(output_dir, bednamehg19))

##### To remove after pull request would be integrated in master
##### Use when rtracklayer::import is not working
BedToGRanges <- function(bedpath) {
    peak_df <- read.table(bedpath)[1:3]
    colnames(peak_df) <- c("seqnames", "start", "end")
    peak <- GRanges(peak_df)
    return(peak)
  }

#####
load_POLR2A_peaks_Myers <- function() {
  peaks_dir <- "output/chip-POLR2A-Myers/peaks"
  ctrl <- "A549_CTRL_POLR2A_ENCFF664KTN.bed"
  dex <- "A549_DEX_POLR2A_ENCFF915LKZ.bed"
    
  peak_pol2_ctrl <- BedToGRanges(file.path(peaks_dir, ctrl))
  peak_pol2_dex <- BedToGRanges(file.path(peaks_dir, dex))
      
  peaks <- GRangesList("POLR2A_CTRL" = peak_pol2_ctrl,
                             "POLR2A_DEX" = peak_pol2_dex)
        
  message("#####################################")
  message("Available set of regions: ")
  print(names(peaks))
  
  return(peaks)
}

# step5: does this regions colocalize with POL2?
POLR2A_peaks <- load_POLR2A_peaks_Myers()
POLR2A_DEX <- POLR2A_peaks[["POLR2A_DEX"]]
ov <- subsetByOverlaps(genomic_regions, POLR2A_DEX)
genomic_regions[!(genomic_regions %in% ov)]
