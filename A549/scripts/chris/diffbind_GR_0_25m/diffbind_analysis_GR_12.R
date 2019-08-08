# setwd("/home/chris/Bureau/sb_cofactor_hr/A549")

library(knitr)
library(DiffBind)
library(GenomicRanges)
library(ComplexHeatmap)
library(GenomicOperations)
source("scripts/chris/metagene2_Reddy.utils.R")
source("scripts/chris/diffbind_GR_0_25m/diffbind.utils.R")
source("scripts/ckn_utils.R")

##### Gather all down GR regions
timepoint <- c("0m", "5m", "10m", "15m", "20m", "25m")
ltp <- length(timepoint)

GR_diffbind_downreg <- list()
for (i in 1:(ltp-1)) {
  for (j in (i+1):ltp) {
    for (TF in c(TRUE)) {
      tp1 <- timepoint[i]
      tp2 <- timepoint[j]
      contrast_name <- paste0(tp2, "VS", tp1)
      report <- open_diffBind("GR", tp1, tp2, pval = TF, reps = "12", output_dir = "output/analyses/GR_diffbind")
      
      if (!is.null(report)) {
        upreg <- report %>% dplyr::filter(Fold > 0)
        downreg <- report %>% dplyr::filter(Fold < 0)
        if (nrow(downreg) != 0) {
          GR_diffbind_downreg[[paste0(contrast_name, "_downreg")]] <- makeGRangesFromDataFrame(downreg)
        }
      }
    }
  }
}

names(GR_diffbind_downreg)
sapply(GR_diffbind_downreg, length)
saveRDS(GR_diffbind_downreg, file = "output/analyses/annotate_peaks_with_hic/GR_diffbind_downreg_pval.rds")

# EP300_DOWN between 0h and 1h
EP300_DOWN <- open_diffBind("EP300", tp1 = "0h", tp2 = "1h", pval = TRUE, reps = "123", output_dir = "output/analyses/EP300_diffbind") %>%
  dplyr::filter(Fold < 0) %>% makeGRangesFromDataFrame

# Overlaps Cofactors_down and GR_DOWN
cofactors <- load_diffbind_cofactors_peaks()
cofactors <- cofactors[grep("_DOWN", names(cofactors))]
cofactors <- append(cofactors, GRangesList("EP300_DOWN" = EP300_DOWN))
names(cofactors)

# List of GRanges regions that will be display in the upsetplot
set_GR_list_with_cofactors <- c(cofactors, GR_diffbind_downreg)
names(set_GR_list_with_cofactors)
sapply(set_GR_list_with_cofactors, length)

# Build the intersection matrix
inter_cofactors <- GenomicOperations::GenomicOverlaps(set_GR_list_with_cofactors)
matrix_cofactors <- inter_cofactors@matrix
sum(matrix_cofactors > 1)
# build_intersect count the number of overlaps, we only need the absence/presence in that purpose
matrix_cofactors[matrix_cofactors > 1] <- 1
sum(matrix_cofactors > 1)
colnames(matrix_cofactors)

# UpSet plot
m4 <- make_comb_mat(matrix_cofactors, remove_empty_comb_set = TRUE)
m4_10_plot <- displayUpSet(m4, threshold = 10)
m4_10_plot

# Extract intersection regions
cs <- comb_size(m4)
de <- comb_degree(m4)
set_GR_list <- list()
for (id in names(cs)){
  i <- extract_comb(m4, id)
  setname <- idToName(id, colnames(matrix_cofactors))
  message("##### ", id, " | degree ", de[id], " | ", setname)
  regions <- inter_cofactors$Regions[i]
  
  set_GR_list[[setname]] <- regions
}

names(set_GR_list)
sapply(set_GR_list, length)

# Metagene plot to control the intersection set
# for (i in c(35, 30, 36, 28, 38, 33, 31, 27, 21)) {
#   contrast <- names(set_GR_list)[i]
#   message("##### ", contrast)
#   df_set_GR_list <- make_df_metagene_Reddy(chip_target = c("GR", "EP300"), peaks = set_GR_list[[i]], merge_replicates = TRUE, reps = "12")
#   title_group <- paste(contrast, paste(length(set_GR_list[[i]]), "regions"), sep = " | ")
#   p <- plot_metagene_Reddy(df_set_GR_list, title = title_group)
#   saveMetagene(metagene_plot = p,
#                output_dir = "output/analyses/GR_diffbind_downreg_setlist_2",
#                output_file = paste("GR_diffbind_downreg_GR_EP300", contrast, "reps12", sep = "_"),
#                width = 20, height = 9)
# }

# Example code to access the coordinates of the regions of the intersection set
# t <- set_GR_list[["NIPBL+BRD4+CDK9+MED1+SMC1A+25mVS10m"]] %>% as.data.frame
# t$seqnames <- as.character(t$seqnames)
# t <- t %>% arrange(seqnames) %>%
#   mutate(coord = paste0(seqnames, ":", start, "-", end),
#          index = 1:length(coord))
# kable(t)
