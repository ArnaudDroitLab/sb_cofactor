# setwd("/home/chris/Bureau/sb_cofactor_hr/A549")

library(knitr)
library(DiffBind)
library(GenomicRanges)
library(ComplexHeatmap)
source("scripts/chris/metagene2_Reddy.utils.R")
source("scripts/chris/diffbind_GR_0_25m/diffbind.utils.R")
source("scripts/ckn_utils.R")

##### Gather all down EP300 regions
timepoint <- c("0m", "5m", "10m", "15m", "20m", "25m")
ltp <- length(timepoint)

# explore EP300 diffbind 
EP300_diffbind_downreg <- list()
for (i in 1:(ltp-1)) {
  for (j in (i+1):ltp) {
    for (TF in c(TRUE)) {
      tp1 <- timepoint[i]
      tp2 <- timepoint[j]
      contrast_name <- paste0(tp2, "VS", tp1)
      report <- open_diffBind("EP300", tp1, tp2, pval = TF, reps = "123", output_dir = "output/analyses/EP300_diffbind")
      
      if (!is.null(report)) {
        upreg <- report %>% dplyr::filter(Fold > 0)
        downreg <- report %>% dplyr::filter(Fold < 0)
        if (nrow(downreg) != 0) {
          EP300_diffbind_downreg[[paste0(contrast_name, "_downreg")]] <- makeGRangesFromDataFrame(downreg)
        }
      }
    }
  }
}

names(EP300_diffbind_downreg)
sapply(EP300_diffbind_downreg, length)

# Overlaps Cofactors_down and GR_DOWN
cofactors <- load_diffbind_cofactors_peaks()
cofactors <- cofactors[grep("_DOWN", names(cofactors))]
# cofactors <- append(cofactors, GRangesList("EP300_DOWN" = EP300_DOWN))
names(cofactors)

set_EP300_list_with_cofactors <- c(cofactors, EP300_diffbind_downreg)
names(set_EP300_list_with_cofactors)
sapply(set_EP300_list_with_cofactors, length)

inter_cofactors <- GenomicOperations::GenomicOverlaps(set_EP300_list_with_cofactors)
matrix_cofactors <- inter_cofactors@matrix
sum(matrix_cofactors > 1)
matrix_cofactors[matrix_cofactors > 1] <- 1
sum(matrix_cofactors > 1)
colnames(matrix_cofactors)

m4 <- make_comb_mat(matrix_cofactors, remove_empty_comb_set = TRUE)
m4_plot <- displayUpSet(combMat = m4, threshold = 50)
m4_plot