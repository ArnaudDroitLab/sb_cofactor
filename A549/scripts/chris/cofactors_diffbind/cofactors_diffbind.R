# setwd("/home/chris/Bureau/sb_cofactor_hr/A549")

library(ComplexHeatmap)
library(GenomicOperations)
source("scripts/ckn_utils.R")

#####
cofactors <- load_diffbind_cofactors_peaks()
cofactors <- cofactors[grep("NIPBL|BRD4", names(cofactors))]
sapply(cofactors, length)

inter_cofactors <- GenomicOperations::GenomicOverlaps(cofactors)
matrix_cofactors <- inter_cofactors@matrix
sum(matrix_cofactors > 1)
matrix_cofactors[matrix_cofactors > 1] <- 1
sum(matrix_cofactors > 1)
colnames(matrix_cofactors)

cofactors_names <- c("BRD4", "NIPBL")
cofactors_set_order <- c(paste0(cofactors_names, "_UP"), paste0(cofactors_names, "_UNBIASED"), paste0(cofactors_names, "_DOWN"))
m4 <- make_comb_mat(matrix_cofactors, remove_empty_comb_set = TRUE)

upset_brd4_nipbl <- displayUpSet(combMat = m4, threshold = 30)
