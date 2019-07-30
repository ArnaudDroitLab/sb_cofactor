# setwd("/home/chris/Bureau/sb_cofactor_hr/A549")

library(ComplexHeatmap)
source("scripts/ckn_utils.R")

#####
displayUpSet <- function(combMat, threshold = 1, customSetOrder = FALSE) {
  combMat <- combMat[comb_size(combMat) >= threshold]
  annot_top <- HeatmapAnnotation("Intersection\nsize" = anno_barplot(comb_size(combMat), 
                                                                     border = FALSE,
                                                                     gp = gpar(fill = "black"),
                                                                     height = unit(3, "cm")), 
                                 "Size" = anno_text(comb_size(combMat),
                                                    rot = 0,
                                                    just = "center",
                                                    location = 0.25),
                                 annotation_name_side = "left", annotation_name_rot = 0)
  annot_right <- rowAnnotation("Set size" = anno_barplot(set_size(combMat), 
                                                         border = FALSE, 
                                                         gp = gpar(fill = "black"), 
                                                         width = unit(2, "cm")),
                               "Size" = anno_text(set_size(combMat)))
  
  if (customSetOrder != FALSE) {
    UpSet(combMat, top_annotation = annot_top, right_annotation = annot_right,
          set_order = customSetOrder)
  } else {
    UpSet(combMat, top_annotation = annot_top, right_annotation = annot_right)
  }
}

#####
cofactors <- load_diffbind_cofactors_peaks()
cofactors <- cofactors[grep("NIPBL|BRD4", names(cofactors))]
sapply(cofactors, length)

inter_cofactors <- build_intersect(cofactors)
matrix_cofactors <- inter_cofactors$Matrix
sum(matrix_cofactors > 1)
matrix_cofactors[matrix_cofactors > 1] <- 1
sum(matrix_cofactors > 1)
colnames(matrix_cofactors)

cofactors_names <- c("BRD4", "NIPBL")
cofactors_set_order <- c(paste0(cofactors_names, "_UP"), paste0(cofactors_names, "_UNBIASED"), paste0(cofactors_names, "_DOWN"))
m4 <- make_comb_mat(matrix_cofactors, remove_empty_comb_set = TRUE)

combMat <- displayUpSet(combMat = m4, threshold = 30)