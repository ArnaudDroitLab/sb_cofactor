# setwd("/home/chris/Bureau/sb_cofactor_hr/A549")

library(ComplexHeatmap)
source("scripts/ckn_utils.R")

cofactors <- load_diffbind_cofactors_peaks()
cofactors <- cofactors[grep("_DOWN|_UP", names(cofactors))]
sapply(cofactors, length)

inter_cofactors <- build_intersect(cofactors)
matrix_cofactors <- inter_cofactors$Matrix
sum(matrix_cofactors > 1)
matrix_cofactors[matrix_cofactors > 1] <- 1
sum(matrix_cofactors > 1)
colnames(matrix_cofactors)

m4 <- make_comb_mat(matrix_cofactors, remove_empty_comb_set = TRUE)
m4 <- m4[comb_size(m4) >= 30]
UpSet(m4)
comb_size(m4)

annot_top <- HeatmapAnnotation("Intersection\nsize" = anno_barplot(comb_size(m4), 
                                                                   border = FALSE, gp = gpar(fill = "black"), height = unit(3, "cm")), 
                               annotation_name_side = "left", annotation_name_rot = 0,
                               "Size" = anno_text(comb_size(m4), rot = 0, just = "center", location = 0.25))
annot_right <- rowAnnotation("Set size" = anno_barplot(set_size(m4), 
                                                       border = FALSE, 
                                                       gp = gpar(fill = "black"), 
                                                       width = unit(2, "cm")),
                             "Size" = anno_text(set_size(m4))
)

UpSet(m4, top_annotation = annot_top, right_annotation = annot_right)