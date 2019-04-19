setwd("/home/chris/Bureau/sb_cofactor_hr/mESC")

library(ComplexHeatmap)
library(circlize)

#
col_fun = colorRamp2(c(0, 0.5, 1), c("#0f4259", "white", "#800020"))
# col_fun = colorRamp2(c(0, 0.5, 1), c("#800020", "white", "#0f4259"))

###
get_nth_element <- function(lst, n) {
        sapply(lst, "[", n)
}
###

data <- read.table("scripts/framptongram/comparison_matrix_20190408.txt", header = TRUE)

mat <- as.matrix(data[, 4:ncol(data)])

sample_names_tmp <- as.character(data$NAME)
sample_names_split <- strsplit(sample_names_tmp, "/")
sample_names_mESC <- get_nth_element(sample_names_split, 4)
sample_names <- gsub("mESC_", "", sample_names_mESC)

colnames(mat) <- sample_names
rownames(mat) <- sample_names

Heatmap(mat,
        row_names_side = "left",
        row_dend_side = "right",
        column_names_side = "top",
        column_names_rot = 45,
        column_dend_side = "bottom",
        row_dend_width = unit(50, "mm"),
        column_dend_height = unit(50, "mm"),
        col = col_fun,
        rect_gp = gpar(col = "white", lwd = 1))

#
# Heatmap(mat,
#         row_names_side = "left",
#         row_dend_side = "right",
#         column_names_side = "top",
#         column_names_rot = 45,
#         column_dend_side = "bottom",
#         row_dend_width = unit(50, "mm"),
#         column_dend_height = unit(50, "mm"),
#         col = col_fun,
#         rect_gp = gpar(col = "white", lwd = 5),
#         cell_fun = function(j, i, x, y, width, height, fill) {
#           grid.text(sprintf("%.2f", mat[i, j]), x, y, gp = gpar(fontsize = 10))
#         })
# ?Heatmap
