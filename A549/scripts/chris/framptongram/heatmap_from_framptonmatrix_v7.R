setwd("/home/chris/Bureau/sb_cofactor_hr/A549")

library(ComplexHeatmap)
library(circlize)

#
col_fun = colorRamp2(c(-1, 0, 1), c("#0f4259", "white", "#800020"))
# col_fun = colorRamp2(c(0, 1), c("#0f4259", "#800020"))

###
get_nth_element <- function(lst, n) {
  sapply(lst, "[", n)
}

###
data <- read.table("scripts/chris/framptongram/comparison_matrix_20190430_v7.txt", header = TRUE)

mat <- as.matrix(data[, 4:ncol(data)])

sample_names_tmp <- as.character(data$NAME)
sample_names_split <- strsplit(sample_names_tmp, "/")
sample_names_mESC <- get_nth_element(sample_names_split, 4)
sample_names <- gsub("mESC_", "", sample_names_mESC)
sample_names_split <- strsplit(sample_names, "_")
sample_names <- get_nth_element(sample_names_split, 1)

colnames(mat) <- sample_names
rownames(mat) <- sample_names

# h <- Heatmap(mat, name = "Correlation",
#         row_names_side = "left",
#         row_names_gp = gpar(fontsize = 11),
#         row_dend_side = "right",
#         column_names_side = "top",
#         column_names_gp = gpar(fontsize = 11),
#         column_names_rot = 45,
#         column_dend_side = "bottom",
#         row_dend_width = unit(50, "mm"),
#         column_dend_height = unit(50, "mm"),
#         col = col_fun,
#         rect_gp = gpar(col = "white", lwd = 1),
#         cell_fun = function(j, i, x, y, width, height, fill) {
#         grid.text(sprintf("%.2f", mat[i, j]), x, y, gp = gpar(fontsize = 7))
#         })

# pdf("L1000_KD_16prot_20190417.pdf")
h2 <- Heatmap(mat,
              row_names_side = "left",
              row_names_gp = gpar(fontsize = 11),
              row_dend_side = "right",
              column_names_side = "top",
              column_names_gp = gpar(fontsize = 11),
              column_names_rot = 45,
              column_dend_side = "bottom",
              show_column_dend = FALSE,
              row_dend_width = unit(50, "mm"),
              column_dend_height = unit(50, "mm"),
              col = col_fun,
              rect_gp = gpar(col = "white", lwd = 1))
h2
# dev.off()


# Add lignée cellulaire sur le cote

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
