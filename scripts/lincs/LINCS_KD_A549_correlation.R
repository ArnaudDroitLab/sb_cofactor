# setwd("/home/chris/Bureau/sb_cofactor_hr")
setwd("/Users/chris/Desktop/sb_cofactor_hr")

library(ComplexHeatmap)
library(circlize)
source("scripts/lincs/lincs.utils.R")

# Load fold change gene expression values
A549_LINCs_KD_matrix <- readRDS("output/analysis/lincs/A549_LINCs_KD_matrix.rds")

# Remove gene name in order to build a matrix
rownames(A549_LINCs_KD_matrix) <- A549_LINCs_KD_matrix$Name_GeneSymbol
A549_LINCs_KD_matrix <- A549_LINCs_KD_matrix %>% select(-Name_GeneSymbol)

# As matrix
mat <- as.matrix(A549_LINCs_KD_matrix)
max(mat)
min(mat)

# Load mutated cofactors names
MUTATED_COFACTORS <- get_mutated_cofactors()
# MUTCOF_A549 is the list of available KD of mutated cofactors in A549 LINCS
MUTCOF_A549 <- MUTATED_COFACTORS[MUTATED_COFACTORS %in% colnames(A549_LINCs_KD_matrix)]

#######
# Pearson correlation matrix between each KD sample
cor.pearson <- cor(mat, method = "pearson")
max(cor.pearson - diag(nrow(cor.pearson)))
min(cor.pearson)

#######
# Which KD are the closer from my mutated cofactors in term of gene signature?
col_fun = colorRamp2(c(-2, 0, 2), c("#0f4259", "white", "#800020")) # Define colors for heatmap
output <- "output/analysis/lincs/closeToMutatedCofactors"
dir.create(output)
date <- "20190510"

for (COF in MUTCOF_A549) {
  closeToCOF <- sort(cor.pearson[COF, ], decreasing = TRUE)[1:20]
  matCOF <-  as.matrix(A549_LINCs_KD_matrix[, names(closeToCOF)])
  
  heatmapCOF <- Heatmap(matCOF, name = "LogDiffExp",
          show_row_names = FALSE,
          show_row_dend = FALSE,
          column_names_side = "top",
          column_names_rot = 45,
          cluster_columns = FALSE,
          col = col_fun,
          bottom_annotation = HeatmapAnnotation(correlation = anno_text(round(closeToCOF,3),
                                                                        rot = 0,
                                                                        just = "center")))
  
  pdf_name <- file.path(output, paste0(COF, "_top20closerKD_", date, ".pdf"))
  pdf(file = pdf_name, width = 10, height = 5)
  print(heatmapCOF)
  dev.off()
}

#######
# Pearson correlation matrix between each KD sample
cor.pearson.mutcof <- cor.pearson[MUTCOF_A549, MUTCOF_A549]
max(cor.pearson.mutcof - diag(nrow(cor.pearson.mutcof)))
min(cor.pearson.mutcof)

col_fun2 = colorRamp2(c(-1, 0, 1), c("#0f4259", "white", "#800020"))
Heatmap(cor.pearson.mutcof, name = "Pearson correlation",
        row_names_side = "left",
        row_dend_side = "right",
        column_names_side = "top",
        column_dend_side = "bottom",
        row_dend_width = unit(30, "mm"),
        column_dend_height = unit(30, "mm"),
        column_dend_reorder = TRUE,
        col = col_fun2,
        rect_gp = gpar(col = "white", lwd = 1),
        cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(sprintf("%.2f", cor.pearson.mutcof[i, j]), x, y, gp = gpar(fontsize = 10))
        })

#######
# Is there a set of genes highly correlated in KD samples?
mat_mutcof <- mat[, MUTCOF_A549]
cor.pearson.mutcof.genes <- cor(t(mat_mutcof), method = "pearson")
max(cor.pearson.mutcof.genes - diag(nrow(cor.pearson.mutcof.genes)))
min(cor.pearson.mutcof.genes)

Heatmap(cor.pearson.mutcof.genes,
        col = col_fun2)