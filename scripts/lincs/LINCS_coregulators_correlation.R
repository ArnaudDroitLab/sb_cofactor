## DEPRECATED

# setwd("/home/chris/Bureau/sb_cofactor_hr")
setwd("/Users/chris/Desktop/sb_cofactor_hr")

library(ComplexHeatmap)
library(circlize)
source("scripts/lincs/lincs.utils.R")

###
get_nth_element <- function(lst, n) {
  sapply(lst, "[", n)
}
###

# Load
A549 <- readRDS("output/analysis/lincs/LINCS_KD_matrix/A549_96h_LINCS_KD_matrix.rds")
A375 <- readRDS("output/analysis/lincs/LINCS_KD_matrix/A375_96h_LINCS_KD_matrix.rds")
ASC <- readRDS("output/analysis/lincs/LINCS_KD_matrix/ASC_96h_LINCS_KD_matrix.rds")
HA1E <- readRDS("output/analysis/lincs/LINCS_KD_matrix/HA1E_96h_LINCS_KD_matrix.rds")
HEKTE <- readRDS("output/analysis/lincs/LINCS_KD_matrix/HEKTE_96h_LINCS_KD_matrix.rds")

colnames(A549)[2:length(colnames(A549))] <- paste0("A549_", colnames(A549)[2:length(colnames(A549))])
colnames(A375)[2:length(colnames(A375))] <- paste0("A375_", colnames(A375)[2:length(colnames(A375))])
colnames(ASC)[2:length(colnames(ASC))] <- paste0("ASC_", colnames(ASC)[2:length(colnames(ASC))])
colnames(HA1E)[2:length(colnames(HA1E))] <- paste0("HA1E_", colnames(HA1E)[2:length(colnames(HA1E))])
colnames(HEKTE)[2:length(colnames(HEKTE))] <- paste0("HEKTE_", colnames(HEKTE)[2:length(colnames(HEKTE))])

allCellLine <- list(A549, A375, ASC, HA1E, HEKTE)
allCellLine_mat <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "Name_GeneSymbol", all.x = TRUE), allCellLine)
allCellLine_mat <- allCellLine_mat %>% select(-Name_GeneSymbol)

# As matrix
mat <- as.matrix(allCellLine_mat)
max(mat)
min(mat)

# Load mutated cofactors names
MUTATED_COFACTORS <- get_mutated_cofactors()
# MUTCOF_A549 is the list of available KD of mutated cofactors in A549 LINCS
MUTCOF_A549 <- paste0("A549_", MUTATED_COFACTORS[paste0("A549_", MUTATED_COFACTORS) %in% colnames(mat)])
# MUTCOF_A375 is the list of available KD of mutated cofactors in A375 LINCS
MUTCOF_A375 <- paste0("A375_", MUTATED_COFACTORS[paste0("A375_", MUTATED_COFACTORS) %in% colnames(mat)])
# MUTCOF_ASC is the list of available KD of mutated cofactors in A375 LINCS
MUTCOF_ASC <- paste0("ASC_", MUTATED_COFACTORS[paste0("ASC_", MUTATED_COFACTORS) %in% colnames(mat)])
# MUTCOF_HA1E is the list of available KD of mutated cofactors in HA1E LINCS
MUTCOF_HA1E <- paste0("HA1E_", MUTATED_COFACTORS[paste0("HA1E_", MUTATED_COFACTORS) %in% colnames(mat)])
# MUTCOF_HEKTE is the list of available KD of mutated cofactors in HEKTE LINCS
MUTCOF_HEKTE <- paste0("HEKTE_", MUTATED_COFACTORS[paste0("HEKTE_", MUTATED_COFACTORS) %in% colnames(mat)])
# All MUTCOF in all cell lines
MUTCOF <- c(MUTCOF_A549, MUTCOF_A375, MUTCOF_ASC, MUTCOF_HA1E, MUTCOF_HEKTE)

# New matrix
mat_mutcof <- mat[, MUTCOF]

#######
# Pearson correlation matrix between each KD sample (mutated cofactors)
cor.pearson.mutcof <- cor(mat_mutcof, method = "pearson")
max(cor.pearson.mutcof - diag(nrow(cor.pearson.mutcof))) # -diag(n) allow to remove the identity matrix
min(cor.pearson.mutcof)

# correlation color scale
col_fun2 = colorRamp2(c(-1, 0, 1), c("#0f4259", "white", "#800020"))

# annotation per cell lines
col_celllines_tmp <- strsplit(colnames(cor.pearson.mutcof), "_")
col_celllines <- get_nth_element(col_celllines_tmp, 1)
ha <- HeatmapAnnotation(CellLine = col_celllines,
                        col = list(CellLine = c("HA1E" = "#A4036F", "A549" = "#048BA8", "ASC" = "#16DB93", "A375" = "#EFEA5A", "HEKTE" = "#F29E4C")))
rowha = rowAnnotation(CellLine = col_celllines,
                      col = list(CellLine = c("HA1E" = "#A4036F", "A549" = "#048BA8", "ASC" = "#16DB93", "A375" = "#EFEA5A", "HEKTE" = "#F29E4C")),
                      show_legend = FALSE)
# heatmap
heatmap_MUTCOF <- Heatmap(cor.pearson.mutcof, name = "Pearson correlation",
                          row_names_side = "left",
                          row_dend_side = "right",
                          column_names_side = "top",
                          column_names_rot = 45,
                          column_dend_side = "bottom",
                          row_dend_width = unit(30, "mm"),
                          column_dend_height = unit(30, "mm"),
                          column_dend_reorder = TRUE,
                          col = col_fun2,
                          rect_gp = gpar(col = "white", lwd = 1),
                          top_annotation = ha,
                          left_annotation = rowha)
pdf(file = "output/analysis/lincs/heatmap_MUTCOF_simplepearson_20190515.pdf", width = 25, height = 22)
print(heatmap_MUTCOF)
dev.off()

#######
# Pearson correlation matrix between each KD sample (mutated cofactors) with absolute values
cor.pearson.mutcof.abs <- cor(abs(mat_mutcof), method = "pearson")
max(cor.pearson.mutcof.abs - diag(nrow(cor.pearson.mutcof.abs))) # -diag(n) allow to remove the identity matrix
min(cor.pearson.mutcof.abs)

col_fun2 = colorRamp2(c(-1, 0, 1), c("#0f4259", "white", "#800020"))
heatmap_MUTCOF <- Heatmap(cor.pearson.mutcof.abs, name = "Pearson correlation",
                          row_names_side = "left",
                          row_dend_side = "right",
                          column_names_side = "top",
                          column_names_rot = 45,
                          column_dend_side = "bottom",
                          row_dend_width = unit(30, "mm"),
                          column_dend_height = unit(30, "mm"),
                          column_dend_reorder = TRUE,
                          col = col_fun2,
                          rect_gp = gpar(col = "white", lwd = 1),
                          top_annotation = ha,
                          left_annotation = rowha)
heatmap_MUTCOF
pdf(file = "output/analysis/lincs/heatmap_MUTCOF_absolute_simplepearson_20190515.pdf", width = 25, height = 22)
print(heatmap_MUTCOF)
dev.off()
