# setwd("/home/chris/Bureau/sb_cofactor_hr")
setwd("/Users/chris/Desktop/sb_cofactor_hr")

library(ComplexHeatmap)
library(circlize)
source("scripts/lincs/lincs.utils.R")

# Load fold change gene expression values
A549_LINCs_KD_matrix <- readRDS("output/analysis/lincs/LINCS_KD_matrix/A549_96h_LINCS_KD_matrix.rds")

# Remove gene name in order to build a matrix
rownames(A549_LINCs_KD_matrix) <- A549_LINCs_KD_matrix$Name_GeneSymbol
A549_LINCs_KD_matrix <- A549_LINCs_KD_matrix %>% select(-Name_GeneSymbol)

# As matrix
mat <- as.matrix(A549_LINCs_KD_matrix)
max(mat)
min(mat)

###
CONTROL_GENES <- c("MED1", "MED12L", "MED15","MED21", "MED26", "MED28", "MED4", "MED6", "MED7",
                   "SMC3", "SMC4",
                   "SMARCA2", "SMARCA4", "SMARCA5", "SMARCAD1", "SMARCB1", "SMARCC1", "SMARCC2", "SMARCD2", "SMARCE1",
                   "ELL3", "AFF4",
                   "POLA1", "POLA2", "POLB", "POLD4", "POLE2", "POLE3", "POLG", "POLQ", "POLR1A", "POLR1C",
                   "POLR2A", "POLR2C", "POLR2D", "POLR2E", "POLR2F", "POLR2H", "POLR2I", "POLR2K",
                   "POLR3B", "POLR3C", "POLR3D", "POLR3E", "POLR3F", "POLR3K")
CONTROL_GENES_A549 <- CONTROL_GENES[CONTROL_GENES %in% colnames(A549_LINCs_KD_matrix)]

### New matrix
mat_control <- mat[, CONTROL_GENES_A549]

# correlation color scale
col_fun2 = colorRamp2(c(-1, 0, 1), c("#0f4259", "white", "#800020"))

# Pearson correlation matrix (signed values)
cor.pearson.control <- cor(mat_control, method = "pearson")
max(cor.pearson.control - diag(nrow(cor.pearson.control))) # -diag(n) allow to remove the identity matrix
min(cor.pearson.control)

# heatmap
heatmap_control <- Heatmap(cor.pearson.control, name = "Pearson correlation",
                          row_names_side = "left",
                          row_dend_side = "right",
                          column_names_side = "top",
                          column_names_rot = 45,
                          column_dend_side = "bottom",
                          row_dend_width = unit(30, "mm"),
                          column_dend_height = unit(30, "mm"),
                          column_dend_reorder = TRUE,
                          col = col_fun2,
                          rect_gp = gpar(col = "white", lwd = 1))
pdf(file = "output/analysis/lincs/heatmap_CONTROLA549_simplepearson_20190522.pdf", width = 25, height = 22)
print(heatmap_control)
dev.off()

### Pearson correlation matrix (absolute values)
cor.pearson.mutcof.abs <- cor(abs(mat_control), method = "pearson")
max(cor.pearson.mutcof.abs - diag(nrow(cor.pearson.mutcof.abs))) # -diag(n) allow to remove the identity matrix
min(cor.pearson.mutcof.abs)

heatmap_control.abs <- Heatmap(cor.pearson.mutcof.abs, name = "Pearson correlation",
                          row_names_side = "left",
                          row_dend_side = "right",
                          column_names_side = "top",
                          column_names_rot = 45,
                          column_dend_side = "bottom",
                          row_dend_width = unit(30, "mm"),
                          column_dend_height = unit(30, "mm"),
                          column_dend_reorder = TRUE,
                          col = col_fun2,
                          rect_gp = gpar(col = "white", lwd = 1))
pdf(file = "output/analysis/lincs/heatmap_CONTROLA549_absolute_simplepearson_20190522.pdf", width = 25, height = 22)
print(heatmap_control.abs)
dev.off()