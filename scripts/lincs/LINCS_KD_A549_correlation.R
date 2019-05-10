# setwd("/home/chris/Bureau/sb_cofactor_hr")
setwd("/Users/chris/Desktop/sb_cofactor_hr")

library(ComplexHeatmap)
library(circlize)
source("scripts/lincs/lincs.utils.R")

### Define colors for heatmap
col_fun = colorRamp2(c(-2, 0, 2), c("#0f4259", "white", "#800020"))

# Load fold change gene expression values
A549_LINCs_KD_matrix <- readRDS("output/analysis/lincs/A549_LINCs_KD_matrix.rds")

# Load mutated cofactors names
MUTATED_COFACTORS <- get_mutated_cofactors()

# Remove gene name in order to build a matrix
rownames(A549_LINCs_KD_matrix) <- A549_LINCs_KD_matrix$Name_GeneSymbol
A549_LINCs_KD_matrix <- A549_LINCs_KD_matrix %>% select(-Name_GeneSymbol)

# As matrix
mat <- as.matrix(A549_LINCs_KD_matrix)
max(mat)
min(mat)

# Pearson correlation matrix between each KD sample
cor.pearson <- cor(mat, method = "pearson")
max(cor.pearson - diag(nrow(cor.pearson)))
min(cor.pearson)

# Which KD protein are close to the one I'm interested
closeToNIPBL <- sort(cor.pearson["NIPBL", ], decreasing = TRUE)[1:20]
matNIPBL <- as.matrix(A549_LINCs_KD_matrix[, names(closeToNIPBL)])
Heatmap(matNIPBL, name = "LogDiffExp",
        top_annotation = HeatmapAnnotation(correlation = anno_text(round(closeToNIPBL,2))),
        show_row_names = FALSE,
        show_row_dend = FALSE,
        column_names_side = "top",
        cluster_columns = FALSE,
        col = col_fun)

for (cof in MUTATED_COFACTORS) {
  
}

