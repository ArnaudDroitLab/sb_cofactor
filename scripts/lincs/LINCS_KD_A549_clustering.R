# setwd("/home/chris/Bureau/sb_cofactor_hr/mESC")
setwd("/Users/chris/Desktop/sb_cofactor_hr/mESC")

library(dplyr)
library(factoextra)
library(ComplexHeatmap)
library(circlize)
source("scripts/ilincs/lincs.utils.R")

col_fun = colorRamp2(c(-2, 0, 2), c("#0f4259", "white", "#800020"))

A549_ILINCs_KD_matrix <- readRDS("output/analysis/ILINCS/A549_ILINCs_KD_matrix.rds")

rownames(A549_ILINCs_KD_matrix) <- A549_ILINCs_KD_matrix$Name_GeneSymbol
A549_ILINCs_KD_matrix <- A549_ILINCs_KD_matrix %>% select(-Name_GeneSymbol)

protKD <- colnames(A549_ILINCs_KD_matrix)
mutated <- ifelse(protKD %in% mutated_cofactors, TRUE, FALSE)
isMutated <- data.frame(protKD, mutated)

mat <- as.matrix(A549_ILINCs_KD_matrix)
max(mat)
min(mat)

# Heatmap
Heatmap(mat, name = "LogDiffExp",
        row_names_side = "left",
        row_names_gp = gpar(fontsize = 0),
        row_dend_side = "right",
        show_row_dend = FALSE,
        show_column_names = FALSE,
        column_names_side = "top",
        column_names_gp = gpar(fontsize = 11),
        column_names_rot = 45,
        column_dend_side = "bottom",
        show_column_dend = FALSE,
        row_dend_width = unit(50, "mm"),
        column_dend_height = unit(50, "mm"),
        col = col_fun,
        rect_gp = gpar(col = "white", lwd = 0))
        
# test method from package factoextra to plot distance matrix
# res.dist <- get_dist(t(A549_ILINCs_KD_matrix), stand = TRUE, method = "pearson")
# fviz_dist(res.dist, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

#
dmat <- dist(t(mat))
clusters <- hclust(dmat)
plot(clusters)

clusterCut3 <- cutree(clusters, 3)
table(clusterCut3, isMutated$mutated)

clusterCut4 <- cutree(clusters, 4)
table(clusterCut4, isMutated$mutated)

clusterCut5 <- cutree(clusters, 5)
table(clusterCut5, isMutated$mutated)

clusterCut7 <- cutree(clusters, 7)
table(clusterCut7, isMutated$mutated)

clusterCut10 <- cutree(clusters, 10)
table(clusterCut10, isMutated$mutated)

clusterCut20 <- cutree(clusters, 20)
table(clusterCut20, isMutated$mutated)

clusterCut30 <- cutree(clusters, 30)
table(clusterCut30, isMutated$mutated)
