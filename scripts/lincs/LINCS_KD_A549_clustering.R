# setwd("/home/chris/Bureau/sb_cofactor_hr")
setwd("/Users/chris/Desktop/sb_cofactor_hr")

library(factoextra)
library(ComplexHeatmap)
library(circlize)
library(ape)
source("scripts/lincs/lincs.utils.R")

###
col_fun = colorRamp2(c(-2, 0, 2), c("#0f4259", "white", "#800020"))

A549_LINCs_KD_matrix <- readRDS("output/analysis/lincs/A549_LINCs_KD_matrix.rds")
MUTATED_COFACTORS <- get_mutated_cofactors()

rownames(A549_LINCs_KD_matrix) <- A549_LINCs_KD_matrix$Name_GeneSymbol
A549_LINCs_KD_matrix <- A549_LINCs_KD_matrix %>% select(-Name_GeneSymbol)

protKD <- colnames(A549_LINCs_KD_matrix)
mutated <- ifelse(protKD %in% MUTATED_COFACTORS, TRUE, FALSE)
isMutated <- data.frame(protKD, mutated)

mat <- as.matrix(A549_LINCs_KD_matrix)
max(mat)
min(mat)

# Heatmap
# column_names_rot does not exist anymore
Heatmap(mat, name = "LogDiffExp",
        row_names_side = "left",
        row_names_gp = gpar(fontsize = 0),
        row_dend_side = "right",
        show_row_dend = FALSE,
        show_column_names = FALSE,
        column_names_side = "top",
        column_names_gp = gpar(fontsize = 11),
        # column_names_rot = 45,
        column_dend_side = "bottom",
        show_column_dend = FALSE,
        row_dend_width = unit(50, "mm"),
        column_dend_height = unit(50, "mm"),
        col = col_fun,
        rect_gp = gpar(col = "white", lwd = 0))
        
# test method from package factoextra to plot distance matrix
# res.dist <- get_dist(t(A549_LINCs_KD_matrix), stand = TRUE, method = "pearson")
# fviz_dist(res.dist, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

# Dissimilarity measure: Euclidian distance
dmat <- dist(t(mat), method = "euclidian")
dmat.mat <- as.matrix(dmat)

sort(dmat.mat["NIPBL", ])[1:40]
sort(dmat.mat["IGF1R", ])[1:40]
sort(dmat.mat["SCP2", ])[1:40]

clusters <- hclust(dmat)

pdf(file = "output/analysis/lincs/clustering_A549_20190506_v1.pdf", width = 500, height = 500)
plot(clusters)
dev.off()

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

clusterCut100 <- cutree(clusters, 100)
table(clusterCut100, isMutated$mutated)

# http://www.sthda.com/english/wiki/print.php?id=237
# Try correlation-based distance: 
res.cor <- cor(mat, method = "pearson")
d.cor <- as.dist(1 - res.cor)

test <- as.matrix(d.cor)
clusters.cor <- hclust(d.cor)

cut3 <- cutree(clusters.cor, 3)
table(cut3, isMutated$mutated)

cut10 <- cutree(clusters.cor, 10)
table(cut10, isMutated$mutated)
