# setwd("/home/chris/Bureau/sb_cofactor_hr/mESC")
setwd("/Users/chris/Desktop/sb_cofactor_hr/mESC")

mutated_cofactors <- c("NIPBL", "SMC1A", "SMC3", "RAD21", "HDAC8", "BRD4", "CREBBP", "EP300", "KAT6A", "KAT6B",
                       "SRCAP", "BRPF1", "RAI1", "MBD5", "EHMT1", "HDAC4", "ARID1A", "ARID1B", "SMARCA4", "SMARCB1",
                       "SMARCE1", "ANKRD11", "KMT2A", "SMARCA2", "AFF4", "KMT2D", "KDM6A", "TAF6", "ESCO2", "CHD7",
                       "DDX11", "STAG1", "STAG2", "SETD5", "MED12", "MED13L", "MED17", "MED23", "MED25")

A549_ILINCs_KD_matrix <- readRDS("output/analysis/ILINCS/A549_ILINCs_KD_matrix.rds")

rownames(A549_ILINCs_KD_matrix) <- A549_ILINCs_KD_matrix$Name_GeneSymbol
A549_ILINCs_KD_matrix <- A549_ILINCs_KD_matrix %>% select(-Name_GeneSymbol)

a549.kd.pca <- prcomp(A549_ILINCs_KD_matrix, center = TRUE, scale = TRUE)
summary(a549.kd.pca)
str(a549.kd.pca)

group.cofactors <- ifelse(rownames(A549_ILINCs_KD_matrix) %in% mutated_cofactors, "YES", "NO")

library(ggbiplot)
ggbiplot(a549.kd.pca)
ggbiplot(a549.kd.pca, labels=rownames(A549_ILINCs_KD_matrix))
ggbiplot(a549.kd.pca, ellipse = TRUE, labels=rownames(A549_ILINCs_KD_matrix), groups = group.cofactors)
ggbiplot(a549.kd.pca, ellipse = TRUE, circle = TRUE, labels=rownames(A549_ILINCs_KD_matrix), groups = group.cofactors)

ggbiplot(a549.kd.pca, ellipse = TRUE, obs.scale = 1, var.scale = 1, labels=rownames(A549_ILINCs_KD_matrix), groups = group.cofactors)

ggbiplot(a549.kd.pca, ellipse = TRUE, var.axes = FALSE, labels=rownames(A549_ILINCs_KD_matrix), groups = group.cofactors)
ggbiplot(a549.kd.pca, ellipse = TRUE, var.axes = FALSE, groups = group.cofactors)
