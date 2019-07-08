setwd("/home/chris/Bureau/sb_cofactor_hr/A549")
# setwd("/Users/chris/Desktop/sb_cofactor_hr/A549")

library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(purrr)

# Load raw data GSE104714
raw <- read_tsv("input/GSE104714/GSE104714_DEX_EtOH_expression.tsv") %>% dplyr::select(-X2)
dim(raw)
names(raw)

# Make matrix for correlation
mat <- raw %>% dplyr::select(-transcript) %>% as.matrix
rownames(mat) <- raw$transcript

mat <- mat[rowSums(mat) > 50,]

# Annotation for heatmap > condition: EtOH or DEX
annot_condition <- strsplit(colnames(mat), split = "_") %>% purrr::map(1) %>% unlist
ha <- HeatmapAnnotation(Condition = annot_condition,
                        col = list(Condition = c("EtOH" = "#16DB93", "DEX" = "#EFEA5A")))
rowha = rowAnnotation(Condition = annot_condition,
                      col = list(Condition = c("EtOH" = "#16DB93", "DEX" = "#EFEA5A")),
                      show_legend = FALSE)

# Correlation analysis : Pearson method
cor_pearson <- cor(mat, method = "pearson")
min(cor_pearson)
col_pearson <- colorRamp2(c(0.8, 0.9, 1), c("#0f4259", "white", "#800020"))

cor_pearson_heatmap <- Heatmap(cor_pearson, name = "Pearson correlation",
                               row_names_side = "right",
                               row_dend_side = "right",
                               row_dend_width = unit(25, "mm"),
                               column_names_side = "top", column_names_rot = 45,
                               show_column_dend = FALSE,
                               top_annotation = ha,
                               right_annotation = rowha,
                               rect_gp = gpar(col = "white", lwd = 0.5),
                               col = col_pearson)
cor_pearson_heatmap

# Correlation analysis : Spearman method
cor_spearman <- cor(mat, method = "spearman")
min(cor_spearman)
col_spearman <- colorRamp2(c(0.9, 0.95, 1), c("#0f4259", "white", "#800020"))

cor_spearman_heatmap <- Heatmap(cor_spearman, name = "Spearman correlation",
                               row_names_side = "right",
                               row_dend_side = "right",
                               row_dend_width = unit(25, "mm"),
                               column_names_side = "top", column_names_rot = 45,
                               show_column_dend = FALSE,
                               top_annotation = ha,
                               right_annotation = rowha,
                               rect_gp = gpar(col = "white", lwd = 0.5),
                               col = col_spearman)
cor_spearman_heatmap
