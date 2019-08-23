# setwd("/home/chris/Bureau/sb_cofactor_hr/A549")
# setwd("/Users/chris/Desktop/sb_cofactor_hr/A549")

library(ComplexHeatmap)
library(circlize)
library(knitr)
library(tidyverse)

deg <- readRDS("output/analyses/deg.rds")
kable(deg$deg_numbers)

 raw <- deg$is_deg$FC1
colnames(raw)

FC1 <- raw %>% dplyr::select(-gene_id)
colnames(FC1)

FC1_mat <- as.matrix(FC1)

deg_genes <- rowSums(FC1_mat) != 0
deg_genes_id <- raw$gene_id[deg_genes]
  
deg_genes_mat <- FC1_mat[deg_genes, ]

# deg_genes_mat[deg_genes_mat == 0] <- NA
rownames(deg_genes_mat) <- deg_genes_id

# heatmap
col_fun = colorRamp2(c(0, 1, 2), c("black", "red", "green"))

Heatmap(deg_genes_mat,
        col = col_fun, na_col = "black",
        row_names_side = "left", show_row_names = FALSE,
        column_names_side = "top", column_names_rot = 45,
        row_dend_side = "right",
        column_dend_side = "bottom",
        column_order = colnames(deg_genes_mat),
        # rect_gp = gpar(col = "white", lwd = 0.01),
        )
