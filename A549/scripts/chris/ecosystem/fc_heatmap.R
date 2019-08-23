# setwd("/home/chris/Bureau/sb_cofactor_hr/A549")
# setwd("/Users/chris/Desktop/sb_cofactor_hr/A549")

library(ComplexHeatmap)
library(circlize)
library(knitr)
library(tidyverse)

deg <- readRDS("output/analyses/deg.rds")
kable(deg$deg_numbers)

reg_genes <- c(deg$gene_list$FC1$upreg, deg$gene_list$FC1$downreg)

raw <- read_delim("results/a549_dex_time_points/FC_mat.csv", delim = ",")
raw2 <- raw %>% dplyr::filter(gene_id %in% reg_genes)
FC <- raw2 %>% dplyr::select(-gene_id)
colnames(FC)

FCmat <- as.matrix(FC)

# heatmap
col_fun = colorRamp2(c(-2, -1, 0, 1, 2), c("red", "red", "white", "green","green"))

Heatmap(FCmat,
        col = col_fun, na_col = "black",
        row_names_side = "left", show_row_names = FALSE,
        column_names_side = "top", column_names_rot = 45,
        row_dend_side = "right",
        column_dend_side = "bottom",
        column_order = colnames(FCmat),
        # rect_gp = gpar(col = "white", lwd = 0.01),
)