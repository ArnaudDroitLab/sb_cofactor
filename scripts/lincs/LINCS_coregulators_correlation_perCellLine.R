# setwd("/home/chris/Bureau/sb_cofactor_hr")
setwd("/Users/chris/Desktop/sb_cofactor_hr")

library(ComplexHeatmap)
library(circlize)
source("scripts/lincs/lincs.utils.R")

# Load mutated_cofactors gene symbol
MUTATED_COFACTORS <- get_mutated_cofactors()

# Load all available gene KD and OE signatures from LINCS database
all_signatures_KD <- read.csv("input/lincs/LINCS_consensus_gene_KD_Signatures_all_37275.xls", sep = "\t")
all_signatures_OE <- read.csv("input/lincs/LINCS_gene_overexpression_Signatures_all_9291.xls", sep = "\t")

### Parameters
# Heatmap: correlation color scale
col_fun2 <- colorRamp2(c(-1, 0, 1), c("#0f4259", "white", "#800020"))
# Date
date <- "20190527"
# Map cellLines and timepoint
cellLines <- c("VCAP", "HEK293T")
timeline <- c("96 h", "48 h")
cellLines <- c(cellLines, "A375", "A549", "ASC", "HA1E", "HCC515", "HEKTE", "HEPG2", "HT29", "MCF7", "NPC", "PC3", "SKL", "SW480")
timeline <- c(timeline, rep.int("96 h", 13))
cellLines <- c(cellLines, "PC3", "SHSY5Y", "VCAP")
timeline <- c(timeline, rep.int("120 h", 3))
cellLines <- c(cellLines, "MCF7", "PC3")
timeline <- c(timeline, rep.int("144 h", 2))
#
param <- data.frame(cellLines, timeline)

# Analysis including all 978 genes
for (i in seq(nrow(param))) {
  cLine <- param[i, ] %>% pull(cellLines) %>% as.character
  time <- param[i, ] %>% pull(timeline) %>% as.character
  time_str <- gsub(" ", "", time, fixed = TRUE)
    
  message("#####\t", cLine, " | ", time)
  
  signMat_wGeneName <- get_signMat_KDOE(cellLine = cLine, time = time, KD_matrix = all_signatures_KD, OE_matrix = all_signatures_OE)
  
  saveSignMat(matrix = signMat_wGeneName, cellLine = cLine,
              output_dir = "output/analysis/lincs/coregulators_correlation_perCellLine_978genes",
              output_file = paste0("signMat_978genes", "_", cLine, "_", time_str, "_", date))
  
  signMat <- signMat_wGeneName %>% select(-Name_GeneSymbol)
  
  #
  mat <- as.matrix(signMat)
  max(mat)
  min(mat)
  
  #
  cor.pearson.mutcof <- cor(mat, method = "pearson")
  max(cor.pearson.mutcof - diag(nrow(cor.pearson.mutcof))) # -diag(n) allow to remove the identity matrix
  min(cor.pearson.mutcof)
  
  # annotation KD/OE
  kd_oe_tmp <- strsplit(colnames(cor.pearson.mutcof), "_")
  col_kd_oe <- get_nth_element(kd_oe_tmp, 1)
  ha <- HeatmapAnnotation(Experiment = col_kd_oe,
                          col = list(Experiment = c("KD" = "#16DB93", "OE" = "#EFEA5A")))
  rowha = rowAnnotation(Experiment = col_kd_oe,
                        col = list(Experiment = c("KD" = "#16DB93", "OE" = "#EFEA5A")),
                        show_legend = FALSE)
  
  # heatmap
  heatmap_MUTCOF <- Heatmap(cor.pearson.mutcof, name = "Pearson correlation",
                            # Title
                            column_title = paste0(cLine, " | ", time) ,
                            column_title_gp = gpar(fill = "black", col = "white", fontface = "bold", fontsize = "20"),
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
                            left_annotation = rowha,
                            cell_fun = function(j, i, x, y, width, height, fill) {
                              grid.text(sprintf("%.2f", cor.pearson.mutcof[i, j]), x, y, gp = gpar(fontsize = 10))
                            })
  
  saveHeatmap(heatmap_obj = heatmap_MUTCOF,
              output_dir = "output/analysis/lincs/coregulators_correlation_perCellLine_978genes",
              output_file = paste0("heatmap_MUTCOF_", cLine, "_", time_str, "_simplepearson_", date),
              format = "pdf",
              width = 15, height = 12)
}

# Analysis only with significant differential genes
# NB: A gene is kept in the analysis if its expression is diffenrentially significant in at least one sample

param <- data.frame("cellLines" = "A549", "timeline" = "96 h")

for (i in seq(nrow(param))) {
  cLine <- param[i, ] %>% pull(cellLines) %>% as.character
  time <- param[i, ] %>% pull(timeline) %>% as.character
  time_str <- gsub(" ", "", time, fixed = TRUE)
  
  message("#####\t", cLine, " | ", time)
  
  signMat_wGeneName <- get_signMat_KDOE(cellLine = cLine, time = time, KD_matrix = all_signatures_KD, OE_matrix = all_signatures_OE)
  signMat_pval_wGeneName <- get_signMat_KDOE(cellLine = cLine, time = time, KD_matrix = all_signatures_KD, OE_matrix = all_signatures_OE, pval = TRUE)
  # saveSignMat(matrix = signMat_wGeneName, cellLine = cLine,
  #             output_dir = "output/analysis/lincs/coregulators_correlation_perCellLine_978genes",
  #             output_file = paste0("signMat_978genes", "_", cLine, "_", time_str, "_", date))
  # 
  # signMat <- signMat_wGeneName %>% select(-Name_GeneSymbol)
  # 
  # #
  # mat <- as.matrix(signMat)
  # max(mat)
  # min(mat)
  # 
  # #
  # cor.pearson.mutcof <- cor(mat, method = "pearson")
  # max(cor.pearson.mutcof - diag(nrow(cor.pearson.mutcof))) # -diag(n) allow to remove the identity matrix
  # min(cor.pearson.mutcof)
  # 
  # # annotation KD/OE
  # kd_oe_tmp <- strsplit(colnames(cor.pearson.mutcof), "_")
  # col_kd_oe <- get_nth_element(kd_oe_tmp, 1)
  # ha <- HeatmapAnnotation(Experiment = col_kd_oe,
  #                         col = list(Experiment = c("KD" = "#16DB93", "OE" = "#EFEA5A")))
  # rowha = rowAnnotation(Experiment = col_kd_oe,
  #                       col = list(Experiment = c("KD" = "#16DB93", "OE" = "#EFEA5A")),
  #                       show_legend = FALSE)
  # 
  # # heatmap
  # heatmap_MUTCOF <- Heatmap(cor.pearson.mutcof, name = "Pearson correlation",
  #                           # Title
  #                           column_title = paste0(cLine, " | ", time) ,
  #                           column_title_gp = gpar(fill = "black", col = "white", fontface = "bold", fontsize = "20"),
  #                           row_names_side = "left",
  #                           row_dend_side = "right",
  #                           column_names_side = "top",
  #                           column_names_rot = 45,
  #                           column_dend_side = "bottom",
  #                           row_dend_width = unit(30, "mm"),
  #                           column_dend_height = unit(30, "mm"),
  #                           column_dend_reorder = TRUE,
  #                           col = col_fun2,
  #                           rect_gp = gpar(col = "white", lwd = 1),
  #                           top_annotation = ha,
  #                           left_annotation = rowha,
  #                           cell_fun = function(j, i, x, y, width, height, fill) {
  #                             grid.text(sprintf("%.2f", cor.pearson.mutcof[i, j]), x, y, gp = gpar(fontsize = 10))
  #                           })
  # 
  # saveHeatmap(heatmap_obj = heatmap_MUTCOF,
  #             output_dir = "output/analysis/lincs/coregulators_correlation_perCellLine_978genes",
  #             output_file = paste0("heatmap_MUTCOF_", cLine, "_", time_str, "_simplepearson_", date),
  #             format = "pdf",
  #             width = 15, height = 12)
}

# signMat_pval_wGeneName < 0.05
# t <- signMat_pval_wGeneName < 0.05
# rowSums(t)
# rowSums(t) > 1
# sum(rowSums(t) > 1)