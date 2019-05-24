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
  
# 
cellLines <- c("A375", "A549", "ASC", "HA1E", "HCC515", "HEKTE", "HEPG2", "HT29", "MCF7", "NPC", "PC3", "SKL", "SW480")

#
for (cLine in cellLines) {
  message("#####\t", cLine)
  # download KD matrix
  dfKD_cLine <- all_signatures_KD %>% filter(CellLine == cLine, TargetGene %in% MUTATED_COFACTORS)
  signIds_KD <- get_signIds(dfKD_cLine, cell_line = cLine, time = "96 h")
  targets_KD <- paste0("KD_", get_target(dfKD_cLine, cell_line = cLine, time = "96 h"))
  KD_mat <- downloadSignatureInBatch(signIds_KD, targets_KD)
  
  # download OE matrix if available and merge with KD matrix
  dfOE_cLine <- all_signatures_OE %>% filter(CellLine == cLine, TargetGene %in% MUTATED_COFACTORS)
  if (nrow(dfOE_cLine) != 0) {
    signIds_OE <- get_signIds(dfOE_cLine, cell_line = cLine, time = "96 h")
    targets_OE <- paste0("OE_", get_target(dfOE_cLine, cell_line = cLine, time = "96 h"))
    OE_mat <- downloadSignatureInBatch(signIds_OE, targets_OE)
    signMat_wGeneName <- merge(KD_mat, OE_mat, by = "Name_GeneSymbol")
  } else {
    signMat_wGeneName <- KD_mat
  }
  
  saveSignMat(matrix = signMat_wGeneName, cellLine = cLine,
              output_dir = "output/analysis/lincs/coregulators_correlation_perCellLine",
              output_file = paste0("signMat", "_", cLine))
  
  signMat <- signMat_wGeneName %>% select(-Name_GeneSymbol)
  
  #
  mat <- as.matrix(signMat)
  max(mat)
  min(mat)
  
  #
  cor.pearson.mutcof <- cor(mat, method = "pearson")
  max(cor.pearson.mutcof - diag(nrow(cor.pearson.mutcof))) # -diag(n) allow to remove the identity matrix
  min(cor.pearson.mutcof)
  
  # correlation color scale
  col_fun2 = colorRamp2(c(-1, 0, 1), c("#0f4259", "white", "#800020"))
  
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
                            column_title = cLine, column_title_gp = gpar(fill = "black", col = "white", fontface = "bold", fontsize = "20"),
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
  
  # saveHeatmap(heatmap_obj = heatmap_MUTCOF,
  #             output_dir = "output/analysis/lincs/coregulators_correlation_perCellLine",
  #             output_file = paste0("heatmap_MUTCOF_", cLine, "_simplepearson_", "20190523"),
  #             format = "pdf",
  #             width = 15, height = 12)
}
